#!/usr/bin/python

#This code is Copyright (c) the Centre for Proteomic and Genomic Research
#(2011-2012) and is distributed under the terms of Gnu General Public License
#version 3 or higher.  For the full terms and conditions please see the file
#COPYING which should have been included in this distribution.

import sys
import re
import copy

#print a usage message and exit
def usage():
	print """Pepfrag produces a list of peptide fragments providing complete coverage of the
input sequences for spotting onto a microarray slide. 

pepfrag [-l n] [-L n] [-o n] [-O n] [-C n] [-p filename] [-t] [sequence file]

    -l number
        Minimum length of output fragments (default 8)
    -L number
        Maximum length of output fragments (default 30)
    -o number
        Minimum amount of overlap between adjacent fragments (default 2)
    -O number
        Maximum amount of overlap between adjacent fragments (default 4)
    -C number
        Approximate chunk size. Input sequences are split into chunks of
        approximately this length to make computation more tractable. Actual
        size of the chunks is heuristically determined by computing the
        frequency of repetition of potential fragments surrounding the putative
        split point(s) and choosing the points that breaks the least repeated 
        fragment. Larger chunk sizes produce better results, but increase
        computation time exponentially. As a rule of thumb make chunk size a 
        multiple of max length - min overlap. The default value for chunk size
        is 4 times this value, i.e. 4*(maxlen-minoverlap)
    -p filename
        Accept previous results and merge them with those from the current
        input sequences. 'filename' can be '-' to indicate standard input for
        piping purposes. The previous results input stream can be the same as
        the sequence stream (i.e. they can be the same file, or both be on 
        standard input) but then previous results must come first. Previous 
        results must provided in a format identical to the output of pepfrag
        (see Understanding the output below).
    -t
        Do not post-process data using naive algorithm, even if it would 
        produce better results, i.e. always use the pepfrag algorithm. Used 
        mainly for testing and algorithm comparision purposes.
    -n
        Only output the results of the naive algorithm. Used mainly for
        testing and algorithm comparison purposes.
    sequence file
        A file containing the sequence(s) to process in FASTA format. Multiple
        sequences are allowed, and their results will be merged. If not 
        specified, the sequences will be read from standard input

Understanding the output

Output is split into two sections. The first section is an indexing of input
sequences and the second section is a list of fragments and their starting
positions within the indexed sequences. The sequence index takes the form of an
index number, followed by a colon and a space, followed by the fasta title of
the sequence.

The second section lists each fragment on a single line followed by a tab,
followed by a comma separated list of sequence index,position pairs. The number
before the comma indicates which sequence the fragment was found in, and the
number after the comma indicates where in the sequence the fragment can be
found (indexed from 0).

As an example:

--- Example Input File ---
>example test sequence
ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ

> bogus test sequence 2
ABCDEFGHIJKLMNOPQRSTuVWXYZABCDEFGHIJKLMNOpQRSTUVWXYZ
--- End ---

--- Example output ---
1: example test sequence
2: bogus test sequence 2
ABCDEFGHIJKLMNO 1:0, 1:26, 2:0, 2:26
NOPQRSTUVWXYZAB 1:13
NOpQRSTUVWXYZ   2:39
NOPQRSTuVWXYZAB 2:13
NOPQRSTUVWXYZ   1:39
--- End ---
"""
	sys.exit()

#The brute force algortithm to analyse a single chunk, fragment set tester. We
#test each set in the space using a depth first traversal starting #rom
#maxlen,minoverlap to minlen,maxoverlap; There will be fewer fragments for
#maxlen,minoverlap so evaluating this set first sets a good basline for
#pruning. Any other sets that are partially evaluated and grow in size past
#this baseline are discarded (pruned)

def bruteforce(pos, seq, basepos):
	global best, current, titleidx
	stack = []
	jump_to_return = False

	while True:
		#START
		jump_to_start = False
		i = maxlen
		while (i > minlen and not jump_to_start) or jump_to_return:
			if not jump_to_return:
				realminoverlap = min(minoverlap,max(0,pos-minoverlap))
				realmaxoverlap = min(maxoverlap,max(0,pos-maxoverlap))
				j = realminoverlap
			while (j <= realmaxoverlap and not jump_to_start) or jump_to_return:
				if not jump_to_return:
					if seqlen - (pos-j+1) < minlen:
						jump_to_return = True
					elif pos-j+i >= seqlen:
						if seq[pos-j:] in current:
							current[seq[pos-j:]].append((titleidx, basepos+pos-j))
						else:
							if best is not None and (len(current)+1 > len(best) or len(current)+1 > nbestlen):
								jump_to_return = True
							else:
								current[seq[pos-j:]] = [(titleidx, basepos+pos-j)]
						if not jump_to_return:
							if best is None or len(current) < len(best):
								best = copy.deepcopy(current)
							current[seq[pos-j:]].pop()
							if len(current[seq[pos-j:]]) == 0:
								del current[seq[pos-j:]]
					else:
						if seq[pos-j:pos-j+i] in current:
							current[seq[pos-j:pos-j+i]].append((titleidx, basepos+pos-j))
						else:
							if best is not None and (len(current)+1 > len(best) or len(current)+1 > nbestlen):
								jump_to_return = True
							else:
								current[seq[pos-j:pos-j+i]] = [(titleidx, basepos+pos-j)]
						if not jump_to_return:
							stack.append((pos,seq,basepos,i,j,realminoverlap,realmaxoverlap))
							pos = pos-j+i
							jump_to_return = False
							jump_to_start = True
				#RETURN
				if jump_to_return:
					if len(stack) == 0:
						return
					jump_to_return = False
					(pos,seq,basepos,i,j,realminoverlap,realmaxoverlap) = stack.pop()
					current[seq[pos-j:pos-j+i]].pop()
					if len(current[seq[pos-j:pos-j+i]]) == 0:
						del current[seq[pos-j:pos-j+i]]
				j += 1
			i -= 1
		if not jump_to_start:
			jump_to_return = True

#Counts the number of occurences of needle in haystack
def find_all_count(needle, haystack):
	return len([m.start() for m in re.finditer(needle, haystack)])

#Setup up variables and default values
sequences = []
best = None
current = {}
minlen = 8
maxlen = 30
minoverlap = 2
maxoverlap = 4
chunklength = 4*(maxlen-minoverlap)+minoverlap
accept_prev = None
usenaive = None
filename = None
apf = None
titleidx = 1

i = 1
while i < len(sys.argv):
	if sys.argv[i] == "-l":
		i += 1
		minlen = int(sys.argv[i])
	elif sys.argv[i] == "-L":
		i += 1
		maxlen = int(sys.argv[i])
	elif sys.argv[i] == "-o":
		i += 1
		minoverlap = int(sys.argv[i])
	elif sys.argv[i] == "-O":
		i += 1
		maxoverlap = int(sys.argv[i])
	elif sys.argv[i] == "-C":
		i += 1
		chunklength = int(sys.argv[i])
	elif sys.argv[i] == "-p":
		i += 1
		accept_prev = sys.argv[i]
	elif sys.argv[i] == "-t":
		usenaive = False
	elif sys.argv[i] == "-n":
		usenaive = True
	else:
		if sys.argv[i].startswith("-"):
			usage()
		filename = sys.argv[i]
	i += 1

if chunklength < minlen:
	sys.stderr.write("Chunk length cannot be less the minimum fragment length\n")
	sys.exit(1)

if maxlen < minlen:
	sys.stderr.write("Maximum fragment length cannot be less than minimum fragment length\n")
	sys.exit(1)
	
if maxoverlap < minoverlap:
	sys.stderr.write("Maximum overlap length cannot be less than minimum overlap length\n")
	sys.exit(1)

if filename is not None:
	f = open(filename, "r")
else:
	f = sys.stdin

if accept_prev is not None:
	if accept_prev == "-":
		if f is sys.stdin:
			apf = f
		else:
			apf = sys.stdin
	else:
		if accept_prev == filename:
			apf = f
		else:
			apf = open(accept_prev, "r")

	line = apf.readline()
	while line.strip() != "" and not line.startswith(">"):
		if re.match("^\d+:", line):
			print line.strip()
			titleidx += 1
		else:
			fragment, poslist = line.split(None, 1)
			current[fragment] = [(int(x.split(":")[0]), int(x.split(":")[1])) for x in poslist.split(',')]
		line = apf.readline()

if apf is not f:
	line = f.readline()
while line != "":
	while line.strip() == "":
		line = f.readline()
	if not line.startswith(">"):
		sys.stderr.write("Sequence input not in valid fasta format\n")
		sys.exit(1)
	title = line[1:].strip()
	seq = ""
	line = f.readline()
	while line.strip() != "" and not line.startswith(">"):
		seq += line.strip()
		line = f.readline()
	sequences.append((seq, title))

#pre calculate baseline using naive filtered method
if usenaive is not False:
	nbest = copy.deepcopy(current)
	ntitleidx = titleidx
	for seq, title in sequences:
		if usenaive is True:
			print "{0}: {1}".format(ntitleidx, title)
		for i in range(0, len(seq), maxlen-minoverlap):
			if seq[i:i+maxlen] in nbest:
				nbest[seq[i:i+maxlen]].append((ntitleidx,i))
			else:
				nbest[seq[i:i+maxlen]] = [(ntitleidx,i)]
		ntitleidx += 1
	nbestlen = len(nbest)
else:
	nbestlen = sys.maxint

if usenaive is not True:
	for seq, title in sequences:
		print "{0}: {1}".format(titleidx, title)
		chunk = []
		seqlen = len(seq)
		if chunklength > seqlen:
			chunk.append(seq)
		else:
			startidx = 0
			while startidx+chunklength < seqlen:
				min_penalty = 0
				penalties = [0]*(4*minlen+1)
				for i in range(4*minlen+1):
					penalty = find_all_count(re.escape(seq[startidx+chunklength-(2*minlen)+i:startidx+chunklength-minlen+i]), seq)
					for j in range(min(minlen, (4*minlen+1)-i)):
						penalties[i+j] += penalty
						if min_penalty < penalties[i+j]:
							min_penalty = penalties[i+j]
				min_penalty_idx = minlen
				for i in range(minlen, 3*minlen+1):
					if penalties[i] < min_penalty:
						min_penalty_idx = i
						min_penalty = penalties[i]
				min_penalty_idx = startidx+chunklength-(2*minlen)+min_penalty_idx
				chunk.append(seq[startidx:min_penalty_idx])
				startidx = min_penalty_idx-(minoverlap-1)
			chunk.append(seq[startidx:])

		if len(chunk[-1]) < minlen:
			chunk[-2] += chunk[-1]
			del chunk[-1]

		basepos = 0
		for c in chunk:
			seqlen = len(c)
			best = None
			bruteforce(0,c,basepos)
			current = copy.deepcopy(best)
			basepos += len(c)-minoverlap+1
		
		titleidx += 1

if usenaive is True or (usenaive is None and (best is None or len(best) > len(nbest))):
	best = nbest

for fragment in best:
	print "{0}\t{1}".format(fragment, ", ".join([str(t)+":"+str(p) for t, p in best[fragment]]))
