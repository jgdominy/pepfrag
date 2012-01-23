#!/usr/bin/python

#This code is Copyright (c) the Centre for Proteomic and Genomic Research
#(2011-2012) and is distributed under the terms of Gnu General Public License
#version 3 or higher.  For the full terms and conditions please see the file
#COPYING which should have been included in this distribution.

import sys
import random
import os

aa_upper = "ARNDCEQGHILKMFPSTWYVUO"
aa_lower = "arndceqghilkmfpstwyvuo"

numsubs = 0
numindels = 0
numboth = 0
numseqs = 0
files = []

def get_first_sequence(fname):
	f = open(fname, "r")
	line = f.readline()
	while line != "" and line == "\n":
		line = f.readline()
	if not line.startswith(">"):
		print >> sys.stderr, "{} not a fasta file".format(fname)
		return None
	else:
		title = line[1:]
		seq = ""
		line = f.readline()
		while line != "" and line != "\n":
			seq += line.strip()
			line = f.readline()
		f.close()
		return seq, title

i = 1
while i < len(sys.argv):
	if sys.argv[i].startswith("-"):
		if sys.argv[i][1] == "s":
			i += 1
			numsubs = int(sys.argv[i])
		elif sys.argv[i][1] == "i":
			i += 1
			numindels = int(sys.argv[i])
		elif sys.argv[i][1] == "m":
			i += 1
			numboth = int(sys.argv[i])
		elif sys.argv[i][1] == "r":
			i += 1
			numseqs = int(sys.argv[i])
		else:
			print >> sys.stderr, "Command line switch {} not recognised".format(sys.argv[i])
	else:
		files.append(sys.argv[i])
	i += 1

for filename in files:
	try:
		seq, title = get_first_sequence(filename)
		basename, ext = os.path.splitext(filename)
		if seq[0].isupper():
			aa_list = aa_upper
		else:
			aa_list = aa_lower

		if numsubs > 0:
			f = open(basename+"_subs"+ext, "w")
			print >> f, ">"+title.rstrip()
			print >> f, seq.rstrip()
			for i in range(numseqs):
				print >> f, ">{}| with {} random substitutions | {}".format(title.rstrip(), numsubs, i+1)
				used = []
				outseq = seq
				for j in range(numsubs):
					pos = random.randint(0, len(seq))
					while pos in used:
						pos = random.randint(0, len(seq))
					outseq = outseq[0:pos]+random.choice(aa_list)+outseq[pos+1:]
				print >> f, outseq
			f.close()
		
		if numindels > 0:
			f = open(basename+"_indels"+ext, "w")
			print >> f, ">"+title.rstrip()
			print >> f, seq.rstrip()
			for i in range(numseqs):
				print >> f, ">{}| with {} random indels | {}".format(title.rstrip(), numsubs, i+1)
				used = []
				outseq = seq
				for j in range(numindels):
					pos = random.randint(0, len(seq))
					while pos in used:
						pos = random.randint(0, len(seq))
					if random.randint(0,1) == 0:
						outseq = outseq[0:pos]+outseq[pos+1:]
					else:
						outseq = outseq[0:pos]+random.choice(aa_list)+outseq[pos:]
				print >> f, outseq
			f.close()

		if numboth > 0:
			f = open(basename+"_muts"+ext, "w")
			print >> f, ">"+title.rstrip()
			print >> f, seq.rstrip()
			for i in range(numseqs):
				print >> f, ">{}| with {} random mutations | {}".format(title.rstrip(), numsubs, i+1)
				used = []
				outseq = seq
				for j in range(numboth):
					pos = random.randint(0, len(seq))
					while pos in used:
						pos = random.randint(0, len(seq))
					action = random.randint(0,2)
					if action == 0:
						outseq = outseq[0:pos]+random.choice(aa_list)+outseq[pos+1:]
					elif action == 1:
						outseq = outseq[0:pos]+outseq[pos+1:]
					else:
						outseq = outseq[0:pos]+random.choice(aa_list)+outseq[pos:]
				print >> f, outseq
			f.close()
	except TypeError:
		pass	
