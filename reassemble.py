#!/usr/bin/python

#This code is Copyright (c) the Centre for Proteomic and Genomic Research
#(2011-2012) and is distributed under the terms of Gnu General Public License
#version 3 or higher.  For the full terms and conditions please see the file
#COPYING which should have been included in this distribution.

import sys
import re

titles = []
sequences = []
current = {}

if len(sys.argv) == 1 or sys.argv[1] == "-":
	apf = sys.stdin
else:
	apf = open(sys.argv[1], "r")

line = apf.readline()
while line != "":
	if re.match("^\d+:", line):
		sequences.append("")
		titles.append(line.split(":",2)[1])
	else:
		fragment, poslist = line.split(None, 1)
		current[fragment] = [(int(x.split(":")[0])-1, int(x.split(":")[1])) for x in poslist.split(',')]
	line = apf.readline()

for f in current:
	for seq, pos in current[f]:
		if pos > len(sequences[seq]):
			sequences[seq] = sequences[seq]+(" "*(pos+len(f)-len(sequences[seq])))
		sequences[seq] = sequences[seq][0:pos]+f+sequences[seq][pos+len(f):]

for i in range(len(titles)):
	print ">"+titles[i].rstrip()
	print sequences[i]
