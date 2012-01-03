#!/usr/bin/python

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
