/* This code is Copyright (c) the Centre for Proteomic and Genomic Research
 * (2011-2012) and is distributed under the terms of Gnu General Public License
 * version 3 or higher.  For the full terms and conditions please see the file
 * COPYING which should have been included in this distribution. */

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>

using namespace std;

//maps of vectors of seqpos together form the main data structure used to
//store a list of peptide fragments found/used and their respective positions

//a type for storing fragment positions within multiple sequences
typedef struct seqpos {
	int idx; //the index of the sequence of this position
	int pos; //the index of the starting position of the fragment within the sequence
} seqpos;

//a type for storing explicit stack frames within the depth first search
typedef struct frame {
	unsigned int pos;
	unsigned int l;
	unsigned int o;
	string F;
	seqpos tspos;
	frame *prev;
} frame;

//a struture to keep a record of sequences read in from file, and their titles
typedef struct fastaseq {
	string title;
	string sequence;
} fastaseq;

//Setup up variables and default values
map <string, vector<seqpos> > best;
map <string, vector<seqpos> > current;
map <string, vector<seqpos> > nbest;
unsigned int minlen = 8;
unsigned int maxlen = 30;
unsigned int minoverlap = 2;
unsigned int maxoverlap = 4;
unsigned int titleidx = 1;
unsigned int nbestlen;
unsigned int seqlen;
unsigned int basepos;

//maths convenience function
int max(int a, int b) {
	return (a > b?a:b);
}
	
int min(int a, int b) {
	return (a < b?a:b);
}

//string convenience functions

//checks whether a string consists entirely of white space
int allspace(string s) {
	unsigned int i;
	i = 0;
	while (s[i] != 0) {
		if (!isspace(s[i])) return 0;
		i++;
	}
	return 1;
}

int isnumeric(string s) {
	unsigned int i;
	i = 0;
	while (s[i] != 0) {
		if (!isdigit(s[i])) return 0;
		i++;
	}
	return 1;
}

//strips white space from either end of the string
char *strip(const char *s) {
	unsigned int i, j, len;

	len = strlen(s);
	i = 0;
	while ((i < len)&&(isspace(s[i]))) i++;
	j = len-1;
	while ((j > 0)&&(isspace(s[j]))) j--;
	return strndup(s+i,j-i+1);
}

//print a usage message and exit
void usage() {
	cerr << "Pepfrag produces a list of peptide fragments providing complete coverage of the\n"
		"input sequences for spotting onto a microarray slide. \n"
		"\n"
		"pepfrag [-l n] [-L n] [-o n] [-O n] [-C n] [-p filename] [-t] [sequence file]\n"
		"\n"
		"    -l number\n"
		"        Minimum length of output fragments (default 8)\n"
		"    -L number\n"
		"        Maximum length of output fragments (default 30)\n"
		"    -o number\n"
		"        Minimum amount of overlap between adjacent fragments (default 2)\n"
		"    -O number\n"
		"        Maximum amount of overlap between adjacent fragments (default 4)\n"
		"    -C number\n"
		"        Approximate chunk size. Input sequences are split into chunks of\n"
		"        approximately this length to make computation more tractable. Actual\n"
		"        size of the chunks is heuristically determined by computing the\n"
		"        frequency of repetition of potential fragments surrounding the putative\n"
		"        split point(s) and choosing the points that breaks the least repeated \n"
		"        fragment. Larger chunk sizes produce better results, but increase\n"
		"        computation time exponentially. As a rule of thumb make chunk size a \n"
		"        multiple of max length - min overlap. The default value for chunk size\n"
		"        is 4 times this value, i.e. 4*(maxlen-minoverlap)\n"
		"    -p filename\n"
		"        Accept previous results and merge them with those from the current\n"
		"        input sequences. 'filename' can be '-' to indicate standard input for\n"
		"        piping purposes. The previous results input stream can be the same as\n"
		"        the sequence stream (i.e. they can be the same file, or both be on \n"
		"        standard input) but then previous results must come first. Previous \n"
		"        results must provided in a format identical to the output of pepfrag\n"
		"        (see Understanding the output below).\n"
		"    -t\n"
		"        Do not post-process data using naive algorithm, even if it would \n"
		"        produce better results, i.e. always use the pepfrag algorithm. Used \n"
		"        mainly for testing and algorithm comparision purposes.\n"
		"    -n\n"
		"        Only output the results of the naive algorithm. Used mainly for\n"
		"        testing and algorithm comparison purposes.\n"
		"    sequence file\n"
		"        A file containing the sequence(s) to process in FASTA format. Multiple\n"
		"        sequences are allowed, and their results will be merged. If not \n"
		"        specified, the sequences will be read from standard input\n"
		"\n"
		"Understanding the output\n"
		"\n"
		"Output is split into two sections. The first section is an indexing of input\n"
		"sequences and the second section is a list of fragments and their starting\n"
		"positions within the indexed sequences. The sequence index takes the form of an\n"
		"index number, followed by a colon and a space, followed by the fasta title of\n"
		"the sequence.\n"
		"\n"
		"The second section lists each fragment on a single line followed by a tab,\n"
		"followed by a comma separated list of sequence index,position pairs. The number\n"
		"before the comma indicates which sequence the fragment was found in, and the\n"
		"number after the comma indicates where in the sequence the fragment can be\n"
		"found (indexed from 0).\n"
		"\n"
		"As an example:\n"
		"\n"
		"--- Example Input File ---\n"
		">example test sequence\n"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ\n"
		"\n"
		"> bogus test sequence 2\n"
		"ABCDEFGHIJKLMNOPQRSTuVWXYZABCDEFGHIJKLMNOpQRSTUVWXYZ\n"
		"--- End ---\n"
		"\n"
		"--- Example output ---\n"
		"1: example test sequence\n"
		"2: bogus test sequence 2\n"
		"ABCDEFGHIJKLMNO 1:0, 1:26, 2:0, 2:26\n"
		"NOPQRSTUVWXYZAB 1:13\n"
		"NOpQRSTUVWXYZ   2:39\n"
		"NOPQRSTuVWXYZAB 2:13\n"
		"NOPQRSTUVWXYZ   1:39\n"
		"--- End ---\n";
	exit(EXIT_FAILURE);
}

//A recursive brute force depth first search of the fragment set space
//implemented using an explicit stack (S). Each set in the space is tested
//using a depth first traversal starting from maxlen,minoverlap to
//minlen,maxoverlap; There will be fewer fragments for maxlen,minoverlap so
//evaluating this set first sets a good basline for pruning. Any other sets
//that are partially evaluated and grow in size past this baseline are
//discarded (pruned)
void dfsnr(const string seq, unsigned int pos = 0) {
	unsigned int l, o;
	string F;
	seqpos tspos;
	frame *S = NULL;
	frame *temp = NULL;

	start:
	tspos.idx = titleidx;
	tspos.pos = pos+basepos;
	if ((best.size() == 0)||(current.size()+1 < best.size())) { //prune based on set size
		for (l = maxlen; l >= minlen; l--) {
			F = seq.substr(pos,l);
			current[F].push_back(tspos);
			for (o = minoverlap; o <= maxoverlap && l > o; o++) {
				if (seq.length() - (pos+l-o) >= minlen) {//is there a potential to recurse, i.e. if the substring starting at pos+l-o still greater than or equal minlen
					temp = new frame; temp->pos = pos; temp->l = l; temp->o = o; temp->F = F; temp->tspos = tspos; temp->prev = S; S = temp; //push stack
					pos += l-o;
					goto start; //recurse
					retpos:
					pos = S->pos; l = S->l; o = S->o; F = S->F; tspos = S->tspos; temp = S; S = S->prev; delete temp; //pop stack
				} else {
					if (seq.length() - pos <= maxlen) {//is the remainer of the sequence small enough to append to the last fragment
						current[F].pop_back(); if (current[F].size() == 0) current.erase(F);
						F = seq.substr(pos);
						current[F].push_back(tspos);
						best.clear();
						best = current;
					}
				}
			}
			current[F].pop_back(); if (current[F].size() == 0) current.erase(F);
		}
	}
	if (S == NULL) return;
	goto retpos;
}

//Counts the number of occurences of needle in haystack
int find_all_count(string needle, string haystack) {
	unsigned int i, c = 0, hlen = haystack.length(), nlen = needle.length();
	for (i = 0; i < hlen-nlen; i++) {
		if (haystack.substr(i,nlen) == needle) c++;
	}
	return c;
}

char *readfile(FILE *f, unsigned int& filesize) {
	char *retbuf = NULL, *tmpbuf = NULL;
	unsigned int bufsize = 4096, bytesread = 0;
	
	filesize = 0;
	retbuf = new char[bufsize];
	bytesread = fread(retbuf, 1, 4096, f);
	do {
		filesize += bytesread;
		if (filesize >= bufsize) {
			bufsize += 4096;
			tmpbuf = new char[bufsize];
			memcpy(tmpbuf, retbuf, bufsize-4096);
			delete retbuf;
			retbuf = tmpbuf;
			tmpbuf = NULL;
		}
		bytesread = fread(retbuf+bufsize-4096, 1, 4096, f);
	} while (bytesread != 0);
	if (feof(f)) {
		return retbuf;
	} else {
		cerr << "Error reading from file\n";
		exit(EXIT_FAILURE);
	}
}

char *getnextline(const char *buffer, const unsigned int bufsize, unsigned int &pos) {
	unsigned int i = pos;
	char *ret = NULL;
	if (pos >= bufsize) return ret;
	while ((pos < bufsize)&&(buffer[pos] != '\n')) pos++;
	if (buffer[pos] == '\n') {
		ret = new char[pos-i+2];
		memcpy(ret, buffer+i, pos-i+1);
		*(ret+pos+-i+1) = 0;
	} else {
		ret = new char[pos-i+1];
		memcpy(ret, buffer+i, pos-i+1);
		*(ret+pos-i) = 0;
	}
	pos++;
	return ret;
}

int main(int argc, char**argv) {
	vector <fastaseq> sequences;
	unsigned int chunklength = 4*(maxlen-minoverlap)+minoverlap;
	char *accept_prev = NULL;
	int usenaive = 2; //2 = default/unset i.e. use pepfrag but return naive if it is better, 1 = force naive only, 3 = force pepfrag only
	char *filename = NULL;
	FILE *f = NULL, *apf = NULL;
	int a;
	unsigned int i, j ,k;
	char *buffer = NULL, *line = NULL;
	unsigned int filesize, fpos;
	char *fragment = NULL, *poslist = NULL, *convstr = NULL, *sline = NULL;
	seqpos tspos;
	unsigned int parse_error;
	int *penalties;
	unsigned int numsequences = 0;
	fastaseq fseq;
	unsigned int ntitleidx;
	string tseq;
	unsigned int numchunks;
	int min_penalty;
	int min_penalty_idx;
	vector <string> chunk;
	unsigned int startidx;
	unsigned int penalty;
	vector <seqpos> positions;
	map <string, vector <seqpos> >::iterator it;

	a = 1;
	while (a < argc) {
		if (strcmp(argv[a],"-l") == 0) {
			a++;
			if (!isnumeric(argv[a])) {usage();}
			minlen = atoi(argv[a]);
		} else if (strcmp(argv[a],"-L") == 0) {
			a++;
			if (!isnumeric(argv[a])) {usage();}
			maxlen = atoi(argv[a]);
		} else if (strcmp(argv[a],"-o") == 0) {
			a++;
			if (!isnumeric(argv[a])) {usage();}
			minoverlap = atoi(argv[a]);
		} else if (strcmp(argv[a],"-O") == 0) {
			a++;
			if (!isnumeric(argv[a])) {usage();}
			maxoverlap = atoi(argv[a]);
		} else if (strcmp(argv[a],"-C") == 0) {
			a++;
			if (!isnumeric(argv[a])) {usage();}
			chunklength = atoi(argv[a]);
		} else if (strcmp(argv[a],"-p") == 0) {
			a++;
			accept_prev = argv[a];
		} else if (strcmp(argv[a],"-t") == 0) {
			usenaive = 1;
		} else if (strcmp(argv[a],"-n") == 0) {
			usenaive = 3;
		} else {
			if (strncmp(argv[a],"-",1) == 0) {
				usage();
			}
			filename = argv[a];
		}
		a++;
	}

	if (chunklength < minlen) { cerr << "Chunk length cannot be less the minimum fragment length\n"; return EXIT_FAILURE;}
	if (maxlen < minlen) {cerr << "Maximum fragment length cannot be less than minimum fragment length\n"; return EXIT_FAILURE;}
	if (maxoverlap < minoverlap) {cerr << "Maximum overlap length cannot be less than minimum overlap length\n"; return EXIT_FAILURE;}
	if (minoverlap >= minlen) {cerr << "Minimum fragment length must ge greater than minimum overlap length\n"; return EXIT_FAILURE;}

	if (filename != NULL) {
		f = fopen(filename, "r");
		if (f == NULL) {
			cerr << "Error opening " << filename << ": " << strerror(errno) << "\n";
			return EXIT_FAILURE;
		}
	} else {
		f = stdin;
	}

	if (accept_prev != NULL) {
		if (strcmp(accept_prev,"-") == 0) {
			if (f == stdin) {
				apf = f;
			} else {
				apf = stdin;
			}
		} else {
			if (strcmp(accept_prev,filename) == 0) {
				apf = f;
			} else {
				apf = fopen(accept_prev, "r");
				if (apf == NULL) {
					cerr << "Error opening " << accept_prev << ": " << strerror(errno) << "\n";
					return EXIT_FAILURE;
				}
			}
		}

		buffer = readfile(apf, filesize);
		fpos = 0;
		do {
			if (line != NULL) delete line;
			line = getnextline(buffer, filesize, fpos);
			if (sscanf(line, "%*u:") != EOF) {
				cout << line;
				titleidx++;
			} else {
				i = 0;
				while ((i < strlen(line))&&(line[i] != '\t')) i++;
				fragment = strndup(line, i-1);
				poslist = strdup(line+i+1);
				i = 0;
				while (i < strlen(poslist)) {
					parse_error = 1;
					while ((i < strlen(poslist))&&(isspace(poslist[i]))) i++; j = i;
					while ((i < strlen(poslist))&&(isdigit(poslist[i]))) i++; 
					if (poslist[i] == ':') {
						convstr = new char[i-j]; strncpy(convstr, poslist+j, i-j-1); tspos.idx = atoi(convstr); delete convstr;
						i++;
						j = i;
					} else {
						cerr << "Parse error reading previous results: " << line << "\n";
						return EXIT_FAILURE;
					}
					while ((i < strlen(poslist))&&(isdigit(poslist[i]))) i++; k = i;
					while ((i < strlen(poslist))&&(isspace(poslist[i]))) i++;
					if (poslist[i] == ',') {
						convstr = new char[k-j]; strncpy(convstr, poslist+j, k-j-1); tspos.pos = atoi(convstr); delete convstr;
						current[fragment].push_back(tspos);
						free(fragment);
						parse_error = 0;
						i++;
					} else {
						cerr << "Parse error reading previous results: " << line << "\n";
						return EXIT_FAILURE;
					}
				}
				free(poslist);
				if (parse_error) {
					cerr << "Parse error reading previous results: " << line << "\n";
					return EXIT_FAILURE;
				}
			}
		} while ((fpos < filesize)&&(line != NULL)&&(line[0] != '>'));
	}

	if (apf != f) {
		if (buffer != NULL) delete buffer;
		buffer = readfile(f, filesize);
		fpos = 0;
		if (line != NULL) delete line;
		line = getnextline(buffer, filesize, fpos);
	}
	while ((fpos < filesize)&&(line != NULL)) {
		while (allspace(line)) {
			if (line != NULL) delete line;
			line = getnextline(buffer, filesize, fpos);
		}
		if (line[0] != '>') {
			cerr << "Sequence input not in valid fasta format\n";
			return EXIT_FAILURE;
		}

		numsequences++;
		fseq.title = line+1;
		if (line != NULL) delete line;
		line = getnextline(buffer, filesize, fpos);
		fseq.sequence = "";
		while ((line != NULL)&&(!allspace(line))&&(line[0] != '>')) {
			sline = strip(line);
			fseq.sequence += sline;
			free(sline);
			delete line;
			line = getnextline(buffer, filesize, fpos);
		}
		sequences.push_back(fseq);
	}

	//pre calculate baseline using naive filtered method
	if (usenaive < 3) {
		nbest = current;
		ntitleidx = titleidx;
		for (i = 0; i < numsequences; i++) {
			if (usenaive == 1) {
				cout << ntitleidx << ": " << sequences[i].title;
			}
			seqlen = sequences[i].sequence.length();
			j = 0;
			do {
				tspos.idx = ntitleidx; tspos.pos = j; nbest[sequences[i].sequence.substr(j, maxlen)].push_back(tspos);
				j += maxlen-minoverlap;
			} while ((j < seqlen)&&(j+maxlen < seqlen));
			tspos.idx = ntitleidx; tspos.pos = j; nbest[sequences[i].sequence.substr(j, maxlen)].push_back(tspos);
			ntitleidx += 1;
		}
		nbestlen = nbest.size();
	} else {
		nbestlen = (int)(pow(2,sizeof(int)*8)-1);
	}

	if (usenaive > 1) {
		penalties = new int[4*minlen];
		for (i = 0; i < numsequences; i++) {
			//cout << titleidx << ": " << sequences[i].title;
			numchunks = 0;
			seqlen = sequences[i].sequence.length();
			if (chunklength > seqlen) {
				numchunks = 1;
				chunk.push_back(sequences[i].sequence);
			} else {
				startidx = 0;
				while (startidx+chunklength+minlen < seqlen) {
					numchunks++;
					min_penalty = -1;
					memset(penalties, 0, (4*minlen)*sizeof(int));
					for (j = 0; j <= 3*minlen; j++) {
						penalty = find_all_count(sequences[i].sequence.substr(startidx+chunklength-(2*minlen)+j, minlen), sequences[i].sequence);
						for (k = 0; k < min(minlen, (4*minlen+1)-j); k++) {
							penalties[j+k] += penalty;
							if ((min_penalty == -1)||(penalties[j+k] < min_penalty)) {
								min_penalty = penalties[j+k];
							}
						}
					}
					min_penalty_idx = minlen;
					for (j = minlen; j <= 3*minlen; j++) {
						if (penalties[j] <= min_penalty) {
							min_penalty_idx = j;
							min_penalty = penalties[j];
						}
					}
					min_penalty_idx = startidx+chunklength-(2*minlen)+min_penalty_idx;
					chunk.push_back(sequences[i].sequence.substr(startidx, min_penalty_idx-startidx));
					startidx = min_penalty_idx-(minoverlap-1);
				}
				numchunks++;
				chunk.push_back(sequences[i].sequence.substr(startidx));
			}

			if (chunk[numchunks-1].length() < minlen) {
				chunk[numchunks-2] += chunk[numchunks-1];
				chunk.pop_back();
				numchunks--;
			}
			//cout << numchunks << "\n";

			basepos = 0;
			for (j = 0; j < numchunks; j++) {
				//cout << chunk[j] << "\n";
				seqlen = chunk[j].length();
				best.clear();
				dfsnr(chunk[j]);
				current = best;
				basepos += seqlen-minoverlap+1;
			}
			titleidx++;
			chunk.clear();
			numchunks = 0;
		}
		delete penalties;
	}

	if ((usenaive == 1)||((usenaive == 2)&&((best.empty())||(best.size() > nbest.size())))) {
		best = nbest;
	}

	for (i = 0; i < numsequences; i++) {
		cout << i+1 << ": " << sequences[i].title;
	}
	for (it = best.begin(); it != best.end(); it++) {
		cout << (*it).first << "\t";
		for (i = 0; i < (*it).second.size(); i++) {
			cout << (*it).second[i].idx << ":" << (*it).second[i].pos;
			if (i == (*it).second.size()-1) {
				cout << "\n";
			} else {
				cout << ", ";
			}
		}
	}
	return EXIT_SUCCESS;
}
