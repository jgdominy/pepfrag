/* This code is Copyright (c) the Centre for Proteomic and Genomic Research
 * (2011-2012) and is distributed under the terms of Gnu General Public License
 * version 3 or higher.  For the full terms and conditions please see the file
 * COPYING which should have been included in this distribution. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>

//seqpos, element and collection together form the main data structure used to
//store a list of peptide fragments found/used and their respective positions
//within all the sequences. It is implemented as an ordered list using a dynamic
//array, to facilitate binary searching and avoid stack usage. 
typedef struct seqpos {
	int seqidx;
	int pos;
	struct seqpos *next;
} seqpos;

typedef struct element {
	char *seq;
	seqpos *positions;
} element;

typedef struct collection {
	int capacity;
	int length;
	element *items;
} collection;

//Setup up variables and default values
collection *best = NULL;
collection *current = NULL;
collection *nbest = NULL;
int minlen = 8;
int maxlen = 30;
int minoverlap = 2;
int maxoverlap = 4;
int titleidx = 1;
int nbestlen;
int seqlen;

//a struture to keep a record of sequences read in from file, and their titles
typedef struct {
	char *title;
	char *sequence;
} sequence_rec;

//maths convenience function
int max(int a, int b) {
	return (a > b?a:b);
}
	
int min(int a, int b) {
	return (a < b?a:b);
}

//string convenience functions

//checks whether a string consists entirely of white space
int empty(const char *s) {
	int i;
	i = 0;
	while (s[i] != 0) {
		if (!isspace(s[i])) return 0;
		i++;
	}
	return 1;
}

int isnumeric(const char *s) {
	int i;
	i = 0;
	while (s[i] != 0) {
		if (!isdigit(s[i])) return 0;
		i++;
	}
	return 1;
}

//substr allocates and returns a new string which is a substring of s from i to
//j exclusive, where i < j; If i or j are negative they refer to distance from
//the end of the s
char *substr(const char *s, int i, int j) {
	if (i < 0) i = strlen(s)+i;
	if (j < 0) j = strlen(s)+j;
	return strndup(s+i,j-i);
}

//strips white space from either end of the string
char *strip(const char *s) {
	int i, j, len;

	len = strlen(s);
	if (len == 0) return strdup(s);
	i = 0;
	while ((isspace(s[i]))&&(i < len)) {
		i++;
	}
	j = len-1;
	while ((isspace(s[j]))&&(j > 0)) {
		j--;
	}
	return strndup(s+i, j-i+1);
}

//elemcmp is a callback function for bsearch
int elemcmp(const void *n, const void *h) {
	return strcmp((char *)n, ((element *)h)->seq);
}

//pushspos inserts a fragment into a collection. If the fragment exists, the
//position information is appended to its position list, otherwise the
//fragement is added and the position list is set to contain only the added
//position.
void pushspos(collection *c, const char *seq, int sidx, int spos) {
	element *ipos;
	seqpos *newpos = NULL;
	int i;

	if (strlen(seq) < minlen) {
		fprintf(stderr, "Push attempt too short: %s @ %i\n", seq, spos);
	}
	newpos = malloc(sizeof(seqpos));
	newpos->seqidx = sidx;
	newpos->pos = spos;
	if (c->length == 0) {
		c->length = 1;
		c->capacity = 64;
		c->items = realloc(c->items, 64*sizeof(element));
		newpos->next = NULL;
		c->items[0].seq = strdup(seq);
		c->items[0].positions = newpos;
	} else {
		ipos = bsearch(seq,c->items,c->length,sizeof(element),elemcmp);
		if (ipos != NULL) {
			//append to sequence's position list
			newpos->next = ipos->positions;
			ipos->positions = newpos;
		} else {
			//insert new sequence at correct position
			newpos->next = NULL;
			if (c->length+1 > c->capacity) {
				c->capacity = (c->capacity==0?64:c->capacity*2);
				c->items = realloc(c->items, c->capacity*sizeof(element));
			}
			i = 0;
			while ((i < c->length)&&(strcmp(seq,c->items[i].seq) == 1)) {
				i++;
			}
			if (i < c->length) memmove(c->items+i+1, c->items+i, ((c->length-i)*sizeof(element)));
			c->items[i].seq = strdup(seq);
			c->items[i].positions = newpos;
			c->length++;
		}
	}
	//printf("DEBUG: push %s\n", seq);
}

//popspos pops the most recent position off a given sequences position stack, in a given collection.
void popspos(collection *c, const char *seq) {
	element *ipos;
	seqpos *t;

	//printf("DEBUG: pop %s\n", seq);
	ipos = bsearch(seq,c->items,c->length,sizeof(element),elemcmp);
	if (ipos != NULL) {
		//pop the latest position off the stack
		t = ipos->positions;
		if (t == NULL) {
			fprintf(stderr, "Sequence '%s' in collection has and empty psotion list\n", seq);
			exit(EXIT_FAILURE);
		}
		ipos->positions = ipos->positions->next;
		free(t);
		if (ipos->positions == NULL) {
			free(ipos->seq);
			memmove(ipos, ipos+1, (c->length-(ipos+1-c->items))*sizeof(element));
			c->length--;
		}
	} else {
		fprintf(stderr, "Could not find sequence '%s' in collection\n", seq);
		exit(EXIT_FAILURE);
	}
}

//copy
collection *copy(collection *c) {
	collection *newcol = NULL;
	seqpos *newpos = NULL, *posptr, *trail = NULL;
	//element *newelem = NULL;
	int i;

	if (c == NULL) return c;
		
	newcol = malloc(sizeof(collection));
	if (c->length == 0) {
		newcol->items = NULL;
		newcol->length = 0;
		newcol->capacity = 0;
	} else {
		newcol->items = malloc(sizeof(element)*c->length);
		newcol->length = c->length;
		newcol->capacity = c->length;
	}
	for (i = 0; i < c->length; i++) {
		newcol->items[i].seq = strdup(c->items[i].seq);
		newcol->items[i].positions = NULL;
		posptr = c->items[i].positions;
		while (posptr != NULL) {
			newpos = malloc(sizeof(seqpos));
			if (trail != NULL) {
				trail->next = newpos;
			}
			newpos->seqidx = posptr->seqidx;
			newpos->pos = posptr->pos;
			newpos->next = NULL;
			trail = newpos;
			if (newcol->items[i].positions == NULL) {
				newcol->items[i].positions = newpos;
			}
			posptr = posptr->next;
		}
	}
	return newcol;
}

//destroy
void destroy(collection *c) {
	int i;
	seqpos *posptr = NULL, *next = NULL;

	if (c == NULL) return;
	for (i = 0; i < c->length; i++) {
		free(c->items[i].seq);
		posptr = c->items[i].positions;
		while (posptr != NULL) {
			next = posptr->next;
			free(posptr);
			posptr = next;
		}
	}
	free(c->items);
	free(c);
}

//print a usage message and exit
void usage() {
	fprintf(stderr, "Pepfrag produces a list of peptide fragments providing complete coverage of the\n"
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
		"--- End ---\n"
	);
	exit(EXIT_FAILURE);
}

//The brute force algortithm to analyse a single chunk. We test each set in the
//space using a depth first traversal starting from maxlen,minoverlap to
//minlen,maxoverlap; There will be fewer fragments for maxlen,minoverlap so
//evaluating this set first sets a good basline for pruning. Any other sets
//that are partially evaluated and grow in size past this baseline are
//discarded (pruned)

typedef struct stack {
	int pos;
	//int basepos;
	int i;
	int j;
	int realminoverlap;
	int realmaxoverlap;
	char *tkey;
	struct stack *prev;
} stack;

stack *push(stack *S, int pos, int i, int j, int realminoverlap, int realmaxoverlap, char *tkey) {
	stack *nS;
	nS = malloc(sizeof(stack));
	nS->pos = pos;
	//nS->seq = strdup(seq);
	//nS->basepos = basepos;
	nS->i = i;
	nS->j = j;
	nS->realminoverlap = realminoverlap;
	nS->realmaxoverlap = realmaxoverlap;
	nS->tkey = tkey;
	nS->prev = S;
	return nS;
}

stack *pop(stack *S, int *pos, int *i, int *j, int *realminoverlap, int *realmaxoverlap, char **tkey) {
	stack *prev = S->prev;
	*pos = S->pos;
	//free(*seq);
	//*seq = S->seq;
	//*basepos = S->basepos;
	*i = S->i;
	*j = S->j;
	*realminoverlap = S->realminoverlap;
	*realmaxoverlap = S->realmaxoverlap;
	*tkey = S->tkey;
	free(S);
	return prev;
}

void bruteforce(int pos, const char *seq, const int basepos) {
	stack *S = NULL;
	int i, j, realminoverlap, realmaxoverlap;
	char *tkey = NULL;

	S = push(S,pos,i,j,realminoverlap,realmaxoverlap,tkey);
	start:
	realminoverlap = min(minoverlap,max(0,pos-minoverlap));
	realmaxoverlap = min(maxoverlap,max(0,pos-maxoverlap));
	i = maxlen;
	while (i >= minlen) {
		j = realminoverlap;
		while (j <= realmaxoverlap) {
			if (seqlen - (pos-j+1) < minlen) {
				S = pop(S,&pos,&i,&j,&realminoverlap,&realmaxoverlap,&tkey);
				if (S == NULL) return; else goto retpos;
			}
			if (pos-j+i >= seqlen) {
				if (bsearch(seq+pos-j, current->items, current->length, sizeof(element), elemcmp)) {
					pushspos(current, seq+pos-j, titleidx, basepos+pos-j);
				} else {
					if ((best != NULL)&&((current->length+1 > best->length)||(current->length+1 > nbestlen))) {
						S = pop(S,&pos,&i,&j,&realminoverlap,&realmaxoverlap,&tkey);
						if (S == NULL) return; else goto retpos;
					} else {
						pushspos(current, seq+pos-j, titleidx, basepos+pos-j);
					}
				}
				if ((best == NULL)||(current->length < best->length)) {
					destroy(best); best = NULL;
					best = copy(current);
				}
				popspos(current, seq+pos-j);
			} else {
				tkey = strndup(seq+pos-j,i);
				if (bsearch(tkey, current->items, current->length, sizeof(element), elemcmp)) {
					pushspos(current, tkey, titleidx, basepos+pos-j);
				} else {
					if ((best != NULL)&&((current->length+1 > best->length)||(current->length+1 > nbestlen))) {
						S = pop(S,&pos,&i,&j,&realminoverlap,&realmaxoverlap,&tkey);
						if (S == NULL) return; else goto retpos;
					} else {
						pushspos(current, tkey, titleidx, basepos+pos-j);
					}
				}
				S = push(S,pos,i,j,realminoverlap,realmaxoverlap,tkey);
				pos = pos-j+i;
				goto start;	
				retpos:
				popspos(current, tkey);
				free(tkey);
			}
			j++;
		}
		i--;
	}
	if (S != NULL) {
		S = pop(S,&pos,&i,&j,&realminoverlap,&realmaxoverlap,&tkey);
		if (S == NULL) return; else goto retpos;
	}
}

//Counts the number of occurences of needle in haystack
int find_all_count(const char *needle, const char *haystack) {
	int i, c = 0;
	for (i = 0; i < strlen(haystack)-strlen(needle); i++) {
		if (strncmp(haystack+i,needle,strlen(needle)) == 0) c++;
	}
	return c;
}

int main(int argc, char**argv) {
	sequence_rec *sequences = NULL;
	int chunklength = 4*(maxlen-minoverlap)+minoverlap;
	char *accept_prev = NULL;
	int usenaive = 2; //2 = default/unset i.e. use pepfrag but return naive if it is better, 1 = force naive only, 3 = force pepfrag only
	char *filename = NULL;
	FILE *apf = NULL;
	FILE *f = NULL;
	int i, j ,k;
	char *line = NULL;
	char *sline = NULL;
	char *fragment = NULL;
	char *poslist = NULL;
	size_t linelen;
	int rcount;
	int seqidx;
	int seqloc;
	int parse_error;
	char *convstr = NULL;
	int *penalties;
	int numsequences = 0;
	int ntitleidx;
	char *tseq;
	int numchunks, maxchunks;
	int min_penalty;
	int min_penalty_idx;
	char **chunk;
	int startidx;
	int penalty;
	int basepos;
	seqpos *tmpspos;	

	i = 1;
	while (i < argc) {
		if (strcmp(argv[i],"-l") == 0) {
			i++;
			if (!isnumeric(argv[i])) {usage();}
			minlen = atoi(argv[i]);
		} else if (strcmp(argv[i],"-L") == 0) {
			i++;
			if (!isnumeric(argv[i])) {usage();}
			maxlen = atoi(argv[i]);
		} else if (strcmp(argv[i],"-o") == 0) {
			i++;
			if (!isnumeric(argv[i])) {usage();}
			minoverlap = atoi(argv[i]);
		} else if (strcmp(argv[i],"-O") == 0) {
			i++;
			if (!isnumeric(argv[i])) {usage();}
			maxoverlap = atoi(argv[i]);
		} else if (strcmp(argv[i],"-C") == 0) {
			i++;
			if (!isnumeric(argv[i])) {usage();}
			chunklength = atoi(argv[i]);
		} else if (strcmp(argv[i],"-p") == 0) {
			i++;
			accept_prev = argv[i];
		} else if (strcmp(argv[i],"-t") == 0) {
			usenaive = 1;
		} else if (strcmp(argv[i],"-n") == 0) {
			usenaive = 3;
		} else {
			if (strncmp(argv[i],"-",1) == 0) {
				usage();
			}
			filename = argv[i];
		}
		i++;
	}

	if (chunklength < minlen) {fprintf(stderr, "Chunk length cannot be less the minimum fragment length\n"); return EXIT_FAILURE;}
	if (maxlen < minlen) {fprintf(stderr, "Maximum fragment length cannot be less than minimum fragment length\n"); return EXIT_FAILURE;}
	if (maxoverlap < minoverlap) {fprintf(stderr, "Maximum overlap length cannot be less than minimum overlap length\n"); return EXIT_FAILURE;}

	if (filename != NULL) {
		f = fopen(filename, "r");
		if (f == NULL) {
			fprintf(stderr, "Error opening %s: %s\n", filename, strerror(errno));
			return EXIT_FAILURE;
		}

	} else {
		f = stdin;
	}

	current = malloc(sizeof(collection)); current->length = 0; current->items = malloc(sizeof(element)*64); current->capacity = 64;
	if (accept_prev != NULL) {
		if (strcmp(accept_prev, "-") == 0) {
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
					fprintf(stderr, "Error opening %s: %s\n", filename, strerror(errno));
					return EXIT_FAILURE;
				}
			}
		}

		rcount = getline(&line, &linelen, apf);
		while ((rcount != 0)&&(line[0] != '>')) {
			if (sscanf(line, "%*u:") != EOF) {
				printf("%s", line);
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
					while ((i < strlen(poslist))&&(isdigit(poslist[i]))) i++; if (poslist[i] == ':') {
						convstr = malloc(i-j); strncpy(convstr, poslist+j, i-j-1); seqidx = atoi(convstr); free(convstr);
						i++;
						j = i;
					} else {
						fprintf(stderr, "Parse error reading previous results: %s\n", line);
						return EXIT_FAILURE;
					}
					while ((i < strlen(poslist))&&(isdigit(poslist[i]))) i++; k = i;
					while ((i < strlen(poslist))&&(isspace(poslist[i]))) i++;
					if (poslist[i] == ',') {
						convstr = malloc(k-j); strncpy(convstr, poslist+j, k-j-1); seqloc = atoi(convstr); free(convstr);
						pushspos(current, fragment, seqidx, seqloc);
						free(fragment);
						parse_error = 0;
						i++;
					} else {
						fprintf(stderr, "Parse error reading previous results: %s\n", line);
						return EXIT_FAILURE;
					}
				}
				if (parse_error) {
					fprintf(stderr, "Parse error reading previous results: %s\n", line);
					return EXIT_FAILURE;
				}
			}
			rcount = getline(&line, &linelen, apf);
		}
	}

	if (apf != f) rcount = getline(&line, &linelen, f); else rcount = 0;
	while (rcount != -1) {
		while (empty(line)) rcount = getline(&line, &linelen, f);
		if (line[0] != '>') {
			fprintf(stderr,"Sequence input not in valid fasta format\n");
			return EXIT_FAILURE;
		}

		numsequences++;
		sequences = realloc(sequences,sizeof(sequence_rec)*numsequences);
		sequences[numsequences-1].title = strdup(line+1); strip(sequences[numsequences-1].title);
		rcount = getline(&line, &linelen, f);
		sequences[numsequences-1].sequence = malloc(1); sequences[numsequences-1].sequence[0] = 0;
		while ((!empty(line))&&(line[0] != '>')&&(rcount != -1)) {
			sline = strip(line);
			sequences[numsequences-1].sequence = realloc(sequences[numsequences-1].sequence, strlen(sequences[numsequences-1].sequence)+strlen(sline)+1);
			strcat(sequences[numsequences-1].sequence,sline);
			free(sline);
			rcount = getline(&line, &linelen, f);
		}
	}

	//pre calculate baseline using naive filtered method
	if (usenaive < 3) {
		nbest = copy(current);
		ntitleidx = titleidx;
		for (i = 0; i < numsequences; i++) {
			if (usenaive == 1) {
				printf("%i: %s", ntitleidx, sequences[i].title);
			}
			seqlen = strlen(sequences[i].sequence);
			j = 0;
			do {
				tseq = substr(sequences[i].sequence, j, j+maxlen);
				pushspos(nbest, tseq, ntitleidx, j);
				free(tseq);
				j += maxlen-minoverlap;
			} while ((j < strlen(sequences[i].sequence))&&(j+maxlen < seqlen));
			ntitleidx += 1;
		}
		nbestlen = nbest->length;
	} else {
		nbestlen = (int)(pow(2,sizeof(int)*8)-1);
	}

	if (usenaive > 1) {
		penalties = malloc((4*minlen)*sizeof(int));
		for (i = 0; i < numsequences; i++) {
			printf("%i: %s", titleidx, sequences[i].title);
			numchunks = 0;
			seqlen = strlen(sequences[i].sequence);
			if (chunklength > seqlen) {
				numchunks = 1;
				maxchunks = 1;
				chunk = malloc(sizeof(char *)*maxchunks);
				chunk[0] = strdup(sequences[i].sequence);
			} else {
				maxchunks = (seqlen/(chunklength-(2*minoverlap)))+2;
				chunk = malloc(sizeof(char *)*maxchunks);
				startidx = 0;
				while (startidx+chunklength < seqlen) {
					numchunks++;
					if (numchunks > maxchunks) {
						maxchunks = maxchunks*3/2;
						chunk = realloc(chunk, sizeof(char *)*maxchunks);
					}
					min_penalty = -1;
					memset(penalties, 0, (4*minlen)*sizeof(int));
					for (j = 0; j <= 3*minlen; j++) {
						tseq = substr(sequences[i].sequence, startidx+chunklength-(2*minlen)+j, startidx+chunklength-minlen+j);
						penalty = find_all_count(tseq, sequences[i].sequence);
						free(tseq);
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
					chunk[numchunks-1] = strndup(sequences[i].sequence+startidx, min_penalty_idx-startidx);
					startidx = min_penalty_idx-(minoverlap-1);
				}
				numchunks++;
				if (numchunks > maxchunks) {
					maxchunks = maxchunks*3/2;
					chunk = realloc(chunk, sizeof(char *)*maxchunks);
				}
				chunk[numchunks-1] = strdup(sequences[i].sequence+startidx);
			}

			if (strlen(chunk[numchunks-1]) < minlen) {
				chunk[numchunks-2] = realloc(chunk[numchunks-2], strlen(chunk[numchunks-2])+strlen(chunk[numchunks-1])+1);
				strcat(chunk[numchunks-2], chunk[numchunks-1]);
				free(chunk[numchunks-1]);
				numchunks--;
			}

			//for (j = 0; j < numchunks; j++) {
			//	printf("%s\n", chunk[j]);
			//}

			basepos = 0;
			for (j = 0; j < numchunks; j++) {
				seqlen = strlen(chunk[j]);
				destroy(best); best = NULL;
				bruteforce(0,chunk[j],basepos);
				current = copy(best);
				basepos += seqlen-minoverlap+1;
			}
			titleidx++;
			while (numchunks > 0) {free(chunk[numchunks-1]); numchunks--;}
		}
		free(penalties);
	}

	if ((usenaive == 1)||((usenaive == 2)&&((best == NULL)||(best->length > nbest->length)))) {
		best = copy(nbest);
	}

	for (i = 0; i < numsequences; i++) {
		printf("%i: %s\n", i, sequences[i].title);
	}
	for (i = 0; i < best->length; i++) {
		printf("%s\t", best->items[i].seq);
		tmpspos = best->items[i].positions;
		while (tmpspos != NULL) {
			printf("%i:%i", tmpspos->seqidx, tmpspos->pos);
			if (tmpspos->next == NULL) {
				printf("\n");
			} else {
				printf(", ");
			}
			tmpspos = tmpspos->next;
		}
	}
	return EXIT_SUCCESS;
}
