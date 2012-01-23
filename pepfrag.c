#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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
	int length;
	element *items;
} collection;

//a struture to keep a record of sequences read in from file, and their titles
typedef struct {
	char *title;
	char *sequence;
} sequence_rec;

//maths convenience function
int max(a, b) {
	return (a > b?a:b);
}
	
int min(a, b) {
	return (a < b?a:b);
}

//string convenience functions

//checks whether a string consists entirely of white space
int empty(const char *s) {
	char *i = s;
	while (*i != 0) {
		if (!isspace(*i)) return 0;
		i++;
	}
	return 1;
}

//substr allocates and returns a new string which is a substring of s from i to
//j exclusive, where i < j; If i or j are negative they refer to distance from
//the end of the s
char *substr(const char *s, int i, int j) {
	char *ret;
	if (i < 0) i = strlen(s)-i;
	if (j < 0) j = strlen(s)-j;
	ret = malloc(j-i+1);
	strncpy(ret,s,j-i);
	return ret;
}

//strips white space from either end of the string
void strip(char **s) {
	int i, j, len;
	char *tmp = *s;
	len = strlen(*s);
	i = 0;
	while (isspace(*s[i]))&&(i < len) i++;
	j = strlen(*s)-1;
	while ((isspace(*s[j]))&&(j > 0)) j--;
	*s = strndup(*s+i, j-i);
	free(tmp);
}

//incollection uses a binary search to locate a particular peptide fragment
//within the given collection returning the items index plus 1 if found, otherwise 0
int incollection(const collection *c, const char *s) {
	int min = 0;
	int max = c->length-1;
	int mid, d;

	do {
		mid = (min+max)/2;
		d = strcmp(s,c->items[mid].seq);
		if (d < 0) {
			max = mid-1;
		} else {
			min = mid+1;
		}
	} while ((d != 0)&&(min < max));
	return (d == 0);
}

//col_insert inserts a fragment into a collection. If the fragment exists, the
//position information is appended to its position list, otherwise the
//fragement is added and the position list is set to contain only the added
//position.
void insert(collection *c, const char *seq, int sidx, int spos) {
	int min = 0;
	int max = c->length-1;
	int mid;
	int d;
	seqpos *newpos = NULL;
	element *newelem = NULL;

	newpos = malloc(sizeof(seqpos));
	newpos->seqidx = sidx;
	newpos->pos = spos;
	if (c->length == 0) {
		c->length = 1;
		realloc(c->items, c->length*sizeof(element));
		newpos->next = NULL;
		newelem = malloc(sizeof(element));
		newelem->seq = strdup(seq);
		newelem->positions = newpos;
		c->items = newelem;
	} else {
		do {
			mid = (min+max)/2;
			d = strcmp(seq,c->items[mid].seq);
			if (d < 0) {
				max = mid-1;
			} else {
				min = mid+1;
			}
		} while ((d != 0)&&(min < max));
		if (d == 0) {
			//append to sequence's position list
			newpos->next = c->items[mid].positions;
			c->items[mid].positions = newpos;
		} else {
			//col_insert new sequence at correct position
			newpos->next = NULL;
			newelem = malloc(sizeof(element));
			newelem->seq = strdup(seq);
			newelem->positions = newpos;
			c->length++;
			realloc(c->items, c->length*sizeof(element));
			if (d > 0) mid--;	
			memmove(c->items+mid+2, c->items+mid+1, c->length-mid-2);
		}
	}
}

//remove

//pushidx

//popidx

//copy
collection *copy(collection *c) {
	collection *newcol = NULL;
	seqpos *newpos = NULL, *posptr;
	//element *newelem = NULL;
	int i;

	if (c == NULL) return c;
	newcol = malloc(sizeof(collection));
	newcol->items = malloc(sizeof(element)*c->length);
	for (i = 0; i < c->length; i++) {
		newcol->items[i].seq = strdup(c->items[i].seq);
		newcol->items[i].positions = NULL;
		posptr = c->items[i].positions;
		while (posptr != NULL) {
			newpos = malloc(sizeof(seqpos));
			newpos->seqidx = posptr->seqidx;
			newpos->pos = posptr->pos;
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
		"followed by a comma separated list of sequence index:position pairs. The number\n"
		"before the colon indicates which sequence the fragment was found in, and the\n"
		"number after the colon indicates where in the sequence the fragment can be\n"
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
	exit(1);
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
	//char *seq;
	struct stack *prev;
} stack;

stack *push(stack *S, int pos, int i, int j, int realminoverlap, int realmaxoverlap) {
	stack *nS;
	nS = malloc(sizeof(stack));
	nS->pos = pos;
	//nS->seq = strdup(seq);
	//nS->basepos = basepos;
	nS->i = i;
	nS->j = j;
	nS->realminoverlap = realminoverlap;
	nS->realmaxoverlap = realmaxoverlap;
	nS->prev = S;
	return nS;
}

stack *pop(stack *S, int *pos, int *i, int *j, int *realminoverlap, int *realmaxoverlap) {
	stack *prev = S->prev;
	*pos = S->pos;
	//free(*seq);
	//*seq = S->seq;
	//*basepos = S->basepos;
	*i = S->i;
	*j = S->j;
	*realminoverlap = S->realminoverlap;
	*realmaxoverlap = S->realmaxoverlap;
	free(S);
	return prev;
}

//Setup up variables and default values
collection *best = NULL;
collection *current = NULL;
int minlen = 8;
int maxlen = 30;
int minoverlap = 2;
int maxoverlap = 4;
int titleidx = 1;
int nbestlen;
int seqlen;

void bruteforce(int pos, const char *seq, const int basepos) {
	stack *S = NULL;
	int i, j, realminoverlap, realmaxoverlap;
	char *tkey;

	start:
	i = maxlen;
	while (i >= minlen) {
		realminoverlap = min(minoverlap,max(0,pos-minoverlap));
		realmaxoverlap = min(maxoverlap,max(0,pos-maxoverlap));
		j = realminoverlap;
		while (j < realmaxoverlap) {
			if (seqlen - (pos-j+1) < minlen) {
				S = pop(S,&pos,&i,&j,&realminoverlap,&realmaxoverlap);
				if (S == NULL) return; else goto retpos;
			}
			if (pos-j+i >= seqlen) {
				if (incollection(current, seq+pos-j)) {
					pushidx(current, seq+pos-j, titleidx, basepos+pos-j);
				} else {
					if ((best != NULL)&&((current->length+1 > best->length)||(current->length+1 > nbestlen))) {
						S = pop(S,&pos,&i,&j,&realminoverlap,&realmaxoverlap);
						if (S == NULL) return; else goto retpos;
					} else {
						col_insert(current, seq+pos-j, titleidx, basepos+pos-j);
					}
				}
				if ((best == NULL)||(current->length < best->length)) {
					destroy(best); best = NULL;
					best = copy(current);
				}
				if (popidx(current, seq+pos-j) == 0) {
					col_remove(current, seq+pos-j);
				}
			} else {
				tkey = strndup(seq+pos-j,i);
				if (incollection(current, tkey)) {
					pushidx(current, tkey, titleidx, basepos+pos-j);
				} else {
					if ((best != NULL)&&((current->length+1 > best->length)||(current->length+1 > nbestlen))) {
						S = pop(S,&pos,&i,&j,&realminoverlap,&realmaxoverlap);
						if (S == NULL) return; else goto retpos;
					} else {
						col_insert(current, tkey, titleidx, basepos+pos-j);
					}
				}
				S = push(S,pos,i,j,realminoverlap,realmaxoverlap);
				pos = pos-j+i;
				goto start;	
				retpos:
				if (popidx(current, tkey) == 0) {
					col_remove(current, tkey);
				}
				free(tkey);
			}
			j++;
		}
		i--;
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
	char **sequences = NULL;
	int chunklength = 4*(maxlen-minoverlap)+minoverlap;
	char *accept_prev = NULL;
	int usenaive = 2;
	char *filename = NULL;
	FILE *apf = NULL;
	FILE *f = NULL;
	int i, j ,k;
	char *line = NULL;
	char *fragment = NULL;
	char *poslist = NULL;
	int linelen;
	int rcount;
	int seqidx;
	int seqpos;
	int parse_error;
	char *convstr = NULL;
	int *penalties;

	i = 1;
	while (i < argc) {
		if (strcmp(argv[i],"-l") == 0) {
			i++;
			minlen = atoi(argv[i]);
		} else if (strcmp(argv[i],"-L") == 0) {
			i++;
			maxlen = atoi(argv[i]);
		} else if (strcmp(argv[i],"-o") == 0) {
			i++;
			minoverlap = atoi(argv[i]);
		} else if (strcmp(argv[i],"-O") == 0) {
			i++;
			maxoverlap = atoi(argv[i]);
		} else if (strcmp(argv[i],"-C") == 0) {
			i++;
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

	if (chunklength < minlen) {fprintf(stderr, "Chunk length cannot be less the minimum fragment length\n"); exit(1);}
	if (maxlen < minlen) {fprintf(stderr, "Maximum fragment length cannot be less than minimum fragment length\n"); exit(1);}
	if (maxoverlap < minoverlap) {fprintf(stderr, "Maximum overlap length cannot be less than minimum overlap length\n"); exit(1);}

	if (filename != NULL) {
		f = fopen(filename, "r");
	} else {
		f = stdin;
	}

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
			}
		}

		rcount = getline(&line, &linelen, apf);
		while ((rcount != 0)&&(line[0] != '>')) {
			if (sscanf(line, "%*u:") != EOF) {
				printf(line);
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
						exit(1);
					}
					while ((i < strlen(poslist))&&(isdigit(poslist[i]))) i++; k = i;
					while ((i < strlen(poslist))&&(isspace(poslist[i]))) i++;
					if (poslist[i] == ',') {
						convstr = malloc(k-j); strncpy(convstr, poslist+j, k-j-1); seqpos = atoi(convstr); free(convstr);
						col_insert(current, fragment, seqidx, seqpos);
						free(fragment);
						parse_error = 0;
						i++;
					} else {
						fprintf(stderr, "Parse error reading previous results: %s\n", line);
						exit(1);
					}
				}
				if (parse_error) {
					fprintf(stderr, "Parse error reading previous results: %s\n", line);
					exit(1);
				}
			}
			rcount = getline(&line, &linelen, apf);
		}
	}

	if (apf != f) rcount = getline(&line, &linelen, f);
	while (rcount != 0) {
		while (empty(line)) rcount = getline(&line, &linelen, f);
		if (line[0] != '>') {
			fprintf(stderr,"Sequence input not in valid fasta format\n");
			exit(1);
		}

		realloc(sequences,sizeof(sequence_rec)*numsequences);
		sequences[numsequences-1].title = malloc(strlen(line)); strncpy(title,line+1); strip(&title);
		rcount = getline(&line, &linelen, f);
		sequences[numsequences-1].sequence = malloc(1); sequences[numsequences-1].sequence[0] = 0;
		while ((!empty(line))&&(line[0] != '>')) {
			strip(&line);
			realloc(sequences[numsequences-1].sequence, strlen(sequences[numsequences-1].sequence)+strlen(line)+1);
			strcat(sequences[numsequences-1].sequence,line);
			rcount = getline(&line, &linelen, f);
		}
	}

	//pre calculate baseline using naive filtered method
	if (usenaive > 1) {
		nbest = copy(current);
		ntitleidx = titleidx;
		for (i = 0; i < numsequences; i++) {
			if (usenaive == 3) {
				printf("%i: %s", ntitleidx, sequences[i].title);
			}
			for (j = 0; j < strlen(sequences[i].sequence); j += maxlen-minoverlap) {
				tseq = substr(sequences[i].sequence, i, i+maxlen);
				col_insert(nbest, tseq, ntitleidx, j);
				free(tseq);
			}
			ntitleidx += 1
		}
		nbestlen = nbest.length;
	} else {
		nbestlen = (int)(pow(2,sizeof(int)*8)-1);
	}

	if (usenaive < 2) {
		penalties = malloc((4*minlen+1)*sizeof(int));
		for (i = 0; i < numsequences; i++) {
			printf("%i: %s", titleidx, sequences[i].title);
			numchunks = 1;
			chunk = malloc(sizeof(char *));
			seqlen = strlen(sequences[i].sequence);
			if (chunklength > seqlen) {
				chunk[0] = malloc(seqlen+1);
				strcpy(chunk[0],sequences[i].sequence);
			} else {
				startidx = 0;
				while (startidx+chunklength < seqlen) {
					min_penalty = 0;
					memset(penalties, 0, (4*minlen+1)*sizeof(int));
					for (j = 0; j < 4*minlen+1; j++) {
						tstr = substr(sequences[i].sequence, startidx+chunklength-(2*minlen)+j, startidx+chunklength-minlen+j);
						penalty = find_all_count(tstr, seq);
						free(tstr);
						for (k = 0; k < min(minlen, (4*minlen+1)-j); k++) {
							penalties[j+k] += penalty;
							if (min_penalty < penalties[j+k]) {
								min_penalty = penalties[j+k];
							}
						}
					}
					min_penalty_idx = minlen;
					for (j = minlen; j < 3*minlen+1; j++) {
						if (penalties[j] < min_penalty) {
							min_penalty_idx = j;
							min_penalty = penalties[j];
						}
					}
					min_penalty_idx = startidx+chunklength-(2*minlen)+min_penalty_idx;
					numchunks++;
					realloc(chunks, sizeof(char *)*numchunks);
					chunk[numchunks-1] = malloc(min_penalty_idx-startidx+1);
					strncpy(chunk[numchunks-1], sequences[i]->sequence+startidx, min_penalty_idx-startidx);
					startidx = min_penalty_idx-(minoverlap-1);
				}
				numchunks++;
				realloc(chunks, sizeof(char *)*numchunks);
				chunk[numchunks-1] = malloc(strlen(sequences[i]->sequence)-startidx+1);
				strcpy(chunk[numchunks-1], sequences[i].sequence+startidx);
			}

			if (strlen(chunks[numchunks-1]) < minlen) {
				realloc(chunks[numchunks-2], strlen(chunks[numchunks-2])+strlen(chunks[numchunks-1])+1);
				strcat(chunks[numchunks-2], chunks[numchunks-1]);
				free(chunks[numchunks-1]);
				numchunks--;
				realloc(chunks, sizeof(char *)*numchunks);
			}

			basepos = 0;
			for (c = 0; c < numchunks; c++) {
				seqlen = strlen(chunks[c]);
				destroy(best); best = NULL;
				bruteforce(0,chunks[c],basepos);
				current = copy(best);
				basepos += seqlen-minoverlap+1;
			}
			titleidx++;
		}
	}

	if ((usenaive == 2)||((usenaive == 1)&&((best == NULL)||(best.length > nbest.length)))) {
		best copy(nbest);
	}

	printout(best);
}
