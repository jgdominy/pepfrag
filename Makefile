pepfrag: pepfrag.cc
	g++ -g3 -o pepfrag -Wall pepfrag.cc

test: sequences/flus.results pepfrag
	./pepfrag -l 8 -L 15 -o 11 -O 13 -C 32 sequences/flus.fasta > sequences/flus.results
	
