pepfrag: pepfrag.cc
	g++ -O3 -o pepfrag -Wall pepfrag.cc

test: pepfrag
	make -C sequences
	
