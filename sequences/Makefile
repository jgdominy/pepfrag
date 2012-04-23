test: flus-N.results flus-C32.results flus-C48.results flus-C64.results

flus-N.results: flus.fasta
	@../pepfrag -l 12 -L 15 -o 11 -O 13 -n $< | grep -v "^[0-9]\+:" | wc -l > $@

flus-C32.results: flus.fasta
	@../pepfrag -l 12 -L 15 -o 11 -O 13 -C 32 $< | grep -v "^[0-9]\+:" | wc -l > $@

flus-C48.results: flus.fasta
	@../pepfrag -l 12 -L 15 -o 11 -O 13 -C 48 $< | grep -v "^[0-9]\+:" | wc -l > $@

flus-C64.results: flus.fasta
	@../pepfrag -l 12 -L 15 -o 11 -O 13 -C 64 $< | grep -v "^[0-9]\+:" | wc -l > $@

