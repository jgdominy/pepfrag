SOURCE_FILES = $(wildcard *.fasta)
RES_FILES = $(SOURCE_FILES:.fasta=_naive.results)
RES_FILES += $(SOURCE_FILES:.fasta=_p56.results)
RES_FILES += $(SOURCE_FILES:.fasta=_p84.results)
RES_FILES += $(SOURCE_FILES:.fasta=_p112.results)
RES_FILES += flus.results

test: $(RES_FILES)

flus.results: flus.fasta
	@python ../pepfrag.py -l 8 -L 15 -o 11 -O 13 -C 32 $< | grep -v "^[0-9]\+:" | wc -l > $@

%_naive.results: %.fasta
	@python ../pepfrag.py -L 30 -o 2 -n $< | grep -v "^[0-9]\+: " | wc -l > $@
	
%_p56.results: %.fasta
	@python ../pepfrag.py -l 8 -L 30 -o 2 -O 4 -C 56 -t $< | grep -v "^[0-9]\+: " | wc -l > $@
	
%_p84.results: %.fasta
	@python ../pepfrag.py -l 8 -L 30 -o 2 -O 4 -C 84 -t $< | grep -v "^[0-9]\+: " | wc -l > $@
	
%_p112.results: %.fasta
	@python ../pepfrag.py -l 8 -L 30 -o 2 -O 4 -C 112 -t $< | grep -v "^[0-9]\+: " | wc -l > $@
