noinst_HEADERS = *.h *.hpp

bin_PROGRAMS = genome_index display_trie associate_regions icisj \
	getdna_simple interspersed_regions make_dnas_file fasta2cisfasta \
	get_uncovered_regions pileup_histogram filter_maq_pileup

getdna_simple_SOURCES = getdna_simple.cc dna.cc dnacol.cc \
	region.cc enum.cc string_tools.cc litestream.cc

get_uncovered_regions_SOURCES = get_uncovered_regions.cc

pileup_histogram_SOURCES = pileup_histogram.cc

make_dnas_file_SOURCES = make_dnas_file.cc dna.cc dnacol.cc litestream.cc string_tools.cc

filter_maq_pileup_SOURCES = filter_maq_pileup.cc

fasta2cisfasta_SOURCES = fasta2cisfasta.cc

associate_regions_SOURCES = \
	associate_regions.cc litestream.cc \
	dnacol.cc dna.cc region.cc hit.cc misc.cc \
	string_tools.cc nested.cc \
	region_association.cc enum.cc 

interspersed_regions_SOURCES = interspersed_regions.cc \
	dna.cc dnacol.cc region.cc nested.cc enum.cc string_tools.cc \
	litestream.cc

genome_index_SOURCES = genome_index.cc index_trie.cc \
	dnastats.cc dna_scanning.cc dna.cc dnacol.cc nested.cc \
	region.cc hit.cc string_tools.cc enum.cc litestream.cc

display_trie_SOURCES = display_trie.cc index_trie_scan.cc index_trie.cc \
	dna_scanning.cc dnacol.cc dna.cc litestream.cc

icisj_SOURCES = \
	icisj.cc search.cc index_trie_scan.cc index_trie.cc dna_scanning.cc \
	dna.cc nested.cc region.cc hit.cc misc.cc \
	cluster.cc diversity_cluster.cc \
	cluster_base.cc pattern.cc pwm.cc collection.cc \
	dnacol.cc enum.cc string_tools.cc litestream.cc

maqgene_binaries = getdna_simple get_uncovered_regions pileup_histogram \
	make_dnas_file filter_maq_pileup fasta2cisfasta associate_regions

maqgene: $(maqgene_binaries)

.PHONY: maqgene
.PHONY: maqgene-install


maqgene-install: maqgene
	$(INSTALL_PROGRAM) $(maqgene_binaries) $(DESTDIR)$(bindir)

# revcomp_SOURCES = revcomp.cc pwm.cc

#comp_SOURCES = comp.cc
#align_filter_SOURCES = align_filter.cc alignment.cc output.cc load.cc structs.cc annotation.cc
#site_score_SOURCES = site_score.cc structs.cc search.cc pwm.cc load.cc

# cis_phylo_SOURCES = module.cc region.cc hit.cc 

INCLUDES = @INCLUDES@
LIBS = @LIBS@
