SHELL = /bin/bash

srcdir = .
SOURCES = $(shell find $(srcdir) -name "*.cc")

#BIN = evaluate_sam_alignment fastq2fastq shuffle_paired_fastq \
# filter_sam_by_score filter_matching_lines get_peak_regions
# test/matrix_test samutil_check sam_eval

BIN = shuffle_paired_fastq fastq2fastq union_regions align_eval			\
	 pretty_plot make_sq_header gtf_annotate_regions					\
	 filter_reads_for_tophat make_dnas_file fasta2cisfasta split_reads	\
	 scatter_smoothing

BIN += sim
#BIN += test_cigar
BIN += samutil


OPT = -O0

CPPFLAGS = -I.
LDFLAGS = -lgsl -lgslcblas


ifeq ($(findstring el5,$(shell uname -r)), el5)
bindir = $(HOME)/usr_el5/bin
OBJDIR = obj_el5
CXXFLAGS = -ggdb3 -Wall $(OPT) $(DEBUG)
else
bindir = $(HOME)/usr/bin
OBJDIR = obj
CXXFLAGS = -ggdb3 -Wall -std=c++0x $(OPT) $(DEBUG)
endif






.PHONY: all clean

all: $(BIN)


#FIRST=$(shell echo $(MAKECMDGOALS) | cut -f 1 -d '_')
#SECOND=$(shell echo $(MAKECMDGOALS) | cut -f 2 -d '_')
#$(foreach v,$(.VARIABLES),$(info $(v)=$($(v))))

#clean:
#	-rm $(BIN) *.o *.d cisortho/*.o cisortho/*.d dep/*.o dep/*.d

#sam_eval_OBJS = $(addprefix $(OBJDIR)/, sam_eval.o sam_helper.o \
#	cigar_ops.o dep/alignment_stats.o dep/tools.o file_utils.o)

#sam_eval : $(sam_eval_OBJS)
#	$(CXX) $(CXXFLAGS) -lz -o $@ $^


fastq2fastq_OBJS = $(addprefix $(OBJDIR)/, fastq2fastq.o fastq_tools.o)
fastq2fastq: $(fastq2fastq_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

file_bytes_OBJS = $(addprefix $(OBJDIR)/, file_bytes.o file_utils.o)
file_bytes: $(file_bytes_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^


shuffle_paired_fastq_OBJS = $(addprefix $(OBJDIR)/, shuffle_paired_fastq.o)
shuffle_paired_fastq: $(shuffle_paired_fastq_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


#filter_sam_by_score: filter_sam_by_score.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

#filter_matching_lines: filter_matching_lines.o cisortho/string_tools.o	\
#	line_tools.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

union_regions_OBJS = $(addprefix $(OBJDIR)/, union_regions.o dep/tools.o)
union_regions: $(union_regions_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

#sam_stats: sam_stats.o sam_stats_raw.o sam_stats_out.o sam_stats_aux.o	\
#	matrix_tools.o dep/tools.o cigar_ops.o sam_helper.o					\
#	dep/pileup_tools.o dep/nucleotide_stats.o dep/stats_tools.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


sim_OBJS = $(addprefix $(OBJDIR)/,\
	sim.o sim_reads.o sim_expression.o transcript_generator.o	\
	readsim_aux.o cigar_ops.o sam_helper.o dep/tools.o			\
	dep/simulation.o dep/nucleotide_stats.o cisortho/dna.o		\
	cisortho/dnacol.o cisortho/litestream.o dep/stats_tools.o	\
	sam_buffer.o fragment_generator.o gtf.o file_utils.o)

sim: $(sim_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lz -o $@ $^



samutil_OBJS = $(addprefix $(OBJDIR)/, samutil.o						\
	sam_transcript_to_genome.o dep/tools.o	\
	gtf.o cisortho/dna.o sam_score_mapq.o sam_raw_score_aux.o			\
	align_eval_raw.o sam_raw_score_dist.o cisortho/region.o				\
	readsim_aux.o file_utils.o sam_buffer.o sam_helper.o cigar_ops.o	\
	cisortho/nested.o dep/nucleotide_stats.o get_spliced_sequence.o		\
	dep/stats_tools.o cisortho/litestream.o cisortho/enum.o				\
	cisortho/dnacol.o)

samutil: $(samutil_OBJS)
	$(CXX) $(CXXFLAGS) -lz -lgsl -lgslcblas -o $@ $^



#get_peak_regions: get_peak_regions.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

pretty_plot_OBJS = $(addprefix $(OBJDIR)/, \
	pretty_plot.o pretty_plot_graph.o pretty_plot_gtf_sam.o	\
	cigar_ops.o nclist.o sam_helper.o gtf.o dep/tools.o meta_gene.o		\
	file_utils.o)

pretty_plot: $(pretty_plot_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^

#test/matrix_test: test/matrix_test.o matrix_tools.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


align_eval_OBJS = $(addprefix $(OBJDIR)/, \
	align_eval.o align_eval_raw.o align_eval_sort.o				\
	align_eval_aux.o align_eval_checksort.o align_eval_mask.o			\
	align_eval_coverage.o align_eval_stats.o cigar_ops.o file_utils.o	\
	sam_helper.o dep/tools.o)

align_eval: $(align_eval_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^


#samutil_check: samutil_check.o sam_helper.o cigar_ops.o file_utils.o
#	$(CXX) $(CXXFLAGS) -lz -o $@ $^


#test_cigar: test_cigar.o cigar_ops.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


make_sq_header_OBJS = $(addprefix $(OBJDIR)/, make_sq_header.o)
make_sq_header : $(make_sq_header_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

gtf_annotate_regions_OBJS = $(addprefix $(OBJDIR)/, \
	gtf_annotate_regions.o nclist.o gtf.o cigar_ops.o dep/tools.o)

gtf_annotate_regions : $(gtf_annotate_regions_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^


filter_reads_for_tophat_OBJS = $(addprefix $(OBJDIR)/, \
	filter_reads_for_tophat.o file_utils.o fastq_tools.o dep/tools.o)

filter_reads_for_tophat : $(filter_reads_for_tophat_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^


split_reads_OBJS = $(addprefix $(OBJDIR)/, split_reads.o cisortho/string_tools.o \
	dep/tools.o file_utils.o)
split_reads : $(split_reads_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^


make_dnas_file_OBJS = $(addprefix $(OBJDIR)/, \
	make_dnas_file.o cisortho/dnacol.o					\
	cisortho/dna_scanning.o cisortho/dna.o cisortho/string_tools.o	\
	cisortho/litestream.o)

make_dnas_file : $(make_dnas_file_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^


fasta2cisfasta_OBJS = $(addprefix $(OBJDIR)/, fasta2cisfasta.o dep/tools.o file_utils.o)

fasta2cisfasta : $(fasta2cisfasta_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^


scatter_smoothing_OBJS = $(addprefix $(OBJDIR)/, \
	scatter_smoothing.o file_utils.o dep/tools.o cisortho/string_tools.o)

scatter_smoothing : $(scatter_smoothing_OBJS)
	$(CXX) $(CXXFLAGS) -lz $(LDFLAGS) -o $@ $^


-include $(subst .cc,$(OBJDIR)/.d,$(SOURCES))

define make-depend
$(CXX) -MM -MF $3 -MP -MT $2 $(CXXFLAGS) $(CPPFLAGS) $1
endef


$(OBJDIR)/%.o: %.cc
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@


.PHONY: install

install: all
	install $(BIN) $(bindir)
