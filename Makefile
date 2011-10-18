SHELL = /bin/bash

srcdir = .
SOURCES = $(shell find $(srcdir) -name "*.cc")

#BIN = evaluate_sam_alignment fastq2fastq \
# filter_sam_by_score filter_matching_lines get_peak_regions
# test/matrix_test samutil_check sam_eval

BIN = fastq2fastq make_sq_header make_dnas_file fasta2cisfasta	\
	 shuffle_paired_fastq fastq_type

BIN += samutil align_eval deal_fastq
# BIN += sim
# BIN += pretty_plot
# BIN += gtf_annotate_regions

# BIN += filter_reads_for_tophat split_reads scatter_smoothing \
# shuffle_paired_fastq union_regions 

# BIN += test_cigar

OPT = -O0

CPPFLAGS = -I.
# LDFLAGS = -lgsl -lgslcblas -static-libgcc -L.
# LDFLAGS = -lgsl -lgslcblas
LDFLAGS = 

bindir = $(HOME)/usr/bin
OBJDIR = obj
CXXFLAGS = -ggdb3 -Wall -Wno-pointer-arith -std=c++0x -fopenmp $(OPT) $(DEBUG)


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


fastq_type_OBJS = $(addprefix $(OBJDIR)/, fastq_type.o fastq_tools.o)
fastq_type: $(fastq_type_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^


file_bytes_OBJS = $(addprefix $(OBJDIR)/, file_bytes.o file_utils.o)
file_bytes: $(file_bytes_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^


shuffle_paired_fastq_OBJS = $(addprefix $(OBJDIR)/, shuffle_paired_fastq.o)
shuffle_paired_fastq: $(shuffle_paired_fastq_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


deal_fastq_OBJS = $(addprefix $(OBJDIR)/, deal_fastq.o file_utils.o	\
dep/tools.o)

deal_fastq: $(deal_fastq_OBJS)
	$(CXX) $(CXXFLAGS) -lz -o $@ $^

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


sim_OBJS = $(addprefix $(OBJDIR)/, sim.o sim_reads.o sim_expression.o	\
	transcript_generator.o readsim_aux.o cigar_ops.o sam_order.o		\
	sam_helper.o dep/tools.o dep/simulation.o dep/nucleotide_stats.o	\
	cisortho/dna.o cisortho/dnacol.o cisortho/litestream.o				\
	dep/stats_tools.o sam_buffer.o fragment_generator.o gtf.o			\
	file_utils.o)

sim: $(sim_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lz -o $@ $^



samutil_OBJS = $(addprefix $(OBJDIR)/, samutil.o sam_tx2genome.o		\
	dep/tools.o sam_score.o sam_truncate.o sam_seqindex.o sam_rejoin.o	\
	sam_sort.o sam_checksort.o align_eval_aux.o sam_score_aux.o			\
	sam_aux.o seq_projection.o gtf.o file_utils.o sam_buffer.o			\
	sam_helper.o sam_order.o cigar_ops.o)


samutil: $(samutil_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lz -lpthread -o $@ $^



#get_peak_regions: get_peak_regions.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

pretty_plot_OBJS = $(addprefix $(OBJDIR)/, pretty_plot.o			\
	pretty_plot_graph.o pretty_plot_gtf_sam.o cigar_ops.o nclist.o	\
	sam_order.o sam_helper.o gtf.o seq_projection.o dep/tools.o		\
	meta_gene.o file_utils.o)

pretty_plot: $(pretty_plot_OBJS)
	$(CXX) $(CXXFLAGS) -lz $(LDFLAGS) -o $@ $^

#test/matrix_test: test/matrix_test.o matrix_tools.o
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


align_eval_OBJS = $(addprefix $(OBJDIR)/, align_eval.o				\
	align_eval_raw.o align_eval_aux.o align_eval_mask.o				\
	align_eval_coverage.o align_eval_stats.o cigar_ops.o			\
	seq_projection.o file_utils.o sam_helper.o sam_order.o gtf.o	\
	dep/tools.o)

align_eval: $(align_eval_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lz -o $@ $^



make_sq_header_OBJS = $(addprefix $(OBJDIR)/, make_sq_header.o)
make_sq_header : $(make_sq_header_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

gtf_annotate_regions_OBJS = $(addprefix $(OBJDIR)/, \
	gtf_annotate_regions.o nclist.o gtf.o cigar_ops.o dep/tools.o)

gtf_annotate_regions : $(gtf_annotate_regions_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


make_dnas_file_OBJS = $(addprefix $(OBJDIR)/, \
	make_dnas_file.o cisortho/dnacol.o					\
	cisortho/dna_scanning.o cisortho/dna.o cisortho/string_tools.o	\
	cisortho/litestream.o)


make_dnas_file : $(make_dnas_file_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^


fasta2cisfasta_OBJS = $(addprefix $(OBJDIR)/, fasta2cisfasta.o dep/tools.o file_utils.o)

fasta2cisfasta : $(fasta2cisfasta_OBJS)
	$(CXX) $(CXXFLAGS) -lz $(LDFLAGS) -o $@ $^


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
