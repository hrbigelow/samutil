SHELL = /bin/bash

#### Start of system configuration section. ####

srcdir = .
prefix = $(HOME)/usr
bindir = $(prefix)/bin

CC = g++
INSTALL = /usr/bin/install -c
INSTALLDATA = /usr/bin/install -c -m 644
CPPFLAGS = -I. -I$(HOME)/usr/include
CXXFLAGS = -ggdb3 -O0 -Wall
LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -lm -lgmp

#LDFLAGS = -L$(HOME)/usr/lib -lgsl -lgslcblas -llevmar -lm -lgmpxx -lgmp -llapack -lblas -lgfortran -lcblas -latlas

SOURCES = $(shell find $(srcdir) -name "*.cc")

EXE = dep test_dirichlet quantile_test filter_sam_by_score

.PHONY : all
all : $(EXE)

dep : dep.o comp.o discomp.o \
	simp.o simc.o bqs.o bqslocus.o bqs2jpd.o metropolis.o sampling.o \
	integrands.o tools.o transformation.o error_estimate.o pileup_tools.o stats_tools.o dirichlet.o \
	slice_sampling.o hilbert.o simulation.o nucleotide_stats.o nucleotide_reader.o \
	base_qual_strand_reader.o usage_strings.o anomaly.o anomaly_tools.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

quantile_test : quantile_test.o sampling.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)


test_dirichlet : test_dirichlet.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

filter_sam_by_score : filter_sam_by_score.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)


-include $(subst .cc,.d,$(SOURCES))

define make-depend
$(CXX) -MM -MF $3 -MP -MT $2 $(CXXFLAGS) $(CPPFLAGS) $1
endef


%.o: %.cc
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@


.PHONY: install
install: $(EXE)
	$(INSTALL) $^ $(bindir)


.PHONY: clean
clean:
	rm -f *.o
