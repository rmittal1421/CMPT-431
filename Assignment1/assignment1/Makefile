
ifdef USE_INT
MACRO = -DUSE_INT
endif

#compiler setup
PCC = g++
PCFLAGS = -std=c++14 -O3 -pthread $(MACRO)

COMMON= utils.h cxxopts.h get_time.h graph.h quick_sort.h
SERIAL= pi_calculation triangle_counting page_rank
PARALLEL= pi_calculation_parallel triangle_counting_parallel page_rank_parallel page_rank_parallel_atomic
ALL= $(SERIAL) $(PARALLEL)


all : $(ALL)

% : %.cpp $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

$(COMMON):
	ln -s core/$@ .

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
	rm $(COMMON)
