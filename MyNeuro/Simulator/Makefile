.SECONDEXPANSION:

mode = 2

ifeq ($(mode),1)
	main=MyNeuro
else ifeq ($(mode),2)
    main=Retransformer
else ifeq ($(mode),3)
	main=Soft
endif

ifeq ($(main),MyNeuro)
    MAINS := main.cpp
	FILTER1 := Retransformer.cpp
	FILTER2 := soft.cpp
else ifeq ($(main),Retransformer)
    MAINS := Retransformer.cpp
	FILTER1 := main.cpp
	FILTER2 := soft.cpp
else
	MAINS := soft.cpp
	FILTER1 := Retransformer.cpp
	FILTER2 := main.cpp
endif

ALLSRC := $(wildcard *.cpp)
SRC := $(filter-out $(FILTER1) $(FILTER2) $(MAINS),$(ALLSRC))
ALLOBJ := $(ALLSRC:.cpp=.o)
OBJ := $(SRC:.cpp=.o)

CXX := g++
CXXFLAGS := -fopenmp -O3 -std=c++0x -w	# -w disables warnings

.PHONY: all clean Mode
# all: Mode
all: Mode $(MAINS:.cpp=)
Mode:
	@echo "The mode is"  $(main)

$(MAINS:.cpp=): $(OBJ) $$@.o
	$(CXX) $(CXXFLAGS) $^ -o $@
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

depend: .depend
.depend: $(ALLSRC)
	@$(RM) .depend
	@$(CXX) -MM $(CXXFLAGS) $^ > .depend;
include .depend

clean:
	$(RM) $(MAINS:.cpp=)
	$(RM) $(ALLOBJ)

