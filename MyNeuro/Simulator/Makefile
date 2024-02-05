.SECONDEXPANSION:

MAINS := main.cpp
ALLSRC := $(wildcard *.cpp)
SRC := $(filter-out $(MAINS),$(ALLSRC))
ALLOBJ := $(ALLSRC:.cpp=.o)
OBJ := $(SRC:.cpp=.o)

CXX := g++
CXXFLAGS := -fopenmp -O3 -std=c++0x -w	# -w disables warnings

.PHONY: all clean
all: $(MAINS:.cpp=)

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

