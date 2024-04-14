.SECONDEXPANSION:

mode = 2

ifeq ($(mode),1)
	main=TransSeg
else ifeq ($(mode),2)
    main=STAR
else ifeq ($(mode),3)
	main=Soft
endif

ifeq ($(main),TransSeg)
    MAINS := TransSeg.cpp
	FILTER1 := STAR.cpp
	FILTER2 := soft.cpp
else ifeq ($(main),STAR)
    MAINS := STAR.cpp
	FILTER1 := TransSeg.cpp
	FILTER2 := soft.cpp
else
	MAINS := soft.cpp
	FILTER1 := STAR.cpp
	FILTER2 := TransSeg.cpp
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
	@echo "==> The mode is"  $(main)

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

