.PHONY: clean
SAMTOOLDIR=libs/samtools-0.1.19/
BOOST=libs/boost_1_45_0_subset/

SRCDIR=src
OBJDIR=obj
BINDIR=bin
DEF= -DSAMTOOLS=\"$(shell pwd)/libs/samtools-0.1.19/samtools\"
CXX=g++

CPPFLAGS=#-DLOGDEBUG -DDEBUGREADS #-DNDEBUG -DMMTEST
CXXFLAGS=-I$(SAMTOOLDIR) -I$(BOOST) -I$(INCLUDEDIR) -I$(SRCDIR) -Wno-deprecated -O3

LDFLAGS= -L$(SAMTOOLDIR) -lbam -lz -lpthread

SOURCES := $(shell find $(SRCDIR) -name "*.cpp")
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)


all: pre-build main-build

pre-build:
	mkdir -p obj bin obj/haplotype_phaser

main-build: ohvarfinder

ohvarfinder: $(OBJECTS)
	$(CXX) -o $(BINDIR)/$@ $(CXXFLAGS) $(CPPFLAGS) $(DEF) $(OBJECTS) $(LDFLAGS)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEF) -c $< -o $@

dependencies:
	cd libs; make

clean:
	rm -f $(OBJECTS)

test:
	cd tests;make test
