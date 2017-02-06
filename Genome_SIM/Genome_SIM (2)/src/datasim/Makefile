# Makefile for datasim library
CC = g++
LIBDIR = ../../lib/
LIB= $(LIBDIR)libdatasim.a

CFLAGS= -O2 -Wall -Wno-deprecated

OBJS = Locus.o Model.o Individual.o \
  DatasimExcept.o Population.o ProbSim.o \
  DataSimPreferences.o Stringmanip.o \
  Sim.o PopSim.o Chromosome.o Chrompool.o \
  IndPopulation.o

  
INTERFACES = Locus.h Model.h Individual.h Ran2.h \
  DatasimExcept.h Population.h ProbSim.h \
  DataSimPreferences.h Stringmanip.h \
  Sim.h PopSim.h Chromosome.h Chrompool.h \
  Enums.h IndPopulation.h


# GNU compiler flags (works with g++ and egcs)
CC          = g++
PRELINK     =
LD          = g++ -w
AR          = ar rv
INSTALL     = install -c
RM          = rm

# INC_DIRS= -I..

%.o: %.cpp $(INTERFACES)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB): $(OBJS) 
	$(PRELINK)
	$(AR) $(LIB) $?
	@echo $(LIB) is now up-to-date
.cpp.a:;

clean:
	rm -f *.o *~ core 
