# -----------------------------------------------------------------------------
# Sample Makefile for data simulation library
# -----------------------------------------------------------------------------


DATASIM_DIR = ../datasim
DATASIM_INC_DIR = -I$(DATASIM_DIR)
DATASIM_LIB_DIR = -L../../lib
BINDIR = ../../bin

# compiler
CC = /usr/bin/g++

# compiler flags
CCFLAGS = -Wall -O2 $(DATASIM_INC_DIR)

OBJS = genomeSIM.o
LIBS = -lm -ldatasim
INTERFACES = 
PROG_NAME = genomeSIM

all: $(PROG_NAME)

$(PROG_NAME): $(OBJS)
	$(CC) -o $(BINDIR)/$(PROG_NAME) $(OBJS) $(DATASIM_LIB_DIR) $(LIBS)

%.o: %.cpp $(INTERFACES)
	$(CC) -c $(CCFLAGS) $<

clean:
	rm -f *.o *~ core 


