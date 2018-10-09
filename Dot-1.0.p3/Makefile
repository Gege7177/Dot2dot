CC = gcc 
CFLAGS = -O3 -finline-functions -fomit-frame-pointer -funroll-loops -pthread -Wall -pedantic -g -D NO_GRAPH 

HEADERS = $(wildcard *.h)
SOURCES = $(wildcard *.c)
OBJECTS= $(SOURCES:.c=.o)

exe = dot

all: .depend $(OBJECTS) $(exe)

dot: $(OBJECTS)
	$(CC) -pthread $(OBJECTS) -o $(exe) 

.PHONY: clean 

clean:
	@-rm $(OBJECTS)
	@-rm -f *~ .depend $(exe) 

.depend:
	$(CC) -MM $(SOURCES) > .depend

ifeq ($(wildcard .depend), .depend)
include .depend
endif
# DO NOT DELETE
