# Makefile for GEDI simulator tools

LIBS = -lm 
INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) 
CFLAGS += -Wall
CFLAGS += -O3
#CFLAGS += -g
ARCH=$(shell uname -m)

CC = gcc

THIS=chooseMetricPrints

$(THIS):	$(THIS).o 
		$(CC) $(CFLAGS) $@.o -o $@ $(LIBS) $(CFLAGS) $(INCLS)

.c.o:		$<
		$(CC) $(CFLAGS) -I. $(INCLS) -D$(ARCH)  -c $<

clean:
		\rm -f *% *~ *.o

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/$(THIS).old
		cp $(THIS) $(HOME)/bin/$(ARCH)/$(THIS)


