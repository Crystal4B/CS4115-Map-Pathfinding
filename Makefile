#
# 
# Makefile  for CS4115 labs


LIBS=
SOURCES= streets.cc


CC=gcc
CFLAGS=-g -Wall

CXX=g++
CXXFLAGS=$(CFLAGS)

LD=g++
LDFLAGS=


OBJECTS=${SOURCES:.cc=.o}
EXECS=${SOURCES:.cc=}


.PHONY: clean
.SUFFIXES: .o .c .cc


all: $(EXECS)

clean:
	$(RM) $(OBJECTS) streets

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<

streets: $(OBJECTS)
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

glut7: glut7.o
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

glut8: glut8.o
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
