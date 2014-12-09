CC = gcc
CFLAGS = -std=gnu99 -g
LDFLAGS = -L.. -lcgds -lm ../synclust.so
INCLUDES = -I..
TARGET = testexec

all: testexec

clean:
	rm -f *.o
	rm -f testexec

.PHONY: all clean

SOURCES = main.c helpers.c t.connexity.c t.convexSolver.c t.dijkstra.c t.kmeansClustering.c t.neighbors.c
OBJECTS = $(SOURCES:%.c=%.o)

testexec: $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
