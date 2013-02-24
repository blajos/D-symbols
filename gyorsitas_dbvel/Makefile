CC=g++
CFLAGS=-c -Wall -g -std=gnu++11
LDFLAGS=-ldb_cxx
SOURCES=backtrack.cxx
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=a.out

all: $(SOURCES) $(EXECUTABLE)
	rm -f data/*

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cxx
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(EXECUTABLE)
