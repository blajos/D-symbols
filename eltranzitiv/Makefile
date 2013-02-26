CC=g++
CFLAGS=-c -Wall -g -O2 -std=gnu++11
LDFLAGS=-ldb_cxx
SOURCES=backtrack.cxx
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=a.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
	objcopy --only-keep-debug $@ $@.dbg
	objcopy --strip-debug $@
	objcopy --add-gnu-debuglink=$@.dbg $@

%.o: %.cxx
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(EXECUTABLE)
