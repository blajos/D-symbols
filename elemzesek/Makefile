CC=g++
CFLAGS=-c -Wall -g -std=gnu++11
LDFLAGS=-ldb_cxx
SOURCES=base.cxx ddiag.cxx dsym.cxx mxfunction.cxx param.cxx simplex.cxx \
  svg.cxx dlist.cxx test.cxx main.cxx
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=a.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cxx
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(EXECUTABLE)
