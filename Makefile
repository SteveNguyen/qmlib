TARGET_LIB = libqm.a

CC = mpicc

CXXFLAGS = -I/usr/lib/openmpi/include/ -O3 -DNDEBUG #-g
LIBS = -lstdc++ -lboost_serialization -ltbb -lboost_system -lboost_mpi -lmpi -lmpi_cxx -lboost_graph_parallel

SOURCES_LIB = qmlib.cc prob_utils.cc
OBJECTS_LIB  := $(SOURCES_LIB:%.cc=%.o)
HEADERS_LIB := $(SOURCES_LIB:%.cc=%.h)

EXAMPLE_TARGET = test_pendule
EXAMPLE_SOURCE = test_pendule.cc
EXAMPLE_OBJ = test_pendule.o

.PHONY: all clean

all: lib example

lib: $(TARGET_LIB)

example: $(EXAMPLE_TARGET)

$(EXAMPLE_TARGET): $(EXAMPLE_OBJ)
	$(CC) -o $(EXAMPLE_TARGET) $(EXAMPLE_OBJ) $(TARGET_LIB) $(LIBS) -lm

$(EXAMPLE_OBJ): $(EXAMPLE_SOURCE)
	$(CC) -c $(EXAMPLE_SOURCE) -o $@ $(CXXFLAGS)

$(TARGET_LIB): $(OBJECTS_LIB)
	$(AR) rcs $@ $(OBJECTS_LIB)

$(OBJECTS_LIB): $(SOURCES_LIB) $(HEADERS_LIB)
	$(CC) -c $*.cc -o $@ $(CXXFLAGS)

clean:
	rm $(TARGET_LIB) $(OBJECTS_LIB)
