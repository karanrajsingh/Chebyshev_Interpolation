CC	=g++
CFLAGS	=-c -Wall -DNDEBUG -O4 -ffast-math -ffinite-math-only -I ~/Dropbox/Eigen/
LDFLAGS	=
SOURCES	=./Chebyshev_Interpolation_1D.cpp ./Test_Chebyshev_1D.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./Cheb1D

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.out ./*.o ./Cheb1D