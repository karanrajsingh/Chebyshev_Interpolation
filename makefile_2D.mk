CC	=g++
CFLAGS	=-c -Wall -DNDEBUG -O4 -ffast-math -ffinite-math-only -I ~/Dropbox/Eigen/
LDFLAGS	=
SOURCES	=./Chebyshev_Interpolation_1D.cpp ./Chebyshev_Interpolation_2D.cpp ./Test_Chebyshev_2D.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./Cheb2D

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.out ./*.o ./Cheb2D