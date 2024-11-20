CPP = g++
CFLAGS = -Wall -Wextra -O3 -std=c++14

SRCS = source/data.cpp source/functions.cpp source/matrix.cpp source/linear_algebra.cpp source/space_scheme.cpp source/time_scheme.cpp source/main.cpp
TARGET = run

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CPP) $(CFLAGS) -o source/$(TARGET) $(SRCS)


exec:
	./source/run "./input/data.dat"

clean:
	rm -f output/*.dat
	rm -f output/Implicite_Euler_cas2/*.vtk

clear:
	rm -f $(TARGET)