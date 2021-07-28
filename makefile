CC=mpic++
CFLAGS=-Wall -std=c++11 -fopenmp
OBJECTS=$(SOURCES:.cpp=.o)

SOURCES=main.cpp savebmp.cpp

EXECUTABLE=project.x

all: 
	 $(CC) $(CFLAGS) -o $(EXECUTABLE) $(SOURCES) 

clean:
	rm *.o *.x *.bmp
