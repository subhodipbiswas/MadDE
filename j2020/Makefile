

HEADER = 
TARGET = main

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
 CC = g++

OPTION = -Ofast -Wall

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(OPTION) -lm 

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) $(OPTION)  -c $<

clean:
	rm -rf *.o
