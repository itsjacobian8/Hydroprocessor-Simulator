# compiler
CC = gcc

# compiler flags
CFLAGS = -g -Wall -O3

# library flags
LIBS = -lyaml -lm

# library  and include paths (needed for homebrew yaml)
LFLAGS = -L/usr/local/lib
INCLUDES = -I/usr/local/include 

# target executable
TARGET = simulator

# source files
SRCS = simulator.c submodels.c 

# object files
OBJS = $(SRCS:.c=.o)

.PHONY: depend clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC)	-c	$< 	-o	$@

depend: $(SRCS)
	makedepend $(INCLUDES) $^

clean:
	$(RM) *.o *~ $(TARGET)