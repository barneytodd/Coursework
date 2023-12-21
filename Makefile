CC = gcc
CFLAGS = -Wall -Wextra -g
DEPS = LLL_Reduction.h Enumeration.h
SOURCE_FILES = LLL_Reduction.c Enumeration.c main.c
OBJECTS = LLL_Reduction.o Enumeration.o main.o

all: $(OBJECTS)
	$(CC) -o runme main.c LLL_Reduction.o -lm

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

#LLL_Reduction.o: LLL_Reduction.c LLL_Reduction.h
#	gcc -c LLL_Reduction.c

test: 

clean: 
	rm -rf runme LLL_Reduction.o

#need a test option
