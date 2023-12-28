CC = gcc
CFLAGS = -Wall -Wextra -g
DEPS = LLL_Reduction.h Enumeration.h
SOURCE_FILES = LLL_Reduction.c Enumeration.c main.c  test.c
OBJECTS = LLL_Reduction.o Enumeration.o main.o #test.o

all: $(OBJECTS)
	$(CC) -o runme $(OBJECTS) -lm $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: test.o $(OBJECTS)
	$(CC) -o testme $(OBJECTS) -lm $(CFLAGS)
	./testme

clean: 
	rm -rf runme testme *.o result.txt
