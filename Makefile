CC = gcc
CFLAGS = -Wall -Wextra -g
DEPS = LLL_Reduction.h Enumeration.h


all: LLL_Reduction.o main.c
	gcc -o runme main.c LLL_Reduction.o -lm

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(FLAGS)

#LLL_Reduction.o: LLL_Reduction.c LLL_Reduction.h
#	gcc -c LLL_Reduction.c

test: 

clean: 
	rm -rf runme LLL_Reduction.o

#need a test option
