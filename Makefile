all: LLL_Reduction.o main.c
	gcc -o runme main.c LLL_Reduction.o -lm

LLL_Reduction.o: LLL_Reduction.c LLL_Reduction.h
	gcc -c LLL_Reduction.c

clean: 
	rm -rf runme LLL_Reduction.o
