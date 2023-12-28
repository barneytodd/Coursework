#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void runTests(int n, ull expr)
{
    ull act = perimeter(n);
    printf("For numbers: %d Expected: %llu Got: %llu\n", n, expr, act);
    if(act != expr)
        printf("Expected %llu, got %llu\n", expr, act);
    assert(act == expr);
}

int main() {
    runTests(5, 80);
    runTests(7, 216);
    runTests(20, 114624);
    runTests(30, 14098308);
    return 0;
