/**
 * Nullity test of Gradient operator
 */

#include <iostream>
#include "mole.h"

int main() {
    int k = 2;
    int m = 2*k;
    float dx = 1;
    float tol = 1e-16;
    
    Gradient G(k, m, dx);
    vec field(m+2);
    
    vec sol = G*field;
    
    norm(sol) < tol ? cout << "\033[1;32mTest PASSED!\033[0m\n" : cout << "\033[1;31mTest FAILED!\033[0m\n";
    
    return 0;
}
