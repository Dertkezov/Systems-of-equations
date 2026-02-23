#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include "full_matrix.hpp"
#include "csr_matrix.hpp"

int main() {
    std::ofstream out("data.txt");
    
    std::size_t N = 1000;
    int m = 10;
    std::vector<double> y;
    
    for (int s = 100; s <= 10000; s += 100) {
        full_matrix A(N, N);
        int k = 0;
        
        for (std::size_t i = 0; i < N; i++) {
            for (std::size_t j = 0; j < N; j++) {
                if ((i * N + j) % s == 0) {
                    A(i, j) = 0.1319;
                    k++;
                }
            }
        }
        
        std::vector<double> x(N, 1.0);

        CSR_matrix A_CSR(A);
        
        clock_t start = clock();
        for (int i = 0; i < m; i++){
            y = A * x;   
        }
        clock_t end = clock();
        double t1 = (double)(end - start) / CLOCKS_PER_SEC / m;

        start = clock();
        for (int i = 0; i < m; i++){
            y = A_CSR * x;   
        }
        end = clock();
        double t2 = (double)(end - start) / CLOCKS_PER_SEC / m;

        out << k << " " << t1 << " " << t2 << std::endl;
    }
    
    out.close();
}
