#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <random>
#include "full_matrix.hpp"
#include "csr_matrix.hpp"

int main() {
    std::ofstream out("data.txt");
    
    std::size_t N = 3000;
    int m = 10;
    std::vector<double> y;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> distrib_double(-1000.0, 1000.0);
    
    std::uniform_int_distribution<int> distrib_int(0, N - 1);

    full_matrix A(N, N);

    int non_zeros = 0;

    int i, j;
    for (int k = 0; k <= 1000; k++) {
        i = distrib_int(gen);
        j = distrib_int(gen);
        while (A(i, j) != 0) {
            i = distrib_int(gen);
            j = distrib_int(gen);
        }
        A(i, j) = distrib_double(gen);
    }
    CSR_matrix A_CSR(A);

    std::vector<double> x(N, 0.0);
    for (std::size_t k = 0; k < N; k++) {
        x[k] = distrib_double(gen);
    }

    clock_t start = clock();
    for (int s = 0; s < m; s++){
        y = A * x;   
    }
    clock_t end = clock();
    double t1 = (double)(end - start) / CLOCKS_PER_SEC / m;

    start = clock();
    for (int s = 0; s < m; s++){
        y = A_CSR * x;   
    }
    end = clock();
    double t2 = (double)(end - start) / CLOCKS_PER_SEC / m;

    out << 1000 << " " << t1 << " " << t2 << std::endl;

    for (int q = 0; q < 80; q++){
        for (int k = 0; k < 100000; k++) {
            i = distrib_int(gen);
            j = distrib_int(gen);
            while (A(i, j) != 0) {
                i = distrib_int(gen);
                j = distrib_int(gen);
            }
            A(i, j) = distrib_double(gen);
        }

        CSR_matrix A_CSR(A);
            
        clock_t start = clock();
        for (int s = 0; s < m; s++){
            y = A * x;   
        }
        clock_t end = clock();
        double t1 = (double)(end - start) / CLOCKS_PER_SEC / m;

        start = clock();
        for (int s = 0; s < m; s++){
            y = A_CSR * x;   
        }
        end = clock();
        double t2 = (double)(end - start) / CLOCKS_PER_SEC / m;

        non_zeros += 100000;
        out << non_zeros << " " << t1 << " " << t2 << std::endl;
        std::cout << q << std::endl;
    }

    out.close();
}
