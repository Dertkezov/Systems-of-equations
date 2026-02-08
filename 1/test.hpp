#pragma once
#include <iostream>
#include "solver.hpp"
#include <vector>

void test1 (){
    int N = 3;

    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    std::vector<double> b(N, 0.0);
    
    for (int i = 1; i < N; i++) {
        A[i][i - 1] = 1;
    }

    for (int i = 0; i < N; i++) {
        A[i][i] = 5;
    }

    for (int i = 1; i < N; i++) {
        A[i - 1][i] = 2;
    }

    b[0] = 9;
    b[1] = 17;
    b[2] = 17;
    
    std::vector<double> x = x_solution(A, b);

    double d = 1e-9;
    if (((1 - d <= x[0]) && (x[0] <= 1 + d)) && 
    ((2 - d <= x[1]) && (x[1] <= 2 + d)) && 
    ((3 - d <= x[2]) && (x[2] <= 3 + d))){
        std::cout << "Test for matrix 3 * 3 is correct!" << std::endl;
    }
    else{
        std::cout << "Incorrect, check solver.cpp and main.cpp!!!" << std::endl;
    }
}

void test2 (){
    int N = 5;

    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    std::vector<double> b(N, 0.0);
    
    for (int i = 1; i < N; i++) {
        A[i][i - 1] = 1;
    }

    for (int i = 0; i < N; i++) {
        A[i][i] = 5;
    }

    for (int i = 1; i < N; i++) {
        A[i - 1][i] = 2;
    }

    b[0] = 9;
    b[1] = 17;
    b[2] = 25;
    b[3] = 33;
    b[4] = 29;
    
    std::vector<double> x = x_solution(A, b);

    double d = 1e-9;
    if (((1 - d <= x[0]) & (x[0] <= 1 + d)) && 
    ((2 - d <= x[1]) && (x[1] <= 2 + d)) && 
    ((3 - d <= x[2]) && (x[2] <= 3 + d)) && 
    ((4 - d <= x[3]) && (x[3] <= 4 + d)) && 
    ((5 - d <= x[4]) && (x[4] <= 5 + d))){
        std::cout << "Test for matrix 5 * 5 is correct!" << std::endl;
    }
    else{
        std::cout << "Incorrect, check solver.cpp and main.cpp!!!" << std::endl;
    }
}

bool matrix_check(std::vector<std::vector<double>>& A){
    int N = A.size();
    std::cout << "Check of |b_i| > |a_i| + |c_i|:" << std::endl;
    if ((std::abs(A[0][0]) > std::abs(A[0][1])) && (std::abs(A[N - 1][N - 1]) > std::abs(A[N - 1][N - 2]))){
        for (int i = 1; i < N - 1; i++){
            if (std::abs(A[i][i]) <= (std::abs(A[i][i - 1]) + std::abs(A[i][i + 1]))){
                return(false);
            }
        }
    }
    else{
        return(false);
    }
    return(true);
}