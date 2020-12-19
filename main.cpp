#include <iostream>
#include <utility>
#include <string>
#include <functional>
#include <vector>
#include <unistd.h>
#include "matrix.hpp"
//#include "testint.hpp"

using namespace sjtu;
int main() {
    Matrix<int> int_Matrix(20,30);
    for (int i = 0; i < int_Matrix.row_max; ++i) {
        for (int j = 0; j < int_Matrix.column_max; ++j) {
            int_Matrix(i, j) = i * int_Matrix.column_max + j;
        }
    }

    int_Matrix.print_matrix();

    int_Matrix.resize(30, 20);
    int_Matrix.print_matrix();
}