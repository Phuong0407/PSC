#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <iostream>
#include <stdexcept>

class matrix{
    public:
    unsigned int row = 0, col = 0;
    double* elements;


    inline unsigned int index(unsigned int i, unsigned int j) const {
        if (i >= row || j >= col)
            throw std::out_of_range("Index out of bounds");
        return i * col + j;
    }


    matrix(unsigned int rows, unsigned int cols) : row(rows), col(cols) {
        elements = new double[row * col];
        for (unsigned int i = 0; i < row * col; ++i)
            elements[i] = 0.0;
    }


    ~matrix() {
        delete[] elements;
    }


    matrix(const matrix& other) : row(other.row), col(other.col) {
        elements = new double[row * col];
        for (unsigned int i = 0; i < row * col; ++i)
            elements[i] = other.elements[i];
    }


    matrix& operator=(const matrix& other) {
        if (this == &other)
            return *this;
        
        delete[] elements;

        row = other.row;
        col = other.col;
        elements = new double[row * col];
        for (unsigned int i = 0; i < row * col; ++i)
            elements[i] = other.elements[i];
        return *this;
    }


    void print() const {
        for (unsigned int i = 0; i < row; ++i) {
            for (unsigned int j = 0; j < col; ++j) {
                std::cout << elements[index(i, j)] << "  ";
            }
            std::cout << std::endl;
        }
    }


    matrix add(const matrix& other) const {
        if (row != other.row || col != other.col)
            throw std::invalid_argument("matrix dimensions must match for addition.");

        matrix result(row, col);
        for (unsigned int i = 0; i < row; ++i) {
            for (unsigned int j = 0; j < col; ++j) {
                result.elements[index(i, j)] = elements[index(i, j)] + other.elements[index(i, j)];
            }
        }
        return result;
    }


    matrix multiply(const matrix& other) const {
        if (col != other.row)
            throw std::invalid_argument("matrix dimensions are incompatible for multiplication.");

        matrix result(row, other.col);
        for (unsigned int i = 0; i < row; ++i) {
            for (unsigned int j = 0; j < other.col; ++j) {
                for (unsigned int k = 0; k < col; ++k) {
                    result.elements[result.index(i, j)] += elements[index(i, k)] * other.elements[other.index(k, j)];
                }
            }
        }
        return result;
    }
};

#endif