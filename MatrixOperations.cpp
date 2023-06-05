#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

class Matrix {
protected:
    vector<vector<double>> matrix;
    int n, m;
public:
    Matrix() {}

    Matrix(int row, int column) {
        this->n = row;
        this->m = column;
        for (int i = 0; i < n; i++)
            matrix[i].resize(m, 0);
    }

    Matrix operator+(Matrix v) {
        Matrix sumMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                sumMatrix.matrix[i][j] = matrix[i][j] + v.matrix[i][j];
            }
        }
        return sumMatrix;
    }

    Matrix operator-(Matrix v) {
        Matrix subMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                subMatrix.matrix[i][j] = matrix[i][j] - v.matrix[i][j];
            }
        }
        return subMatrix;
    }

    Matrix operator*(Matrix v) {
        Matrix mulMatrix(n, v.m);
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < v.m; i++) {
                double sum = 0.0;
                for (int j = 0; j < m; j++) {
                    sum += matrix[k][j] * v.matrix[j][i];
                }
                mulMatrix.matrix[k][i] = sum;
            }
        }
        return mulMatrix;
    }

    friend istream &operator>>(istream &input, Matrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                input >> v.matrix[i][j];
            }
        }
        return input;
    }

    friend ostream &operator<<(ostream &output, Matrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                output << v.matrix[i][j] << ' ';
            }
            output << '\n';
        }
        return output;
    }

    Matrix &operator=(Matrix v) {

        Matrix assignMatrix(n, m);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                assignMatrix.matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }

    Matrix Transpose() {

        Matrix transposeMatrix(m, n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposeMatrix.matrix[i][j] = matrix[j][i];
            }
        }
        return transposeMatrix;
    }
};


class Vector : public Matrix {
public:
    Vector(int n) : Matrix(n, 1) {}

    double Norm() {
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            norm += pow(matrix[i][0], 2);
        }
        return sqrt(norm);
    }
};

int main() {
    int n, m;
    cin >> n >> m;
    Matrix A(n, m);
    Vector b(n);
    cin >> A;
    cin >> b;
    Vector x = A * b;
    cout << x;
    return 0;
}
