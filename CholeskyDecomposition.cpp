#include<iomanip>
#include<iostream>
#include<cfloat>
#include<cmath>

using namespace std;

class SymmetricMatrix {
public:
    int n, m;
    int d = 1.0;
    double matrix[100][100];

    SymmetricMatrix(int row, int column) {
        this->n = row;
        this->m = column;
    }

    SymmetricMatrix operator+(SymmetricMatrix v) {
        SymmetricMatrix sumMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                sumMatrix.matrix[i][j] = matrix[i][j] + v.matrix[i][j];
            }
        }
        return sumMatrix;
    }

    SymmetricMatrix operator-(SymmetricMatrix v) {
        SymmetricMatrix subMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                subMatrix.matrix[i][j] = matrix[i][j] - v.matrix[i][j];
            }
        }
        return subMatrix;
    }

    SymmetricMatrix operator*(SymmetricMatrix v) {
        SymmetricMatrix mulMatrix(n, v.m);
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

    bool operator==(SymmetricMatrix x) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (matrix[i][j] != x.matrix[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    friend istream &operator>>(istream &input, SymmetricMatrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                input >> v.matrix[i][j];
            }
        }
        return input;
    }

    friend ostream &operator<<(ostream &output, SymmetricMatrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                output << v.matrix[i][j] << ' ';
            }
            output << '\n';
        }
        return output;
    }

    SymmetricMatrix &operator=(SymmetricMatrix v) {

        SymmetricMatrix assignMatrix(n, m);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                assignMatrix.matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }

    SymmetricMatrix Transpose() {

        SymmetricMatrix transposeMatrix(m, n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposeMatrix.matrix[i][j] = matrix[j][i];
            }
        }
        return transposeMatrix;
    }

    void CholeskyDecomposition() {
        SymmetricMatrix L(n,n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = 0;
                if (j == i) {
                    for (int k = 0; k < j; k++)
                        sum += pow(L.matrix[j][k], 2);
                    L.matrix[j][j] = sqrt(matrix[j][j] - sum);
                } else {

                    for (int k = 0; k < j; k++)
                        sum += (L.matrix[i][k] * L.matrix[j][k]);
                    L.matrix[i][j] = (matrix[i][j] - sum) / L.matrix[j][j];
                }
            }
        }

        cout << "L:\n";
        cout << L;
        cout << "L transposed:\n";
        SymmetricMatrix L_T = L.Transpose();
        cout << L_T;
    }

};


int main() {

    int n;
    cin >> n;
    SymmetricMatrix S(n, n);
    cin >> S;
    S.CholeskyDecomposition();
    return 0;
}
