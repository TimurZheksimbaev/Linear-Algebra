#include<iomanip>
#include<iostream>
#include<math.h>

using namespace std;

class Matrix {
public:

    int n, m;
    double matrix[100][100];

    Matrix(int r, int c) {
        this->n = r;
        this->m = c;
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

    Matrix Transpose() {
        Matrix transposeMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposeMatrix.matrix[i][j] = matrix[j][i];
            }
        }
        return transposeMatrix;
    }

    Matrix Inverse() {
        Matrix id(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) id.matrix[i][j] = 1;
                else id.matrix[i][j] = 0;
            }
        }

        for (int t = 0; t < n - 1; t++) {
            /// --------------  Pivoting ---------------
            int I = t;
            double mx = matrix[t][t];
            for (int i = t + 1; i < n; i++) {
                if (abs(matrix[i][t]) > mx) {
                    mx = abs(matrix[i][t]);
                    I = i;
                }
            }
            if (I != t) {
                for (int j = 0; j < n; j++) {
                    swap(matrix[t][j], matrix[I][j]);
                    swap(id.matrix[t][j], id.matrix[I][j]);
                }
            }

            /// --------- Forward Elimination ----------
            for (int i = t + 1; i < n; i++) {
                double T = -matrix[i][t] / matrix[t][t];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] += T * matrix[t][j];
                    id.matrix[i][j] += T * id.matrix[t][j];
                }
            }
        }

        /// -------------- Way Back ---------------
        for (int t = n - 1; t >= 0; t--) {
            for (int i = t - 1; i >= 0; i--) {
                double T = -matrix[i][t] / matrix[t][t];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] += T * matrix[t][j];
                    id.matrix[i][j] += T * id.matrix[t][j];
                }
            }
        }
        /// ------- Diagonal Normalization --------
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                id.matrix[i][j] /= matrix[i][i];
            }
            matrix[i][i] = 1.0;
        }

        return id;
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
        int r;
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                output << fixed << setprecision(4) << v.matrix[i][j] << ' ';
            }
            output << '\n';
        }
        return output;
    }
};

int main() {
    int n, m;
    cin >> m;
    int t[m], b[m];
    for (int i = 0; i < m; i++) {
        cin >> t[i] >> b[i];
    }
    cin >> n;
    Matrix A(m, n + 1);

    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < A.m; j++) {
            A.matrix[i][j] = double(pow(t[i], j));
        }
    }

    Matrix B(m, 1);
    for (int i = 0; i < m; i++) {
        B.matrix[i][0] = b[i];
    }

    Matrix R1 = A.Transpose() * A;
    Matrix R2 = A.Transpose() * B;

    cout << "A:\n" << A << "A_T*A:\n" << R1;

    R1 = R1.Inverse();
    Matrix R3 = R1 * R2;

    cout << "(A_T*A)^-1:\n" << R1;
    cout << "A_T*b:\n" << R2;
    cout << "x~:\n" << R3;
}


