#include<iomanip>
#include<iostream>
#include<cfloat>
#include<cmath>

using namespace std;

class Matrix {
public:
    int n, m;
    int d = 1.0;
    double matrix[100][100];

    Matrix(int row, int column) {
        this->n = row;
        this->m = column;
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
                output << fixed << setprecision(4) << v.matrix[i][j] << ' ';
            }
            output << '\n';
        }
        return output;
    }

//    Matrix &operator=(Matrix v) {
//
//        Matrix assignMatrix(n, m);
//        for (int i = 0; i < v.n; i++) {
//            for (int j = 0; j < v.m; j++) {
//                assignMatrix.matrix[i][j] = v.matrix[i][j];
//            }
//        }
//        return assignMatrix;
//    }

    Matrix Transpose() {

        Matrix transposeMatrix(m, n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposeMatrix.matrix[i][j] = matrix[j][i];
            }
        }
        return transposeMatrix;
    }

    double Norm() {
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            norm += pow(matrix[i][0], 2);
        }
        return sqrt(norm);
    }
};

class SquareMatrix : public Matrix {
public:

    SquareMatrix(int n) : Matrix(n, n) {}

    SquareMatrix &operator=(Matrix v) {

        SquareMatrix assignMatrix(n);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.n; j++) {
                matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }

    int pivot(int pi, int pj) {
        int pivot = pi;
        double mx = matrix[pi][pj];
        for (int i = pi + 1; i < n; i++) {
            if (abs(matrix[i][pj]) > mx) {
                mx = abs(matrix[i][pj]);
                pivot = i;
            }
        }
        return pivot;
    }


    /// (((((((((((((((((((((((((((((((((((((( DETERMINANT )))))))))))))))))))))))))))))))))))))))))))))))))))))))))

//    void pivot(int row, int column) {
//        int piv = row;
//        double maxnum = matrix[row][column];
//
//        for (int i = row + 1; i < n; i++) {
//            if (abs(matrix[i][column]) > maxnum) {
//                maxnum = abs(matrix[i][column]);
//                piv = i;
//            }
//        }
//
//        if (piv >= row) {
//            cout << "step #" << step++ << ": permutation\n";
//            for (int j = 0; j < n; j++)
//                swap(matrix[row][j], matrix[piv][j]);
//            d *= -1;
//            print();
//        }
//    }

    void determinant() {
        int step = 1;
        double determinant = 1.0;
        for (int k = 0; k < n - 1; k++) {

            pivot(k, k);

            if (matrix[k][k] == 0.0) {
                cout << "result:\n0.00\n";
                exit(0);
            }
            determinant *= matrix[k][k];
            for (int i = k + 1; i < n; i++) {

                if (matrix[i][k] == 0.0) continue;

                cout << "step #" << step++ << ": elimination\n";

                double multiplier = -matrix[i][k] / matrix[k][k];

                for (int j = k; j < n; j++)
                    matrix[i][j] += multiplier * matrix[k][j];
            }
        }
        determinant *= matrix[n - 1][n - 1];
        cout << "result:\n" << determinant * d << '\n';
    }

    /// (((((((((((((((((((((((((((((((((((((((((( INVERSE )))))))))))))))))))))))))))))))))))))))))))))))))))))))))

    SquareMatrix Inverse() {
        SquareMatrix id(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) id.matrix[i][j] = 1.0;
                else id.matrix[i][j] = 0.0;
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


    /// (((((((((((((((((((((((((((((((((((((((( LINEAR SYSTEMS ))))))))))))))))))))))))))))))))))))))))))))))))))))

    void printLinSys(Matrix b) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (abs(matrix[i][j]) < 0.01) matrix[i][j] = 0.0;
                cout << fixed << setprecision(2) << matrix[i][j] << " ";
            }
            cout << "\n";
        }
        for (int i = 0; i < n; i++) {
            if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0;
            cout << fixed << setprecision(2) << b.matrix[i][0] << "\n";
        }
    }


    void linearSystems(Matrix b) {

        int step = 1; // step of eliminations and permutations

        cout << "step #0:\n";
        printLinSys(b);
        for (int k = 0; k < n; k++) {

            int p = pivot(k, k);

            if (p != k) {
                cout << "step #" << step++ << ": permutation\n";
                for (int j = 0; j < n; j++) {
                    swap(matrix[k][j], matrix[p][j]);
                    swap(b.matrix[k][j], b.matrix[p][j]);
                }
                printLinSys(b);
            }

            for (int i = k + 1; i < n; i++) {
                if (matrix[i][k] == 0.0) continue; // If we don't need elimination.
                cout << "step #" << step++ << ": elimination\n";

                double multiplier = matrix[i][k] / matrix[k][k];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] -= multiplier * matrix[k][j];
                    if (!j) b.matrix[i][j] -= multiplier * b.matrix[k][j];


                    if (signbit(matrix[i][j])) {
                        if (j == i) matrix[i][i] = 1.0;
                        else matrix[i][j] = 0.0;
                        if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0;
                    }
//                    } else {
//                        if (j == i) matrix[i][i] = 1.0;
//                        else matrix[i][j] = 0.0;
//                        if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0;
//                    }

                }
                printLinSys(b);
            }
        }


        for (int t = n - 1; t >= 0; t--) {
            for (int i = t - 1; i >= 0; i--) {
                cout << "step #" << step++ << ": elimination\n";
                double multiplier = matrix[i][t] / matrix[t][t];
                for (int j = 0; j < n; j++) {
                    matrix[i][j] -= multiplier * matrix[t][j];
                    if (!j) b.matrix[i][j] -= multiplier * b.matrix[t][j];


                    if (signbit(matrix[i][j])) {
                        if (j == i) matrix[i][i] = 1.0;
                        else matrix[i][j] = 0.0;
                        if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0;
                    }
//                    } else {
//                        if (j == i) matrix[i][i] = 1.0;
//                        else matrix[i][j] = 0.0;
//                        if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0;
//                    }


                }
                printLinSys(b);
            }
        }

        cout << "Diagonal normalization:\n";
        for (int i = 0; i < n; i++) {
            b.matrix[i][0] /= matrix[i][i];


            if (signbit(matrix[i][i])) {
                for (int j = 0; j < n; j++) {
                    if (j == i) matrix[i][i] = 1.0;
                    else matrix[i][j] = 0.0;
                }
                if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0;
            } else {
                for (int j = 0; j < n; j++) {
                    if (j == i) matrix[i][i] = 1.0;
                    else matrix[i][j] = 0.0;
                }
                if (abs(b.matrix[i][0]) < 0.01) b.matrix[i][0] = 0.0;
            }
        }
        printLinSys(b);

        cout << "result:\n";
        cout << b;
    }

};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i == j)
                    matrix[i][j] = 1;
                else
                    matrix[i][j] = 0;
    }

};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(int p1, int p2, int n) : IdentityMatrix(n) {
        p1--;
        p2--;
        for (int i = 0; i < n; i++) {
            swap(matrix[p1][i], matrix[p2][i]);
        }
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix(int e1, int e2, int n, SquareMatrix v) : IdentityMatrix(n) {
        e1--;
        e2--;
        matrix[e1][e2] -= (v.matrix[e1][e2] / v.matrix[e1 - 1][e2]);
    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector(int n) : Matrix(n, 1) {}

    ColumnVector operator=(Matrix v) {

        ColumnVector assignMatrix(n);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.n; j++) {
                matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }

};

int main() {
    int n, m;
    cin >> n;
    SquareMatrix A(n), alpha(n);
    Matrix *pA = &A;
    Matrix *palpha = &alpha;
    cin >> A;


    cin >> m;
    ColumnVector b(m), beta(m);
    Matrix *pb = &b;
    Matrix *pbeta = &beta;
    cin >> b;

    double e;
    cin >> e;


    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) sum += A.matrix[i][j];
        }
        if (sum > A.matrix[i][i]) {
            cout << "The method is not applicable!";
            exit(0);
        }
    }

    SquareMatrix D(n);
    Matrix *pD = &D;
    IdentityMatrix I(n);
    Matrix *pI = &I;


    // inverse of D and Identity matrix
    for (int i = 0; i < n; i++) {
        D.matrix[i][i] = 1.0 / A.matrix[i][i];
    }

    alpha = I - (D * A);
    beta = D * b;

    cout << "alpha:\n" << alpha << "beta:\n" << beta;

    int index = 1;
    double error = DBL_MAX, new_error = DBL_MAX;
    Matrix x = beta;

    cout << "x(0):\n" << beta;

    while (error >= e) {
        Matrix new_x = x;
        x = alpha * new_x + beta;
        new_x = new_x - x;
        error = new_x.Norm();
        cout << "e: " << error << "\n";
        cout << "x(" << index++ << "):\n" << x;
    }
    return 0;
}


