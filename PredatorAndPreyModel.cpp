
#include<iostream>
#include<vector>
#include<iomanip>
#include<cmath>

using namespace std;

class MyMatrix {

public:
    double matrix[100][100];
    double d = 1.0;
    int n, m;

    MyMatrix(int row, int column) {
        this->n = row;
        this->m = column;
    }


    MyMatrix operator+(MyMatrix v) {

        MyMatrix sumMatrix(n, m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                sumMatrix.matrix[i][j] = matrix[i][j] + v.matrix[i][j];
            }
        }
        return sumMatrix;
    }

    MyMatrix operator-(MyMatrix v) {

        MyMatrix subMatrix(n, m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                subMatrix.matrix[i][j] = matrix[i][j] - v.matrix[i][j];
            }
        }
        return subMatrix;
    }

    MyMatrix operator*(MyMatrix v) {

        MyMatrix mulMatrix(n, v.m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < v.m; j++) {
                for (int k = 0; k < m; k++) {
                    mulMatrix.matrix[i][j] += matrix[i][k] * v.matrix[k][j];
                }
            }
        }
        return mulMatrix;
    }

    MyMatrix Transpose() {

        MyMatrix transposeMatrix(m, n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                transposeMatrix.matrix[i][j] = matrix[j][i];
            }
        }
        return transposeMatrix;
    }

    MyMatrix &operator=(MyMatrix v) {

        MyMatrix assignMatrix(n, m);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                assignMatrix.matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }


    friend ostream &operator<<(ostream &output, MyMatrix &v) {

        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                output << fixed << setprecision(4) << v.matrix[i][j] << ' ';
            }
            output << '\n';
        }
        return output;
    }

    friend istream &operator>>(istream &input, MyMatrix &v) {
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.m; j++) {
                input >> v.matrix[i][j];
            }
        }
        return input;
    }

    double magnitude() {
        double res = 0.0;
        for (int i = 0; i < n; i++) {
            res += matrix[i][0] * matrix[i][0];
        }
        return sqrt(res);
    }


    MyMatrix Inverse() {
        MyMatrix id(n, n);
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
};

class SquareMatrix : public MyMatrix {
public:

    SquareMatrix(int n) : MyMatrix(n, n) {}

    SquareMatrix &operator=(MyMatrix v) {

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

    void printLinSys(MyMatrix b) {
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


    MyMatrix linearSystems(MyMatrix b) {

        int step = 1; // step of eliminations and permutations

        for (int k = 0; k < n; k++) {

            int p = pivot(k, k);

            if (p != k) {
                for (int j = 0; j < n; j++) {
                    swap(matrix[k][j], matrix[p][j]);
                    swap(b.matrix[k][j], b.matrix[p][j]);
                }
            }

            for (int i = k + 1; i < n; i++) {
                if (matrix[i][k] == 0.0) continue; // If we don't need elimination.

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
            }
        }


        for (int t = n - 1; t >= 0; t--) {
            for (int i = t - 1; i >= 0; i--) {
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
            }
        }

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


        return b;
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

class ColumnVector : public MyMatrix {
public:
    ColumnVector(int n) : MyMatrix(n, 1) {}

    ColumnVector operator=(MyMatrix v) {

        ColumnVector assignMatrix(n);
        for (int i = 0; i < v.n; i++) {
            for (int j = 0; j < v.n; j++) {
                matrix[i][j] = v.matrix[i][j];
            }
        }
        return assignMatrix;
    }
};

double victim(double t, double a1, double a2, double b1, double b2, double k0, double v0) {
    if (t == 0) return v0 + a2/b2;
    return v0 * cos(sqrt(a1 * a2) * t) - k0 * (sqrt(a2) * b1) / (b2 * sqrt(a1)) * sin(sqrt(a1 * a2) * t) + a2 / b2;
}

double killer(double t, double a1, double a2, double b1, double b2, double k0, double v0) {
    if (t == 0) return k0 + a1/b1;
    return v0 * (sqrt(a1) * b2) / (b1 * sqrt(a2)) * sin(sqrt(a1 * a2) * t) + k0 * cos(sqrt(a1 * a2) * t) + a1 / b1;
}


int main() {
    double k0, v0, T, N, a1, a2, b1, b2;
    cin >> v0 >> k0 >> a1 >> b1 >> a2 >> b2 >> T >> N;

    v0 -= a2 / b2;
    k0 -= a1 / b1;

    cout << "t:\n";
    for (double t = 0.0; t <= T; t += T / N) {
        cout << fixed << setprecision(2) << t << ' ';
    }

    cout << "\nv:\n";
    for (double t = 0.0; t <= T; t += T / N) {
        cout << fixed << setprecision(2) << victim(t, a1, a2, b1, b2, k0, v0) << ' ';
    }

    cout << "\nk:\n";
    for (double t = 0.0; t <= T; t += T / N) {
        cout << fixed << setprecision(2) << killer(t, a1, a2, b1, b2, k0, v0) << ' ';
    }

    return 0;
}

/*
 * input
110
40
0.4
0.01
0.3
0.005
50
200
*/