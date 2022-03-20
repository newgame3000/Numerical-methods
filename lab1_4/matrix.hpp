#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>

using namespace std;

template <class T>
class matrix {
private: 
    vector<vector<T>> data;
public:

    matrix(int n, int m) {
        data = vector<vector<T>>(n, vector<T>(m));
    }

    matrix(const matrix<T> &a) {
        data = a.data;
    }

    vector<T>& operator[](int i) {
        return data[i];
    }

    void Read() {
        for (uint i = 0; i < data.size(); ++i) {
            for (uint j = 0; j < data[i].size(); ++j) {
                cin >> data[i][j];
            }
        }
    }

    void Print() {
        cout.precision(3);
        cout << fixed;
        for (uint i = 0; i < data.size(); ++i) {
            for (uint j = 0; j < data[i].size(); ++j) {
                cout << setw(6) << data[i][j] << " ";
            }
            cout << endl;
        }
    }

    uint Rows_count() const{
        return data.size();
    }

    uint Columns_count() const{
        return data[0].size();
    }

    void Swap_rows(int i, int j) {
        vector<T> s = data[i];
        data[i] = data[j];
        data[j] = s;
    }

    void Swap_columns(int i, int j) {
        vector<T> s(data.size());
        for (uint k = 0; k < data.size(); ++k) {
            s[k] = data[k][i];
        }

        for (uint k = 0; k < data.size(); ++k) {
            data[k][i] = data[k][j];
        }

        for (uint k = 0; k < data.size(); ++k) {
            data[k][j] = s[k];
        }    
    }

    void Fill_in_the_column(vector<T> a, int k) {
        for (uint i = 0; i < a.size(); ++i) {
            data[i][k] = a[i];
        }
    }

    double Norm() {
        double n = 0;
        for (uint j = 0; j < data[0].size(); ++j) {
            n += abs(data[0][j]); 
        }

        for (uint i = 0; i < data.size(); ++i) {
            double k = 0;
            for (uint j = 0; j < data[i].size(); ++j) {
                k += abs(data[i][j]); 
            }
            if (k > n) {
                n = k;
            }
        }

        return n;
    }

    double SqrtDiagonal() {
        double res = 0;
        for (uint i = 0; i < data.size(); ++i) {
            for (uint j = i + 1; j < data[i].size(); ++j) {
                res += data[i][j] * data[i][j];
            }
        }
        res = sqrt(res);
        return res;
    }

    matrix<T> Transpose() {
        matrix<T> res(data[0].size(), data.size()) ;
        for (uint i = 0; i < data.size(); ++i) {
            for (uint j = 0; j < data[i].size(); ++j) {
               res[i][j] = data[j][i];
            }
        }
        return res;
    }

    matrix<T>& operator=(const matrix<T> &l) {
        for (uint i = 0; i < data.size(); ++i) {
            data[i] = l.data[i];
        }
        return *this;
    }

    friend matrix<T> operator*(const matrix<T> &l, const matrix<T> &r) {
        
        matrix<T> res(l.Rows_count(), r.Columns_count());

        for (uint i = 0; i < l.Rows_count(); ++i) {
            for (uint j = 0; j < r.Columns_count(); ++j) {
                for (uint t = 0; t < l.Columns_count(); ++t) {
                    res.data[i][j] += l.data[i][t] * r.data[t][j];
                }
            }
        }
        return res;
    }

    friend vector<T> operator*(const matrix<T> &l, const vector<T> &r) {
        vector<T> res(l.Rows_count());

        for (uint i = 0; i < l.Rows_count(); ++i) {
                for (uint t = 0; t < l.Columns_count(); ++t) {
                    res[i] += l.data[i][t] * r[t];
            }
        }
        return res;
    }
};

#endif
