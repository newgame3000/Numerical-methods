#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>
#include <vector>
#include <string>
#include <iomanip>

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

    matrix Inverse() {
        
    }

    void Print() {
        int max = 10000000;
        for (uint i = 0; i < data.size(); ++i) {
            for (uint j = 0; j < data[i].size(); ++j) {
                if (abs(data[i][j] > max)) {

                }
            }
        }

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
};

#endif