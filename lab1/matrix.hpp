#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>
#include <vector>
#include <string>

using namespace std;

template <class T>
class matrix {
private: 
    vector<vector<T>> data;
public:

    matrix(int n, int m) {
        data = vector<vector<T>>(n, vector<T>(m));
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
        for (uint i = 0; i < data.size(); ++i) {
            for (uint j = 0; j < data[i].size(); ++j) {
                cout << data[i][j] << fixed << " ";
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