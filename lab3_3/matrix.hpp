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

    int Sign(T a){
        if (a > 0) {
            return 1;
        } else {
            if (a < 0) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    double Norm(vector<T> b, int k) {
        double ans = 0;

        for (uint i = k; i < b.size(); ++i) {
            ans += b[i] * b[i];
        }

        return (sqrt(ans));
    }

    matrix<T> VVT(vector<T> v, vector<T> vt) {
        matrix<T> ans(v.size(), v.size());

        for (uint i = 0; i < v.size(); ++i) {
            for (uint j = 0; j < vt.size(); ++j) {
                ans[i][j] = v[i] * vt[j];
            }
        }

        return ans;
    }

    T VTV(vector<T> v, vector<T> vt) {
        T ans = 0;
        for (uint i = 0; i < v.size(); ++i) {
                ans += v[i] * vt[i];
        }

        return ans;
    }

    void Maximum(matrix<T> &u, matrix<T> &l,vector<pair<int, int>> &sw, int &swaps, int n, int j) {
        double max = u[j][j];
        int index = j;
        for (int i = j; i < n; ++i) {
            if (u[i][j] > u[j][j] && u[i][j] > max) {
                index = i;
                max = u[i][j];
            }
        }

        if (index != j) {
            sw.push_back(make_pair(j, index));
            swaps += 1;
            l.Swap_rows(j, index);
            l.Swap_columns(j, index);
            u.Swap_rows(j, index);
        }
    }

public:

    matrix(int n, int m) {
        data = vector<vector<T>>(n, vector<T>(m));
    }

    matrix(const matrix<T> &a) {
        data = a.data;
    }

    matrix(initializer_list<vector<T>> list): data(list) {}

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
        //cout << data.size() << " " << data[0].size() << endl;
        matrix<T> res(data[0].size(), data.size()) ;
        for (uint i = 0; i < data.size(); ++i) {
            for (uint j = 0; j < data[i].size(); ++j) {
                res[j][i] = data[i][j];
            }
        }
        return res;
    }

    void Identity() {
        for (uint i = 0; i < data.size(); ++i) {
            data[i][i] = 1;
        }
    }

    void QR(matrix<T> &q, matrix<T> &r) {
        r.data = data;
        q = matrix<T>(r.Rows_count(), r.Columns_count());
        q.Identity();

        matrix<T> E(r.Rows_count(), r.Columns_count());
        matrix<T> H(r.Rows_count(), r.Columns_count());

        E.Identity();

        for (uint k = 0; k < r.Columns_count() - 1; ++k) {
            vector<T> b(r.Rows_count());
            for (uint i = 0; i < r.Rows_count(); ++i) {
                b[i] = r[i][k];
            }
            vector<T> v(r.Rows_count());

            for (uint i = 0; i < r.Rows_count(); ++i) {
                
                v[i] = b[i];

                if (i < k) {
                    v[i] = 0;
                }

                if (i == k) {
                    v[i] += Sign(b[i]) * Norm(b, k);
                }
            }


            H = E  - 2 * (VVT(v, v) / VTV(v, v));

            q = q * H;
            r = H * r;
        }
    }

    void LU(matrix<T> &l, matrix<T> &m, vector<pair<int, int>> &sw, int &swaps, int n) {
        l[0][0] = 1;
        for (int j = 0; j < n; ++j) {
            l[j][j] = 1;
        }
        int k = 1;
        for (int j = 0; j < n; ++j) {
            Maximum(m, l, sw, swaps, n, j);
            for (int i = k; i < n; ++i) {
                l[i][j] = m[i][j] / m[j][j];
                for (int t = 0; t < n; ++t) {
                    m[i][t] -= m[j][t] * l[i][j];
                }
            }
            k += 1;
        }
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

    friend matrix<T> operator/(const matrix<T> &l, const T &r) {
        matrix<T> res(l.Rows_count(), l.Columns_count());
        for (uint i = 0; i < l.Rows_count(); ++i) {
            for (uint j = 0; j < l.Columns_count(); ++j) {
                res.data[i][j] = l.data[i][j] / r;
            }
        }

        return res;
    }

    friend matrix<T> operator*(const T &r, const matrix<T> &l) {
        matrix<T> res(l.Rows_count(), l.Columns_count());
        for (uint i = 0; i < l.Rows_count(); ++i) {
            for (uint j = 0; j < l.Columns_count(); ++j) {
                res.data[i][j] = l.data[i][j] * r;
            }
        }
        return res;
    }

    friend matrix<T> operator-(const matrix<T> &l, const matrix<T> &r) {
        
        matrix<T> res(l.Rows_count(), l.Columns_count());

        for (uint i = 0; i < l.Rows_count(); ++i) {
            for (uint j = 0; j < r.Columns_count(); ++j) {
                res.data[i][j] = l.data[i][j] - r.data[i][j];
            }
        }
        return res;
    }

};

#endif
