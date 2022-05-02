#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP 

#include <iostream>
#include <vector>


using namespace std;

template <class T>
vector<T> operator* (const vector<T>& l, const T r) {
    vector<T> res(l.size());
    for (uint i = 0; i < l.size(); ++i) {
        res[i] = l[i] * r;
    }

    return res;
}

template <class T>
vector<T> operator+ (const vector<T>& l, const vector<T>& r) { 

    uint minlen = min(l.size(), r.size());
    vector<T> res(max(l.size(), r.size()));

    for (uint i = 0; i < minlen; ++i) {
        res[i] = l[i] + r[i];
    }

    if (l.size() < r.size()) {
        for (uint i = minlen; i < r.size(); ++i) {
            res[i] = r[i];
        }
    } 

    if (l.size() > r.size()) {
        for (uint i = minlen; i < l.size(); ++i) {
            res[i] = l[i];
        }
    } 

    while(res.size() > 0 && res[res.size() - 1] == 0) {
        res.pop_back();
    }

    return res;
}

template <class T>
vector<T> operator* (const vector<T>& l, const vector<T>& r) { 

    vector<T> res(l.size() + r.size());

    for (uint i = 0; i < l.size(); ++i) {
        for (uint j = 0; j < r.size(); ++j) {
            res[i + j] += l[i] * r[j];
        }
    }

    while(res.size() > 0 && res[res.size() - 1] == 0) {
        res.pop_back();
    }

    return res;
}

template <class T>
class polynomial {
private:

    //a0 + a1x + a2x^2
    //[a0, a1, a2]
    vector<T> data;
    vector <vector<double>> fuctions;
    vector <vector<bool>> calculate;

    double f(int l, int r, vector<T> x, vector<T> y) {
        if (calculate[l][r]) {
            return fuctions[l][r];
        }
        calculate[l][r] = true;
        double res;
        if (l + 1 == r) {
            res = (y[l] - y[r]) / (x[l] - x[r]);
        } else {
            res = (f(l, r - 1, x, y) - f(l + 1, r, x ,y)) / (x[l] - x[r]);
        }
        return fuctions[l][r] = res;
    }

public:

    polynomial(int n) {
        //cout << n << endl;
        data = vector<T>(n);
        calculate = vector<vector<bool>>(n, vector<bool>(n));
        fuctions = vector<vector<double>>(n, vector<double>(n));
    }

    polynomial(const polynomial<T> &a) {
        data = a.data;
    }

    polynomial(const vector<T> &a) {
        data = a;
    }

    polynomial() {
        data = vector<T>();
    }

    T& operator[](int i) {
        return data[i];
    }

    void operator= (const polynomial<T> &r) {
        data = r.data;
    }

    friend polynomial<T> operator*(const polynomial<T>& l, const polynomial<T>& r) {
        polynomial<T> res(l.data.size());
        res.data = l.data * r.data;
        return res;
    }


    friend polynomial<T> operator+(const polynomial<T>& l, const polynomial<T>& r) {
        polynomial<T> res(l.data.size());
        res.data = l.data + r.data;
        return res;
    }

    friend polynomial<T> operator+(const polynomial<T>& l, T r) {
        polynomial<T> res(l);
        res.data[0] += r;
        return res;
    }

    friend polynomial<T> operator*(const polynomial<T>& l, const T& r) {
        polynomial<T> res(l.data.size());
        res.data = l.data * r;
        return res;
    }

    T operator()(T x) {
        T res = 0;
        for (uint i = 0; i < data.size(); ++i) {
            res += pow(x, i) * data[i];
        }
        return res;
    }

    void LagrangePolynomial(vector<T> x, vector<T> y) {
        if (x.size() != y.size()) {
            throw runtime_error("Different length of arrays");
        }

        data = vector<T>(x.size() - 1);

        for (uint i = 0; i < x.size(); ++i) {

            vector<T> l_i(1);
            l_i[0] = 1;

            for (uint j = 0; j < x.size(); ++j) {
                if (j == i) {
                    continue;
                }

                vector<T> step(2);
                step[0] = - x[j] / (x[i] - x[j]);
                step[1] = 1.0 / (x[i] - x[j]);

                l_i = l_i * step;
            }
            data = data + l_i * y[i];
        }

        while(data.size() > 0 && data[data.size() - 1] == 0) {
            data.pop_back();
        }
    }

    void NewtonPolynomial(vector<T> x, vector<T> y) {
        if (x.size() != y.size()) {
            throw runtime_error("Different length of arrays");
        }

        data = vector<T>(x.size() - 1);

        data[0] = y[0];

        vector<T> step(1);
        step[0] = 1;
        for(uint i = 0; i < x.size() - 1; ++i) {

            vector<T> x_i(2);
            x_i[0] = -x[i];
            x_i[1] = 1;
            step = step * x_i;

            data = data + (step * f(0, i + 1, x, y));
        }

        while(data.size() > 0 && data[data.size() - 1] == 0) {
            data.pop_back();
        }
    }


    void Debug_Print() {
        for (uint i = 0; i < data.size(); ++i) {
            cout << data[i] << " ";
        }

        cout << endl;
    }

    void Print() {

        if (data.size() == 0) {
            cout << "0\n";
            return;
        }

        for (int i = (int)data.size() - 1; i >= 0 ; --i) {
            if (i == 0) {
                if (data[0] < 0) {
                    cout << " - " << -data[0] << endl;
                } else {
                    if (data[0] > 0) {
                        cout << " + " << data[0] << endl;
                    }
                }
                continue;
            } 

            if (i == (int)data.size() - 1 && data[i] >= 0) {
                cout << "x^" << i << " * " << data[i];
                continue;
            }

            if (data[i] >= 0) {
                cout << " + x^" << i << " * " << data[i];
            } else {
                cout << " - x^" << i << " * " << -data[i];
            }
        }
    }


};


#endif