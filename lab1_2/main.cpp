#include <iostream>
#include <vector>

using namespace std;

int main() {
    int n;
    cin >> n;
    vector<double> x(n);
    vector<double> m((n - 2) * 3 + 4);
    vector<double> q(n);
    vector<double> p(n);
    vector<double> d(n);

    for (uint i = 0; i < m.size(); ++i) {
        cin >> m[i];
    }

    for (int i = 0; i < n; ++i) {
        cin >>d[i];
    }    

    p[0] = -m[1] / m[0];
    q[0] = d[0] / m[0];
    int k = 1;
    for (uint i = 2; i < m.size() - 2; i += 3) { //m[i] - a, m[i + 1] - b, m[i + 2] - c
        p[k] = - m[i + 2] / (m[i + 1] + m[i] * p[k - 1]);
        q[k] = (d[k] - m[i] * q[k - 1])/(m[i + 1] + m[i] * p[k - 1]);
        k += 1;
    }
    cout << m[n - 2] << m[n - 1] << endl;
    p[k] = 0;
    q[k] = (d[k] - m[m.size() - 2] * q[k - 1]) / (m[m.size() - 1] + m[m.size() - 2] * p[m.size() - 1]);

    cout << q[k] << endl;

    x[n - 1] = q[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    for (int i = 0; i < n; ++i) {
        cout << x[i] << " ";
    }

    cout << endl;
}