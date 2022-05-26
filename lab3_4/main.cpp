#include <iostream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <functional>


using namespace std;

double FirstDerivative(vector<double> x, vector<double> y, double tx) {
    if (x.size() != y.size()) {
        return 0;
    }

    int j = 0;

    for (uint i = 0; i < x.size(); ++i) {

        if (tx >= x[i] && tx <= x[i + 1]) {
            j = i - 1;
        }
    }

    return(y[j + 1] - y[j]) / (x[j + 1] - x[j]) + ((y[j + 2] - y[j + 1]) / (x[j + 2] - x[j + 1]) - (y[j + 1] - y[j]) / (x[j + 1] - x[j])) / (x[j + 2] - x[j]) * (2 * tx - x[j] - x[j + 1]);
}

double SecondDerivative(vector<double> x, vector<double> y, double tx) {
    if (x.size() != y.size()) {
        return 0;
    }

    int j = 0;

    for (uint i = 0; i < x.size(); ++i) {

        if (tx >= x[i] && tx <= x[i + 1]) {
            j = i - 1;
        }
    }

    return 2 * ((y[j + 2] - y[j + 1]) / (x[j + 2] - x[j + 1]) - (y[j + 1] - y[j]) / (x[j + 1] - x[j])) / (x[j + 2] - x[j]);
}

int main() {
    int n;
    cin >> n;

    vector<double> x(n);
    vector<double> y(n);

    for (int i = 0; i < n; ++i) {
        cin >> x[i];
        cin >> y[i];
    }

    double tx;
    cin >> tx;

    cout << "Первая производная: " << FirstDerivative(x, y, tx)  << endl;

    cout << "Вторая производная: " << SecondDerivative(x, y, tx)  << endl;

}
