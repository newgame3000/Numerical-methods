#include <iostream>
#include <vector>
#include <math.h>
#include <unistd.h>

using namespace std;

// x = [0, l]
// t = [0, t]
// h - шаг по  x
// tau - шаг по t
// sigma = alpha * tau / h ^ 2 (<= 1/2)
// n - количество точек по x
// k - количество точек по  t

const double PI = 3.141592653589793;

double f(double x, double t) {
    return 0;
}

double alpha = 1;
double betta = 1;

double gamma0(double a, double b, double c, double t) {
    return exp((c - a) * t) * (cos(b * t) + sin(b * t));
}

double gamma1(double a, double b, double c,  double t) {
    return -exp((c - a) * t)* (cos(b * t) + sin(b * t));
}

vector<double> threeDiag(int n, vector<double> m, vector<double> d) {
    vector<double> x(n);
    vector<double> q(n);
    vector<double> p(n);

    p[0] = -m[1] / m[0];
    q[0] = d[0] / m[0];
    
    int k = 1;
    for (uint i = 2; i < m.size() - 2; i += 3) { //m[i] - a, m[i + 1] - b, m[i + 2] - c
        p[k] = - m[i + 2] / (m[i + 1] + m[i] * p[k - 1]);
        q[k] = (d[k] - m[i] * q[k - 1])/(m[i + 1] + m[i] * p[k - 1]);
        k += 1;
    }

    p[k] = 0;
    q[k] = (d[k] - m[m.size() - 2] * q[k - 1]) / (m[m.size() - 1] + m[m.size() - 2] * p[k - 1]);

    x[n - 1] = q[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    return x;
}

vector<vector<double>> threeDiagSolve(vector<vector<double>> u, int n, int k, double a, double c, double b, double tau, double h, int type = 0) {
    

    for (int j = 1; j < k; ++j) {

        vector<double> m((n - 2) * 3 + 4);
        vector<double> d(n);

        if (type == 0) {
            m[0] = -alpha / h + betta;
            m[1] = alpha / h;

            d[0] = gamma0(a, b, c, tau * (j + 1));
            d[n - 1] = gamma1(a, b, c, tau * (j + 1));

            m[m.size() - 2] = -alpha / h;
            m[m.size() - 1] = alpha / h + betta;
        }

        int mi = 2;


        for (int i = 1; i < n - 1; ++i) {
            m[mi] = (a * tau / (h * h)  - b * tau / (2 * h));
            m[mi + 1] = -1 - 2 * a * tau / (h * h) + c * tau;
            m[mi + 2] = (a * tau / (h * h)  + b * tau / (2 * h));
            mi += 3;
            d[i] = -u[i][j - 1] - tau * f(h * i, j * tau);
            
        }

        if (type == 1) {
            double coeff = (- alpha / (2 * h)) / m[4];
            m[0] = (-3 * alpha / (2 * h)) + betta - coeff * m[2];
            m[1] = alpha * 2 / h - coeff * m[3];
            d[0] = gamma0(a, b, c, (j + 1) * tau) - coeff * d[1];

            coeff = (alpha / (2 * h))  / m[m.size() - 5];
            m[m.size() - 2] = -alpha * 2 / h - coeff * m[m.size() - 4];
            m[m.size() - 1] = (3 * alpha / (2 * h)) + betta - m[m.size() - 3] * coeff;
            d[n - 1] = gamma1(a, b, c, (j + 1) * tau) -  d[n - 2] * coeff;

        }

        if (type == 2) {

            double coeff = h - b * h * h / (2 * a);
            m[0] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / coeff * alpha + betta;
            m[1] = alpha / coeff;

            d[0] = gamma0(a, b, c, (j + 1) * tau) - (u[0][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / coeff;


            coeff = -h - b * h * h / (2 * a);
            m[m.size() - 2] = alpha / coeff;
            m[m.size() - 1] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / coeff * alpha + betta;

            d[n - 1] = gamma1(a, b, c, (j + 1) * tau) - (u[n - 1][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / coeff;
        }

        d = threeDiag(n, m, d);

        for (int i = 0; i < n; ++i) {
            u[i][j] = d[i];
        }
    }

    return u;

}

vector<vector<double>> finiteDifferenceMethod(vector<vector<double>> u, int n, int k, double a, double c, double b, double tau, double h, int type) {

    for (int j = 0; j < k - 1; ++j) {
        for (int i = 1; i < n - 1; ++i) {
            u[i][j + 1] = tau * (a * (u[i - 1][j] - 2 * u[i][j] + u[i + 1][j]) / (h * h) + 
                          b * (u[i + 1][j] - u[i - 1][j]) / (2 * h) + 
                          c * u[i][j] + f(h * i, j * tau)) + u[i][j];

            if (type == 0) {
                u[0][j + 1] = (gamma0(a, b, c, (j + 1) * tau) - alpha / h * u[1][j + 1]) / (-alpha / h + betta);
                u[n - 1][j + 1] = (gamma1(a, b, c, (j + 1) * tau) + alpha / h * u[n - 2][j + 1]) / (alpha / h + betta);
            }

            if (type == 1) {
                u[0][j + 1] = (gamma0(a, b, c, (j + 1) * tau) - alpha * 2 / h * u[1][j + 1] + alpha / (2 * h) * u[2][j + 1]) / (-3 * alpha / ( 2 * h) + betta);
                u[n - 1][j + 1] = (gamma1(a, b, c, (j + 1) * tau) - alpha / (h * 2) * u[n - 3][j + 1] + 2 * alpha / h * u[n - 2][j + 1]) / (3 * alpha / (h * 2) + betta);
            }

            if (type == 2) {
                double coeff = h - b * h * h / (2 * a);
                u[0][j + 1] = (gamma0(a, b, c, (j + 1) * tau) - alpha / coeff * u[1][j + 1] - (alpha * h * h / (2 * tau * a) * u[0][j] + alpha *  f(0, j * tau) * h * h / (2 * a)) / coeff) / 
                              (-(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / coeff * alpha + betta);

                coeff = -h - b * h * h / (2 * a);
                u[n - 1][j + 1] = (gamma1(a, b, c, (j + 1) * tau) - alpha / coeff * u[n - 2][j + 1] - (alpha * h * h / (2 * tau * a) * u[n - 1][j] + alpha *  f(0, j * tau) * h * h / (2 * a)) / coeff) / 
                              (-(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / coeff * alpha + betta);
            }
        }
    }

    return u;

}

vector<vector<double>> combinedMethod(vector<vector<double>> u, int n, int k, double a, double c, double b, double tau, double h, double theta, int type) {
    for (int j = 1; j < k; ++j) {

        vector<double> m((n - 2) * 3 + 4);
        vector<double> d(n);

        if (type == 0) {
            m[0] = -alpha / h + betta;
            m[1] = alpha / h;

            d[0] = gamma0(a, b, c, tau * (j + 1));
            d[n - 1] = gamma1(a, b, c, tau * (j + 1));

            m[m.size() - 2] = -alpha / h;
            m[m.size() - 1] = alpha / h + betta;
        }

        int mi = 2;


        for (int i = 1; i < n - 1; ++i) {
            m[mi] = theta * (a * tau / (h * h)  - b * tau / (2 * h));
            m[mi + 1] = -1 + theta * (- 2 * a * tau / (h * h) + c * tau);
            m[mi + 2] = theta * (a * tau / (h * h)  + b * tau / (2 * h));
            mi += 3;


            d[i] = -u[i][j - 1]  - tau * (theta * f(h * i, j * tau) + (1 - theta) * ((u[i - 1][j - 1] - 2 * u[i][j - 1] + u[i + 1][j - 1]) / (h * h) * a  + 
                (u[i + 1][j - 1] - u[i - 1][j - 1]) / (2 * h) * b + 
                c * u[i][j - 1]));
        }

          if (type == 1) {
            double coeff = (- alpha / (2 * h)) / m[4];
            m[0] = (-3 * alpha / (2 * h)) + betta - coeff * m[2];
            m[1] = alpha * 2 / h - coeff * m[3];
            d[0] = gamma0(a, b, c, (j + 1) * tau) - coeff * d[1];

            coeff = (alpha / (2 * h))  / m[m.size() - 5];
            m[m.size() - 2] = -alpha * 2 / h - coeff * m[m.size() - 4];
            m[m.size() - 1] = (3 * alpha / (2 * h)) + betta - m[m.size() - 3] * coeff;
            d[n - 1] = gamma1(a, b, c, (j + 1) * tau) -  d[n - 2] * coeff;
        }


        if (type == 2) {

            double coeff = h - b * h * h / (2 * a);
            m[0] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / coeff * alpha + betta;
            m[1] = alpha / coeff;

            d[0] = gamma0(a, b, c, (j + 1) * tau) - (u[0][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / coeff;


            coeff = -h - b * h * h / (2 * a);
            m[m.size() - 2] = alpha / coeff;
            m[m.size() - 1] = -(1 + h * h / (2 * a * tau) - c * h * h / (2 * a)) / coeff * alpha + betta;

            d[n - 1] = gamma1(a, b, c, (j + 1) * tau) - (u[n - 1][j - 1] * alpha * h * h / (2 * a * tau) + alpha * f(0, j * tau) * h * h / (2 * a)) / coeff;
        }

        d = threeDiag(n, m, d);

        for (int i = 0; i < n; ++i) {
            u[i][j] = d[i];
        }

    }

    return u;

}

void begin(vector<vector<double>> &u, int n, double h){
    for (int i = 0; i < n; ++i) {
        u[i][0] = sin(i * h);

    }
}

int main() {

    double a, b, c, sigma, l, T, tau, h, theta;

    int n, k;
    cin >> l >> n >> T >> sigma >> a >> b >> c >>theta;

    h = l / (n - 1);
    tau = sigma * h * h / a;
    k = T / tau + 1;


    vector<vector<double>> u(n, vector<double>(k));
    begin(u, n, h);
    vector<vector<double>> sol(n, vector<double>(k));


    FILE * mainData;

    if ((mainData = fopen("data.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    fprintf(mainData,"%s", (to_string(a) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(b) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(c) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(tau) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(h) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(k) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(n) + "\n").c_str());



    //2Т1П

    sol = finiteDifferenceMethod(u, n, k, a, c, b, tau, h, 0);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }


    sol = threeDiagSolve(u, n, k, a, c, b, tau, h, 0);



    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }

    sol = combinedMethod(u, n, k, a, c, b, tau, h, theta, 0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }




    //3Т2П

    sol = finiteDifferenceMethod(u, n, k, a, c, b, tau, h, 1);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }


    sol = threeDiagSolve(u, n, k, a, c, b, tau, h, 1);



    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }

    sol = combinedMethod(u, n, k, a, c, b, tau, h, theta, 1);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }



    //2Т2П

    sol = finiteDifferenceMethod(u, n, k, a, c, b, tau, h, 2);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }


    sol = threeDiagSolve(u, n, k, a, c, b, tau, h, 2);



    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }

    sol = combinedMethod(u, n, k, a, c, b, tau, h, theta, 2);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            fprintf(mainData, "%s", (to_string(sol[i][j]) + "\n").c_str());
        }
    }

    fclose(mainData);
}