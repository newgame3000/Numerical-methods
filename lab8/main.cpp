#include <iostream>
#include <vector>
#include <math.h>
#include <unistd.h>

using namespace std;

// x = [0, lx]
// y = [0, ly]
//  t = [0, T]
// hx - шаг по  x
// hy - шаг по y
// tau - шаг по t
// nx - количество точек по x
// ny - количество точек по y
// k - количество точек по t

const double PI = 3.141592653589793;

double alpha = 0;
double betta = 1;

double f(double x, double y, double t) {
    return - x * y * sin(t);
}

double f0yt(double y, double t) {
    return 0;
}

double flyt(double y, double t) {
    return y * cos(t);
}

double fx0t(double x, double t) {
    return 0;
}

double fxlt(double x, double t) {
    return x * cos(t);
}

double fxy0(double x, double y) {
    return x * y;
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

void begin(vector<vector<vector<double>>> &u, int nx, int ny, double hx, double hy){
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            u[i][j][0] = fxy0(i * hx, j * hy);
        }
    }
}

void alternatingDirectionMethod(vector<vector<vector<double>>> &u, double a, double b1, double b2, double c, int nx, int ny, int k, double hx, double hy, double tau) {

    for (int kk = 0; kk < k - 1; ++kk) {

        vector<vector<double>> layer(nx, vector<double>(ny));
        vector<double> m((nx - 2) * 3 + 4);
        vector<double> d(nx);

        m[0] = -alpha / hx + betta;;
        m[1] = alpha / hx;
        m[m.size() - 2] = -alpha / hx;
        m[m.size() - 1] = alpha / hx + betta;

        for (int j = 1; j < ny - 1; ++j) {

           
            d[0] = f0yt(j * hy, kk * tau + 0.5 * tau);
            d[nx - 1] = flyt(j * hy, kk * tau + 0.5 * tau);


            int mi = 2;
            for (int i = 1; i < nx - 1; ++i) {
                m[mi] = -a / (hx * hx) + b1 / (2 * hx);
                m[mi + 1] = 2.0 / tau + 2 * a / (hx * hx) - c;
                m[mi + 2] = -a / (hx * hx) - b1 / (2 * hx); 
                mi += 3;

                d[i] = 2.0 / tau * u[i][j][kk] + a * (u[i][j + 1][kk] - 2 * u[i][j][kk] + u[i][j - 1][kk]) / (hy * hy) + 
                       f(i * hx, j * hy, kk * tau + 0.5 * tau) + b2 / (2 * hy) * (u[i][j + 1][kk] - u[i][j - 1][kk]);
            }

            d = threeDiag(nx, m, d);

            for (int i = 0; i < nx; ++i) {
                layer[i][j] = d[i];
            }
        }

        for (int i = 0; i < nx; ++i) {
            layer[i][0] = fx0t(hx * i, kk * tau + 0.5 * tau);
            layer[i][ny - 1] = fxlt(hx * i, kk * tau + 0.5 * tau);
        }

        d = vector<double>(ny);
        m = vector<double> ((ny - 2) * 3 + 4);

        m[0] = -alpha / hy + betta;;
        m[1] = alpha / hy;
        m[m.size() - 2] = -alpha / hy;
        m[m.size() - 1] = alpha / hy + betta;

        for (int i = 1; i < nx - 1; ++i) {
 
            d[0] = fx0t(i * hx, (kk + 1) * tau);
            d[ny - 1] = fxlt(i * hx, (kk + 1) * tau);

            int mi = 2;
            for (int j = 1; j < ny - 1; ++j) {
                m[mi] = -a / (hy * hy) + b2 / (2 * hy);
                m[mi + 1] = 2.0 / tau + 2 * a / (hy * hy) - c;
                m[mi + 2] = -a / (hy * hy) - b2 / (2 * hy);
                mi += 3;

                d[j] = 2.0 / tau * layer[i][j] + a * (layer[i + 1][j] - 2 * layer[i][j] + layer[i - 1][j]) / (hx * hx) + 
                       f(i * hx, j * hy, (kk + 1) * tau) + b1 / (2 * hx) * (layer[i + 1][j] - layer[i - 1][j]);
            } 

            d = threeDiag(ny, m, d);

            for (int j = 0; j < ny; ++j) {
                u[i][j][kk + 1] = d[j];
            } 
        }

        for (int j = 0; j < ny; ++j) {
            u[0][j][kk + 1] = f0yt(hy * j, (kk + 1) * tau);
            u[nx - 1][j][kk + 1] = flyt(hy * j, (kk + 1) * tau);
        }
    }
}

void fractionalStepMethod(vector<vector<vector<double>>> &u, double a, double b1, double b2, double c, int nx, int ny, int k, double hx, double hy, double tau) {
     for (int kk = 0; kk < k - 1; ++kk) {

        vector<vector<double>> layer(nx, vector<double>(ny));
        vector<double> m((nx - 2) * 3 + 4);
        vector<double> d(nx);

        m[0] = -alpha / hx + betta;;
        m[1] = alpha / hx;
        m[m.size() - 2] = -alpha / hx;
        m[m.size() - 1] = alpha / hx + betta;

        for (int j = 1; j < ny - 1; ++j) {

            d[0] = f0yt(j * hy, kk * tau + 0.5 * tau);
            d[nx - 1] = flyt(j * hy, kk * tau + 0.5 * tau);

            int mi = 2;
            for (int i = 1; i < nx - 1; ++i) {
                m[mi] = -a / (hx * hx) + b1 / (2 * hx);
                m[mi + 1] = 1.0 / tau + 2.0 / (hx * hx) - c;
                m[mi + 2] = -a / (hx * hx) - b1 / (2 * hx); 
                mi += 3;

                d[i] = u[i][j][kk] / tau + f(i * hx, j * hy, kk * tau + 0.5 * tau) / 2;
            }   

            d = threeDiag(nx, m, d);

            for (int i = 0; i < nx; ++i) {
                layer[i][j] = d[i];
            }
        }

        for (int i = 0; i < nx; ++i) {
            layer[i][0] = fx0t(hx * i, kk * tau + 0.5 * tau);
            layer[i][ny - 1] = fxlt(hx * i, kk * tau + 0.5 * tau);
        }

        d = vector<double>(ny);
        m = vector<double> ((ny - 2) * 3 + 4);

        m[0] = -alpha / hy + betta;;
        m[1] = alpha / hy;
        m[m.size() - 2] = -alpha / hy;
        m[m.size() - 1] = alpha / hy + betta;

        for (int i = 1; i < nx - 1; ++i) {
 
            d[0] = fx0t(i * hx, (kk + 1) * tau);
            d[ny - 1] = fxlt(i * hx, (kk + 1) * tau);

            int mi = 2;
            for (int j = 1; j < ny - 1; ++j) {
                m[mi] = -a / (hy * hy) + b2 / (2 * hy);
                m[mi + 1] = 1.0 / tau + 2 / (hy * hy) - c;
                m[mi + 2] = -a / (hy * hy) - b2 / (2 * hy);
                mi += 3;

                d[j] = layer[i][j] / tau + f(i * hx, j * hy, kk * tau + 0.5 * tau) / 2;
            } 

            d = threeDiag(ny, m, d);

            for (int j = 0; j < ny; ++j) {
                u[i][j][kk + 1] = d[j];
            } 
        }

        for (int j = 0; j < ny; ++j) {
            u[0][j][kk + 1] = f0yt(hy * j, kk * tau + 0.5 * tau);
            u[nx - 1][j][kk + 1] = flyt(hy * j, kk * tau + 0.5 * tau);
        }
    }
}


int main() {
    double a, b1, b2, c, lx, ly, T, hx, hy, tau ;

    int nx, ny, k;
    cin >> lx >> ly >> T >> nx >> ny >> k >> a >> b1 >> b2 >> c;

    hx = lx / (nx - 1);
    hy = ly / (ny - 1);
    tau = T / (k - 1);

    vector<vector<vector<double>>> u(nx, vector<vector<double>>(ny, vector<double> (k)));
    begin(u, nx, ny, hx, hy);

    vector<vector<vector<double>>> sol(u);

    FILE * mainData;

    if ((mainData = fopen("data.txt", "w")) == NULL) {
        cout << "Ошибка открытия файла\n";
        return -1;
    }

    fprintf(mainData,"%s", (to_string(a) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(b1) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(b2) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(c) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(hx) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(hy) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(tau) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(nx) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(ny) + "\n").c_str());
    fprintf(mainData, "%s", (to_string(k) + "\n").c_str());


    alternatingDirectionMethod(sol, a, b1, b2, c, nx, ny, k, hx, hy, tau);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int kk = 0; kk < k; ++kk) {
                fprintf(mainData, "%s", (to_string(sol[i][j][kk]) + " ").c_str());
            }
        }
    }

     fprintf(mainData, "\n");

    fractionalStepMethod(u, a, b1, b2, c, nx, ny, k, hx, hy, tau);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int kk = 0; kk < k; ++kk) {
                fprintf(mainData, "%s", (to_string(u[i][j][kk]) + " ").c_str());
            }
        }
    }

    fclose(mainData);
}