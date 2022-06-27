#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

constexpr double l = 10;
constexpr double T0 = 300;
constexpr double R = 0.5;

constexpr double k0 = 1;

constexpr double a1 = 0.0134;
constexpr double b1 = 1;
constexpr double c1 = 4.35e-4;
constexpr int m1 = 1;

constexpr double a2 = 2.049;
constexpr double b2 = 0.563e-3;
constexpr double c2 = 0.528e5;
constexpr int m2 = 1;

constexpr double alpha0 = 0.05;
constexpr double alphaN = 0.01;

double lambda(double T)
{
    return a1 * (b1 + c1 * pow(T, m1));
}

double c(double T)
{
    return a2 + b2 * pow(T, m2) - c2 / (T * T);
}

double alpha(double x)
{
    return 0.05 * l / (4 * x + l);
}

constexpr double Fmax = 50;
constexpr double tmax = 5;


double F0(double t)
{
    return (t / tmax) * Fmax * exp(-(t / tmax - 1));
}

struct Range
{
    double start = 0.0;
    double stop;
    double step;

    Range(double _stop, double _step) : stop(_stop), step(_step)
    {

    }
};

struct DoubleRange
{
    Range xrange;
    Range trange;

    DoubleRange(const Range &_xrange, const Range &_trange) : xrange(_xrange), trange(_trange) {}
};

struct triplet
{
    double f, s, t;

    triplet(double x, double y, double z) : f(x), s(y), t(z) {}

    double operator [](int i) const
    {
        switch (i)
        {
            case 0:
                return f;
            case 1:
                return s;
            case 2:
                return t;
        }

    }
};

vector<double> solve_matrix(const std::vector<triplet>& mat, const vector<double>& vec)
{
    vector<double> ksi, eta;

    int n = mat.size() - 1;

    ksi.push_back(0.0);
    eta.push_back(0.0);

    ksi.push_back(-mat[0][1] / mat[0][0]);
    eta.push_back(vec[0] / mat[0][0]);


    for (int i = 1; i < mat.size() - 1; i++)
    {
        double a = mat[i][0];
        double b = -mat[i][1];
        double c = mat[i][2];

        ksi.push_back(c / (b - a * ksi[i]));
        eta.push_back((-vec[i] + a * eta[i]) / (b - a * ksi[i]));
    }


    vector<double> res;
    res.push_back((vec[n] - mat[n][1] * eta[n]) / (mat[n][1] * ksi[n] + mat[n][2]));
    for (int i = 0; i < n; i++)
    {
        res.push_back(res[i] * ksi[n - i] + eta[n - i]);
    }

    for (int i = 0; i < res.size() / 2; i++)
    {
        std::swap(res[i], res[res.size() - 1 - i]);
    }

    return res;
}

constexpr double EPS = 1e-2;

bool check(const vector<double>& a, const vector<double>& b)
{
    for (int i = 0; i < a.size(); i++)
        if (fabs(a[i] - b[i]) > EPS)
            return false;

    return true;
}

vector<double> boundary_condition(const Range& range)
{
    vector<double> res;
    for (double x = range.start; x < range.stop + EPS; x += range.step)
    {
        res.push_back(T0);
    }

    return res;
}

double k(double T)
{
    return lambda(T);
}

double p(double x)
{
    return 2. / R * alpha(x);
}

double f(double x, double t, double T)
{
    double tmp = T / 300.;
    double tmp2 = k0 * tmp * tmp;
    return 2 * T0 / R * alpha(x) + tmp2 * F0(t) * exp(-tmp2 * x);
}

std::pair<vector<triplet>, vector<double>> dif_scheme(const vector<double>& prev, const vector<double>& prev_layer, double t, const Range& xrange, const Range& trange)
{
    vector<triplet> tres;
    vector<double> cres;

    double h = xrange.step;
    double tau = trange.step;

    tres.push_back(triplet(
            (c(prev[0]) + (c(prev[0]) + c(prev[1])) / 4.) * h / 4. - alpha0 * tau + ((p(0) + p(h)) / 4. + p(0)) * tau * h / 4. + (k(prev[0]) + k(prev[1])) / (2. * h) * tau,
            (c(prev[0]) + c(prev[1])) * h / 16. - (k(prev[0]) + k(prev[1])) / (2. * h) * tau + (p(prev[0]) + p(prev[1])) * tau * h / 16.,
            0.)
    );

    cres.push_back(((c(prev[0]) * prev_layer[0] + (c(prev[0]) + c(prev[1])) * (prev_layer[0] + prev_layer[1]) / 4.) * h / 4.
                    + alpha0 * T0 * tau +
                    (3. * f(0, t, prev[0]) + f(h, t, prev[1])) * h * tau / 8.));

    int i = 1;

    for (double x = xrange.start + xrange.step; i < prev.size() - 1; x += xrange.step, i++)
    {
        double A = (k(prev[i - 1]) + k(prev[i])) * tau / h / 2.;
        double D = (k(prev[i + 1]) + k(prev[i])) * tau / h / 2.;
        tres.push_back(triplet(
                A,
                -(A + D + c(prev[i]) * h + p(x) * tau * h),
                D
        ));

        cres.push_back(-(f(x, t, prev[i]) * tau * h + c(prev[i]) * prev_layer[i] * h));
    }


    int n = prev.size() - 1;

    tres.push_back(triplet(
            0.,
            ((c(prev[n - 1]) + c(prev[n])) / 4. + c(prev[n])) * h / 4. + (k(prev[n]) + k(prev[n - 1])) * tau / h / 2 + alphaN * tau + p(l) * tau * h / 4. + (p(l) + p(l - h)) * tau * h / 16.,
            (c(prev[n-1]) + c(prev[n])) * h / 16. - (k(prev[n]) + k(prev[n-1])) * tau / h / 2 + (p(l) + p(l-h)) * h * tau / 16.

    ));

    cres.push_back((((c(prev[n - 1]) + c(prev[n])) * (prev_layer[n-1] + prev_layer[n]) / 4. + c(prev[n]) * prev_layer[n]) * h / 4. + alphaN * T0 * tau + (3. * f(l, t, prev[n]) + f(l - h, t, prev[n-1])) * h * tau / 8.));


    return make_pair(tres, cres);
}

vector<double> simple_iter(const vector<double> &prev_layer, double t, const DoubleRange &range)
{
    vector<double> cur = prev_layer;
    vector<double> prev;
    do
    {
        prev = cur;
        std::pair<vector<triplet>, vector<double>> coefs(dif_scheme(cur, prev_layer, t, range.trange, range.xrange));
        cur = solve_matrix(coefs.first, coefs.second);
    } while (!check(cur, prev));

    return cur;
}

vector<vector<double>> create_mesh(const DoubleRange &range)
{
    vector<vector<double>> res;

    res.push_back(boundary_condition(range.xrange));

    const Range& trange = range.trange;
    const Range& xrange = range.xrange;
    int i = 0;

    for (double t = trange.start + trange.step; t < trange.stop + EPS; t += trange.step, i++)
    {
        res.push_back(simple_iter(res[i], t, range));
    }

    return res;
}


void print_mesh(ostream& stream, const vector<vector<double>>& mesh, const DoubleRange & range)
{
    const Range& trange = range.trange;
    const Range& xrange = range.xrange;

    for (double x = xrange.start; x < xrange.stop + EPS; x += xrange.step)
        stream << x << " ";
    stream << "\n";

    for (double t = trange.start; t < trange.stop + EPS; t += trange.step)
        stream << t << " ";
    stream << "\n";

    for (int i = 0; i < mesh.size(); i++)
    {
        for (int j = 0; j < mesh[i].size(); j++)
        {
            stream << mesh[i][j] << " ";
        }
        stream << "\n";
    }
}

int main()
{
    DoubleRange range(Range(10., 1e-2), Range(10, 1e-3));

    {
        ofstream str("result.txt");
        print_mesh(str, create_mesh(range), range);
    }
    
}
