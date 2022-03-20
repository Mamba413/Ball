#include "bf.h"
#include <vector>
#include <iostream>

#ifdef R_BUILD
#include <Rcpp.h>
using namespace Rcpp;
#endif

double sbd(std::vector<std::vector<double> > &xx, std::vector<std::vector<double> > &xy,
           std::vector<std::vector<double> > &yx, std::vector<std::vector<double> > &yy,
           int n1, int n1_total, int n2, int n2_total)
{
    BF bf_x = BF(xx, n1, n1_total);
    BF bf_y = BF(yy, n2, n2_total);
    // std::cout << "initial pass\n";

    std::vector<std::vector<double> > pxx(n1, std::vector<double>(n1_total));
    std::vector<std::vector<double> > pxy(n1, std::vector<double>(n1_total));
    std::vector<std::vector<double> > pyx(n2, std::vector<double>(n2_total));
    std::vector<std::vector<double> > pyy(n2, std::vector<double>(n2_total));

    bf_x.get_fitted(pxx);
    bf_y.get_fitted(pyy);
    // std::cout << "get fitted pass\n";

    bf_x.predict(pxy, xy, n2_total);
    // std::cout << "\n predict bf_x pass\n";
    bf_y.predict(pyx, yx, n1_total);
    // std::cout << "\n predict bf_y pass\n";
    // std::cout << "get predict pass\n");

    double tmp;
    double denominator;
    double A = 0.0, C = 0.0, bd;

    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n1_total; j++)
        {
            tmp = pxx[i][j] - pxy[i][j];
            A += tmp * tmp;
        }
    }
    denominator = 1.0 / ((double)n1 * n1_total);
    A *= denominator;

    for (int i = 0; i < n2; i++)
    {
        for (int j = 0; j < n2_total; j++)
        {
            tmp = pyy[i][j] - pyx[i][j];
            C += tmp * tmp;
        }
    }
    denominator = 1.0 / ((double)n2 * n2_total);
    C *= denominator;

    bd = A + C;

    bf_x.free_BF();
    bf_y.free_BF();
    
    return bd;
}

// [[Rcpp::export]]
double sbd_cpp(Rcpp::NumericMatrix &x, int n1, int n1_total, int n2, int n2_total)
{
    int num = n1_total + n2_total;
    int fix_center_num = n1 + n2;

    std::vector<std::vector<double> > xx(n1, std::vector<double>(n1_total));
    std::vector<std::vector<double> > xy(n1, std::vector<double>(n2_total));
    std::vector<std::vector<double> > yx(n2, std::vector<double>(n1_total));
    std::vector<std::vector<double> > yy(n2, std::vector<double>(n2_total));

    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n1_total; j++)
        {
            xx[i][j] = x(i, j);
        }
    }

    for (int i = 0; i < n1; i++)
    {
        for (int j = n1_total; j < num; j++)
        {
            xy[i][j - n1_total] = x(i, j);
        }
    }

    for (int i = n1; i < fix_center_num; i++)
    {
        for (int j = 0; j < n1_total; j++)
        {
            yx[i - n1][j] = x(i, j);
        }
    }

    for (int i = n1; i < fix_center_num; i++)
    {
        for (int j = n1_total; j < num; j++)
        {
            yy[i - n1][j - n1_total] = x(i, j);
        }
    }

    // std::cout << "xx matrix: " << std::endl;
    // for (int i = 0; i < n1; i++)
    // {
    //     for (int j = 0; j < n1_total; j++)
    //     {
    //         std::cout << xx[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    double sbd_value = sbd(xx, xy, yx, yy, n1, n1_total, n2, n2_total);

    return sbd_value;
}
