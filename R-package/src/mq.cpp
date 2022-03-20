#include "bf.h"
#include <vector>
#include <iostream>

#define R_BUILD
#ifdef R_BUILD
#include <Rcpp.h>
using namespace Rcpp;
#endif

std::vector<double> mq(std::vector<std::vector<double> > &xx, int n1, int n1_total)
{
    BF bf_x = BF(xx, n1, n1_total);
    // std::cout << "initial pass\n";
    std::vector<std::vector<double> > pxx(n1, std::vector<double>(n1_total));
    bf_x.get_fitted(pxx);
    bf_x.free_BF();

    double mq_value;
    std::vector<double> mq_vec(n1_total);
    for (int j = 0; j < n1_total; j++)
    {
        mq_value = 0.0;
        for (int i = 0; i < n1; i++)
        {
            mq_value += pxx[i][j];
        }
        mq_value /= (double (n1));
        mq_vec[j] = mq_value;
    }
    
    return mq_vec;
}

// [[Rcpp::export]]
void mq_cpp(Rcpp::NumericMatrix &x, int n1, int n1_total, Rcpp::NumericVector &mq_res)
{
    n1 = n1_total;  // enforces a square matrix
    std::vector<std::vector<double> > xx(n1, std::vector<double>(n1_total));

    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n1_total; j++)
        {
            xx[i][j] = x(i, j);
        }
    }

    std::vector<double> mq_vec = mq(xx, n1, n1_total);
    for (int j = 0; j < n1_total; j++)
    {
        mq_res(j) = mq_vec[j];
    }
    return;
}
