#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
float calculate_ziw_statistic(NumericVector x, NumericVector y) {
  // total observations in each vector
  int Nx = x.size();
  int Ny = y.size();
  NumericVector N = NumericVector::create(Nx,Ny);
  
  // number of non-zero observations in each vector
  int nx = sum(x != 0);
  int ny = sum(y != 0);
  NumericVector n = NumericVector::create(nx,ny);
  
  // non-zero proportion
  NumericVector prop = n / N ;
  float propmax  = max(prop) ;
  float propmean = mean(prop);
  
  // keep only round(pmax * N) observations in each group
  int Ntrun = round(propmax * N);
  // print(prop);
  //print(propmax);
  //print(propmean);
  return 1.0 * propmax;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
x <- c(0,1,0,2,3)
y <- c(1,3,0,4)
calculate_ziw_statistic(x,y)
*/
