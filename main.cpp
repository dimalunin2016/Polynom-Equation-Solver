#include <iostream>
#include <complex>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>
#include "PolynomEquationSolver.h"

/*
#ifdef Test
#define TEST true
#else
#define TEST false
#endif
*/

const double EPS = 1e-8;

bool IsZero(double x) {
  
  return fabs(x) < EPS;
}


bool IsZero(std::complex<double> x) {
  
  return fabs(x.real()) < EPS & fabs(x.imag()) < EPS;
}


template<typename Field>
bool IsAnswer( 
    const std::vector<Field>& params, 
    const std::vector<Field>& roots) {

  for (const Field& x: roots) {
    const Field ZERO = x - x;
    const Field ONE = x / x;
    Field ans = ZERO;
    Field x_in_degree = ONE; 
    for (int i = params.size() - 1; i >= 0; --i) {
      ans += x_in_degree * params[i];
      x_in_degree *= x;
    }
    if (!IsZero(ans)) {
      return false;
    }
  }
  return true;
}


void RealCheckAns(
    
    const RealEquationSolutions& ans, 
    const std::vector<double>& param,
    int exp_num_of_roots) {

    assert(ans.GetNumOfSolutions() == exp_num_of_roots);
    assert(IsAnswer<double>(param, ans.GetSolutions()));
}


void RealSquareEquationWithoutErrorsTest() {
  
  const int N = 3;
  std::vector<double> param[N];
  param[0] = {0, 0, 1, -1, -1}; 
  param[1] = {4, 8, 4};
  param[2] = {0, 3, 3, 3};
  RealEquationSolutions ans_vec[N];
  
  RealEquationSolutions ans_var[N];
  ans_var[0] = PolynomialSolver(0, 0, 0, 0, 0, 1, -1, -1); 
  ans_var[1] = PolynomialSolver(0, 0, 4, 8, 4);
  ans_var[2] = PolynomialSolver(3, 3, 3);
  
  int num_of_roots = 2;
  for (int i = 0; i < N; ++i) {
    ans_vec[i] = PolynomialSolver(param[i]);
    RealCheckAns(ans_vec[i], param[i], num_of_roots);
    RealCheckAns(ans_var[i], param[i], num_of_roots);
    num_of_roots -= 1;
  }
}


void RealLinearAndConstEquationWithoutErrorsTest() {
  
  const int N = 3;
  std::vector<double> param[N];
  param[0] = {4};
  param[1] = {3, 7};
  param[2] = {0};
  
  RealEquationSolutions ans_var[N];
  ans_var[0] = PolynomialSolver(0, 0, 4);
  ans_var[1] = PolynomialSolver(3, 7);
  ans_var[2] = PolynomialSolver(0, 0, 0);
  
  RealCheckAns(ans_var[0], param[0], 0);
  RealCheckAns(ans_var[1], param[1], 1);
  RealCheckAns(ans_var[2], param[2], PS_INF_ROOTS); 
}


void RealEquationErrorTest() {

  bool is_mistake_found = false;
  try {
    PolynomialSolver(1, -1, 3, 4);
  } catch (...) {
    is_mistake_found = true;
  }
  assert(is_mistake_found);
  
  is_mistake_found = false;
  try {
    PolynomialSolver(1, 3, 4.0/0);
  } catch (...) {
    is_mistake_found = true;
  }
  assert(is_mistake_found);
}


void ComplexCheckAns(
    
    const ComplexEquationSolutions& ans, 
    const std::vector<std::complex<double> >& param,
    int exp_num_of_roots) {

    assert(ans.GetNumOfSolutions() == exp_num_of_roots);
    assert(IsAnswer<std::complex<double> >(param, ans.GetSolutions()));
}


void ComplexSquareEquationWithoutErrorsTest() {
  
  const int N = 2;
  std::vector<std::complex<double> > param[N];
  std::complex<double> 
    c1(4, 0), c2(8, 0), c3(4, 0), 
    c4(7.78, 8.75), c5(6773.0, 377.91), c6(123.93, 3.8);
  param[0] = {c1, c2, c3};
  param[1] = {c4, c5, c6};
  
  ComplexEquationSolutions ans_var[N];
  ans_var[0] = ComplexPolynomialSolver(c1, c2, c3); 
  ans_var[1] = ComplexPolynomialSolver(c4, c5, c6);
  
  ComplexCheckAns(ans_var[0], param[0], 1);
  ComplexCheckAns(ans_var[1], param[1], 2);
}


void ComplexLinearAndConstEquationWithoutErrorsTest() {
  
  const int N = 3;
  std::vector<std::complex<double> > param[N];
  std::complex<double> 
    c1(4.3, 32), c2(8, 34), c3(4, 43.9), c4(0, 0);
  param[0] = {c1};
  param[1] = {c2, c3};
  param[1] = {c4};
  
  ComplexEquationSolutions ans_var[N];
  ans_var[0] = ComplexPolynomialSolver(c1); 
  ans_var[1] = ComplexPolynomialSolver(c2, c3);
  ans_var[2] = ComplexPolynomialSolver(c4);
  
  ComplexCheckAns(ans_var[0], param[0], 0);
  ComplexCheckAns(ans_var[1], param[1], 1);
  ComplexCheckAns(ans_var[2], param[2], PS_INF_ROOTS);
}


void ComplexEquationErrorTest() {
  
  std::complex<double> 
    c1(4.3, 32), c2(8, 34), c3(4, 43.9), c4(0, 0);
  bool is_mistake_found = false;
  try {
    ComplexPolynomialSolver(c1, c2, c3, c4);
  } catch (...) {
    is_mistake_found = true;
  }
  assert(is_mistake_found);
}


void RunTest(bool test) {
  
  if (test) {
    ComplexSquareEquationWithoutErrorsTest();
    ComplexLinearAndConstEquationWithoutErrorsTest();
    ComplexEquationErrorTest();
    RealSquareEquationWithoutErrorsTest();
    RealLinearAndConstEquationWithoutErrorsTest();
    RealEquationErrorTest();
    std::cout<<"ALL TESTS HAVE BEEN PASSED!\n";
  }
}


int main() {
  RunTest(true);
  return 0;
}

