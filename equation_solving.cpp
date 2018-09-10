/*! \file */ 
#include <iostream>
#include <functional>
#include <complex>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <cassert>

/*
#ifdef Test
#define TEST true
#else
#define TEST false
#endif
*/

const int PS_INF_ROOTS = std::numeric_limits<int>::max();


template<typename Field>
class FieldEquationSolutions;


template<typename Field> 
FieldEquationSolutions<Field> FieldPolynomialSolverHelper(
    std::vector<Field>&& polynom_params); 

/**
 * @brief result-class of FieldPolynomialSolver
 * @tparam Field - field in which the answer lies
 */
template<typename Field>
class FieldEquationSolutions {
 private:

  friend FieldEquationSolutions 
    FieldPolynomialSolverHelper<Field>(
        std::vector<Field>&& polynom_params);
  
  int num_of_sol_ = -1;
  std::vector<Field> solutions_;

 public:
  
	/*!
   * @return number of roots;
   *  PS_INF_ROOTS if real number_of_roots is infinity
   */ 
  const int GetNumOfSolutions() const {
    
    return num_of_sol_;
  }


  /*!
   * @returnstd::vector of roots
   */
  const std::vector<Field>& GetSolutions() const {
    
    return solutions_;
  }
};


template<typename Field>
void FindingEquationError(const std::vector<Field>& polynom_params) {

  if(polynom_params.size() >= 4) {
    throw std::runtime_error("Sorry, too much degree of equation");
  }
}


template<>
void FindingEquationError<double>(const std::vector<double>& polynom_params) {

  if(polynom_params.size() >= 4) {
    throw std::runtime_error("Sorry, too much degree of equation");
  }

  for (const double&  param : polynom_params) {
    if(!std::isfinite(param)) {
      throw std::runtime_error("Invalid equation");
    }
  }
}


template<typename Field>
std::vector<Field> ParamsWithoutLeadingZeros(
    std::vector<Field>&& enter_polynom_params) {

  std::vector<Field> polynom_params = std::move(enter_polynom_params);
  if(polynom_params.size() == 0) {
    return polynom_params;
  }

  const Field ZERO = polynom_params[0] - polynom_params[0];

  typename std::vector<Field>::iterator first_not_zero = find_if(
      polynom_params.begin(), 
      polynom_params.end(),
      [& ZERO](Field elem){return elem != ZERO;});
  
  if(first_not_zero == polynom_params.end()) {
    return std::vector<Field>{ZERO};
  }
  polynom_params.erase(
        polynom_params.begin(), first_not_zero);
  return polynom_params;
}


template<typename Field>
std::vector<Field> VectorWithoutRepeatElements(
    const std::vector<Field>& vec) {

  std::vector<Field> ans;
  for (const Field& x: vec) {
    if(ans.empty() || x != ans.back() ) {
      ans.push_back(x);
    }
  }
  return ans;
}


template<typename Field>
std::vector<Field> CreatePolynomParams( 
    std::vector<Field>&& enter_polynom_params) {
   
  std::vector<Field> polynom_params = 
    ParamsWithoutLeadingZeros(
        std::forward<std::vector<Field> >(enter_polynom_params));

  return polynom_params;
}


template<typename Field>
std::vector<Field> CreateVectorOfSqrt(const Field&  x) {
  
  try {

    return {-sqrt(x), sqrt(x)};
  } catch(...) {

    return {};
  }
}


template<>
std::vector<double> CreateVectorOfSqrt(const double&  x) {

  if(x < 0) {
    return {};
  }
  return {-sqrt(x), sqrt(x)};
}


template<typename Field>
void ConstantEquationSolver(
    const std::vector<Field>& polynom_params, 
    int *num_of_sol) {
  
  const Field ZERO = polynom_params[0] - polynom_params[0];
  if(polynom_params[0] == ZERO) {
    *num_of_sol =  PS_INF_ROOTS;
  } else {
    *num_of_sol = 0;
  }
}


template<typename Field>
void LinearEquationSolver(
    const std::vector<Field>& polynom_params,
    int *num_of_sol, std::vector<Field> *solutions) {
  
  const Field ZERO = polynom_params[0] - polynom_params[0];
  *num_of_sol = 1;
  solutions->push_back(polynom_params[1] == ZERO ? 
      ZERO :
      -polynom_params[1] / polynom_params[0]);
}


template<typename Field>
void SquareEquationSolver(
    const std::vector<Field>& polynom_params, 
    int *num_of_sol, std::vector<Field> *solutions) {

  const Field&  a = polynom_params[0];
  const Field&  b = polynom_params[1];
  const Field&  c = polynom_params[2];
  const Field&  D = b * b - a * c - a * c - a * c - a * c;
  std::vector<Field>sol = VectorWithoutRepeatElements<Field>(
      CreateVectorOfSqrt(D));

  *num_of_sol = sol.size();
  for (Field& x : sol) {
    x =(-b + x) /(a + a);
  }
  *solutions = sol;
}


template<typename Field>
FieldEquationSolutions<Field> FieldPolynomialSolverHelper(    
    std::vector<Field>&& enter_polynom_params) {

  FieldEquationSolutions<Field> ans;
  std::vector<Field> polynom_params = CreatePolynomParams(
      std::forward<std::vector<Field> >(enter_polynom_params));

  FindingEquationError<Field>(polynom_params);

  int equation_degree = polynom_params.size() - 1;
  switch(equation_degree) {
    case 2:
      SquareEquationSolver<Field>(
          polynom_params, 
          & ans.num_of_sol_,& ans.solutions_);
      break;
    case 1:
      LinearEquationSolver<Field>(  
          polynom_params, 
          & ans.num_of_sol_,& ans.solutions_);
      break;
    default:
      ConstantEquationSolver<Field>(
          polynom_params, 
          & ans.num_of_sol_);
      break;

  }
  return ans;
}


template<typename Field, typename... Args>
FieldEquationSolutions<Field> FieldPolynomialSolverHelper(
    std::vector<Field>&& previous_params,
    const Field& next_param, 
    Args... args) {
  
  previous_params.push_back(next_param);
  return FieldPolynomialSolverHelper<Field>(
      std::forward<std::vector<Field> >(previous_params),
      args...);
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0
 * @tparam Field - field, in which the equation must be solved
 * @warning degree of Polynom (without leading zeros) must be <= 2
 * @param [in] enter_polynom_params -std::vector of parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return FieldEquationSolutions<Field>(class of equation answers)
 */
template<typename Field>
FieldEquationSolutions<Field> FieldPolynomialSolver(
    std::vector<Field>&& enter_polynom_params) {

  return FieldPolynomialSolverHelper<Field>(
      std::forward<std::vector<Field> >(enter_polynom_params));
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0
 * @tparam Field - field, in which the equation must be solved
 * @warning degree of Polynom (without leading zeros) must be <= 2
 * @param [in] enter_polynom_params -std::vector of parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return FieldEquationSolutions<Field>(class of equation answers)
 */
template<typename Field>
FieldEquationSolutions<Field> FieldPolynomialSolver(    
    const std::vector<Field>& enter_polynom_params) {

  std::vector<Field> polynom_params = enter_polynom_params;
  return FieldPolynomialSolverHelper<Field>(std::move(polynom_params));
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0
 * @tparam Field - field, in which the equation must be solved
 * @warning degree of Polynom (without leading zeros) must be <= 2
 * @param [in] a_0,a_1,...,a_n - parameteres of equation
 * @return FieldEquationSolutions<Field>(class of equation answers)
 */
template<typename Field, typename... Args>
FieldEquationSolutions<Field> FieldPolynomialSolver( 
    const Field& first_param, 
    Args... args) {
  
  std::vector<Field> previous_params = {first_param};
  return FieldPolynomialSolverHelper<Field>(
      std::move(previous_params), args...);
}

/**
 * @breif specialization for FieldEquationSolutions, where answer lies in field of real numbers
 * 
 */
typedef FieldEquationSolutions<double> RealEquationSolutions;
/**
 * @breif specialization for FieldEquationSolutions, where answer lies in field of complex numbers
 */
typedef FieldEquationSolutions<std::complex<double> > ComplexEquationSolutions;


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0 in real numbers
 * @warning degree of Polynom (without leading zeros) must be <= 2  
 * @param [in] a_0,a_1,...,a_n - parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return RealEquationSolutions(class of equation answers)
 */
template<typename... Args>
RealEquationSolutions PolynomialSolver(Args... args) {
  
  return FieldPolynomialSolver<double>(args...);
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0 in real numbers
 * @warning degree of Polynom (without leading zeros) must be <= 2 
 * @param [in] enter_polynom_params -std::vector of parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return RealEquationSolutions(class of equation answers)
 */
RealEquationSolutions PolynomialSolver(
    const std::vector<double>& enter_polynom_params) {
  
  return FieldPolynomialSolver<double>(enter_polynom_params);
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0 in real numbers
 * @warning degree of Polynom (without leading zeros) must be <= 2 
 * @param [in] enter_polynom_params -std::vector of parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return RealEquationSolutions(class of equation answers)
 */
RealEquationSolutions PolynomialSolver(
    std::vector<double>&& enter_polynom_params) {
  
  return FieldPolynomialSolver<double>( 
      std::forward<std::vector<double> >(enter_polynom_params));
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0 in complex numbers
 * @warning degree of Polynom (without leading zeros) must be <= 2  
 * @param [in] a_0,a_1,...,a_n - parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return ComplexEquationSolutions(class of equation answers)
 */
template<typename... Args>
ComplexEquationSolutions ComplexPolynomialSolver(Args... args) {
  
  return FieldPolynomialSolver<std::complex<double> >(args...);
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0 in complex numbers
 * @warning degree of Polynom (without leading zeros) must be <= 2 
 * @param [in] enter_polynom_params -std::vector of parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return ComplexEquationSolutions(class of equation answers)
 */
ComplexEquationSolutions ComplexPolynomialSolver(
    const std::vector<std::complex<double> >& enter_polynom_params) {

  return FieldPolynomialSolver<std::complex<double> >(enter_polynom_params);
}


/**
 * @brief Solves equation like a_0 * x^n + ... + a_(n-1) * x + a_n = 0 in complex numbers
 * @warning degree of Polynom (without leading zeros) must be <= 2 
 * @param [in] enter_polynom_params -std::vector of parameteres of equation
 * @note enter_polynom_params[0] = a_0
 * @return ComplexEquationSolutions(class of equation answers)
 */
ComplexEquationSolutions ComplexPolynomialSolver(
    std::vector<std::complex<double> >&& enter_polynom_params) {

  return FieldPolynomialSolver<std::complex<double> >( 
      std::forward<std::vector<std::complex<double> > >(enter_polynom_params));
}


//----------------------------------------------------------------------------

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

