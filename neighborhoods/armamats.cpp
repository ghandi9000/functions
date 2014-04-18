#define  ARMA_DONT_USE_WRAPPER
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
	mat A = randu<mat>(4,5);
	mat B = randu<mat>(4,5);
	colvec v(10, fill::zeros);
	A(1, 0) = 10;
	v(0)++;
	v(1)++; v(1)++;
	cout << v << endl;
	cout << A << endl;
	cout << A*B.t() << endl;
  
	return 0;
}
