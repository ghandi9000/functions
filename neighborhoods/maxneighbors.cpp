// Wont't compile, just worksheet for inline implementation

// User includes
#include <cstring>
#include <algorithm>

int index(const int val, Rcpp::IntegerVector vec) {
	int ind = std::find(vec.begin(), vec.end(), val) - vec.begin();
	return ind;
}


bool mymax(int i, int j) { return i<j; }

// Function body
int maxnebs_disc() {
	Rcpp::IntegerVector targs(targets);
	Rcpp::IntegerVector nebs(neighbors);
	Rcpp::IntegerVector plot(plots);
	Rcpp::CharacterVector status(stat);
	Rcpp::IntegerVector xx(bqudx);
	Rcpp::IntegerVector yy(bqudy);
	Rcpp::IntegerVector tt(time);
	int rad = Rcpp::as<int>(sr);
	Rcpp::IntegerVector counts(targs.size());
	int max_neb = 0;

	for (int i = 0, n = targs.size(); i < n; i++) {
		int targ = index(targs[i], nebs);
		int num_nebs= 0;
		for (int neb = 0; neb < nebs.size(); neb++) {
			if (targ != neb && plot[targ] == plot[neb] &&
				(std::strcmp("ALIVE", status[neb]) == 0) &&
				(xx[targ] + rad > xx[neb]) && (xx[targ] - rad < xx[neb]) &&
				(yy[targ] + rad > yy[neb]) && (yy[targ] - rad < yy[neb]) &&
				tt[targ] == tt[neb])
				num_nebs++;
		}
		counts[i] = num_nebs;
	}

	max_neb = *std::max_element(counts.begin(), counts.end(), mymax);
	return Rcpp::wrap( max_neb );
}
