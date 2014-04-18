library(inline)

funcs1 <- '
    int index(const int val, Rcpp::IntegerVector vec) {
        int ind = std::find(vec.begin(), vec.end(), val) - vec.begin();
        return ind;
    }

    bool mymax(int i, int j) { return i<j; }
'

inclstxt <- '
#include <cstring>
#include <algorithm>
'

src <- '
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
'
maxnebs_disc <- cxxfunction(signature(targets = "numeric",
                                      neighbors = "integer",
                                      plots = "int",
                                      stat = "char",
                                      bqudx = "int",
                                      bqudy = "int",
                                      time = "int",
                                      sr = "int"),
                            includes = c(inclstxt, funcs1),
                            body = src,
                            plugin = "Rcpp", verbose = FALSE)


## Testing
## maxnebs_disc(targets$tag, as.integer(neighbors$tag),
##              as.integer(neighbors$pplot), neighbors$stat,
##              neighbors$bqudx, neighbors$bqudy, neighbors$time, sr=2)

####################################################################
## Timings
## library(rbenchmark)
## benchmark(maxnebs_disc(targets$tag, neighbors$tag, neighbors$pplot, neighbors$stat,
##              neighbors$bqudx, neighbors$bqudy, neighbors$time, sr=2),
##           maxneighbors_disc(targets, neighbors, sr = 2),
##           columns = c("test", "replications", "elapsed", "relative"),
##           order = "relative", replications = 10)

## (1) maxnebs_disc, (2) maxneighbors_disc
##   replications elapsed relative
## 1           10   0.375    1.000
## 2           10  33.909   90.424
