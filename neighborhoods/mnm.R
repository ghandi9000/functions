source("~/work/functions/neighborhoods/maxneighbors.R")
library(inline)

funcs <- '
    int index(const int val, Rcpp::IntegerVector vec) {
        int ind = std::find(vec.begin(), vec.end(), val) - vec.begin();
        return ind;
    }

    bool mymax(int i, int j) { return i<j; }

double ndist(double tx, double ty, double nx, double ny) {
	double dd = std::sqrt( std::pow(tx - nx, 2) + std::pow(ty - ny, 2) );
	if (dd == 0)
		return 0.5;
	else
		return dd;
}
'

inclstxt <- '
#include <cstring>
#include <algorithm>
#include <cmath>
'

src <- '
	Rcpp::IntegerVector targs(targets);
	Rcpp::IntegerVector nebs(neighbors);
	Rcpp::IntegerVector plot(plots);
	Rcpp::IntegerVector spec(spp);
	Rcpp::CharacterVector status(stat);
	Rcpp::IntegerVector xx(bqudx);
	Rcpp::IntegerVector yy(bqudy);
	Rcpp::IntegerVector tt(time);
	Rcpp::NumericVector indvar(variable);

	int n		= Rcpp::as<int> (maxn);
	int rad		= Rcpp::as<int> (sr);
	int m		= targs.size();
	int nsize	= nebs.size();

	// Initialize matrices
	arma::Col<int> fillpoint(m, arma::fill::zeros);
	arma::mat distances(m, n, arma::fill::zeros);
	arma::mat var(m, n, arma::fill::zeros);
	arma::mat specs(m, n, arma::fill::zeros);
	arma::mat direction_x(m, n, arma::fill::zeros);
	arma::mat direction_y(m, n, arma::fill::zeros);
	arma::mat nebtag(m, n, arma::fill::zeros);

	// Populate matrices
	for (int row = 0; row < m; row++) {
		int targ = index(targs[row], nebs);
		for (int neb = 0; neb < nsize; neb++) {
			if ( (targ != neb) && (plot[targ] == plot[neb]) &&
				 (std::strcmp("ALIVE", status[neb]) == 0) &&
				 ((xx[targ] + rad) > xx[neb]) && ((xx[targ] - rad) < xx[neb]) &&
				 ((yy[targ] + rad) > yy[neb]) && ((yy[targ] - rad) < yy[neb]) &&
				 (tt[targ] == tt[neb]) && (indvar[neb] >= indvar[targ]) ) {
				int current               = fillpoint[row];
				distances(row, current)	  = ndist((double) xx[targ], (double) yy[targ],
													(double) xx[neb], (double) yy[neb]);
				var(row, current)		  = indvar[neb];
				specs(row, current)		  = spec[neb];
				direction_x(row, current) = xx[targ] - xx[neb];
				direction_y(row, current) = yy[targ] - yy[neb];
				nebtag(row, current)	  = nebs[neb];
				fillpoint(row)++;
			}
		}
	}

	// Return list of matrices: list(distances = distances,
	//  variable = bas, species = species, direction_x = direction_x,
	//  direction_y = direction_y)
	return Rcpp::List::create(Rcpp::Named("distances")	 = distances,
							  Rcpp::Named("variable")	 = var,
							  Rcpp::Named("species")	 = specs,
							  Rcpp::Named("direction_x") = direction_x,
							  Rcpp::Named("direction_y") = direction_y,
							  Rcpp::Named("radius")		 = rad,
							  Rcpp::Named("fillpoint")	 = fillpoint,
							  Rcpp::Named("id")			 = targs,
							  Rcpp::Named("neighbor_id") = nebtag);
'

mnmRcpp <- cxxfunction(signature(targets = "int",
                                 neighbors = "int",
                                 plots = "int",
                                 spp = "int",
                                 stat = "char",
                                 bqudx = "int",
                                 bqudy = "int",
                                 time = "int",
                                 variable = "numeric",
                                 sr = "int",
                                 maxn = "int"),
                       includes = c(inclstxt, funcs),
                       body = src,
                       plugin = "RcppArmadillo", verbose = TRUE)

## Testing
## maxn <- maxnebs_disc(targets = targets$id, neighbors = neighbors$id,
##                      plots = neighbors$pplot, stat = neighbors$stat, bqudx = neighbors$bqudx,
##                      bqudy = neighbors$bqudy, time = neighbors$time, sr = 2)

## neighbors$spec <- as.integer(neighbors$spec)
## tst1 <- mnmRcpp(targets = targets$id, neighbors = neighbors$id,
##                 plots = neighbors$pplot, spp = neighbors$spec,
##                 stat = neighbors$stat, bqudx = neighbors$bqudx,
##                 bqudy = neighbors$bqudy, time = neighbors$time,
##                 variable = neighbors$ba, sr = 2, maxn = maxn)

## tst2 <- make.neighbor.matrices(targets, neighbors, sr = 2)

## R wrapper for whole process
mnm <- function(targets, neighbors, sr, ind.var = "ba") {
    maxn <- maxnebs_disc(targets = targets$id, neighbors = neighbors$id,
                         plots = neighbors$pplot, stat = neighbors$stat,
                         bqudx = neighbors$bqudx, bqudy = neighbors$bqudy,
                         time = neighbors$time, sr = sr)

    neighbors$spec <- as.integer(neighbors$spec)
    matrices <- mnmRcpp(targets = targets$id, neighbors = neighbors$id,
                plots = neighbors$pplot, spp = neighbors$spec,
                stat = neighbors$stat, bqudx = neighbors$bqudx,
                bqudy = neighbors$bqudy, time = neighbors$time,
                variable = neighbors[,ind.var], sr = sr, maxn = maxn)

    return ( matrices )
}

####################################################################
## Timings
## library(rbenchmark)
## benchmark(mnm(targets, neighbors, sr = 2),
##           make.neighbor.matrices(targets, neighbors, sr = 2),
##           columns = c("test", "replications", "elapsed", "relative"),
##           order = "relative", replications = 5)

