source("~/work/functions/neighborhoods/maxneighbors.R")
library(inline)

funcs <- '
int locate_by_idyr(const int id, const int yr, const int plt, Rcpp::IntegerVector ids,
				   Rcpp::IntegerVector yrs, Rcpp::IntegerVector plts) {
	for (int i = 0, N = ids.size(); i < N; i++) {
		if (ids[i] == id && yrs[i] == yr && plts[i] == plt)
			return i;
	}
	return -1;
}


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
#include <cstdio>
'

src <- '
	Rcpp::IntegerVector targs(targets);
	Rcpp::IntegerVector targyrs(targtimes);
	Rcpp::IntegerVector targplts(targplots);
	Rcpp::IntegerVector nebs(neighbors);
	Rcpp::IntegerVector plot(plots);
	Rcpp::IntegerVector spec(spp);
	Rcpp::CharacterVector status(stat);
	Rcpp::IntegerVector xx(bqudx);
	Rcpp::IntegerVector yy(bqudy);
	Rcpp::IntegerVector tt(time);
	Rcpp::NumericVector indvar(variable);
	Rcpp::NumericVector zz(z);
	Rcpp::NumericVector slope(slp);
	Rcpp::NumericVector aspect(asp);


	int n        = Rcpp::as<int> (maxn);
	int rad      = Rcpp::as<int> (sr);
	int m        = targs.size();
	int nsize    = nebs.size();

	// Initialize matrices
	arma::Col<int> fillpoint(m, arma::fill::zeros);

	arma::mat distances(m, n, arma::fill::zeros);
	arma::mat var(m, n, arma::fill::zeros);
	arma::mat specs(m, n, arma::fill::zeros);
	arma::mat direction_x(m, n, arma::fill::zeros);
	arma::mat direction_y(m, n, arma::fill::zeros);
	arma::mat nebtag(m, n, arma::fill::zeros);
	arma::mat direction_z(m, n, arma::fill::zeros);

	// Populate matrices
	for (int row = 0; row < m; row++) {
		int targ = locate_by_idyr(targs[row], targyrs[row], targplts[row], nebs, tt, plot);\
        if (targ < 0) printf("Target not found");
		for (int neb = 0; neb < nsize; neb++) {
			if ( (targ != neb) && (plot[targ] == plot[neb]) &&
				 (std::strcmp("ALIVE", status[neb]) == 0) &&
				 ((xx[targ] + rad) > xx[neb]) && ((xx[targ] - rad) < xx[neb]) &&
				 ((yy[targ] + rad) > yy[neb]) && ((yy[targ] - rad) < yy[neb]) &&
				 (tt[targ] == tt[neb]) && (indvar[neb] >= indvar[targ]) ) {
				int current						= fillpoint[row];
				distances(row, current)			= ndist((double) xx[targ], (double) yy[targ],
													(double) xx[neb], (double) yy[neb]);
				var(row, current)				= indvar[neb];
				specs(row, current)				= spec[neb];
				direction_x(row, current)		= xx[targ] - xx[neb];
				direction_y(row, current)		= yy[targ] - yy[neb];
				direction_z(row, current)		= zz[targ] - zz[neb];
				nebtag(row, current)			= nebs[neb];
				fillpoint(row)++;
			}
		}
	}

	// Return list of matrices: list(distances = distances,
	//  variable = bas, species = species, direction_x = direction_x,
	//  direction_y = direction_y)
	return Rcpp::List::create(Rcpp::Named("distances")		  = distances,
							  Rcpp::Named("variable")		  = var,
							  Rcpp::Named("species")		  = specs,
							  Rcpp::Named("direction_x")	  = direction_x,
							  Rcpp::Named("direction_y")	  = direction_y,
							  Rcpp::Named("direction_z")	  = direction_z,
							  Rcpp::Named("radius")			  = rad,
							  Rcpp::Named("number_neighbors") = fillpoint,
							  Rcpp::Named("id")				  = targs,
							  Rcpp::Named("neighbor_id")	  = nebtag,
							  Rcpp::Named("slope")            = slope,
							  Rcpp::Named("aspect")           = aspect,
							  Rcpp::Named("yr")				  = targyrs,
                              Rcpp::Named("plot")             = targplts);
'

## slp, asp, z, direction_z
mnmRcpp <- cxxfunction(signature(targets = "int",
                                 targtimes = "int",
                                 targplots = "int",
                                 neighbors = "int",
                                 plots = "int",
                                 spp = "int",
                                 stat = "char",
                                 bqudx = "int",
                                 bqudy = "int",
                                 z = "numeric",
                                 slp = "numeric",
                                 asp = "int",
                                 time = "int",
                                 variable = "numeric",
                                 sr = "int",
                                 maxn = "int"),
                       includes = c(inclstxt, funcs),
                       body = src,
                       plugin = "RcppArmadillo", verbose = FALSE)

## Testing
## maxn <- maxnebs_disc(targets = targets$id, neighbors = neighbors$id,
##                      plots = neighbors$pplot, stat = neighbors$stat, bqudx = neighbors$bqudx,
##                      bqudy = neighbors$bqudy, time = neighbors$time, sr = 2)

## neighbors$spec <- as.integer(neighbors$spec)
## tst1 <- mnmRcpp(targets = targets$id, targtimes = targets$time,
##                 neighbors = neighbors$id,
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
    matrices <- mnmRcpp(targets = targets$id, targtimes = targets$time,
                        targplots = targets$pplot, neighbors = neighbors$id,
                        plots = neighbors$pplot, spp = neighbors$spec,
                        stat = neighbors$stat, bqudx = neighbors$bqudx,
                        bqudy = neighbors$bqudy, time = neighbors$time,
                        variable = neighbors[,ind.var], sr = sr, maxn = maxn,
                        z = neighbors$z, slp = targets$slope, asp = targets$asp)

    return ( matrices )
}

####################################################################
## Timings
## library(rbenchmark)
## benchmark(mnm(targets, neighbors, sr = 2),
##           make.neighbor.matrices(targets, neighbors, sr = 2),
##           columns = c("test", "replications", "elapsed", "relative"),
##           order = "relative", replications = 5)

