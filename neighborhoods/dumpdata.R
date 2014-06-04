findTargets <-
function(sr) {

}
funcs <-
"\nint locate_by_idyr(const int id, const int yr, Rcpp::IntegerVector ids,\n\t\t\t\t   Rcpp::IntegerVector yrs) {\n\tfor (int i = 0, N = ids.size(); i < N; i++) {\n\t\tif (ids[i] == id && yrs[i] == yr)\n\t\t\treturn i;\n\t}\n\treturn -1;\n}\n\n\ndouble ndist(double tx, double ty, double nx, double ny) {\n\tdouble dd = std::sqrt( std::pow(tx - nx, 2) + std::pow(ty - ny, 2) );\n\tif (dd == 0)\n\t\treturn 0.5;\n\telse\n\t\treturn dd;\n}\n"
funcs1 <-
"\n    int index(const int val, Rcpp::IntegerVector vec) {\n        int ind = std::find(vec.begin(), vec.end(), val) - vec.begin();\n        return ind;\n    }\n\n    bool mymax(int i, int j) { return i<j; }\n"
inclstxt <-
"\n#include <cstring>\n#include <algorithm>\n#include <cmath>\n"
make.neighbor.matrices <-
function(targets, neighbors, sr, ind.var="ba", range=12,
                                   realdist=FALSE, bigger=FALSE) {
    ## define targets and neighbors
    ## targets <- subset(dat, bqudx < (12-sr) & bqudx > (-1 + sr) & bqudy < (12 - sr) &
    ##                   bqudy > (-1 + sr) & stat=="ALIVE")
    ## neighbors <- subset(dat, bqudx < 11 & bqudx > 0 & bqudy < 11 &
    ##                     bqudy > 0 & stat=="ALIVE")
    ifelse(realdist,
           { max.neighbors <- maxneighbors_real(targets, neighbors, sr) },
           { max.neighbors <- maxneighbors_disc(targets, neighbors, sr) })
    ## initialize matrices
    distances <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    bas <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    species <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    direction_x <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    direction_y <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    ## populate matrices
    if(realdist==FALSE) {
        for(i in 1:nrow(targets)) {
            ifelse(bigger==TRUE,
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] & bqudx < targets$bqudx[i]+sr &
                          bqudx > targets$bqudx[i]-sr & bqudy < targets$bqudy[i]+sr &
                          bqudy > targets$bqudy[i]-sr &
                          get(ind.var) >= targets[,ind.var][i] &
                          time==targets$time[i]),
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] & bqudx < targets$bqudx[i]+sr &
                          bqudx > targets$bqudx[i]-sr & bqudy < targets$bqudy[i]+sr &
                          bqudy > targets$bqudy[i]-sr & time==targets$time[i]))
            if (nrow(nebs) > 0) {
                distances[i,1:nrow(nebs)] <-
                neighdist(targets$bqudx[i]*2,targets$bqudy[i]*2,
                          nebs$bqudx*2, nebs$bqudy*2, addifsame = TRUE)
                bas[i,1:nrow(nebs)] <- nebs[,ind.var]
                species[i,1:nrow(nebs)] <- nebs$spec
                direction_x[i,1:nrow(nebs)] <- targets$bqudx[i] - nebs$bqudx
                direction_y[i,1:nrow(nebs)] <- targets$bqudy[i] - nebs$bqudy
            }
        }
    }
    if(realdist==TRUE) {
        for(i in 1:nrow(targets)) {
            ifelse(bigger==TRUE,
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] &
                          neighdist(x,y,targets$x[i], targets$y[i]) <= sr &
                          get(ind.var) >= targets[,ind.var][i] &
                          time == targets$time[i]),
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] &
                          neighdist(x,y,targets$x[i], targets$y[i]) <= sr &
                          get(ind.var) >= targets[,ind.var][i] &
                          time == targets$time[i]))
            if (nrow(nebs) > 0) {
                distances[i,1:nrow(nebs)] <-
                    neighdist(targets$x[i],targets$y[i],
                              nebs$x, nebs$y, addifsame = FALSE)
                bas[i,1:nrow(nebs)] <- nebs[,ind.var]
                species[i,1:nrow(nebs)] <- nebs$spec
                direction_x[i,1:nrow(nebs)] <- targets$x[i] - nebs$x
                direction_y[i,1:nrow(nebs)] <- targets$y[i] - nebs$y
            }
        }
    }
    return( list(distances = distances, variable = bas, species = species,
                 direction_x = direction_x, direction_y = direction_y) )
}
maxnebs_disc <-
structure(function (targets, neighbors, plots, stat, bqudx, bqudy, 
    time, sr) 
.Primitive(".Call")(<pointer: 0x654c14e0>, targets, neighbors, 
    plots, stat, bqudx, bqudy, time, sr), code = "\n// includes from the plugin\n\n#include <Rcpp.h>\n\n\n#ifndef BEGIN_RCPP\n#define BEGIN_RCPP\n#endif\n\n#ifndef END_RCPP\n#define END_RCPP\n#endif\n\nusing namespace Rcpp;\n\n\n// user includes\n\n#include <cstring>\n#include <algorithm>\n\n\n    int index(const int val, Rcpp::IntegerVector vec) {\n        int ind = std::find(vec.begin(), vec.end(), val) - vec.begin();\n        return ind;\n    }\n\n    bool mymax(int i, int j) { return i<j; }\n\n\n// declarations\nextern \"C\" {\nSEXP file14585cff9( SEXP targets, SEXP neighbors, SEXP plots, SEXP stat, SEXP bqudx, SEXP bqudy, SEXP time, SEXP sr) ;\n}\n\n// definition\n\nSEXP file14585cff9( SEXP targets, SEXP neighbors, SEXP plots, SEXP stat, SEXP bqudx, SEXP bqudy, SEXP time, SEXP sr ){\nBEGIN_RCPP\n\n    Rcpp::IntegerVector targs(targets);\n\tRcpp::IntegerVector nebs(neighbors);\n\tRcpp::IntegerVector plot(plots);\n\tRcpp::CharacterVector status(stat);\n\tRcpp::IntegerVector xx(bqudx);\n\tRcpp::IntegerVector yy(bqudy);\n\tRcpp::IntegerVector tt(time);\n\tint rad = Rcpp::as<int>(sr);\n\tRcpp::IntegerVector counts(targs.size());\n\tint max_neb = 0;\n\n\tfor (int i = 0, n = targs.size(); i < n; i++) {\n\t\tint targ = index(targs[i], nebs);\n\t\tint num_nebs= 0;\n\t\tfor (int neb = 0; neb < nebs.size(); neb++) {\n\t\t\tif (targ != neb && plot[targ] == plot[neb] &&\n               (std::strcmp(\"ALIVE\", status[neb]) == 0) &&\n\t\t\t\t(xx[targ] + rad > xx[neb]) && (xx[targ] - rad < xx[neb]) &&\n\t\t\t\t(yy[targ] + rad > yy[neb]) && (yy[targ] - rad < yy[neb]) &&\n\t\t\t\ttt[targ] == tt[neb])\n\t\t\t\tnum_nebs++;\n\t\t}\n\t\tcounts[i] = num_nebs;\n\t}\n\n\tmax_neb = *std::max_element(counts.begin(), counts.end(), mymax);\n    return Rcpp::wrap( max_neb );\n\nEND_RCPP\n}\n\n\n", class = structure("CFunc", package = "inline"))
maxneighbors_disc <-
function(targets, neighbors, sr) {
    max.neighbors <- 0
    for(i in 1:nrow(targets)) {
        nebs <- subset(neighbors, pplot==targets$pplot[[i]] & tag!=targets$tag[[i]] &
                       stat=="ALIVE" & bqudx < targets$bqudx[[i]]+sr & bqudx >
                       targets$bqudx[[i]]-sr & bqudy < targets$bqudy[[i]]+sr & bqudy >
                       targets$bqudy[[i]]-sr & time==targets$time[[i]])
        max.neighbors <- max(max.neighbors, nrow(nebs), na.rm = TRUE)
    }
    max.neighbors
}
maxneighbors_real <-
function(targets, neighbors, sr) {
    max.neighbors <- 0
    for(i in 1:nrow(targets)) {
        nebs <- subset(neighbors, pplot == targets$pplot[i] & tag!=targets$tag[i] &
                       neighdist(x,y,targets$x[i], targets$y[i]) <= sr &
                       get(ind.var) >= targets[,ind.var][i] & stat == "ALIVE" &
                       time==targets$time[i])
        max.neighbors <- max(max.neighbors, nrow(nebs), na.rm = TRUE)
    }
    max.neighbors
}
mnm <-
function(targets, neighbors, sr, ind.var = "ba") {
    maxn <- maxnebs_disc(targets = targets$id, neighbors = neighbors$id,
                         plots = neighbors$pplot, stat = neighbors$stat,
                         bqudx = neighbors$bqudx, bqudy = neighbors$bqudy,
                         time = neighbors$time, sr = sr)

    neighbors$spec <- as.integer(neighbors$spec)
    slp <- unique(neighbors[!is.na(neighbors$slope), "slope"])
    asp <- unique(neighbors[!is.na(neighbors$asp), "asp"])
    matrices <- mnmRcpp(targets = targets$id, targtimes = targets$time,
                        neighbors = neighbors$id,
                        plots = neighbors$pplot, spp = neighbors$spec,
                        stat = neighbors$stat, bqudx = neighbors$bqudx,
                        bqudy = neighbors$bqudy, time = neighbors$time,
                        variable = neighbors[,ind.var], sr = sr, maxn = maxn,
                        z = neighbors$z, slp = slp, asp = asp)

    return ( matrices )
}
mnmRcpp <-
structure(function (targets, targtimes, neighbors, plots, spp, 
    stat, bqudx, bqudy, z, slp, asp, time, variable, sr, maxn) 
.Primitive(".Call")(<pointer: 0x70881450>, targets, targtimes, 
    neighbors, plots, spp, stat, bqudx, bqudy, z, slp, asp, time, 
    variable, sr, maxn), code = "\n// includes from the plugin\n#include <RcppArmadillo.h>\n#include <Rcpp.h>\n\n\n#ifndef BEGIN_RCPP\n#define BEGIN_RCPP\n#endif\n\n#ifndef END_RCPP\n#define END_RCPP\n#endif\n\nusing namespace Rcpp;\n\n\n// user includes\n\n#include <cstring>\n#include <algorithm>\n#include <cmath>\n\n\nint locate_by_idyr(const int id, const int yr, Rcpp::IntegerVector ids,\n\t\t\t\t   Rcpp::IntegerVector yrs) {\n\tfor (int i = 0, N = ids.size(); i < N; i++) {\n\t\tif (ids[i] == id && yrs[i] == yr)\n\t\t\treturn i;\n\t}\n\treturn -1;\n}\n\n\ndouble ndist(double tx, double ty, double nx, double ny) {\n\tdouble dd = std::sqrt( std::pow(tx - nx, 2) + std::pow(ty - ny, 2) );\n\tif (dd == 0)\n\t\treturn 0.5;\n\telse\n\t\treturn dd;\n}\n\n\n// declarations\nextern \"C\" {\nSEXP file14587aeb295e( SEXP targets, SEXP targtimes, SEXP neighbors, SEXP plots, SEXP spp, SEXP stat, SEXP bqudx, SEXP bqudy, SEXP z, SEXP slp, SEXP asp, SEXP time, SEXP variable, SEXP sr, SEXP maxn) ;\n}\n\n// definition\n\nSEXP file14587aeb295e( SEXP targets, SEXP targtimes, SEXP neighbors, SEXP plots, SEXP spp, SEXP stat, SEXP bqudx, SEXP bqudy, SEXP z, SEXP slp, SEXP asp, SEXP time, SEXP variable, SEXP sr, SEXP maxn ){\nBEGIN_RCPP\n\n\tRcpp::IntegerVector targs(targets);\n\tRcpp::IntegerVector targyrs(targtimes);\n\tRcpp::IntegerVector nebs(neighbors);\n\tRcpp::IntegerVector plot(plots);\n\tRcpp::IntegerVector spec(spp);\n\tRcpp::CharacterVector status(stat);\n\tRcpp::IntegerVector xx(bqudx);\n\tRcpp::IntegerVector yy(bqudy);\n\tRcpp::IntegerVector tt(time);\n\tRcpp::NumericVector indvar(variable);\n\tRcpp::NumericVector zz(z);\n\n\tint n        = Rcpp::as<int> (maxn);\n\tint rad      = Rcpp::as<int> (sr);\n\tint m        = targs.size();\n\tint nsize    = nebs.size();\n\tdouble slope = Rcpp::as<double> (slp);\n\tint aspect   = Rcpp::as<int> (asp);\n\n\t// Initialize matrices\n\tarma::Col<int> fillpoint(m, arma::fill::zeros);\n\n\tarma::mat distances(m, n, arma::fill::zeros);\n\tarma::mat var(m, n, arma::fill::zeros);\n\tarma::mat specs(m, n, arma::fill::zeros);\n\tarma::mat direction_x(m, n, arma::fill::zeros);\n\tarma::mat direction_y(m, n, arma::fill::zeros);\n\tarma::mat nebtag(m, n, arma::fill::zeros);\n\tarma::mat direction_z(m, n, arma::fill::zeros);\n\n\t// Populate matrices\n\tfor (int row = 0; row < m; row++) {\n\t\tint targ = locate_by_idyr(targs[row], targyrs[row], nebs, tt);\n\t\tfor (int neb = 0; neb < nsize; neb++) {\n\t\t\tif ( (targ != neb) && (plot[targ] == plot[neb]) &&\n\t\t\t\t (std::strcmp(\"ALIVE\", status[neb]) == 0) &&\n\t\t\t\t ((xx[targ] + rad) > xx[neb]) && ((xx[targ] - rad) < xx[neb]) &&\n\t\t\t\t ((yy[targ] + rad) > yy[neb]) && ((yy[targ] - rad) < yy[neb]) &&\n\t\t\t\t (tt[targ] == tt[neb]) && (indvar[neb] >= indvar[targ]) ) {\n\t\t\t\tint current\t\t\t\t\t\t= fillpoint[row];\n\t\t\t\tdistances(row, current)\t\t\t= ndist((double) xx[targ], (double) yy[targ],\n\t\t\t\t\t\t\t\t\t\t\t\t\t(double) xx[neb], (double) yy[neb]);\n\t\t\t\tvar(row, current)\t\t\t\t= indvar[neb];\n\t\t\t\tspecs(row, current)\t\t\t\t= spec[neb];\n\t\t\t\tdirection_x(row, current)\t\t= xx[targ] - xx[neb];\n\t\t\t\tdirection_y(row, current)\t\t= yy[targ] - yy[neb];\n\t\t\t\tdirection_z(row, current)\t\t= zz[targ] - zz[neb];\n\t\t\t\tnebtag(row, current)\t\t\t= nebs[neb];\n\t\t\t\tfillpoint(row)++;\n\t\t\t}\n\t\t}\n\t}\n\n\t// Return list of matrices: list(distances = distances,\n\t//  variable = bas, species = species, direction_x = direction_x,\n\t//  direction_y = direction_y)\n\treturn Rcpp::List::create(Rcpp::Named(\"distances\")\t\t  = distances,\n\t\t\t\t\t\t\t  Rcpp::Named(\"variable\")\t\t  = var,\n\t\t\t\t\t\t\t  Rcpp::Named(\"species\")\t\t  = specs,\n\t\t\t\t\t\t\t  Rcpp::Named(\"direction_x\")\t  = direction_x,\n\t\t\t\t\t\t\t  Rcpp::Named(\"direction_y\")\t  = direction_y,\n\t\t\t\t\t\t\t  Rcpp::Named(\"direction_z\")\t  = direction_z,\n\t\t\t\t\t\t\t  Rcpp::Named(\"radius\")\t\t\t  = rad,\n\t\t\t\t\t\t\t  Rcpp::Named(\"number_neighbors\") = fillpoint,\n\t\t\t\t\t\t\t  Rcpp::Named(\"id\")\t\t\t\t  = targs,\n\t\t\t\t\t\t\t  Rcpp::Named(\"neighbor_id\")\t  = nebtag,\n\t\t\t\t\t\t\t  Rcpp::Named(\"slope\")            = slope,\n\t\t\t\t\t\t\t  Rcpp::Named(\"aspect\")           = aspect,\n\t\t\t\t\t\t\t  Rcpp::Named(\"yr\")\t\t\t\t  = targyrs);\n\nEND_RCPP\n}\n\n\n", class = structure("CFunc", package = "inline"))
neighdist <-
function(targetx, targety, neighborx, neighbory, addifsame=FALSE) {
  if(addifsame==TRUE) {
      ifelse(sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)==0,
             return(0.5),
             return(sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)))
  }
  if(addifsame==FALSE) sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)
}
percSurround <-
function(sr, targets, nmatrices) {
    nebsize <- surround(sr)
    crowd <- sapply(1:nrow(targets), FUN = function(i) {
        numnebs <- nmatrices$number_neighbors[i]
        rows <- unique(data.frame(
            row_x = nmatrices[["direction_x"]][i, 0:numnebs],
            row_y = nmatrices[["direction_y"]][i, 0:numnebs]))
        nrow(rows) / nebsize
    })
}
src <-
"\n\tRcpp::IntegerVector targs(targets);\n\tRcpp::IntegerVector targyrs(targtimes);\n\tRcpp::IntegerVector nebs(neighbors);\n\tRcpp::IntegerVector plot(plots);\n\tRcpp::IntegerVector spec(spp);\n\tRcpp::CharacterVector status(stat);\n\tRcpp::IntegerVector xx(bqudx);\n\tRcpp::IntegerVector yy(bqudy);\n\tRcpp::IntegerVector tt(time);\n\tRcpp::NumericVector indvar(variable);\n\tRcpp::NumericVector zz(z);\n\n\tint n        = Rcpp::as<int> (maxn);\n\tint rad      = Rcpp::as<int> (sr);\n\tint m        = targs.size();\n\tint nsize    = nebs.size();\n\tdouble slope = Rcpp::as<double> (slp);\n\tint aspect   = Rcpp::as<int> (asp);\n\n\t// Initialize matrices\n\tarma::Col<int> fillpoint(m, arma::fill::zeros);\n\n\tarma::mat distances(m, n, arma::fill::zeros);\n\tarma::mat var(m, n, arma::fill::zeros);\n\tarma::mat specs(m, n, arma::fill::zeros);\n\tarma::mat direction_x(m, n, arma::fill::zeros);\n\tarma::mat direction_y(m, n, arma::fill::zeros);\n\tarma::mat nebtag(m, n, arma::fill::zeros);\n\tarma::mat direction_z(m, n, arma::fill::zeros);\n\n\t// Populate matrices\n\tfor (int row = 0; row < m; row++) {\n\t\tint targ = locate_by_idyr(targs[row], targyrs[row], nebs, tt);\n\t\tfor (int neb = 0; neb < nsize; neb++) {\n\t\t\tif ( (targ != neb) && (plot[targ] == plot[neb]) &&\n\t\t\t\t (std::strcmp(\"ALIVE\", status[neb]) == 0) &&\n\t\t\t\t ((xx[targ] + rad) > xx[neb]) && ((xx[targ] - rad) < xx[neb]) &&\n\t\t\t\t ((yy[targ] + rad) > yy[neb]) && ((yy[targ] - rad) < yy[neb]) &&\n\t\t\t\t (tt[targ] == tt[neb]) && (indvar[neb] >= indvar[targ]) ) {\n\t\t\t\tint current\t\t\t\t\t\t= fillpoint[row];\n\t\t\t\tdistances(row, current)\t\t\t= ndist((double) xx[targ], (double) yy[targ],\n\t\t\t\t\t\t\t\t\t\t\t\t\t(double) xx[neb], (double) yy[neb]);\n\t\t\t\tvar(row, current)\t\t\t\t= indvar[neb];\n\t\t\t\tspecs(row, current)\t\t\t\t= spec[neb];\n\t\t\t\tdirection_x(row, current)\t\t= xx[targ] - xx[neb];\n\t\t\t\tdirection_y(row, current)\t\t= yy[targ] - yy[neb];\n\t\t\t\tdirection_z(row, current)\t\t= zz[targ] - zz[neb];\n\t\t\t\tnebtag(row, current)\t\t\t= nebs[neb];\n\t\t\t\tfillpoint(row)++;\n\t\t\t}\n\t\t}\n\t}\n\n\t// Return list of matrices: list(distances = distances,\n\t//  variable = bas, species = species, direction_x = direction_x,\n\t//  direction_y = direction_y)\n\treturn Rcpp::List::create(Rcpp::Named(\"distances\")\t\t  = distances,\n\t\t\t\t\t\t\t  Rcpp::Named(\"variable\")\t\t  = var,\n\t\t\t\t\t\t\t  Rcpp::Named(\"species\")\t\t  = specs,\n\t\t\t\t\t\t\t  Rcpp::Named(\"direction_x\")\t  = direction_x,\n\t\t\t\t\t\t\t  Rcpp::Named(\"direction_y\")\t  = direction_y,\n\t\t\t\t\t\t\t  Rcpp::Named(\"direction_z\")\t  = direction_z,\n\t\t\t\t\t\t\t  Rcpp::Named(\"radius\")\t\t\t  = rad,\n\t\t\t\t\t\t\t  Rcpp::Named(\"number_neighbors\") = fillpoint,\n\t\t\t\t\t\t\t  Rcpp::Named(\"id\")\t\t\t\t  = targs,\n\t\t\t\t\t\t\t  Rcpp::Named(\"neighbor_id\")\t  = nebtag,\n\t\t\t\t\t\t\t  Rcpp::Named(\"slope\")            = slope,\n\t\t\t\t\t\t\t  Rcpp::Named(\"aspect\")           = aspect,\n\t\t\t\t\t\t\t  Rcpp::Named(\"yr\")\t\t\t\t  = targyrs);\n"
surround <-
function(sr, incl = TRUE) {
    ifelse(incl, (2*sr-1)^2, (2*sr-1)^2-1)
}
