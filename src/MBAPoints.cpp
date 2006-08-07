#include "include/MBA.h"
#include <boost/shared_ptr.hpp>

#include <R.h>
#include <Rinternals.h>
#include <iostream>

extern "C" {

  SEXP MBAPoints(SEXP xyz, SEXP xyzEst, SEXP m, SEXP n, SEXP h) {

    int i,j;
    int nProtect = 0;
    SEXP xyzPts, Z;

    PROTECT(xyzPts = getAttrib(xyz, R_DimSymbol)); nProtect++;      
    int nPts = INTEGER(xyzPts)[0];

    typedef std::vector<double> dVec;
    boost::shared_ptr<dVec> x_arr(new std::vector<double>);
    boost::shared_ptr<dVec> y_arr(new std::vector<double>);
    boost::shared_ptr<dVec> z_arr(new std::vector<double>);

    for(i = 0; i < nPts; i++){
      x_arr->push_back(REAL(xyz)[i]);
      y_arr->push_back(REAL(xyz)[nPts+i]);
      z_arr->push_back(REAL(xyz)[2*nPts+i]);
    }

    //init
    MBA mba(x_arr, y_arr, z_arr);
    
    mba.MBAalg(INTEGER(m)[0], INTEGER(n)[0], INTEGER(h)[0]);

    //retrieve the spline surface and evaluate
    UCBspl::SplineSurface surface = mba.getSplineSurface(); 

    //get the pts
    xyzPts = getAttrib(xyzEst, R_DimSymbol);      
    int nEstPts = INTEGER(xyzPts)[0];

    PROTECT(Z = allocVector(REALSXP, nEstPts)); nProtect++;
    
    for (i = 0; i < nEstPts; i++){
      REAL(Z)[i] = surface.f(REAL(xyzEst)[i], REAL(xyzEst)[nEstPts+i]);
    }

    //clean-up
    mba.cleanup(2);

    //just to be sure
    x_arr.reset();
    y_arr.reset();
    z_arr.reset();

    UNPROTECT(nProtect);
    return(Z);
  }

}
