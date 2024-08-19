#include "include/MBA.h"
#include "include/RMBA.h"
#include <boost/shared_ptr.hpp>
//#include <R.h>
//#include <Rinternals.h>
#include <string>

extern "C" {

  SEXP MBAPoints(SEXP xyz, SEXP xyzEst, SEXP m, SEXP n, SEXP h, SEXP extend, SEXP verbose) {

    int i,j;
    int nProtect = 0;
    SEXP xyzPts, Z;

    //get the surface points
    PROTECT(xyzPts = Rf_getAttrib(xyz, R_DimSymbol)); nProtect++;      
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

    double maxXSurf = *std::max_element((*x_arr).begin(), (*x_arr).end());
    double minXSurf = *std::min_element((*x_arr).begin(), (*x_arr).end());
    double maxYSurf = *std::max_element((*y_arr).begin(), (*y_arr).end());
    double minYSurf = *std::min_element((*y_arr).begin(), (*y_arr).end());

    //get the points
    xyzPts = Rf_getAttrib(xyzEst, R_DimSymbol);      
    int nEstPts = INTEGER(xyzPts)[0];

    double maxXPt = REAL(xyzEst)[0];
    double minXPt = REAL(xyzEst)[0];
    double maxYPt = REAL(xyzEst)[nEstPts];
    double minYPt = REAL(xyzEst)[nEstPts];

    for(i = 0; i < nEstPts; i++){
      if(REAL(xyzEst)[i] > maxXPt) maxXPt = REAL(xyzEst)[i];
      if(REAL(xyzEst)[i] < minXPt) minXPt = REAL(xyzEst)[i];
      if(REAL(xyzEst)[nEstPts+i] > maxYPt) maxYPt = REAL(xyzEst)[nEstPts+i];
      if(REAL(xyzEst)[nEstPts+i] < minYPt) minYPt = REAL(xyzEst)[nEstPts+i];
    }

    std::string extendDirection = "";
    if(INTEGER(extend)[0]){
      if(maxXSurf < maxXPt){
	maxXSurf = maxXPt;
	extendDirection="+x ";
      }
      if(minXSurf > minXPt){
	minXSurf = minXPt;
	extendDirection+="-x ";
      }
      if(maxYSurf < maxYPt){
	maxYSurf = maxYPt;
	extendDirection+="+y ";
      }
      if(minYSurf > minYPt){ 
	minYSurf = minYPt;
	extendDirection+="-y ";
      }
    }

    //init
    MBA mba(x_arr, y_arr, z_arr);
    
    mba.setDomain(minXSurf, minYSurf, maxXSurf, maxYSurf);

    mba.MBAalg(INTEGER(m)[0], INTEGER(n)[0], INTEGER(h)[0]);

    //retrieve the spline surface and evaluate
    UCBspl::SplineSurface surface = mba.getSplineSurface(); 

    double umin, vmin, umax, vmax;
    surface.getDomain(umin, vmin, umax, vmax);
    
    PROTECT(Z = Rf_allocVector(REALSXP, nEstPts)); nProtect++;
    
    int ptsOutSide = 0;
    for (i = 0; i < nEstPts; i++){
      if(REAL(xyzEst)[i] < umin || REAL(xyzEst)[i] > umax ||
	 REAL(xyzEst)[nEstPts+i] < vmin || REAL(xyzEst)[nEstPts+i] > vmax){
	REAL(Z)[i] = NA_REAL;
	ptsOutSide++;
      }else{ 
	REAL(Z)[i] = surface.f(REAL(xyzEst)[i], REAL(xyzEst)[nEstPts+i]);
      }
    }

    if(INTEGER(verbose)[0]){
      if(extendDirection != "")
	Rf_warning("domain extended in the %sdirection(s)\n", extendDirection.c_str());
      
      if(ptsOutSide)
	Rf_warning("%i point(s) fell outside the domain and were set to NA\n", ptsOutSide);
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
