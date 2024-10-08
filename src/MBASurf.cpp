#include "include/MBA.h"
#include "include/RMBA.h"
#include <boost/shared_ptr.hpp>
//#include <R.h>
//#include <Rinternals.h>
#include <iostream>

extern "C" {

  SEXP MBASurf(SEXP xyz, SEXP noX, SEXP noY, SEXP m, SEXP n, SEXP h, SEXP extend, SEXP hpts,
	       SEXP xMinDom, SEXP xMaxDom, SEXP yMinDom, SEXP yMaxDom){

    int i,j,k;
    int nProtect = 0;
    SEXP xyzPts, X, Y, Z;

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

    //init
    MBA mba(x_arr, y_arr, z_arr);
    
    mba.setDomain(REAL(xMinDom)[0],REAL(yMinDom)[0],REAL(xMaxDom)[0],REAL(yMaxDom)[0]);
  	 
    mba.MBAalg(INTEGER(m)[0], INTEGER(n)[0], INTEGER(h)[0]);

    //retrieve the spline surface and evaluate
    UCBspl::SplineSurface surface = mba.getSplineSurface(); 

    //get the surface pts
    int noU = INTEGER(noX)[0];
    int noV = INTEGER(noY)[0];

    PROTECT(X = Rf_allocVector(REALSXP, noU)); nProtect++;
    PROTECT(Y = Rf_allocVector(REALSXP, noV)); nProtect++;
    PROTECT(Z = Rf_allocMatrix(REALSXP, noV, noU)); nProtect++;

    //U=X=col and V=Y=row
    double u, du, v, dv;
    du = (surface.umax() - surface.umin())/(double)(noU-1);
    dv = (surface.vmax() - surface.vmin())/(double)(noV-1);  
    
    for (j = 0; j < noV; j++){      
      v = surface.vmin() + j*dv;
      REAL(Y)[j] = v;
      for (i = 0; i < noU; i++){
	u = surface.umin() + i*du;
	REAL(Z)[j*noU+i] = surface.f(u,v);
      }
    }
    
    for (i = 0; i < noU; i++){
      u = surface.umin() + i*du;
      REAL(X)[i] = u;
    }

    //set NA outside the convex hull
    if(!INTEGER(extend)[0]){
      for(i = 0; i < noU; i++){
	for(j = 0; j < noV; j++){
	  for(k = 0; k < Rf_length(hpts)-1; k++){
	    if(((REAL(Y)[j]-REAL(xyz)[nPts+(INTEGER(hpts)[k])-1])*(REAL(xyz)[(INTEGER(hpts)[k+1])-1]-REAL(xyz)[(INTEGER(hpts)[k])-1]) - (REAL(X)[i]-REAL(xyz)[(INTEGER(hpts)[k])-1])*(REAL(xyz)[nPts+(INTEGER(hpts)[k+1])-1]-REAL(xyz)[nPts+(INTEGER(hpts)[k])-1])) > 0.0){
	      REAL(Z)[j*noU+i] = NA_REAL;
	      break;
	    }
 	  }
	}
      }
    }

    //return obj
    SEXP retList, retListNames;
    
    PROTECT(retList = Rf_allocVector(VECSXP, 3)); nProtect++;
    PROTECT(retListNames = Rf_allocVector(STRSXP, 3)); nProtect++;

    SET_VECTOR_ELT(retList, 0, X);
    SET_VECTOR_ELT(retList, 1, Y);
    SET_VECTOR_ELT(retList, 2, Z);
    
    SET_STRING_ELT(retListNames, 0, Rf_mkChar("x"));
    SET_STRING_ELT(retListNames, 1, Rf_mkChar("y"));
    SET_STRING_ELT(retListNames, 2, Rf_mkChar("z"));

    Rf_namesgets(retList, retListNames);

    //clean-up
    mba.cleanup(2);

    //just to be sure
    x_arr.reset();
    y_arr.reset();
    z_arr.reset();

    UNPROTECT(nProtect);
    return(retList);
  }

}
