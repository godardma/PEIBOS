// Author: Maël GODARD

#include "peibos/PEIBOS.h"


int main()
{
  capd::IMap vectorField("par:sigma,rho,beta;var:x1,x2,x3;fun:10*(x2-x1),28*x1-x2-x1*x3,-2.6*x3+x1*x2;");
  vectorField.setParameter("sigma", 10.);
  vectorField.setParameter("rho", 28.);
  vectorField.setParameter("beta", 8/3);
  
  double tf=0.05;

  VectorVar X(2);
  AnalyticFunction psi0 ({X},{1/sqrt(1+sqr(X[0])+sqr(X[1])),X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))});

  Matrix generators = Matrix({{1,2,3},
                          {2,-1,3},
                          {-3,2,1}});

  double epsilon = 0.1;
  
  PEIBOS(vectorField, 0., psi0, generators, epsilon, "Atlas");
  PEIBOS(vectorField, tf, psi0, generators, epsilon, "Lorenz");

}

  
