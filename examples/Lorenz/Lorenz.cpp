// Author: Maël GODARD

#include "peibos/PEIBOS.h"


int main()
{
  capd::IMap vectorField("par:sigma,rho,beta;var:x1,x2,x3;fun:sigma*(x2-x1),rho*x1-x2-x1*x3,-beta*x3+x1*x2;");
  vectorField.setParameter("sigma", 10.);
  vectorField.setParameter("rho", 28.);
  vectorField.setParameter("beta", 8/3);
  
  double tf=0.1;

  VectorVar X(2);
  AnalyticFunction psi0 ({X},{1/sqrt(1+sqr(X[0])+sqr(X[1])),X[0]/sqrt(1+sqr(X[0])+sqr(X[1])),X[1]/sqrt(1+sqr(X[0])+sqr(X[1]))});

  vector<vector<int>> generators ({{1,2,3},
                                            {-2,1,3},
                                            {3,2,-1}});

  double epsilon = 0.1;
  
  PEIBOS(vectorField, {0.}, psi0, generators, epsilon, "Atlas");

  Figure3D output ("Lorenz");

  PEIBOS(vectorField, {tf}, psi0, generators, epsilon, output);

  // PEIBOS(vectorField, {tf}, psi0, generators, epsilon, {0.,1.,0.}, output);
}

  
