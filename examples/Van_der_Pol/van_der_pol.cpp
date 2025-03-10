// Author: Maël GODARD

#include "peibos/PEIBOS.h"


int main()
{
  capd::IMap vectorField("var:x1,x2;fun:x2,(1-sqr(x1))*x2-x1;");

  double tf = 2.0;
  
  VectorVar X(1);
  AnalyticFunction psi0 ({X},{cos(X[0]*PI/4.),sin(X[0]*PI/4.)});

  vector<initializer_list<int>> generators ({{1,2},
                                            {-2,1}});

  double epsilon = 0.05;
  
  PEIBOS(vectorField, 0., psi0, generators, epsilon, "Atlas");
  PEIBOS(vectorField, tf, psi0, generators, epsilon, "Van der Pol");
}

  
