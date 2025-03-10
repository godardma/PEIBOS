// Author: Maël GODARD

#include "peibos/PEIBOS.h"


int main()
{
  VectorVar y(2);
  double a = 1.4;
  double b = 0.3;
  AnalyticFunction f({y},{y[1]+1-a*sqr(y[0]),b*y[0]});
  AnalyticFunction f_id({y},{y[0],y[1]});
  
  VectorVar X(1);
  AnalyticFunction psi0 ({X},{cos(X[0]*PI/4.),sin(X[0]*PI/4.)});

  vector<vector<int>> generators ({{1,2},
                                  {-2,1}});

  double epsilon = 0.1;
  
  PEIBOS(f_id, psi0, generators, epsilon, "Atlas");
  PEIBOS(f, psi0, generators, epsilon, "Henon");
}

  
