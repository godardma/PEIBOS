#pragma once

// Author: Maël GODARD
// Parallelepipedic Enclosure of the Image of a BOundary of a Set

#include <codac>
#include <codac-capd.h>
#include <capd/capdlib.h>
#include <cmath>
using namespace std;
using namespace codac2;

static ColorMap peibos_cmap()
{
  ColorMap cmap( Model::HSV );
  int i = 0;
        for(int h = 300 ; h > 0 ; h-=10)
        {
          cmap[i]=Color({(float)h,50.,100.,50.},Model::HSV);
          i++;
        }
  return cmap;
}


double distance_from_line_to_origin(Eigen::Matrix<double,3,1> a, Eigen::Matrix<double,3,1> b)
{
  return (a.cross(b)).norm()/((b-a).norm());
}


Matrix from_cauchy (Vector p)
{
  Matrix perm = Matrix::Zero(p.size(),p.size());
  for (int i = 0; i < p.size(); i++)
  {
    perm(abs(p[i])-1,i)=sign(p[i]);
  }
  return perm;
}

template <typename T>
bool contains (vector<Matrix> symmetries, Matrix symmetry, AnalyticFunction<T>& psi_0)
{
  IntervalVector test_box =  Interval(-1.,1.)*IntervalVector::Ones(symmetry.cols()-1);
  IntervalVector psi_0_x = psi_0.eval(EvalMode::NATURAL,test_box);
  
  for (Matrix s : symmetries)
  {
    if ((s*psi_0_x) == (symmetry*psi_0_x))
    {
      return true;
    }
  }
  return false;
}

template <typename T>
vector<Matrix> generate_symmetries (Matrix generators, AnalyticFunction<T>& psi_0)
{
  vector<Matrix> symmetries;

  // Add the generators
  for (int i = 0; i < generators.rows(); i++)
  {
    Matrix symmetry = from_cauchy(generators.row(i));
    symmetries.push_back(symmetry);
  }

  // Add the inverses
  for (int i = 0; i < generators.rows(); i++)
  {
    Matrix symmetry = from_cauchy(generators.row(i));
    if (!contains(symmetries, symmetry.inverse(), psi_0))
    {
      symmetries.push_back(symmetry.inverse());
    }
  }

  // Add the squares
  for (int i = 0; i < generators.rows(); i++)
  {
    Matrix symmetry = from_cauchy(generators.row(i));
    if (!contains(symmetries, symmetry*symmetry, psi_0))
    {
      symmetries.push_back(symmetry*symmetry);
    }
  }

  // Add the products
  for (int i = 0; i < generators.rows(); i++)
  {
    for (int j = 0; j < generators.rows(); j++)
    {
      if (i != j)
      {
        Matrix symmetry1 = from_cauchy(generators.row(i));
        Matrix symmetry2 = from_cauchy(generators.row(j));
        if (!contains(symmetries, symmetry1*symmetry2, psi_0))
        {
          symmetries.push_back(symmetry1*symmetry2);
        }
      }
    }
  }
  return symmetries;
}

template <typename T>
double error(IntervalMatrix JJf, IntervalMatrix JJf_punc, AnalyticFunction<T>& psi_0, IntervalMatrix symmetry, IntervalVector X)
{
  auto xc = X.mid();

  IntervalVector dX=X-xc;

  IntervalMatrix JJg_punc=JJf_punc*symmetry*psi_0.diff(xc);
  IntervalMatrix JJg=JJf*symmetry*psi_0.diff(X);

  IntervalVector E = (JJg - JJg_punc)*dX;
  Interval N = sqr(E[0]) + sqr(E[1]);

  return sqrt(N.ub());
}

Matrix inflate_flat_parallelepiped (IntervalMatrix Jz, double epsilon, double rho)
{
  Eigen::Matrix<double,3,1> a1 ((Jz.col(0) * 0.5 * epsilon).mid());
  Eigen::Matrix<double,3,1> a2 ((Jz.col(1) * 0.5 * epsilon).mid());
  Eigen::Matrix<double,3,1> a3 = a1.cross(a2);

  Eigen::Matrix<double,3,1> v1=a1-a2;
  Eigen::Matrix<double,3,1> v2=a1+a2;
  Eigen::Matrix<double,3,1> v3=-a1+a2;

  double d12=distance_from_line_to_origin(v1, v2);
  double d13=distance_from_line_to_origin(v2, v3);

  double d=max(d12,d13);
  
  a1*=(1+rho/d);
  a2*=(1+rho/d);
  a3*=(rho/a3.norm());

  Matrix A = Matrix({{a1[0], a2[0], a3[0]}, {a1[1], a2[1], a3[1]}, {a1[2], a2[2], a3[2]}});
  
  return A;
}

template <typename T>
void PEIBOS(capd::IMap& gamma, double tf, AnalyticFunction<T>& psi_0, Matrix generators , double epsilon, string output_name)
{
  Figure3D output(output_name);
  output.draw_axes();
  
  capd::IOdeSolver solver(gamma, 20);
  solver.setAbsoluteTolerance(1e-20);
  solver.setRelativeTolerance(1e-20);

  
  capd::ITimeMap timeMap(solver);
  capd::ITimeMap timeMap_punc(solver);

  capd::interval initialTime(0.);
  capd::interval finalTime(tf);

  vector<Matrix> symmetries = generate_symmetries(generators, psi_0);
  for (int i = 0; i < symmetries.size(); i++)
  {
    IntervalMatrix symmetry = symmetries[i];
    for (double t1 = -1; t1 < 1; t1 += epsilon)
    {
      for (double t2 = -1;t2 < 1; t2+=epsilon)
      {

      IntervalVector X({{t1,t1+epsilon},{t2,t2+epsilon}});
      IntervalVector Y = symmetry*psi_0.eval(X);

      capd::IMatrix monodromyMatrix(3,3);
      capd::ITimeMap::SolutionCurve solution(initialTime); 
      capd::IVector c(3);
      c[0] = to_capd(Y[0]);
      c[1] = to_capd(Y[1]);
      c[2] = to_capd(Y[2]);
      capd::C1Rect2Set s(c);
      timeMap(finalTime, s, solution);
      capd::IVector result = timeMap(finalTime, s, monodromyMatrix);
      IntervalMatrix JJf=to_codac(monodromyMatrix);


      auto xc = X.mid();
      auto yc = (symmetry*psi_0.eval(xc)).mid();

      capd::IMatrix monodromyMatrix_punc(3,3);
      capd::ITimeMap::SolutionCurve solution_punct(initialTime);
      capd::IVector c_punct(3);

      c_punct[0] = to_capd(yc[0]);
      c_punct[1] = to_capd(yc[1]);
      c_punct[2] = to_capd(yc[2]);
      capd::C1Rect2Set s_punct(c_punct);
      timeMap_punc(finalTime, s_punct, solution_punct);      
      capd::IVector result_punct = timeMap_punc(finalTime, s_punct, monodromyMatrix_punc);
      IntervalMatrix JJf_punc=to_codac(monodromyMatrix_punc);

      Vector z = Vector(to_codac(result).mid());
      
      double rho = error( JJf, JJf_punc, psi_0, symmetry, X);

      IntervalMatrix Jz = (JJf_punc * symmetry*psi_0.diff(X)).mid();

      Matrix A = inflate_flat_parallelepiped(Jz, epsilon, rho);

      output.draw_parallelepiped(z, A, peibos_cmap().color(((double)i)/((double)symmetries.size()-1.0)));

      }
      
    }
  }

}

template <typename T>
void PEIBOS(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, Matrix generators , double epsilon, string output_name)
{
  Figure3D output(output_name);
  output.draw_axes();

  vector<Matrix> symmetries = generate_symmetries(generators, psi_0);
  for (int i = 0; i < symmetries.size(); i++)
  {
    IntervalMatrix symmetry = symmetries[i];
    for (double t1 = -1; t1 < 1; t1 += epsilon)
    {
      for (double t2 = -1;t2 < 1; t2+=epsilon)
      {

      IntervalVector X({{t1,t1+epsilon},{t2,t2+epsilon}});
      IntervalVector Y = symmetry*psi_0.eval(X);

      IntervalMatrix JJf=f.diff(Y);

      auto xc = X.mid();
      auto yc = (symmetry*psi_0.eval(xc)).mid();

      IntervalMatrix JJf_punc=f.diff(yc).mid();

      Vector z = f.eval(yc).mid();
      
      double rho = error( JJf, JJf_punc, psi_0, symmetry, X);

      IntervalMatrix Jz = (JJf_punc * symmetry*psi_0.diff(X)).mid();

      Matrix A = inflate_flat_parallelepiped(Jz, epsilon, rho);

      output.draw_parallelepiped(z, A, peibos_cmap().color(((double)i)/((double)symmetries.size()-1.0)));

      }
      
    }
  }

}