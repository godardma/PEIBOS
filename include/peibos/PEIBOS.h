#pragma once

// Author: Maël GODARD
// Parallelepipedic Enclosure of the Image of a BOundary of a Set

#include <codac>
#include <codac-capd.h>
#include <capd/capdlib.h>
#include <cmath>
using namespace std;
using namespace codac2;

bool is_valid (Matrix A)
{
  for (int i = 0; i < A.rows(); i++)
  {
    for (int j = 0; j < A.cols(); j++)
    {
      if (std::isnan(A(i,j)) || A(i,j) == oo || A(i,j) == -oo)
      {
        return false;
      }
    }
  }
  return true;
}

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

template <typename T>
bool contains (vector<OctaSym> symmetries, OctaSym symmetry, AnalyticFunction<T>& psi_0)
{
  IntervalVector test_box =  Interval(-1.,1.)*IntervalVector::Ones(symmetry.size()-1);
  IntervalVector psi_0_x = psi_0.eval(EvalMode::NATURAL,test_box);
  
  for (OctaSym s : symmetries)
  {
    if ((s(psi_0_x)) == (symmetry(psi_0_x)))
    {
      return true;
    }
  }
  return false;
}

template <typename T>
vector<OctaSym> generate_symmetries (vector<vector<int>> generators, AnalyticFunction<T>& psi_0)
{
  vector<OctaSym> symmetries;

  // Add the generators
  for (int i = 0; i < generators.size(); i++)
  {
    OctaSym symmetry = OctaSym(generators[i]);
    symmetries.push_back(symmetry);
  }

  // Add the inverses
  for (int i = 0; i < generators.size(); i++)
  {
    OctaSym symmetry = OctaSym(generators[i]);
    if (!contains(symmetries, symmetry.invert(), psi_0))
    {
      symmetries.push_back(symmetry.invert());
    }
  }

  // Add the squares
  for (int i = 0; i < generators.size(); i++)
  {
    OctaSym symmetry = OctaSym(generators[i]);
    if (!contains(symmetries, symmetry*symmetry, psi_0))
    {
      symmetries.push_back(symmetry*symmetry);
    }
  }

  // Add the products
  for (int i = 0; i < generators.size(); i++)
  {
    for (int j = 0; j < generators.size(); j++)
    {
      if (i != j)
      {
        OctaSym symmetry1 = OctaSym(generators[i]);
        OctaSym symmetry2 = OctaSym(generators[j]);
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
double error(IntervalMatrix JJf, IntervalMatrix JJf_punc, AnalyticFunction<T>& psi_0, OctaSym symmetry, IntervalVector X)
{
  auto xc = X.mid();

  IntervalVector dX=X-xc;
  IntervalMatrix JJg_punc=JJf_punc*IntervalMatrix(symmetry.permutation_matrix())*psi_0.diff(xc);

  IntervalMatrix JJg=JJf*IntervalMatrix(symmetry.permutation_matrix())*psi_0.diff(X);

  IntervalVector E = (JJg - JJg_punc)*dX;
  Interval N = sqr(E[0]) + sqr(E[1]);

  return sqrt(N.ub());
}

Matrix inflate_flat_parallelepiped_3D (IntervalMatrix Jz, double epsilon, double rho)
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

Matrix inflate_flat_parallelepiped_2D (IntervalMatrix Jz, double epsilon, double rho)
{
  Vector a1 = Vector((Jz * 0.5 * epsilon).mid());
  Vector a2 ({ -a1[1], a1[0] });

  double norm_a1 = sqrt(sqr(a1[0]) + sqr(a1[1]));

  Matrix A ({{(1+rho/norm_a1)*a1[0], (rho/norm_a1)*a2[0]}, {(1+rho/norm_a1)*a1[1], (rho/norm_a1)*a2[1]}});

  return A;
}

Matrix inflate_flat_parallelepiped (IntervalMatrix Jz, double epsilon, double rho)
{
  if (Jz.rows() == 2)
  {
    return inflate_flat_parallelepiped_2D(Jz, epsilon, rho);
  }
  else if (Jz.rows() == 3)
  {
    return inflate_flat_parallelepiped_3D(Jz, epsilon, rho);
  }
}

// CAPD ODE ///////////////////////////////////////////////////////////////////////////////////////////////////////////

// 3D

template <typename T>
void PEIBOS3D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure3D& figure_3d)
{
  ColorMap cmap = peibos_cmap();
  
  // CAPD solver setup
  capd::IOdeSolver solver(gamma, 20);
  solver.setAbsoluteTolerance(1e-20);
  solver.setRelativeTolerance(1e-20);

  for (double tf : tfs)
  {
    capd::ITimeMap timeMap(solver);
    capd::ITimeMap timeMap_punc(solver);

    capd::interval initialTime(0.);
    capd::interval finalTime(tf);

    // Generate the symmetries from the generators
    vector<OctaSym> symmetries = generate_symmetries(generators, psi_0);
    for (int i = 0; i < symmetries.size(); i++)
    {
      OctaSym symmetry = symmetries[i];
      for (double t1 = -1; t1 < 1; t1 += epsilon)
      {
        for (double t2 = -1;t2 < 1; t2+=epsilon)
        {

          // To get the flow function and its Jacobian (monodromy matrix) for [x]
          IntervalVector X({{t1,t1+epsilon},{t2,t2+epsilon}});
          IntervalVector Y = symmetry(psi_0.eval(X))+offset;

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

          // To get the flow function and its Jacobian (monodromy matrix) for x_hat
          auto xc = X.mid();
          auto yc = (symmetry(psi_0.eval(xc))+offset).mid();

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

          // Center of the parallelepiped
          Vector z = Vector(to_codac(result).mid());
          
          // Maximum error computation
          double rho = error( JJf, JJf_punc, psi_0, symmetry, X);

          IntervalMatrix Jz = (JJf_punc * IntervalMatrix(symmetry.permutation_matrix()) * psi_0.diff(xc)).mid();

          // Inflation of the parallelepiped
          Matrix A = inflate_flat_parallelepiped(Jz, epsilon, rho);

          figure_3d.draw_parallelepiped(z, A, peibos_cmap().color(((double)i)/((double)symmetries.size()-1.0)));

        }
      }
    }
  }
}

template <typename T>
void PEIBOS3D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure3D& figure_3d)
{
  PEIBOS3D(gamma, tfs, psi_0, generators, epsilon, Vector::zero(psi_0.output_size()), figure_3d);
}

template <typename T>
void PEIBOS3D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, string output_name)
{
  // Graphical output
  Figure3D output(output_name);
  output.draw_axes();
  PEIBOS3D(gamma, tfs, psi_0, generators, epsilon, offset, output);
}

template <typename T>
void PEIBOS3D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, string output_name)
{
  // Graphical output
  Figure3D output(output_name);
  output.draw_axes();
  PEIBOS3D(gamma, tfs, psi_0, generators, epsilon, output);
}

// 2D

template <typename T>
void PEIBOS2D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure2D& figure_2d)
{
  ColorMap cmap = ColorMap::rainbow();
  ColorMap cmap_peibos = peibos_cmap();
  
  // CAPD solver setup
  capd::IOdeSolver solver(gamma, 20);
  solver.setAbsoluteTolerance(1e-20);
  solver.setRelativeTolerance(1e-20);

  for (double tf : tfs)
  {
    capd::ITimeMap timeMap(solver);
    capd::ITimeMap timeMap_punc(solver);

    capd::interval initialTime(0.);
    capd::interval finalTime(tf);

    // Generate the symmetries from the generators
    vector<OctaSym> symmetries = generate_symmetries(generators, psi_0);
    for (int i = 0; i < symmetries.size(); i++)
    {
      OctaSym symmetry = symmetries[i];
      for (double t = -1; t < 1; t += epsilon)
      {

        // To get the flow function and its Jacobian (monodromy matrix) for [x]
        IntervalVector X({{t,t+epsilon}});
        IntervalVector Y = symmetry(psi_0.eval(X))+offset;

        capd::IMatrix monodromyMatrix(2,2);
        capd::ITimeMap::SolutionCurve solution(initialTime); 
        capd::IVector c(2);
        c[0] = to_capd(Y[0]);
        c[1] = to_capd(Y[1]);
        capd::C1Rect2Set s(c);
        timeMap(finalTime, s, solution);
        capd::IVector result = timeMap(finalTime, s, monodromyMatrix);
        IntervalMatrix JJf=to_codac(monodromyMatrix);

        // To get the flow function and its Jacobian (monodromy matrix) for x_hat
        auto xc = X.mid();
        auto yc = (symmetry(psi_0.eval(xc))+offset).mid();

        capd::IMatrix monodromyMatrix_punc(2,2);
        capd::ITimeMap::SolutionCurve solution_punct(initialTime);
        capd::IVector c_punct(2);

        c_punct[0] = to_capd(yc[0]);
        c_punct[1] = to_capd(yc[1]);
        capd::C1Rect2Set s_punct(c_punct);
        timeMap_punc(finalTime, s_punct, solution_punct);      
        capd::IVector result_punct = timeMap_punc(finalTime, s_punct, monodromyMatrix_punc);
        IntervalMatrix JJf_punc=to_codac(monodromyMatrix_punc);

        // Center of the parallelepiped
        Vector z = Vector(to_codac(result).mid());
        
        // Maximum error computation
        double rho = error( JJf, JJf_punc, psi_0, symmetry, X);

        IntervalMatrix Jz = (JJf_punc * IntervalMatrix(symmetry.permutation_matrix()) * psi_0.diff(xc)).mid();

        // Inflation of the parallelepiped
        Matrix A = inflate_flat_parallelepiped(Jz, epsilon, rho);

        figure_2d.draw_parallelepiped(z, A, {cmap.color(((double)i)/((double)symmetries.size()-1.0)),cmap_peibos.color(((double)i)/((double)symmetries.size()-1.0))});

      }
    }
  }
}

template <typename T>
void PEIBOS2D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure2D& figure_2d)
{
  PEIBOS2D(gamma, tfs, psi_0, generators, epsilon, Vector::zero(psi_0.output_size()), figure_2d);
}

template <typename T>
void PEIBOS2D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, string output_name)
{
  // Graphical output
  Figure2D output(output_name,GraphicOutput::VIBES | GraphicOutput::IPE);
  output.set_window_properties({50,100},{800,800});
  PEIBOS2D(gamma, tfs, psi_0, generators, epsilon, offset, output);
}

template <typename T>
void PEIBOS2D(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, string output_name)
{
  // Graphical output
  Figure2D output(output_name,GraphicOutput::VIBES | GraphicOutput::IPE);
  output.set_window_properties({50,100},{800,800});
  PEIBOS2D(gamma, tfs, psi_0, generators, epsilon, output);
}

// Generic function, calls according to the output size of the function

template <typename T>
void PEIBOS(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, string output_name)
{
  if (gamma.dimension() == 3)
    PEIBOS3D(gamma, tfs, psi_0, generators, epsilon, offset, output_name);
  else if (gamma.dimension() == 2)
    PEIBOS2D(gamma, tfs, psi_0, generators, epsilon, offset, output_name);
}

template <typename T>
void PEIBOS(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, string output_name)
{
  if (gamma.dimension() == 3)
    PEIBOS3D(gamma, tfs, psi_0, generators, epsilon, output_name);
  else if (gamma.dimension() == 2)
    PEIBOS2D(gamma, tfs, psi_0, generators, epsilon, output_name);
}

// with 2D figure

template <typename T>
void PEIBOS(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure2D& figure_2d)
{
  assert(gamma.dimension() == 2);
  PEIBOS2D(gamma, tfs, psi_0, generators, epsilon, offset, figure_2d);
}

template <typename T>
void PEIBOS(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure2D& figure_2d)
{
  assert(gamma.dimension() == 2);
  PEIBOS2D(gamma, tfs, psi_0, generators, epsilon, figure_2d);
}

// with 3D figure

template <typename T>
void PEIBOS(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure3D& figure_3d)
{
  assert(gamma.dimension() == 3);
  PEIBOS3D(gamma, tfs, psi_0, generators, epsilon, offset, figure_3d);
}

template <typename T>
void PEIBOS(capd::IMap& gamma, vector<double> tfs, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure3D& figure_3d)
{
  assert(gamma.dimension() == 3);
  PEIBOS3D(gamma, tfs, psi_0, generators, epsilon, figure_3d);
}


// ANALYTIC FUNCTION //////////////////////////////////////////////////////////////////////////////////////////////////

// 3D

template <typename T>
void PEIBOS3D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure3D& figure_3d)
{
  ColorMap cmap = peibos_cmap();

  // Generate the symmetries from the generators
  vector<OctaSym> symmetries = generate_symmetries(generators, psi_0);
  for (int i = 0; i < symmetries.size(); i++)
  {
    OctaSym symmetry = symmetries[i];
    
    for (double t1 = -1; t1 < 1; t1 += epsilon)
    {
      for (double t2 = -1;t2 < 1; t2+=epsilon)
      {

        IntervalVector X({{t1,t1+epsilon},{t2,t2+epsilon}});
        IntervalVector Y = symmetry(psi_0.eval(X))+offset;

        IntervalMatrix JJf=f.diff(Y);

        auto xc = X.mid();
        auto yc = (symmetry(psi_0.eval(xc))+offset).mid();

        IntervalMatrix JJf_punc=f.diff(yc).mid();

        // Center of the parallelepiped
        Vector z = f.eval(yc).mid();

        // Maximum error computation
        double rho = error( JJf, JJf_punc, psi_0, symmetry, X);

        IntervalMatrix Jz = (JJf_punc * IntervalMatrix(symmetry.permutation_matrix()) * psi_0.diff(xc)).mid();

        // Inflation of the parallelepiped

        Matrix A = inflate_flat_parallelepiped(Jz, epsilon, rho);
        auto angle = acos((Jz.col(0)/Jz.col(0).norm()).dot(Jz.col(1)/Jz.col(1).norm()));

        if (Jz.col(0)==Jz.col(1) || !is_valid(A) || (abs(angle).ub())<1e-3) // handle degenerated case (and almost degenerated cases)
          {
            A = (f.eval(Y) - z).ub();
          }

        figure_3d.draw_parallelepiped(z, A, cmap.color(((double)i)/((double)symmetries.size()-1.0)));

      }
    }
  }
}

template <typename T>
void PEIBOS3D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure3D& figure_3d)
{
  PEIBOS3D(f, psi_0, generators, epsilon, Vector::zero(f.output_size()), figure_3d);
}

template <typename T>
void PEIBOS3D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, string output_name)
{
  // Graphical output
  Figure3D output(output_name);
  output.draw_axes();
  PEIBOS3D(f, psi_0, generators, epsilon, offset, output);
}

template <typename T>
void PEIBOS3D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, string output_name)
{
  // Graphical output
  Figure3D output(output_name);
  output.draw_axes();
  PEIBOS3D(f, psi_0, generators, epsilon, output);
}

// 2D

template <typename T>
void PEIBOS2D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure2D& figure_2d)
{
  ColorMap cmap = ColorMap::rainbow();
  ColorMap cmap_peibos = peibos_cmap();

  // Generate the symmetries from the generators
  vector<OctaSym> symmetries = generate_symmetries(generators, psi_0);
  for (int i = 0; i < symmetries.size(); i++)
  {
    OctaSym symmetry = symmetries[i];
    
    for (double t = -1; t < 1; t += epsilon)
    {

      IntervalVector X({{t,t+epsilon}});
      IntervalVector Y = symmetry(psi_0.eval(X))+offset;

      IntervalMatrix JJf=f.diff(Y);

      auto xc = X.mid();
      auto yc = (symmetry(psi_0.eval(xc))+offset).mid();

      IntervalMatrix JJf_punc=f.diff(yc).mid();

      // Center of the parallelepiped
      Vector z = f.eval(yc).mid();

      // Maximum error computation
      double rho = error( JJf, JJf_punc, psi_0, symmetry, X);

      IntervalMatrix Jz = (JJf_punc * IntervalMatrix(symmetry.permutation_matrix()) * psi_0.diff(xc)).mid();

      // Inflation of the parallelepiped
      Matrix A = inflate_flat_parallelepiped(Jz, epsilon, rho);

      figure_2d.draw_parallelepiped(z, A, {cmap.color(((double)i)/((double)symmetries.size()-1.0)),cmap_peibos.color(((double)i)/((double)symmetries.size()-1.0))});

    }
  }
}

template <typename T>
void PEIBOS2D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure2D& figure_2d)
{
  PEIBOS2D(f, psi_0, generators, epsilon, Vector::zero(f.output_size()), figure_2d);
}

template <typename T>
void PEIBOS2D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, string output_name)
{
  // Graphical output
  Figure2D output(output_name,GraphicOutput::VIBES | GraphicOutput::IPE);
  output.set_window_properties({50,100},{800,800});
  PEIBOS2D(f, psi_0, generators, epsilon, offset, output);
}

template <typename T>
void PEIBOS2D(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, string output_name)
{
  PEIBOS2D(f, psi_0, generators, epsilon, Vector::zero(f.output_size()), output_name);
}

// Generic function, calls according to the output size of the function

template <typename T>
void PEIBOS(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, string output_name)
{
  if (f.output_size() == 3)
  {
    PEIBOS3D(f, psi_0, generators, epsilon, offset, output_name);
  }
  else if (f.output_size() == 2)
  {
    PEIBOS2D(f, psi_0, generators, epsilon, offset, output_name);
  }
}

template <typename T>
void PEIBOS(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, string output_name)
{
  PEIBOS(f, psi_0, generators, epsilon, Vector::zero(f.output_size()), output_name);
}

// with 2D figure

template <typename T>
void PEIBOS(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure2D& figure_2d)
{
  assert(f.output_size() == 2);
  PEIBOS2D(f, psi_0, generators, epsilon, offset, figure_2d);
}

template <typename T>
void PEIBOS(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure2D& figure_2d)
{
  assert(f.output_size() == 2);
  PEIBOS2D(f, psi_0, generators, epsilon, figure_2d);
}

// with 3D figure

template <typename T>
void PEIBOS(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Vector offset, Figure3D& figure_3d)
{
  assert(f.output_size() == 3);
  PEIBOS3D(f, psi_0, generators, epsilon, offset, figure_3d);
}

template <typename T>
void PEIBOS(AnalyticFunction<T> f, AnalyticFunction<T>& psi_0, vector<vector<int>> generators , double epsilon, Figure3D& figure_3d)
{
  assert(f.output_size() == 3);
  PEIBOS3D(f, psi_0, generators, epsilon, figure_3d);
}