# PEIBOS

## Installation

PEIBOS is a library that provides a guaranteed enclosure of an image set using parallelepipeds.

It requires the installation of both the [CAPD](http://capd.ii.uj.edu.pl/html/index.html) and the [CODAC](https://github.com/codac-team/codac) libraries (The documentation of the version 2 of CODAC is currently being written).

To install PEIBOS in a terminal:

```
git clone git@github.com:godardma/PEIBOS.git
mkdir build && cd build
```

A local installation is highly suggested. To do so you can use the parameter "-DCMAKE_INSTALL_PREFIX" as follows:

```
cmake -DCMAKE_INSTALL_PREFIX=$HOME/peibos/build_install ..
```

Finally :

```
make install
```

If the installation was made locally, you can permanently add it to your CMake path by adding the following line in your .bashrc file (or equivalent):

```
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$HOME/peibos/build_install"
```

## Examples

Two example are available in the "example" folder. 

- Conform : The image of the unit sphere by a function in the form $y=f(x)$
- Lorenz : The integration of the unit sphere through an ODE of the form $\dot{x}=\gamma(x)$

To compile these example in a terminal in the desired folder (Lorenz or Conform):

```
mkdir build && cd build
cmake ..
make && ./peibos_example
```
 Two .obj will be created, one for the Atlas and one for the result. They can be viewed [on this website](https://3dviewer.net/index.html)