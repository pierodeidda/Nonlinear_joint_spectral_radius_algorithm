This repository includes the codes described in the manuscript "Nonlinear Joint Spectral Radius" coauthored by P.Deidda, F.Tudsico, N.Guglielmi (2025) and available at https://arxiv.org/abs/2507.11314

The function polytope_prenorm implements the algorithm described in the cited manuscript and computes the joint spectral radius of a family of continuous, homogeneous and order-preserving maps on a cone by computing an extremal monotone prenorm for the family. 

The function polytope norm is used to study differences in the extremal prenorm and extremal norm in the particular case of a family of linear and order-preserving maps. It assumes prior knowledge of the JSR.

The funcitons order-check is used to study wheter a point x is dominated by any point in a list of points Y with respect to the order induced by R_+^n.

The function pw_method implements the nonlinear power method for computing the spectral radius and dominant eigenvetor of a continuous, homogeneous and order-preserving map on a cone

The function reorder_hull is used to cyclically reorder the vertices of a convex hull so that 0 is always the first point in the list.

The scripts test_2maps_nonlinear and test_2maps_linear include the codes used to produce the figures in the manuscript, testing the polytope_prenorm algorithm on two simple families of continuous, homogeneous and order-preserving maps on the cone R_+^2.   


## Citation

If you use this code in your research, please cite:

P. Deidda, F. Tudsico, N. Guglielmi, "Nonlinear Joint Spectral Radius" (2025). https://arxiv.org/abs/2507.11314

## Requirements

The codes have been tested using MATLAB R2025b.