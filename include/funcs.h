/*!
 * \file funcs.h
 * \author - Original code: SD++ developed by Patrice Castonguay, Antony Jameson,
 *                          Peter Vincent, David Williams (alphabetical by surname).
 *         - Current development: Aerospace Computing Laboratory (ACL)
 *                                Aero/Astro Department. Stanford University.
 * \version 0.1.0
 *
 * High Fidelity Large Eddy Simulation (HiFiLES) Code.
 * Copyright (C) 2014 Aerospace Computing Laboratory (ACL).
 *
 * HiFiLES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HiFiLES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HiFiLES.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include "hf_array.h"

#if defined _GPU
#include "cuda_runtime_api.h"
#include "cusparse_v2.h"
#endif

/*! evaluate lagrange basis */
double eval_lagrange(double in_r, int in_mode, hf_array<double>& in_loc_pts);

/*! evaluate derivative of lagrange basis */
double eval_d_lagrange(double in_r, int in_mode, hf_array<double>& in_loc_pts);

/*! evaluate second derivative of lagrange basis */
double eval_dd_lagrange(double in_r, int in_mode, hf_array<double>& in_loc_pts);

/*! evaluate legendre basis */
double eval_legendre(double in_r, int in_mode);

/*! evaluate derivative of legendre basis */
double eval_d_legendre(double in_r, int in_mode);

/*! evaluate derivative of vcjh basis */
double eval_d_vcjh_1d(double in_r, int in_mode, int in_order, double in_eta);

/*! evaluate derivative of OESFR basis */
double eval_d_oesfr_1d(double in_r, int in_mode, int in_order);

/*! evaluate derivative of Optimized Flux Reconstruction (OFR) basis */
double eval_d_ofr_1d(double in_r, int in_mode, int in_order);

void get_opp_3_tri(hf_array<double>& opp_3, hf_array<double>& loc_upts_tri, hf_array<double>& loc_fpts_tri, hf_array<double>& vandermonde_tri, hf_array<double>& inv_vandermonde_tri, int n_upts_per_tri, int order, double c_tri, int vcjh_scheme_tri);

void get_opp_3_dg(hf_array<double>& opp_3_dg, hf_array<double>& loc_upts_tri, hf_array<double>& loc_fpts_tri, int n_upts_per_tri, int order);

void compute_modal_filter_1d(hf_array <double>& filter_upts, hf_array<double>& vandermonde, hf_array<double>& inv_vandermonde, int N, int order);

void compute_modal_filter_tri(hf_array <double>& filter_upts, hf_array<double>& vandermonde, hf_array<double>& inv_vandermonde, int N, int order);

void compute_modal_filter_tet(hf_array <double>& filter_upts, hf_array<double>& vandermonde, hf_array<double>& inv_vandermonde, int N, int order);

void compute_filt_matrix_tri(hf_array<double>& Filt, hf_array<double>& vandermonde_tri, hf_array<double>& inv_vandermonde_tri, int n_upts_per_ele, int order, double c_tri, int vcjh_scheme_tri, hf_array<double>& loc_upts_tri);

/*! evaluate divergenge of vcjh basis on triangle */
double eval_div_dg_tri(hf_array<double> &in_loc , int in_edge, int in_edge_fpt, int in_order, hf_array<double> &in_loc_fpts_1d);

/*! get intel mkl csr 4 hf_array format (1 indexed column major) */
void hf_array_to_mklcsr(hf_array<double>& in_hf_array, hf_array<double>& out_data, hf_array<int>& out_cols, hf_array<int>& out_b, hf_array<int>& out_e);

void hf_array_to_ellpack(hf_array<double>& in_hf_array, hf_array<double>& out_data, hf_array<int>& out_cols, int& nnz_per_row);

/*! map a square to triangle element */
hf_array<double> rs_to_ab(double in_r, double in_s);

hf_array<double> rst_to_abc(double in_r, double in_s, double in_t);

/*!  helper method to evaluate the gamma function for positive integers */
double eval_gamma(int in_n);

/*!  helper method to evaluate a normalized jacobi polynomial */
double eval_jacobi(double in_r, int in_alpha, int in_beta, int in_mode);

/*!  helper method to evaluate the gradient of a normalized jacobi polynomial */
double eval_grad_jacobi(double in_r, int in_alpha, int in_beta, int in_mode);

/*! evaluate the triangle dubiner basis */
double eval_dubiner_basis_2d(double in_r, double in_s, int in_mode, int in_basis_order);

/*! helper method to evaluate d/dr of triangle dubiner basis */
double eval_dr_dubiner_basis_2d(double in_r, double in_s, int in_mode, int in_basis_order);

/*! helper method to evaluate d/ds of triangle dubiner basis */
double eval_ds_dubiner_basis_2d(double in_r, double in_s, int in_mode, int in_basis_order);

/*! evaluate the tet dubiner basis */
double eval_dubiner_basis_3d(double in_r, double in_s, double in_t, int in_mode, int in_basis_order);

/*! helper method to evaluate gradient of scalar dubiner basis*/
double eval_grad_dubiner_basis_3d(double in_r, double in_s, double in_t, int in_mode, int in_basis_order, int component);

/*! helper method to compute eta for vcjh schemes */
double compute_eta(int vjch_scheme, int order);

/*! helper method to check if number is a perfect square */
bool is_perfect_square(int in_a);

/*! helper method to check if number is a perfect cube */
bool is_perfect_cube(int in_a);

int compare_ints(const void * a, const void *b);

int index_locate_int(int value, int* hf_array, int size);

void eval_isentropic_vortex(hf_array<double>& pos, double time, double& rho, double& vx, double& vy, double& vz, double& p, int n_dims);

void eval_sine_wave_single(hf_array<double>& pos, hf_array<double>& wave_speed, double diff_coeff, double time, double& rho, hf_array<double>& grad_rho, int n_dims);

void eval_sine_wave_group(hf_array<double>& pos, hf_array<double>& wave_speed, double diff_coeff, double time, double& rho, hf_array<double>& grad_rho, int n_dims);

void eval_sphere_wave(hf_array<double>& pos, hf_array<double>& wave_speed, double time, double& rho, int n_dims);

void eval_couette_flow(hf_array<double>& pos, double in_gamma, double in_R_ref, double in_u_wall, double in_T_wall, double in_p_bound, double in_prandtl, double time, double& ene, hf_array<double>& grad_ene, int n_dims);

void eval_poly_ic(hf_array<double>& pos, double rho, hf_array<double>& ics, int n_dims);

int factorial(int in_n);

/*! Functions used in evaluation of shape functions and its 1st and 2nd derivatives
BEGIN:*/
// Convolution function: returns hf_array representation of polynomial that is result of multiplication of polynomial1 and polynomial2
hf_array<double> convol(hf_array<double> & polynomial1, hf_array<double> & polynomial2);

// LagrangeP function: returns lagrange polynomial of order "order", with value of unity at given node number "node", after substituting
// polynomial "subs" where independent variable in polynomial goes
/*
 Example: LagrangeP([0.1 0.3],2,[1 1]) returns:
 lagrange polynomial l^{1} _{2} (x + 1) = ( (x+1) - 0.1)/(0.3 - 0.1)
 note that l^{1} _{2} (x) = 0 at x = xi(1) = 0.1
 l^{1} _{2} (x) = 1 at x = xi(2) = 0.3
*/
hf_array<double> LagrangeP(int order, int node, hf_array<double> & subs);

// shapePoly4Tri function: the shape function T_I(r) in the polynomial format
// it is computed as described in Hughes pp 166
/*
% Array values are coefficients of monomials of increasing order
% I : index of node in triangle along lines (not global indeces)
% following node ordering from Hughes, pp 169
% nNodesSide: number of nodes in each side
*/
hf_array<double> shapePoly4Tri(int in_index, int nNodesSide);


// multPoly function: multiplies polynomials symbolically by stacking them (puts them in different rows)
template <typename T>
hf_array<T> multPoly(hf_array<T> & p1, hf_array<T> & p2);

// nodeFunctionTri function: returns the complete shape function of triangles at a specific global node in_index
// given in_n_spts, the total number of nodes  in the triangle, // and the index_location_hf_array:
// first row in index_location_hf_array contains indeces of r arranged in ascending global node number;
// second row contains indeces of s arranged in ascending global node number;
// third row contains indeces of t arranged in ascending global node number;
// refer to Hughes pp 169 to see link between r/s indeces ordering and global indeces ordering

hf_array<double> nodeFunctionTri(int in_index, int in_n_spts, hf_array<int> & index_location_hf_array);


// linkTriangleNodes function: returns hf_array with three rows as described above;
// output from this function is eventually passed as a parameter to the nodeFunctionTri function

hf_array<int> linkTriangleNodes(int in_n_spts);


// addPoly function: returns a 3D matrix, the stacking in the 3rd dimension represents polynomial addition
// adds polynomials by placing them in different layers

hf_array<double> addPoly(hf_array<double> & p1, hf_array<double> & p2);

// diffPoly function: returns a 3D matrix that represents a polynomial differentiated with respect to a
// dimension.

hf_array<double> diffPoly(hf_array<double> & p, hf_array<int> & term2Diff);

// evalPoly function: returns a double, which is the value of a polynomial "p" evaluated at coordinates coords;
// the height of matrix representing polynomial p must equal the number of elements (columns) of hf_array coords

double evalPoly(hf_array<double> p, hf_array<double> coords);

// createEquispacedArray: returns hf_array with nPoints values equispaced between a and b (in that order)
hf_array<double> createEquispacedArray(double a, double b, int nPoints);

// Check if all contents of polynomial are zero
template <typename T>
bool iszero(hf_array<T> & poly);

// Calculate the number of sides given the number of nodes in triangle
inline int calcNumSides(int nNodes)
{
  return int ( 0.5*(-1 + sqrt( 1 + 8*double(nNodes) ) ) ) ;
}


// eval_dd_nodal_s_basis_new function: new implementation of function that finds nth derivatives with
// respect to r or s at each of the triangle nodes
void eval_dn_nodal_s_basis(hf_array<double> &dd_nodal_s_basis,
                           hf_array<double> in_loc, int in_n_spts, int n_deriv);

/*! Linear equation solution by Gauss-Jordan elimination from Numerical Recipes (http://www.nr.com/) */
void gaussj(int n, hf_array<double>& A, hf_array<double>& b);

/*! Filter resolution function used with Gaussian filter*/
double flt_res(int N, hf_array<double>& wf, hf_array<double>& B, double k_0, double k_c, int ctype);

/*! Set an hf_array to zero*/
void zero_hf_array(hf_array <double>& in_hf_array);

/*! method to add together two hf_arrays M1 and M2*/
hf_array <double> add_hf_arrays(hf_array <double>& M1, hf_array <double>& M2);

/*! method to multiply together two 2-dimensional hf_arrays M1(L*M) by M2(M*N)*/
hf_array <double> mult_hf_arrays(hf_array <double>& M1, hf_array <double>& M2);

/*! method to get inverse of a square matrix*/
hf_array <double> inv_hf_array(hf_array <double>& input);

/*! method to get transpose of a square hf_array*/
hf_array <double> transpose_hf_array(hf_array <double>& in_hf_array);

/*! END */


