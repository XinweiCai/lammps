/*-----------------------------------------------------------------------
   MCF_LAMMPS: Multiscale Complex Fluids simulation with LAMMPS

   Added by Xin Bian (xin_bian@brown.edu),
   the CRUNCH group, Division of Applied Mathematics, Brown University

   Revision: v0, Dec. 22, 2013, Lucy and B-Spline kernels
             v1, Dec. 10, 2014, adding in Wendland kernels

   This is a user-class implementing the SPH kernels and their 
   relevant quantites.
   Reference: SPH  Price, D. J.Comput.Phys. 2012,
                   Liu, G.R., and Liu, M.B., 2003, World Scientific Publishing,
                   Morris, J.F., J. Comput.Phys., 1997,
		   Dehnen, W. and Aly, Hossam, MNRAS, 2012.
------------------------------------------------------------------------- */

#ifndef LMP_SDPD_KERNEL_H
#define LMP_SDPD_KERNEL_H

#define KERNEL_GAUSSIAN    1
#define KERNEL_LUCY        2
#define KERNEL_CUBIC       3
#define KERNEL_QUARTIC     4
#define KERNEL_QUINTIC     5
#define KERNEL_WENDLAND_C2 6
#define KERNEL_WENDLAND_C4 7
#define KERNEL_WENDLAND_C6 8

#define sdpd_pi  3.1415926535897932
#define sdpd_pi_sqrt 1.772453851
//sqrt(sdpd_pi);
 
//#define DEBUG_KERNEL

#include <math.h>

class SDPDKernel {
  
  public:
  SDPDKernel(){};
  SDPDKernel(int, int, double);
  ~SDPDKernel();
  void   init(int, int, double);

  // value of kernels
  double kernel(double);
  double kernel_Gaussian(double);
  double kernel_Lucy(double);
  double kernel_cubic_spline(double);
  double kernel_quartic_spline(double);
  double kernel_quintic_spline(double);
  
  double kernel_Wendland_C2(double);
  double kernel_Wendland_C4(double);
  double kernel_Wendland_C6(double);

  // gradient of kernels
  double kernel_gradient(double);
  double kernel_Gaussian_gradient(double);
  double kernel_Lucy_gradient(double);
  double kernel_cubic_spline_gradient(double);
  double kernel_quartic_spline_gradient(double);
  double kernel_quintic_spline_gradient(double);

  double kernel_Wendland_C2_gradient(double);
  double kernel_Wendland_C4_gradient(double);
  double kernel_Wendland_C6_gradient(double);
    
  // integrals of kernel relevant values
  // They are used for planar wall boundary,
  // which is represented by imaginary particles
  // note the integrals are performed within
  // region smaller than half-sphere (inside wall)
  double kernel_integral_boundary(double);
  double kernel_gradient_e_integral_boundary(double);
  double kernel_gradient_v_integral_boundary(double);
  double kernel_gradient_eve_tangential_integral_boundary(double);
  double kernel_gradient_eve_normal_integral_boundary(double);
  
  // return parameters of the kernel
  double get_sdpd_pi();
  double get_sdpd_pi_sqrt();
  int    get_kernel_type();
  double get_cut_off();
  double get_smoothing_length();
  double get_kappa();
  double get_normalization_sigma();  
  double get_normalization_coef();
  double get_normalization_coef_g();
  
protected:
  
  // num_dim    : number of dimension
  // kernel_type: indicate which kernel to use
  // rc         : cut off radius of the kernel
  // h          : smoothing length of the kernel
  // kappa      : rc/h
  // sigma      : normalization coefficient for the kernel
  // coef       : dimensional normalization coefficient for the kernel
  // coef_g     : dimensional normalization coefficient for the kernel gradient

  int num_dim, kernel_type;

  double rc, h, kappa;
  double sigma;
 //coef, coef_g means differently for B-Spline and Wendland functions
  double coef, coef_g;

//parameters for kernel integral of boundary contribution
  double erf_kappa, exp_kappa2;
  double kib_coef;  
 
};

#endif
