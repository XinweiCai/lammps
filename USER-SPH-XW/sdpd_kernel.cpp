/*-----------------------------------------------------------------------
   MCF_LAMMPS: Multiscale Complex Fluids simulation with LAMMPS

   Added by Xin Bian (xin_bian@brown.edu),
   the CRUNCH group, Division of Applied Mathematics, Brown University

   Revision: v0, Dec. 22, 2013, Lucy and B-Spline kernels
             v1, Dec. 10, 2014, adding in Wendland kernels

   This is a user-class implementing the SPH kernels and their 
   relevant quantites.
   Reference: SPH  Price, D. J. Comput. Phys. 2012,
                   Liu, G.R., and Liu, M.B., 2003, World Scientific Publishing,
                   Morris, J.F., J. Comput. Phys., 1997,
		   Dehnen, W. and Aly, Hossam, MNRAS, 2012.
------------------------------------------------------------------------- */

#include "sdpd_kernel.h"
#include "stdio.h"

SDPDKernel::SDPDKernel(int kdim, int ktype, double kcut)
{
  
  if ( kdim < 1 || kdim > 3 )
    printf("%s\n","SDPDKernel::SDPDKernel() No such dimension in kernel!");
  
  //initialize kernel
  init(kdim, ktype, kcut);
  
#ifdef DEBUG_KERNEL
  //test print onto screen
 printf("SDPDKernel::SDPDKernel() kernel_type: %i\n", kernel_type);
 printf("SDPDKernel::SDPDKernel() rc         : %f\n", rc);
 printf("SDPDKernel::SDPDKernel() h          : %f\n", h);
 printf("SDPDKernel::SDPDKernel() sigma      : %f\n", sigma); 
 printf("SDPDKernel::SDPDKernel() coef       : %f\n", coef);
 printf("SDPDKernel::SDPDKernel() coef_g     : %f\n", coef_g);
#endif
 
}

SDPDKernel::~SDPDKernel()
{
}

void SDPDKernel::init(int kdim, int ktype, double kcut)
{
  num_dim     = kdim;
  kernel_type = ktype;
  rc          = kcut;
  
  // according to different kernels
  // we set smoothing length, normalization coefficients
  switch ( kernel_type )
    {
    case KERNEL_GAUSSIAN:
      // Gaussian kernel has one pieces and 
      // we take rc = 3h
      
      h = rc/3.0;
      
      switch ( num_dim )
	{
	case 1:    
	  
	  sigma = 1.0/sdpd_pi_sqrt;
	  break;
	  
	case 2:
	  
	  sigma = 1.0/sdpd_pi;
	  break;
	  
	case 3:
	  
	  sigma = 1.0/sdpd_pi/sdpd_pi_sqrt;
	  break;	  
	}
      
      break;
      
    case KERNEL_LUCY:
      // Lucy kernel has one pieces and 
      // rc = smoothing length.
      
      h = rc;
      
      switch ( num_dim )
	{
	case 1:    
	  
	  sigma = 5.0/4.0;
	  break;
	  
	case 2:
	  
	  sigma = 5.0/sdpd_pi;
	  break;
	  
	case 3:
	  
	  sigma = 105.0/16.0/sdpd_pi;
	  break;	  
	}
      
      break;
      
    case  KERNEL_CUBIC:
      // cubic spline has two pieces and 
      // rc = two smoothing lengthes.
      
      h = rc/2.0;
      
      switch (num_dim)
	{
	case 1: 
	  
	  sigma = 2.0/3.0;
	  break;
	  
	case 2: 
	  
	  sigma = 10.0/7.0/sdpd_pi;
	  break;
	  
	case 3:
	  
	  sigma = 1.0/sdpd_pi;
	  break;	  
	}
      
      break;
      
    case  KERNEL_QUARTIC:
      // quartic spline has two and half pieces and 
      // rc = two smoothing lengthes.
      
      h = rc/2.5;
      
      switch (num_dim)
	{
	case 1: 
	  
	  sigma = 1.0/24.0;
	  break;
	  
	case 2: 
	  
	  sigma = 96.0/1199.0/sdpd_pi;
	  break;
	  
	case 3:
	  
	  sigma = 1.0/20.0/sdpd_pi;
	  break;
	}
      
      break;
      
    case  KERNEL_QUINTIC:
      // quintic spline has three pieces and 
      // rc = three smoothing lengthes.
      
      h = rc/3.0;
      
      switch (num_dim)
	{
	case 1:    
	  
	  sigma = 1.0/120.0; //two references have conflicting numbers
	  break;
	  
	case 2:
	  
	  sigma = 7.0/478.0/sdpd_pi;
	  break;
	  
	case 3:
	  
	  sigma = 1.0/120.0/sdpd_pi; //two references have conflicting numbers
	  break;
	}
      
      break;
      
    case KERNEL_WENDLAND_C2:

      h = rc;
      
      switch (num_dim)
	{

	case 2:
	  
	  sigma = 7.0/sdpd_pi;
	  break;

	case 3:
	  
	  sigma = 21.0/2.0/sdpd_pi;
	  break;

	default:
	  printf("%s\n", "SDPDKernel::init() Wendland C2 has no such dimension!\n");

	}

      break;

    case KERNEL_WENDLAND_C4:

      h = rc;
      
      switch (num_dim)
	{
	  
	case 2:
	  
	  sigma = 9.0/sdpd_pi;
	  break;

	case 3:
	  
	  sigma = 495.0/32.0/sdpd_pi;
	  break;

	default:
	  printf("%s\n", "SDPDKernel::init() Wendland C4 has no such dimension!\n");

	}

      break;

    case KERNEL_WENDLAND_C6:

      h = rc;

      switch (num_dim)
	{
	  
	case 2:
	  
	  sigma = 78.0/7.0/sdpd_pi;
	  break;
	  
	case 3:
	  
	  sigma = 1365.0/64.0/sdpd_pi;
	  break;
	  
	default:
	  printf("%s\n", "SDPDKernel::init() Wendland C6 has no such dimension!\n");

	}

      break;
      
    default:
      printf("%s\n", "SDPDKernel::init() No such kernel type!\n");
      
    }// switch(kernel_type) 
  
  
  //calculate dimensional normalization coefficient for the kernel
  switch ( num_dim )
    {
    case 1:    
      
      coef  = sigma/h;
      break;
      
    case 2:
      
      coef  = sigma/h/h;
      break;
      
    case 3:
      
      coef  = sigma/h/h/h;
      break;
      
    default:
      printf("%s\n","SDPDKernel::init() No such dimension in kernel!");

    }
  
  
  //dimensional normalization coefficient for the kernel gradient
  coef_g = coef / h;
  
  //number of pieces
  kappa = rc / h;
  
#ifdef DEBUG_KERNEL
  printf("SDPDKernel::SDPDKernel()\n");
  printf("rc    : %g\n", rc);
  printf("h     : %g\n", h);
  printf("coef  : %g\n", coef);
  printf("coef_g: %g\n", coef_g);
  printf("kappa : %g\n", kappa);
#endif

  //parameters for kernel integral of boundary contribution
  erf_kappa  = erf(kappa);
  exp_kappa2 = exp(-kappa*kappa);
  kib_coef   = sigma * sdpd_pi;
  
  return;
}

int SDPDKernel::get_kernel_type()
{
  return kernel_type;
}

double SDPDKernel::get_cut_off()
{
  return rc;
}

double SDPDKernel::get_smoothing_length()
{
  return h;
}

double SDPDKernel::get_kappa()
{
  return kappa;
}

double SDPDKernel::get_sdpd_pi()
{
  return sdpd_pi;
}

double SDPDKernel::get_sdpd_pi_sqrt()
{
  return sdpd_pi_sqrt;
}

double SDPDKernel::get_normalization_sigma()
{
  return sigma;
}

double SDPDKernel::get_normalization_coef()
{
  return coef;
}

double SDPDKernel::get_normalization_coef_g()
{
  return coef_g;
}

double SDPDKernel::kernel(double r)
{
  double w;

   switch ( kernel_type )
     {

     case KERNEL_GAUSSIAN:
       
       w = kernel_Gaussian(r);
       break;
       
     case  KERNEL_LUCY:
       
       w = kernel_Lucy(r);
       break;
       
     case  KERNEL_CUBIC:
       
       w = kernel_cubic_spline(r);
       break;
       
     case  KERNEL_QUARTIC:
       
       w = kernel_quartic_spline(r);
       break;
       
     case KERNEL_QUINTIC:
       
       w = kernel_quintic_spline(r);
       break;

     case KERNEL_WENDLAND_C2:

       w = kernel_Wendland_C2(r);
       break;

     case KERNEL_WENDLAND_C4:

       w = kernel_Wendland_C4(r);
       break;

     case KERNEL_WENDLAND_C6:

       w = kernel_Wendland_C6(r);
       break;
       
     default: 
       printf("%s\n", "SDPDKernel::kernel() No such kernel type!\n");
     }

   return w;
}


double SDPDKernel::kernel_gradient(double r)
{
  double w_g;
  
  switch ( kernel_type )
    {
     
    case KERNEL_GAUSSIAN:
      
      w_g = kernel_Gaussian_gradient(r);
      break;
      
     case  KERNEL_LUCY:
       
      w_g = kernel_Lucy_gradient(r);
      break;
      
    case  KERNEL_CUBIC:
      
      w_g = kernel_cubic_spline_gradient(r);
      break;
      
    case  KERNEL_QUARTIC:
       
      w_g = kernel_quartic_spline_gradient(r);
       break;
  
    case KERNEL_QUINTIC:
      
      w_g = kernel_quintic_spline_gradient(r);
      break;

    case KERNEL_WENDLAND_C2:
       
      w_g = kernel_Wendland_C2_gradient(r);
      break;

    case KERNEL_WENDLAND_C4:
       
      w_g = kernel_Wendland_C4_gradient(r);
      break;

    case KERNEL_WENDLAND_C6:
       
      w_g = kernel_Wendland_C6_gradient(r);
      break;
      
    default: 
      printf("%s\n", "SDPDKernel::kernel_gradient() No such kernel type!\n");
      
    }
  
  return w_g;
}

double SDPDKernel::kernel_Gaussian(double r)
{
  //*********************************************************
  // Calculate Gaussian kernel value
  //
  // r        : inter-particle distance
  //
  // reference: Gingold, R.A.  Monaghan, J.J. MNRAS, 1977.
  //           
  //*********************************************************

  double s, w;
  
  s = r/h;
  
  w = exp(-s*s);
  
  w *= coef;
  
  return w;
}

double SDPDKernel::kernel_Lucy(double r)
{
  //*********************************************************
  // Calculate Lucy kernel value
  //
  // r        : inter-particle distance
  //
  // reference: Lucy,L. Astron.J. 1977
  //            Liu, G.R., and Liu, M.B., 2003, 
  //            World Scientific Publishing
  //*********************************************************

  double s, w;
  
  s = r/h;
  
  if ( s < 1.0 )
    { w = (1.0 + 3.0*s) * (1.0 - s) * (1.0 - s) * (1.0 - s);}
  else
    { w = 0.0;}
  
  w *= coef;
  
  return w;
}

double SDPDKernel::kernel_cubic_spline(double r)
{
  //*********************************************************
  // Calculate cubic spline kernel value
  //
  // r        : inter-particle distance
  //
  // reference: Price, D.J. J.Comput.Phys., 2012.
  //            Liu, G.R., and Liu, M.B., 2003, 
  //            World Scientific Publishing
  //*********************************************************
  double s, w;
  
  s = r/h;
  
  if ( s < 1.0 )
    { w = 1 - 3*s*s/2.0 + 3.0*s*s*s/4.0; }
  else if ( s < 2.0 )
    { w = (2-s) * (2-s) * (2-s)/4.0;}
  else
    { w = 0.0; }
  
  w *= coef;
  
  return w;
}

double SDPDKernel::kernel_quartic_spline(double r)
{
  //*********************************************************
  // Calculate quartic spline kernel value
  //
  // r        : inter-particle distance
  //
  // reference: Price, D.J. J.Comput.Phys., 2012.
  //            Liu, G.R., and Liu, M.B., 2003, 
  //            World Scientific Publishing
  //*********************************************************
  double s, w;
  double s1, s2, s3;
  double s1_4, s2_4, s3_4;
  
  s = r/h;
  
  s1 = 0.5 - s;
  s2 = 1.5 - s;
  s3 = 2.5 - s;
  
  s1_4 = 10.0*s1*s1*s1*s1;
  s2_4 = -5.0*s2*s2*s2*s2;
  s3_4 = s3*s3*s3*s3;
  
  if ( s < 0.5 ) 
    { w = s3_4 + s2_4 + s1_4; }
  else if ( s < 1.5 ) 
    { w = s3_4 + s2_4; }
  else if ( s < 2.5 ) 
    { w = s3_4; }
  else
    { w = 0.0; }
  
  w *= coef;
  
  return w;
}

double SDPDKernel::kernel_quintic_spline(double r)
{
  //*********************************************************
  // Calculate quintic kernel value
  //
  // num_dim  : number of dimension
  // r        : inter-particle distance
  // r_cut    : radius cut off
  //
  // reference: Morris, J.F., J. Comput.Phys., 1997,
  //            Price, D.J., J.Comput.Phys., 2012,
  //            Liu, G.R., and Liu, M.B., 2003, 
  //            World Scientific Publishing
  //*********************************************************
  double s, w;
  double s1, s2, s3;
  double s1_5, s2_5, s3_5;
 
  s  = r/h;
  s1 = 1.0 - s;
  s2 = 2.0 - s;
  s3 = 3.0 - s;
  
  s1_5 = 15.0*s1*s1*s1*s1*s1;
  s2_5 = -6.0*s2*s2*s2*s2*s2;
  s3_5 = s3*s3*s3*s3*s3;
  
  if ( s < 1.0 ) 
    { w = s3_5 + s2_5 + s1_5; }
  else if ( s < 2.0 ) 
    { w = s3_5 + s2_5; }
  else if ( s < 3.0 ) 
    { w = s3_5; }
  else
    { w = 0.0; }
  
  w *= coef;
  
  return w;
}


double SDPDKernel::kernel_Wendland_C2(double r)
{
  //*********************************************************
  // Calculate Wendland C2 kernel value
  //
  // num_dim  : number of dimension
  // r        : inter-particle distance
  // r_cut    : radius cut off
  //
  // reference: Dehnen, W. and Aly, H., MNRAS, 2012
  //*********************************************************

  double rr, w;
  rr = r / rc;
  
  if ( rr < 1.0 ) {
    
    w  = (1.0 - rr) * (1.0 - rr) * (1.0 - rr) * (1.0 - rr) * (1.0  + 4.0 * rr);
    w *= coef;
  }
  else
    { 
      w = 0.0;
    }
  
  return w;
}

double SDPDKernel::kernel_Wendland_C4(double r)
{
  //*********************************************************
  // Calculate Wendland C4 kernel value
  //
  // num_dim  : number of dimension
  // r        : inter-particle distance
  // r_cut    : radius cut off
  //
  // reference: Dehnen, W. and Aly, H., MNRAS, 2012
  //*********************************************************

  double rr, w, temp;
  rr = r / rc;
  
  if ( rr < 1.0 ) {
    temp = (1.0 - rr) * (1.0 - rr) * (1.0 - rr);
    w    = temp * temp * ( 1.0 + 6.0*rr + 35.0*rr*rr/3.0);
    w   *= coef;
  }
  else 
    {
      w = 0.0;
    }
  
  return w;
}

double SDPDKernel::kernel_Wendland_C6(double r)
{
  //*********************************************************
  // Calculate Wendland C6 kernel value
  //
  // num_dim  : number of dimension
  // r        : inter-particle distance
  // r_cut    : radius cut off
  //
  // reference: Dehnen, W. and Aly, H., MNRAS, 2012
  //*********************************************************

  double rr, w, temp;
  rr = r / rc;

  if ( rr < 1.0 ) {
    temp  = (1.0 - rr) * (1.0 - rr);
    temp *= temp;
    w     = temp * temp * ( 1.0 + 8.0*rr + 25.0*rr*rr + 32.0*rr*rr*rr);
    w    *= coef;
  }
  return w;
}


double SDPDKernel::kernel_Gaussian_gradient(double r)
{
  //*********************************************************
  // Calculate Gaussian kernel gradient value
  //
  // r        : inter-particle distance
  //
  // reference: Gingold, R.A.  Monaghan, J.J. MNRAS, 1977.
  //           
  //*********************************************************

  double s, w_g;
  
  s = r/h;
  
  w_g = -2.0 * s * exp(-s*s);
  
  w_g *= coef_g;
  
  return w_g;
}


double SDPDKernel::kernel_Lucy_gradient(double r)
{
  //*********************************************************
  // Calculate Lucy kernel gradient value
  //
  // r        : inter-particle distance
  //
  // reference: Lucy,L. Astron.J. 1977
  //           Liu, G.R., and Liu, M.B., 2003, 
  //           World Scientific Publishing
  //*********************************************************
  double s, w_g;
  
  s = r/h;
  
  if ( s < 1.0 ) 
    { w_g =  -12.0 * s * (1.0 - s) * (1.0 - s); }
  else 
    { w_g = 0.0; }
  
  w_g *= coef_g;
  
  return w_g;
}


double SDPDKernel::kernel_cubic_spline_gradient(double r)
{
  //*********************************************************
  // Calculate cubic spline kernel gradient value
  //
  // r        : inter-particle distance
  //
  // reference: Price, D.J. J.Comput.Phys., 2012.
  //           Liu, G.R., and Liu, M.B., 2003, 
  //           World Scientific Publishing
  //*********************************************************
  
  double s, w_g;
  
  s  = r/h;
  
  if ( s < 1.0 )
    { w_g = -3*s + 9.0*s*s/4.0; }
  else if ( s < 2.0 )
    { w_g = -3.0* (2-s) * (2-s)/4.0;}
  else
    { w_g = 0.0; }
  
  w_g *= coef_g;
  
  return w_g;
}

double SDPDKernel::kernel_quartic_spline_gradient(double r)
{
  //*********************************************************
  // Calculate quartic spline kernel gradient value
  //
  // r        : inter-particle distance
  //
  // reference: Price, D.J. J.Comput.Phys., 2012.
  //            Liu, G.R., and Liu, M.B., 2003, 
  //            World Scientific Publishing
  //*********************************************************
  double s, w;
  double s1, s2, s3;
  double s1_3, s2_3, s3_3;
  
  s = r/h;
  
  s1 = 0.5 - s;
  s2 = 1.5 - s;
  s3 = 2.5 - s;
  
  s1_3 = -40.0*s1*s1*s1;
  s2_3 = 20.0*s2*s2*s2;
  s3_3 = -4*s3*s3*s3;
  
  if ( s < 0.5 ) 
    { w = s3_3 + s2_3 + s1_3; }
  else if ( s < 1.5 ) 
    { w = s3_3 + s2_3; }
  else if ( s < 2.5 ) 
    { w = s3_3; }
  else
    { w = 0.0; }
  
  w *= coef_g;
  
  return w;
}

double SDPDKernel::kernel_quintic_spline_gradient(double r)
{
  //*********************************************************
  // Calculate quintic kernel gradient value
  //
  // r        : inter-particle distance
  //
  // reference: Morris, J.F., J. Comput.Phys., 1997
  //           Price, D.J. J.Comput.Phys., 2012.  
  //           Liu, G.R., and Liu, M.B., 2003, 
  //           World Scientific Publishing
  //*********************************************************

  double s, w_g;
  double s1, s2, s3;
  double s1_4, s2_4, s3_4;
  
  s  = r/h;
  
  s1 = 1.0 - s;
  s2 = 2.0 - s;
  s3 = 3.0 - s;
  
  s1_4 = -75.0  * s1*s1*s1*s1;
  s2_4 = 30.0 * s2*s2*s2*s2;
  s3_4 = -5.0 *   s3*s3*s3*s3;


  if ( s < 1.0 ) 
    { w_g = s3_4 + s2_4 + s1_4; }
  else if ( s < 2.0 ) 
    { w_g = s3_4 + s2_4; }
  else if ( s < 3.0 ) 
    { w_g = s3_4; }
  else
    { w_g = 0.0; }
  
  w_g *= coef_g;
  
  return w_g;
}

double SDPDKernel::kernel_Wendland_C2_gradient(double r)
{
  //*********************************************************
  // Calculate Wendland C2 kernel value
  //
  // num_dim  : number of dimension
  // r        : inter-particle distance
  // r_cut    : radius cut off
  //
  // reference: Dehnen, W. and Aly, H., MNRAS, 2012
  //*********************************************************

  double rr, w_g, temp;
  rr = r / rc;
  
  if ( rr < 1.0 ) {
    
    temp  = (1.0 - rr)*(1.0 - rr)*(1.0 - rr);
    w_g   = -4.0*temp*(1.0  + 4.0*rr) + 4.0*temp*(1.0 - rr);
    w_g  *= coef_g;
  }
  else
    { 
      w_g = 0.0;
    }
  
  return w_g;
}

double SDPDKernel::kernel_Wendland_C4_gradient(double r)
{
  //*********************************************************
  // Calculate Wendland C2 kernel value
  //
  // num_dim  : number of dimension
  // r        : inter-particle distance
  // r_cut    : radius cut off
  //
  // reference: Dehnen, W. and Aly, H., MNRAS, 2012
  //*********************************************************

  double rr, w_g, temp1, temp2;
  rr = r / rc;
  
  if ( rr < 1.0 ) {
    temp1  = (1.0 - rr)*(1.0 - rr);
    temp2 = temp1*(1.0 - rr) ;
    w_g   = -6.0*temp1*temp2*(1.0 + 6.0*rr + 35.0*rr*rr/3.0) + temp2*temp2*(6.0 + 70.0*rr/3.0);
    w_g  *= coef_g;
  }
  else
    { 
      w_g = 0.0;
    }
  
  return w_g;
}

double SDPDKernel::kernel_Wendland_C6_gradient(double r)
{
  //*********************************************************
  // Calculate Wendland C2 kernel value
  //
  // num_dim  : number of dimension
  // r        : inter-particle distance
  // r_cut    : radius cut off
  //
  // reference: Dehnen, W. and Aly, H., MNRAS, 2012
  //*********************************************************

  double rr, w_g, temp1, temp2;
  rr = r / rc;
  
  if ( rr < 1.0 ) {
    
    temp1 = (1.0 - rr)*(1.0 - rr)*(1.0 - rr);
    temp2 = temp1*(1.0 - rr);
    w_g   = -8.0*temp1*temp2*(1.0 + 8.0*rr + 25.0*rr*rr + 32.0*rr*rr*rr) + temp2*temp2*(8.0 + 50.0*rr + 96.0*rr*rr);
    w_g  *= coef_g;
  }
  else
    { 
      w_g = 0.0;
    }
  
  return w_g;
}


double SDPDKernel::kernel_integral_boundary(double r0)
{
  //account for the loss of particle contribution inside wall
  //represented by imaginary particles
  //
  //r0: distance of fluid particle i to the wall

  double kib;

  switch (kernel_type)
    { 
    case KERNEL_GAUSSIAN:
      
      //cut off = rc
      //kib  = kib_coef*(sdpd_pi_sqrt*(erf_kappa - erf(r0/h))/2.0) - kib_coef* (rc-r0) * exp_kappa2/h;
      
      //cut off = inf
      kib   = kib_coef * sdpd_pi_sqrt * (1.0 - erf(r0/h)) / 2.0;
      break;
      
    case KERNEL_LUCY:
      kib =  - (90*pow(rc,7.0) + (-105*r0 - 280*h)*pow(rc,6)  + (336*h*r0 + 252*h*h)*pow(rc,5)
		- 315*h*h *r0*pow(rc,4) - 70*pow(h,4)*pow(rc,3) + 105*pow(h,4)*r0*rc*rc  + 15*pow(r0,7)
		-56*h*pow(r0,6) + 63*h*h*pow(r0,5)
		- 35*pow(h,4)* pow(r0,3) )/(210*pow(h,4));
      kib *= 2.0* kib_coef /h/h/h;					      
      
      break;
      
    default:
      
      printf("%s\n", "SDPDKernel::kernel_integral_boundary() no such kernel type!");
   
    }
  
  return kib;
   
}

double SDPDKernel::kernel_gradient_e_integral_boundary(double r0)
{
  //account for the loss of particle contribution inside wall
  //represented by imaginary particles
  //
  //r0: distance of fluid particle i to the wall

  double kgeib;

  switch (kernel_type)
    { 
    case KERNEL_GAUSSIAN:
      
      kgeib = -kib_coef * exp(-r0*r0/h/h) / h;
	
      break;
     
    case KERNEL_LUCY:
      
      
    default:
      
      printf("%s\n", "SDPDKernel::kernel_gradient_e_integral_boundary() no such kernel type!");
   
    }
  
  return kgeib;
   
}

double SDPDKernel::kernel_gradient_v_integral_boundary(double r0)
{
  //account for the loss of particle contribution inside wall
  //represented by imaginary particles
  //
  //r0: distance of fluid particle i to the wall

  double kgvib;

  switch (kernel_type)
    { 
    case KERNEL_GAUSSIAN:
      
      //kgvib = -kib_coef * h * exp(-r0*r0/h/h) / r0 / h;
      kgvib = - (2.0*r0*r0/h/h+1) * kib_coef * h * exp(-r0*r0/h/h) / r0 / h;
	
      break;
     
      
    default:
      
      printf("%s\n", "SDPDKernel::kernel_gradient_eve_integral_boundary() no such kernel type!");
   
    }
  
  return kgvib;
   
}

double SDPDKernel::kernel_gradient_eve_tangential_integral_boundary(double r0)
{
  //account for the loss of particle contribution inside wall
  //represented by imaginary particles
  //
  //r0: distance of fluid particle i to the wall

  double kgeib;

  switch (kernel_type)
    { 
    case KERNEL_GAUSSIAN:
      
      kgeib = -kib_coef * h * exp(-r0*r0/h/h) / r0 / 2.0;
	
      break;
     
      
    default:
      
      printf("%s\n", "SDPDKernel::kernel_gradient_eve_tangential_integral_boundary() no such kernel type!");
   
    }
  
  return kgeib;
   
}

double SDPDKernel::kernel_gradient_eve_normal_integral_boundary(double r0)
{
  //account for the loss of particle contribution inside wall
  //represented by imaginary particles
  //
  //r0: distance of fluid particle i to the wall

  double kgeib;

  switch (kernel_type)
    { 
    case KERNEL_GAUSSIAN:
      
      kgeib = -kib_coef * (r0*r0 + h*h ) * exp(-r0*r0/h/h) / r0 / h;
	
      break;
     
      
    default:
      
      printf("%s\n", "SDPDKernel::kernel_gradient_eve_normal_integral_boundary() no such kernel type!");
   
    }
  
  return kgeib;
   
}

