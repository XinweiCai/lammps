/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   Reference:Adami2012
   Date:2022-04-01 
   by xinwei
   This pair style only compute the force between boundary and fluid, 
   so it just containthe pressure force 
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph/adami/noslip/boundary,PairSPHAdamiNoslipBoundary)

#else

#ifndef LMP_PAIR_ADAMI_NOSLIP_BOUNDARY_H
#define LMP_PAIR_ADAMI_NOSLIP_BOUNDARY_H

#include <map>
#include <vector>
#include "pair.h"

namespace LAMMPS_NS {

class PairSPHAdamiNoslipBoundary : public Pair {
 public:
	 PairSPHAdamiNoslipBoundary(class LAMMPS *);
  virtual ~PairSPHAdamiNoslipBoundary();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double *soundspeed, *B, *rho0, *dp;
  double **cut,**viscosity;
  int first, groupbit;
  void allocate();
};

}

#endif
#endif
