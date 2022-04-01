/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   Reference:Adami2012
   Date:2022-04-01
   by xinwei
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph/freesurface,PairSPHFreeSurface)

#else

#ifndef LMP_PAIR_FREESURFACE_H
#define LMP_PAIR_FREESURFACE_H

#include <map>
#include <vector>
#include "pair.h"

namespace LAMMPS_NS {

class PairSPHFreeSurface : public Pair {
 public:
	 PairSPHFreeSurface(class LAMMPS *);
  virtual ~PairSPHFreeSurface();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double *soundspeed, *B, *rho0;
  double **cut,**viscosity;
  int first;
  void allocate();
};

}

#endif
#endif
