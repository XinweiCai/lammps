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
   Date:2021-05
   by xinwei
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph/taitwater/adami/noslip,PairSPHTaitwaterAdamiNoslip)

#else

#ifndef LMP_PAIR_TAITWATER_ADAMI_NOSLIP_H
#define LMP_PAIR_TAITWATER_ADAMI_NOSLIP_H

#include <map>
#include <vector>
#include "pair.h"
#include "sdpd_kernel.h"

namespace LAMMPS_NS {

class PairSPHTaitwaterAdamiNoslip : public Pair {
 public:
	 PairSPHTaitwaterAdamiNoslip(class LAMMPS *);
  virtual ~PairSPHTaitwaterAdamiNoslip();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double *soundspeed, *B, *rho0, *dp;
  double **cut,**viscosity;
  std::vector<std::vector<SDPDKernel> > kernel;
  int first, groupbit;
  void allocate();
};

}

#endif
#endif
