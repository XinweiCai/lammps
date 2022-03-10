/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph/rhosum/adami,PairSPHRhoSumAdami)

#else

#ifndef LMP_PAIR_SPH_RHOSUMADAMI_H
#define LMP_PAIR_SPH_RHOSUMADAMI_H

#include <vector>
#include "pair.h"
#include "sdpd_kernel.h"

namespace LAMMPS_NS {

class PairSPHRhoSumAdami : public Pair {
 public:
  PairSPHRhoSumAdami(class LAMMPS *);
  virtual ~PairSPHRhoSumAdami();
  void init_style();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

 protected:
  double **cut;
  std::vector<std::vector<SDPDKernel> > kernel;
  int groupbit, first;
  double gx, gy, gz, dp, rho0, soundspeed, rp;
  void allocate();
};

}

#endif
#endif
