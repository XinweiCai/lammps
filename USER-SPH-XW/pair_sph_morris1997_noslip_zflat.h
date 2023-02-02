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

PairStyle(sph/morris1997/noslip/zflat,PairSPHMorris1997NoslipZFlat)

#else

#ifndef LMP_PAIR_MORRIS_NOSLIP_ZFLAT_H
#define LMP_PAIR_MORRIS_NOSLIP_ZFLAT_H

#include <map>
#include <vector>
#include "pair.h"
#include "sdpd_kernel.h"

namespace LAMMPS_NS {

class PairSPHMorris1997NoslipZFlat : public Pair {
 public:
  PairSPHMorris1997NoslipZFlat(class LAMMPS *);
  virtual ~PairSPHMorris1997NoslipZFlat();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double *soundspeed, *B;
  double **cut,**viscosity;
  std::vector<std::vector<SDPDKernel> > kernel;
  int **rigid_group;
  double **zflat, **beta, **vsurface;
  int first;

  void allocate();
};

}

#endif
#endif
