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

PairStyle(interface/hu,PairInterfaceHu)

#else

#ifndef LMP_PAIR_INTERFACE_HU_H
#define LMP_PAIR_INTERFACE_HU_H

#include "pair.h"

namespace LAMMPS_NS {

class PairInterfaceHu : public Pair {
 public:
  PairInterfaceHu(class LAMMPS *);
  virtual ~PairInterfaceHu();
  void init_style();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
 // virtual double single(int, int, int, int, double, double, double, double &);
  //int pack_forward_comm(int, int *, double *, int, int *);
  //void unpack_forward_comm(int, int, double *);

 protected:
  double **cut;
  int nstep, first;
  
  double typefor_one, typefor_two;
  double tencoeff_one, tencoeff_two;
  double **tencoeff;

  void allocate();
};

}

#endif
#endif
