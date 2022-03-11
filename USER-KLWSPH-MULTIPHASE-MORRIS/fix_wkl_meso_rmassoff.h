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

#ifdef FIX_CLASS

FixStyle(wkl/meso/rmassoff,FixWklMesoRmassoff)

#else

#ifndef LMP_FIX_WKL_MESO_RMASSOFF_H
#define LMP_FIX_WKL_MESO_RMASSOFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWklMesoRmassoff : public Fix {
 public:
  FixWklMesoRmassoff(class LAMMPS *, int, char **);
  int setmask();
  virtual void init();
  virtual void setup_pre_force(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void reset_dt();

 private:
  class NeighList *list;
 protected:
  double dtv,dtf;
  double *step_respa;
  int mass_require;

  class Pair *pair;
};

}

#endif
#endif
