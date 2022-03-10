/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "fix_meso_adami2013.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMesoAdami2013::FixMesoAdami2013(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix meso command requires atom_style with both energy and density");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix meso command");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixMesoAdami2013::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoAdami2013::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}


void FixMesoAdami2013::setup_pre_force(int /*vflag*/)

{
  double **vmom = atom->vest;
  double **vtr = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
     vmom[i][0] = vtr[i][0];
     vmom[i][1] = vtr[i][1];
     vmom[i][2] = vtr[i][2];
     vtr[i][0] =  vtr[i][1] = vtr[i][2] = 0;
     //vtr[i][0] = vmom[i][0];
     //vtr[i][1] = vmom[i][1];
     //vtr[i][2] = vmom[i][2];
    }
  }
}


/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixMesoAdami2013::initial_integrate(int /*vflag*/) {
  // update v and x and rho and e of atoms in group

  double **x = atom->x;
  double **vmom = atom->vest; //动量速度，更新动量
  double **f = atom->f;
  double **fpb = atom->fpb; // background pressure fpb is used to updata vtr(更新粒子的移动速度)
  double **vtr = atom->v; // tansport velocity
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;
  double dtfm;
  double imass;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }
      e[i] += dtf * de[i]; // half-step update of particle internal energy
      rho[i] += dtf * drho[i]; // ... and density
      
      vmom[i][0] += dtfm * f[i][0];
      vmom[i][1] += dtfm * f[i][1];
      vmom[i][2] += dtfm * f[i][2];
      
      vtr[i][0] = vmom[i][0] + dtfm * fpb[i][0];
      vtr[i][1] = vmom[i][1] + dtfm * fpb[i][1];
      vtr[i][2] = vmom[i][2] + dtfm * fpb[i][2];

      x[i][0] += dtv * vtr[i][0];
      x[i][1] += dtv * vtr[i][1];
      x[i][2] += dtv * vtr[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoAdami2013::final_integrate() {

  // update v, rho, and e of atoms in group

  double **vmom = atom->vest;
  double **f = atom->f;
  double **fpb = atom->fpb;
  double **vtr = atom->v;
  double *e = atom->e;
  double *de = atom->de;
  double *rho = atom->rho;
  double *drho = atom->drho;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  double dtfm;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }
      vmom[i][0] += dtfm * f[i][0];
      vmom[i][1] += dtfm * f[i][1];
      vmom[i][2] += dtfm * f[i][2];
      
    // vtr[i][0] = vmom[i][0] + dtfm * fpb[i][0];
    // vtr[i][1] = vmom[i][1] + dtfm * fpb[i][1];
     //vtr[i][2] = vmom[i][2] + dtfm * fpb[i][2];

      e[i] += dtf * de[i];
      rho[i] += dtf * drho[i];

    }
  }
}

/* ------------------------------------------------------------------\---- */

void FixMesoAdami2013::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
