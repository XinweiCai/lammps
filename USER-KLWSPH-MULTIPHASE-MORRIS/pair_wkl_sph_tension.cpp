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

#include "pair_wkl_sph_tension.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairWklSPHTension::PairWklSPHTension(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair
  comm_forward = 1; // curvature test
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairWklSPHTension::~PairWklSPHTension() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairWklSPHTension::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHTension::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double sigma, tmp;
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double *e = atom->e;
  double *rho = atom->rho;
  int *type = atom->type;
  double *mass = atom->mass;
  double **f = atom->f;

  // use atom_vec_line 
  // surfmark indicates whether a calculation need to to be performed at this particle;
  // surfmark[i][0] for curvature; [1] for normal vector; [2] for color smooth; 
  double **surfmark = atom->omega;
  double **surfnorm = atom->torque;
  double *surfmagni = atom->q;
  double *curvature = atom->radius;

  double numb1, sum1;
  numb1 = 0;
  sum1 = 0;

  // check consistency of pair coefficients

  inum = list->inum;
  ilist = list->ilist;

  // recompute density
  // we use a full neighborlist here

//  if (nstep != 0) {
//    if ((update->ntimestep % nstep) == 0) {
    int nnnn = 0;  
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
	itype = type[i];
	nnnn += 1;
        if (surfmark[i][0] > 2 && surfmark[i][1] > 2) {
          tmp = tencoeff[itype] * mass[itype] / rho[i] * curvature[i] * surfmagni[i];
          f[i][0] -= tmp * surfnorm[i][0];
          f[i][1] -= tmp * surfnorm[i][1];
          f[i][2] -= tmp * surfnorm[i][2];

//	  e[i] = tmp;

//	  numb1 = surfnorm[i][0] * surfnorm[i][0] + surfnorm[i][1] * surfnorm[i][1] + surfnorm[i][2] * surfnorm[i][2];
//	  printf("%s %f \n", "coeff = ", tencoeff[itype]);
//          printf("%s %f \n", "mass = ", mass[itype]);
//          printf("%s %f \n", "rho = ", rho[i]);
//          printf("%s %f \n", "curv = ", numb1);

//	  numb1 += 1;
//	  sum1 += surfmagni[i];
        } 

      } 

//      printf("%s %f \n", "sum = ", sum1);
//      printf("%s %f \n", "num = ", numb1);
//      sum1 /= numb1;
//      printf("%s %f \n", "ave = ", sum1);


//    }
//  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    //e[i] = curvature[i];
  }

  // communicate curvature
  comm->forward_comm_pair(this);

}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairWklSPHTension::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");

  memory->create(tencoeff, n + 1, "pair:tencoeff");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairWklSPHTension::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style wkl/sph/tension");
  //nstep = force->inumeric(FLERR,arg[0]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairWklSPHTension::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for wkl/sph/tension coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR,arg[2]);
  double tencoeff_one = force->numeric(FLERR,arg[3]);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    tencoeff[i] = tencoeff_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairWklSPHTension::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair wkl/sph/tension coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairWklSPHTension::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

int PairWklSPHTension::pack_forward_comm(int n, int *list, double *buf,
                                     int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  double *curvature = atom->radius;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = curvature[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHTension::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *curvature = atom->radius;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    curvature[i] = buf[m++];
  }

}

