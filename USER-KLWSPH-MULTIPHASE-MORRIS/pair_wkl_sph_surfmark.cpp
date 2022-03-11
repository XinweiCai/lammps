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

#include "pair_wkl_sph_surfmark.h"
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

PairWklSPHSurfMark::PairWklSPHSurfMark(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair
  comm_forward = 3; //surfmak(3)
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairWklSPHSurfMark::~PairWklSPHSurfMark() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairWklSPHSurfMark::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHSurfMark::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, h, ih, ihsq;
  double lenzerosq, lenonesq, lentwosq;
  int *jlist;
  double wf;
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double *e = atom->e;
  int *type = atom->type;
  double *mass = atom->mass;

  // use atom_vec_line 
  // surfmark indicates whether a calculation need to to be performed at this particle;
  // surfmark[i][0] for curvature; [1] for normal vector; [2] for color smooth; 
  double **surfmark = atom->omega; 

  // check consistency of pair coefficients

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

//  if (nstep != 0) {
//    if ((update->ntimestep % nstep) == 0) {
      
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
	surfmark[i][0] = 0;
	surfmark[i][1] = 0;
	surfmark[i][2] = 0;
      }
      
      lenzerosq = lenzero * lenzero;
      lenonesq = lenone * lenone;
      lentwosq = lentwo * lentwo;
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

	for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;

          jtype = type[j];
	  if (itype != jtype) {

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;  
          if (rsq < lentwosq) {
            surfmark[i][2] = 4;
            if (rsq < lenonesq) {
	      surfmark[i][1] = 6;
              if (rsq < lenzerosq) {
	        surfmark[i][0] = 8;
	      }
	    }
	  }

	  }
	}
      }

      //for (ii = 0; ii < inum; ii++) {
        //i = ilist[ii];
       // e[i] = surfmark[i][0] + surfmark[i][1] + surfmark[i][2];
      //}


//    }
//  }

  // communicate surfmark
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairWklSPHSurfMark::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairWklSPHSurfMark::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style sph/rhosum");
 // nstep = force->inumeric(FLERR,arg[0]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairWklSPHSurfMark::coeff(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR,"Incorrect number of args for sph/rhosum coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR,arg[2]);
  
  lenzero = force->numeric(FLERR,arg[3]);
  lenone = force->numeric(FLERR,arg[4]);
  lentwo = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
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

double PairWklSPHSurfMark::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/rhosum coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairWklSPHSurfMark::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairWklSPHSurfMark::pack_forward_comm(int n, int *list, double *buf,
                                     int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  double **surfmark = atom->omega;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = surfmark[j][0];
    buf[m++] = surfmark[j][1];
    buf[m++] = surfmark[j][2]; 
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHSurfMark::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double **surfmark = atom->omega;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    surfmark[i][0] = buf[m++];
    surfmark[i][1] = buf[m++];
    surfmark[i][2] = buf[m++];
}
