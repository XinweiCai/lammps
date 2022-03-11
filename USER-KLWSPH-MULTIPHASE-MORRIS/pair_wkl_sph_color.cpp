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

#include "pair_wkl_sph_color.h"
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

PairWklSPHColor::PairWklSPHColor(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair
  comm_forward = 1; //color(1)
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairWklSPHColor::~PairWklSPHColor() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairWklSPHColor::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHColor::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, h, ih, ihsq;
  int *jlist;
  double wf;
  double r0, rr1, rr2, sigma;
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double *rho = atom->rho;
  double *e = atom->e;
  int *type = atom->type;
  double *mass = atom->mass;

  // use atom_vec_line 
  // surfmark indicates whether a calculation need to to be performed at this particle;
  // surfmark[i][0] for curvature; [1] for normal vector; [2] for color smooth; 
  // ato->radius is for both color and curvature
  double **surfmark = atom->omega; 
  double *color = atom->radius; // 

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // smooth the color index, use atom types
  
//  if (nstep != 0) {
//    if ((update->ntimestep % nstep) == 0) {
    
      // initialize color with self-contribution,
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        itype = type[i];

	if (surfmark[i][2] > 2) {
          // quintic kernal
          h = 0.333333333333333333 * cut[itype][itype];
          if (domain->dimension == 3) {
            // 3d
            wf = 0.17555809878660322 / (h * h *h);
          }else{
            // 2d
            wf = 0.3076551619600655 / (h * h);
          }
  
          color[i] = mass[itype] / rho[i] * itype * wf;
  	
	} else {
        
	color[i] = itype;

	}

      }

      // add color by atom type at each atom via kernel function overlap
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
	if (surfmark[i][2] > 2) {
          xtmp = x[i][0];
          ytmp = x[i][1];
          ztmp = x[i][2];
          itype = type[i];
          jlist = firstneigh[i];
          jnum = numneigh[i];
  
          sigma = 0;
  
          for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
  
            jtype = type[j];
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx * delx + dely * dely + delz * delz;
  
            if (rsq < cutsq[itype][jtype]) {
              h = 0.333333333333333333 * cut[itype][jtype];
              ih = 1.0 / h;
              ihsq = ih * ih;
  
              r0 = sqrt(rsq) / h;
              wf = (3 - r0) * (3 - r0);
              wf = wf * wf * (3-r0);
  
              if (r0 < 2){
                rr1 = 2 - r0;
                rr2 = 6 * rr1 * rr1 * rr1 * rr1;
                wf -= rr2 * rr1;
                if (r0 < 1){
                  rr1 = 1 - r0;
                  rr2 = 15 * rr1 * rr1 * rr1 * rr1;
                  wf += rr2 * rr1;
                }
              }
              if (domain->dimension == 3) {
                // 3d
                wf *= 0.00265997119373641245 * ihsq * ih;
              }else{
                // 2d
                wf *= 0.00466144184787977995 * ihsq;
              }
              color[i] += mass[jtype] / rho[j] * jtype * wf;
              sigma += mass[jtype] / rho[j] *  wf;
            }
  
          }
  
          // if other kernal is used, here need to be changed
          h = 0.5 * 0.3333333333333333333 * cut[itype][itype];
          if (domain->dimension == 3) {
            // 3d
            wf = 0.17555809878660322 / (h * h *h);
          }else{
            // 2d
            wf = 0.3076551619600655 / (h * h);
          }
          sigma += mass[itype] / rho[i] *  wf;
          //color[i] /= sigma;
	}
      }

      //for (ii = 0; ii < inum; ii++) {
        //i = ilist[ii];
       // e[i] = color[i];
      //}


//    }
//  }

  // communicate color
  comm->forward_comm_pair(this);

}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairWklSPHColor::allocate() {
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

void PairWklSPHColor::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style wkl/sph/color");
  //nstep = force->inumeric(FLERR,arg[0]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairWklSPHColor::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for wkl/sph/color coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR,arg[2]);
  
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

double PairWklSPHColor::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair wkl/sph/color coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairWklSPHColor::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairWklSPHColor::pack_forward_comm(int n, int *list, double *buf,
                                     int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  double *color = atom->radius;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = color[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHColor::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *color = atom->radius;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    color[i] = buf[m++];
  }
}
