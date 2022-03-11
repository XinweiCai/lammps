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

#include "pair_wkl_sph_curvature.h"
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

PairWklSPHCurvature::PairWklSPHCurvature(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair
  comm_forward = 1; //curvature(1)
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairWklSPHCurvature::~PairWklSPHCurvature() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairWklSPHCurvature::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHCurvature::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, jmass, h, ih, ihsq;
  double lenzerosq, lenonesq, lentwosq;
  int *jlist;
  double wf, wfd;
  double r0, rr1, rr2, sigma, tmp;
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double *e = atom->e;
  double *rho = atom->rho;
  int *type = atom->type;
  double *mass = atom->mass;

  // use atom_vec_line 
  // surfmark indicates whether a calculation need to to be performed at this particle;
  // surfmark[i][0] for curvature; [1] for normal vector; [2] for color smooth; 
  double **surfmark = atom->omega;
  double **surfnorm = atom->torque;
  double *surfmagni = atom->q;
  double *curvature = atom->radius;

  // check consistency of pair coefficients

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

//  if (nstep != 0) {
//    if ((update->ntimestep % nstep) == 0) {
      // reset interface normal vector
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        curvature[i] = 0;
      }

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        if (surfmark[i][0] > 2 && surfmark[i][1] > 2) {
          xtmp = x[i][0];
          ytmp = x[i][1];
          ztmp = x[i][2];
          itype = type[i];
          jlist = firstneigh[i];
          jnum = numneigh[i];
          
	  // 这里需要注意一下
	  h = cut[itype][itype];
        if (domain->dimension == 3) {
          /*
          // Lucy kernel, 3d
          wf = 2.0889086280811262819e0 / (h * h * h);
          */

          // quadric kernel, 3d
          wf = 2.1541870227086614782 / (h * h * h);
        } else {
          /*
          // Lucy kernel, 2d
          wf = 1.5915494309189533576e0 / (h * h);
          */

          // quadric kernel, 2d
          wf = 1.5915494309189533576e0 / (h * h);
        }
          //sigma = 0;
	  sigma = mass[itype] / rho[i] *  wf;

          for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
	    if (surfmark[j][1] > 2) {
              j &= NEIGHMASK;
              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx * delx + dely * dely + delz * delz;
              jtype = type[j];
              //jmass = mass[jtype];
  
              ////// quintic kernal //////  support domain is 3h
              if (rsq < cutsq[itype][jtype]) {
                h = 0.3333333333333333333 * cut[itype][jtype];
                ih = 1.0 / h;
                ihsq = ih * ih;
  
                r0 = sqrt(rsq) / h; // normalize r0
                wf = (3 - r0) * (3 - r0);
                wfd = - wf * wf;
                wf = - wfd * (3 - r0);
  
                if (r0 < 2){
                  rr1 = 2 - r0;
                  rr2 = 6 * rr1 * rr1 * rr1 * rr1;
                  wf -= rr2 * rr1;
                  wfd += rr2;
                  if (r0 < 1){
                    rr1 = 1 - r0;
                    rr2 = 15 * rr1 * rr1 * rr1 * rr1;
                    wf += rr2 * rr1;
                    wfd -= rr2;
                  }
                }
  	      if (domain->dimension == 3) {
                  // 3d
                  wf *= 0.00265997119373641245 * ihsq * ih;
                  wfd *= 0.0132998559686820627 * ihsq * ihsq;
                }else{
                  // 2d
                  wf *= 0.00466144184787977995 * ihsq;
                  wfd *= 0.02330720923939889888 * ihsq * ih;
                }
  
                r0 *= h; //give the inital r0


                if (r0 > 1e-32){ // 这个后面考虑一下需不需要，因为好像不包括本身
                  curvature[i] +=  mass[jtype] / rho[j] * ((surfnorm[j][0] - surfnorm[i][0]) * delx + (surfnorm[j][1]
                       - surfnorm[i][1]) * dely + (surfnorm[j][2] - surfnorm[i][2]) * delz) * wfd / r0 ;
                  //curvature[i] +=  mass[jtype] / rho[j] * (surfnorm[j][0] * delx + surfnorm[j][1]
                    //       * dely + surfnorm[j][2] * delz) * wfd / r0 ;
                  //if rhosum is used, the above formula is not right
                  //curvature[i] += ((surfnorm[j][0] - surfnorm[i][0]) * delx + (surfnorm[j][1]
                    //           - surfnorm[i][1]) * dely + (surfnorm[j][2] - surfnorm[i][2]) * delz) * wfd / r0 ;
                   sigma += mass[jtype] / rho[j] *  wf;
                }
       

              } // jrsq

	    } // jmark
          } // jj

	  curvature[i] /= sigma;

        } // imark

      } // ii

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        //e[i] = 0;
      }

//    }
//  }

  // communicate curvature
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairWklSPHCurvature::allocate() {
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

void PairWklSPHCurvature::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style wkl/sph/curvature");
  //nstep = force->inumeric(FLERR,arg[0]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairWklSPHCurvature::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for wkl/sph/curvature coefficients");
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

double PairWklSPHCurvature::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair wkl/sph/curvature coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairWklSPHCurvature::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairWklSPHCurvature::pack_forward_comm(int n, int *list, double *buf,
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

void PairWklSPHCurvature::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *curvature = atom->radius;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
    curvature[i] = buf[m++];
  }
    
}
