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

#include "pair_wkl_sph_normvec.h"
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

PairWklSPHNormVec::PairWklSPHNormVec(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair
  comm_forward = 5; //surfnorm(3) + surfmagni(1) + surfmark[1]
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairWklSPHNormVec::~PairWklSPHNormVec() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairWklSPHNormVec::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHNormVec::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, h, ih, ihsq;
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
  double *color = atom->radius;
  double *surfmagni = atom->q;

  double sum1, numb1;
  sum1 = 0;
  numb1 = 0;

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
        surfnorm[i][0] = 0;
        surfnorm[i][1] = 0;
        surfnorm[i][2] = 0;
      }

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        if (surfmark[i][1] > 2) {
          xtmp = x[i][0];
          ytmp = x[i][1];
          ztmp = x[i][2];
          itype = type[i];
          jlist = firstneigh[i];
          jnum = numneigh[i];

          sigma = 0;


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
            j &= NEIGHMASK;
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx * delx + dely * dely + delz * delz;
            jtype = type[j];
            // jmass = mass[jtype];

            ////// quintic kernal //////  support domain is 3h
            if (rsq < cutsq[itype][jtype]) {
              h = 0.3333333333333333333 * cut[itype][jtype];
              ih = 1.0 / h;
              ihsq = ih * ih;

              r0 = sqrt(rsq) / h; // normalize r0
              wf = (3.0 - r0) * (3.0 - r0);
              wfd = - wf * wf;
              wf = - wfd * (3.0 - r0);

              if (r0 < 2.0){
                rr1 = 2.0 - r0;
                rr2 = 6.0 * rr1 * rr1 * rr1 * rr1;
                wf -= rr2 * rr1;
                wfd += rr2;
                if (r0 < 1.0){
                  rr1 = 1.0 - r0;
                  rr2 = 15.0 * rr1 * rr1 * rr1 * rr1;
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

              if (r0 > 1e-32) {
              //tmp = jmass / rho[j] * (color2 - color1) * wfd / r0;
              //if use rhosum, the above formula is not right
              tmp = mass[jtype] / rho[j] * (color[j] - color[i]) * wfd / r0;
              //tmp = mass[jtype] / rho[j] * color[j] * wfd / r0;
             // tmp = (jtype - itype) * wfd / r0;
              surfnorm[i][0] += tmp * delx;
              surfnorm[i][1] += tmp * dely;
              surfnorm[i][2] += tmp * delz;

              sigma += mass[jtype] / rho[j] * wf;

              }


            } // jrsq

          } // jj

          surfmagni[i] = sqrt(surfnorm[i][0] * surfnorm[i][0] + surfnorm[i][1] * surfnorm[i][1] + surfnorm[i][2] * surfnorm[i][2]);
	  if (surfmagni[i] > 0.01 * ih){  // 后面优化的时候这里要改成 0.01 * ih
            surfnorm[i][0] /= surfmagni[i];
            surfnorm[i][1] /= surfmagni[i];
            surfnorm[i][2] /= surfmagni[i];

	    sum1 += surfmagni[i];
	    numb1 += 1;

	   // surfmagni[i] /= sigma;

          }else{
            surfnorm[i][0] = 0;
            surfnorm[i][1] = 0;
            surfnorm[i][2] = 0;

            surfmark[i][1] = 0;
          }

	 // sum1 += surfmagni[i];
         // numb1 += 1;

        } // imark

      } // ii

//      printf("%s %f \n", "sum = ", sum1);
//      printf("%s %f \n", "num = ", numb1);
//      sum1 /= numb1;
//      printf("%s %f \n", "ave = ", sum1);

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
       // e[i] = surfmagni[i];
      }

//    }
//  }

  // communicate surfnorm
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairWklSPHNormVec::allocate() {
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

void PairWklSPHNormVec::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style wkl/sph/nromvec");
  //nstep = force->inumeric(FLERR,arg[0]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairWklSPHNormVec::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for wkl/sph/nromvec coefficients");
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

double PairWklSPHNormVec::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair wkl/sph/nromvec coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairWklSPHNormVec::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairWklSPHNormVec::pack_forward_comm(int n, int *list, double *buf,
                                     int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  double **surfnorm = atom->torque;
  double *surfmagni = atom->q;
  double **surfmark = atom->omega;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = surfnorm[j][0];
    buf[m++] = surfnorm[j][1];
    buf[m++] = surfnorm[j][2]; 
    buf[m++] = surfmagni[j];
    buf[m++] = surfmark[j][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairWklSPHNormVec::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double **surfnorm = atom->torque;
  double *surfmagni = atom->q;
  double **surfmark = atom->omega;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    surfnorm[i][0] = buf[m++];
    surfnorm[i][1] = buf[m++];
    surfnorm[i][2] = buf[m++];
    surfmagni[i] = buf[m++];
    surfmark[i][1] = buf[m++];
  }
}
