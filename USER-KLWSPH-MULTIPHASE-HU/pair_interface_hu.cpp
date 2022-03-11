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

#include "pair_interface_hu.h"
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

PairInterfaceHu::PairInterfaceHu(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair
  comm_forward = 5; //surfnorm(3) + surfmagni(1) + surfmark[1]
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairInterfaceHu::~PairInterfaceHu() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);

   // memory->destroy(tencoeff); // 这个不知道是干什么的
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairInterfaceHu::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairInterfaceHu::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype, i1, j1;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, jmass, h, ih, ihsq;
  int *jlist;
  double wf, wfd;
  double r0, rr1, rr2, tmp;
  
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  double sigmai, sigmaj, vsphi, vsphj;
  double gradc[3];
  double gcmagni;
  double surfstress[3][3];
  double tencoeffi, coeffcount;

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
  double *color = atom->radius;
  double *surfmagni = atom->q;
//  double **surfnorm = atom->february;
//  double *color = atom->monday;
//  double *surfmagni = atom->tuesday;

  // check consistency of pair coefficients

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

  // reset interface normal vector
  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      e[i] = 0;}

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];

      if (surfmark[i][0] > 2) {
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        imass = mass[itype];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        sigmai = rho[i] / imass;
//	vsphi = imass / rho[i]; //这个在文章中老换来换去，还是用体积来讲好

        for (i1 = 0; i1 < 3; i1++) {
          gradc[i1] = 0;
          for (j1 = 0; j1 < 3; j1++)
            surfstress[i1][j1] = 0;
        }
        tencoeffi = 0;
        coeffcount = 0;

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;

          if (itype != type[j]) {
            delx = x[j][0] - xtmp;
            dely = x[j][1] - ytmp;
            delz = x[j][2] - ztmp;
            rsq = delx * delx + dely * dely + delz * delz;
            jtype = type[j];
            jmass = mass[jtype];
            sigmaj = rho[j] / jmass;

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

              if (r0 > 1e-32) {
                tmp = sigmai / (sigmaj * sigmaj) * wfd / r0;
                gradc[0] += tmp * delx;
                gradc[1] += tmp * dely;
                gradc[2] += tmp * delz;
               // printf("%f \n", tencoeff[itype][jtype]);
                // mean value for coefficient if there are more than 2 types.
                tencoeffi += tencoeff[itype][jtype] * wf / sigmaj;
                coeffcount += wf / sigmaj;
              }

              //sigma += wf;

            } // jrsq

          } // jtype

        } // jj

        //// 这里有个隐患，文章中是将每对流体产生的应力张量先算好再求和的 
        gcmagni = sqrt(gradc[0] * gradc[0] + gradc[1] * gradc[1] + gradc[2] * gradc[2]);
        if (gcmagni > 0.00000001 * ih){
          for (i1 = 0; i1 < 3; i1++) {
            for (j1 = 0; j1 < 3; j1++) {
              surfstress[i1][j1] = - gradc[i1] * gradc[j1];
              if (i1 == j1) surfstress[i1][j1] += gcmagni * gcmagni / (domain->dimension);
              surfstress[i1][j1] *= tencoeffi / coeffcount / gcmagni;
            }
          }

        }


        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          jtype = type[j];

          delx = x[j][0] - xtmp;
          dely = x[j][1] - ytmp;
          delz = x[j][2] - ztmp ;
          rsq = delx * delx + dely * dely + delz * delz;
          jtype = type[j];
          jmass = mass[jtype];
          sigmaj = rho[j] / jmass;

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

            if (r0 > 1e-32) {
              tmp = wfd / r0 / (sigmai * sigmai);
      	for (j1 = 0; j1 < 3; j1++){
//	  if (itype == 1){	
            f[i][j1] -= tmp * (delx * surfstress[0][j1] + dely * surfstress[1][j1] + delz * surfstress[2][j1]);
//	  }
	}
      	  
      	for (j1 = 0; j1 < 3; j1++){
//	  if (jtype == 1){
      	    f[j][j1] += tmp * (delx * surfstress[0][j1] + dely * surfstress[1][j1] + delz * surfstress[2][j1]);
//	  }
	}
            }

          } // jrsq

        } // jj

         e[i] = sqrt(surfstress[0][0] * surfstress[0][0] + surfstress[1][1] * surfstress[1][1] + surfstress[2][2] * surfstress[2][2]);
      } // imark

    } // ii

  // communicate surfnorm
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairInterfaceHu::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");

  memory->create(tencoeff, n + 1, n + 1, "pair:tencoeff");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairInterfaceHu::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style wkl/sph/nromvec");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairInterfaceHu::coeff(int narg, char **arg) {
  if (narg != 5)
    error->all(FLERR,"Incorrect number of args for wkl/sph/nromvec coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR,arg[2]);

  double tencoeff_one = force->numeric(FLERR,arg[3]);
  double tencoeff_two = force->numeric(FLERR,arg[4]);

/*  if (ilo != jlo || ihi != jhi) {
    error->all(FLERR,"Now this pair-style can only be set for two specific atom type");
    printf("%s","No hurry! set them one by one!");
  }
  typefor_one = ilo;
  typefor_two = jlo;
  */

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;

      tencoeff[i][j] = tencoeff_one;
      tencoeff[j][i] = tencoeff_two;

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

double PairInterfaceHu::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair wkl/sph/nromvec coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

//double PairInterfaceHu::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
//    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
//  fforce = 0.0;
//
//  return 0.0;
//}
//
///* ---------------------------------------------------------------------- */
//
//int PairInterfaceHu::pack_forward_comm(int n, int *list, double *buf,
//                                     int /*pbc_flag*/, int * /*pbc*/) {
//  int i, j, m;
//  double **surfnorm = atom->torque;
//  double *surfmagni = atom->q;
//  double **surfmark = atom->omega;
//
//  m = 0;
//  for (i = 0; i < n; i++) {
//    j = list[i];
//    buf[m++] = surfnorm[j][0];
//    buf[m++] = surfnorm[j][1];
//    buf[m++] = surfnorm[j][2]; 
//    buf[m++] = surfmagni[j];
//    buf[m++] = surfmark[j][1];
//  }
//  return m;
//}
//
///* ---------------------------------------------------------------------- */
//
//void PairInterfaceHu::unpack_forward_comm(int n, int first, double *buf) {
//  int i, m, last;
//  double **surfnorm = atom->torque;
//  double *surfmagni = atom->color;
//  double **surfmark = atom->omega;
//  
//  m = 0;
//  last = first + n;
//  for (i = first; i < last; i++) {
//    surfnorm[i][0] = buf[m++];
//    surfnorm[i][1] = buf[m++];
//    surfnorm[i][2] = buf[m++];
//    surfmagni[i] = buf[m++];
//    surfmark[i][1] = buf[m++];
//  }
//}
