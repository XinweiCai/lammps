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

#include "pair_wkl_sph_rhosum_taitwater_hu.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairWklSPHRhosumTaitwaterHu::PairWklSPHRhosumTaitwaterHu(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  first = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairWklSPHRhosumTaitwaterHu::~PairWklSPHRhosumTaitwaterHu() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
  }
}

/* ---------------------------------------------------------------------- */

void PairWklSPHRhosumTaitwaterHu::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, h, ih, ihsq, velx, vely, velz;
  double rsq, tmp, wfd, delVdotDelR, deltaE;

  double wf, r0, rr1, rr2;

  double viscij;
  double sigmai, sigmaj;


  ev_init(eflag, vflag);

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *de = atom->de;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 1.e-32) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
                  i, j, sqrt(cutsq[i][j]));
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    sigmai = rho[i] / imass;

    // compute pressure of atom i with Tait EOS
    //tmp = rho[i] / rho0[itype];
    //fi = tmp * tmp * tmp;
    //fi = B[itype] * (fi * fi * tmp - 1.0) / (rho[i] * rho[i]);

    //fi = soundspeed[itype] * soundspeed[itype] * rho[i];
    fi = soundspeed[itype] * soundspeed[itype] * (rho[i] - rho0[itype]); 


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        
        jmass = mass[jtype];
	sigmaj = rho[j] / jmass;

        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;


	////// quintic kernal //////  support domain is 3h
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

        r0 *= h; //give the real r0
	wfd /= r0; // then wfd misses a factor of r



//        wfd = h - sqrt(rsq);
//        if (domain->dimension == 3) {
//          // Lucy Kernel, 3d
//          // Note that wfd, the derivative of the weight function with respect to r,
//          // is lacking a factor of r.
//          // The missing factor of r is recovered by
//          // (1) using delV . delX instead of delV . (delX/r) and
//          // (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
//          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
//        } else {
//          // Lucy Kernel, 2d
//          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
//        }

        //// compute pressure  of atom j with Tait EOS
        //tmp = rho[j] / rho0[jtype];
        //fj = tmp * tmp * tmp;
        //fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);
	
	//fj = soundspeed[jtype] * soundspeed[jtype] * rho[j]; 
	fj = soundspeed[jtype] * soundspeed[jtype] * (rho[j] - rho0[jtype]);
	
        velx = vxtmp - v[j][0];
        vely = vytmp - v[j][1];
        velz = vztmp - v[j][2];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        // Morris Viscosity (Morris, 1996)
//        fvisc = 2 * viscosity[itype][jtype] / (rho[i] * rho[j]);
//        fvisc *= imass * jmass * wfd;
        
        if (itype == jtype){
          viscij = viscosity[itype][itype];
        } else {
          if (viscosity[itype][jtype] < 1E-15) {
            viscij = 0;
	  } else {
            viscij = 2 * viscosity[itype][itype] * viscosity[jtype][jtype] / (viscosity[itype][itype] + viscosity[jtype][jtype]);
	  }
        }

	fvisc = viscij * (1.0 / (sigmai * sigmai) + 1.0 /(sigmaj * sigmaj)) * wfd;
	//fvisc = viscosity[itype][jtype] * (1.0 / (sigmai * sigmai) + 1.0 /(sigmaj * sigmaj)) * wfd;


        // total pair force & thermal energy increment
//        fpair = -imass * jmass * (fi + fj) * wfd / (rho[i] * rho[j]);
        fpair = - (fi / (sigmai * sigmai) + fj /(sigmaj * sigmaj)) * wfd;


        deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely*vely + velz*velz));

       // printf("testvar= %f, %f \n", delx, dely);

        f[i][0] += delx * fpair + velx * fvisc;
        f[i][1] += dely * fpair + vely * fvisc;
        f[i][2] += delz * fpair + velz * fvisc;

        // and change in density
        //drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        de[i] += deltaE;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair + velx * fvisc;
          f[j][1] -= dely * fpair + vely * fvisc;
          f[j][2] -= delz * fpair + velz * fvisc;
          de[j] += deltaE;
         // drho[j] += imass * delVdotDelR * wfd;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairWklSPHRhosumTaitwaterHu::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairWklSPHRhosumTaitwaterHu::settings(int narg, char **/*arg*/) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style sph/taitwater/morris");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairWklSPHRhosumTaitwaterHu::coeff(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR,
        "Incorrect args for pair_style sph/taitwater/morris coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double rho0_one = force->numeric(FLERR,arg[2]);
  double soundspeed_one = force->numeric(FLERR,arg[3]);
  double viscosity_one = force->numeric(FLERR,arg[4]);
  double cut_one = force->numeric(FLERR,arg[5]);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
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

double PairWklSPHRhosumTaitwaterHu::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/taitwater/morris coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];

  return cut[i][j];
}

