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

#include <math.h>
#include <stdlib.h>
#include "pair_sph_taitwater_morris_noslip_flat.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "group.h"
#include "math_vector2.h"
#include <iostream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHTaitwaterMorrisNoslipFlat::PairSPHTaitwaterMorrisNoslipFlat(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHTaitwaterMorrisNoslipFlat::~PairSPHTaitwaterMorrisNoslipFlat() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
    memory->destroy(rigid_group);
    memory->destroy(zflat);
    memory->destroy(beta);
  }
}

/* ---------------------------------------------------------------------- */

inline vector<double, 3> noslip(const double z, const vector<double, 3> &xb, const vector<double, 3> &xf,
                                const vector<double, 3> &vb, const vector<double, 3> &vf, const double beta) {
  const double df = xf[2] - z;
  const double db = xb[2] - z;
  double scale =  - db / df;
  if (scale < 0)
    return vb;
  if (scale > beta) scale = beta;
  return - scale * vf;
}

void PairSPHTaitwaterMorrisNoslipFlat::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, velx, vely, velz;
  double rsq, tmp, wfd, delVdotDelR, deltaE;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

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
  int *mask = atom->mask;

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

    // constant density for the boundary
    if (itype > 2) rho[i] = rho0[itype];

    // compute pressure of atom i with Tait EOS
    tmp = rho[i] / rho0[itype];
    fi = tmp * tmp * tmp;
    fi = B[itype] * (fi * fi * tmp - 1.0) / (rho[i] * rho[i]);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      // constant density for the boundary
      if (jtype > 2) rho[j] = rho0[jtype];

      if (rsq < cutsq[itype][jtype]) {
        // no slip

        vector<double, 3> vi = {vxtmp, vytmp, vztmp};
        vector<double, 3> vj = {v[j][0], v[j][1], v[j][2]};
        const int groupid = rigid_group[itype][jtype];
        const int groupbit = group->bitmask[groupid];
        if ((mask[i] & groupbit) && !(mask[j] & groupbit))
          vi = noslip(zflat[itype][jtype], {xtmp, ytmp, ztmp}, {x[j][0], x[j][1], x[j][2]},
                      vi, vj, beta[itype][jtype]);
        else if (!(mask[i] & groupbit) && (mask[j] & groupbit))
          vj = noslip(zflat[itype][jtype], {x[j][0], x[j][1], x[j][2]}, {xtmp, ytmp, ztmp},
                      vj, vi, beta[itype][jtype]);
        const vector<double, 3> vel = vi - vj;

        // kernel

        const double r = sqrt(rsq);
        wfd = kernel[itype][jtype].kernel_gradient(r) / r;

        // compute pressure  of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        fj = tmp * tmp * tmp;
        fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);

        // Neumann boundary condition for pressure
        if (itype > 2)
          fi = fj * (rho[j] * rho[j]) / (rho[i] * rho[i]);
        else if (jtype > 2)
          fj = fi * (rho[i] * rho[i]) / (rho[j] * rho[j]);

        velx = vel[0];
        vely = vel[1];
        velz = vel[2];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        // Morris Viscosity (Morris, 1996)

        /*
        fvisc = 2 * viscosity[itype][jtype] / (rho[i] * rho[j]);

        fvisc *= imass * jmass * wfd;
        */
        fvisc = imass * jmass / (rho[i] * rho[j]) * wfd * viscosity[itype][jtype] * 5.0/3.0;

        // total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj) * wfd;
        fpair += fvisc * delVdotDelR / rsq;
        deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely*vely + velz*velz));

       // printf("testvar= %f, %f \n", delx, dely);

        f[i][0] += delx * fpair + velx * fvisc;
        f[i][1] += dely * fpair + vely * fvisc;
        f[i][2] += delz * fpair + velz * fvisc;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        de[i] += deltaE;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair + velx * fvisc;
          f[j][1] -= dely * fpair + vely * fvisc;
          f[j][2] -= delz * fpair + velz * fvisc;
          de[j] += deltaE;
          drho[j] += imass * delVdotDelR * wfd;
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

void PairSPHTaitwaterMorrisNoslipFlat::allocate() {
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
  memory->create(rigid_group, n + 1, n + 1, "pair:rigid_group");
  memory->create(zflat, n + 1, n + 1, "pair:zflat");
  memory->create(beta, n + 1, n + 1, "pair:beta");

  kernel.resize(n + 1);
  for (int i = 1; i <= n; ++i)
    kernel[i].resize(n + 1);
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHTaitwaterMorrisNoslipFlat::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/taitwater/morris");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHTaitwaterMorrisNoslipFlat::coeff(int narg, char **arg) {
  if (narg != 10)
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
  int kernel_one = force->inumeric(FLERR,arg[6]);
  int rigid_group_one = group->find(arg[7]);
  double zflat_one = force->numeric(FLERR,arg[8]);
  double beta_one = force->numeric(FLERR,arg[9]);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      rho0[j] = rho0_one;
      soundspeed[j] = soundspeed_one;
      B[j] = B_one;

      viscosity[i][j] = viscosity_one;
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      kernel[i][j].init(domain->dimension, kernel_one, cut_one);
      rigid_group[i][j] = rigid_group_one;
      zflat[i][j] = zflat_one;
      beta[i][j] = beta_one;

      setflag[i][j] = 1;
      setflag[i][i] = 1;
      setflag[j][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHTaitwaterMorrisNoslipFlat::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/taitwater/morris coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];
  kernel[j][i].init(domain->dimension, kernel[i][j].get_kernel_type(), cut[j][i]);
  rigid_group[j][i] = rigid_group[i][j];
  zflat[j][i] = zflat[i][j];
  beta[j][i] = beta[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHTaitwaterMorrisNoslipFlat::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
