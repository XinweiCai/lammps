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
#include "pair_sph_taitwater_angular_noslip_sphere.h"
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

PairSPHTaitwaterAngularNoslipSphere::PairSPHTaitwaterAngularNoslipSphere(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHTaitwaterAngularNoslipSphere::~PairSPHTaitwaterAngularNoslipSphere() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
    memory->destroy(rigid_group);
    memory->destroy(rsphere);
    memory->destroy(beta);
  }
}

/* ---------------------------------------------------------------------- */

inline vector<double, 3> noslip(const vector<double, 3> &xc, const double R, const vector<double, 3> &xb, const vector<double, 3> &xf,
                                const vector<double, 3> &vc, const vector<double, 3> &omega, const vector<double, 3> &vb, const vector<double, 3> &vf, const double beta) {
  const vector<double, 3> dfc = xf - xc;
  const double fc = norm(dfc);
  const double df = fc - R;
  if (df < 0)
    return vb;
  const vector<double, 3> dbc = xb - xc;
  const double db = R - dot(dfc, dbc) / fc;
  const vector<double ,3> xs = xc + dfc / fc * R;
  const vector<double, 3> vs = vc + cross(omega, xs - xc);
  double scale = db / df;
  if (scale > beta) scale = beta;
  return - scale * (vf - vs) + vs;
}

void PairSPHTaitwaterAngularNoslipSphere::compute(int eflag, int vflag) {
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

  // no slip: update center of mass, v, omega

  std::map<int, vector<double, 3> > xcm;
  std::map<int, vector<double, 3> > vcm;
  std::map<int, vector<double, 3> > omega;
  for (std::map<int, double>::iterator it = masstotal.begin(); it != masstotal.end(); ++it) {
    vector<double, 3> x_one;
    group->xcm(it->first, it->second, x_one.data());
    xcm[it->first] = x_one;

    vector<double, 3> v_one;
    group->vcm(it->first, it->second, v_one.data());
    vcm[it->first] = v_one;

    double lmom[3];
    group->angmom(it->first, x_one.data(), lmom);
    double inertia[3][3];
    group->inertia(it->first, x_one.data(), inertia);
    vector<double, 3> omega_one;
    group->omega(lmom, inertia, omega_one.data());
    omega[it->first] = omega_one;
  }

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

      if (rsq < cutsq[itype][jtype]) {
        // no slip

        vector<double, 3> vi = {vxtmp, vytmp, vztmp};
        vector<double, 3> vj = {v[j][0], v[j][1], v[j][2]};
        const int groupid = rigid_group[itype][jtype];
        const int groupbit = group->bitmask[groupid];
        if ((mask[i] & groupbit) && !(mask[j] & groupbit))
          vi = noslip(xcm[groupid], rsphere[itype][jtype], {xtmp, ytmp, ztmp}, {x[j][0], x[j][1], x[j][2]},
                      vcm[groupid], omega[groupid], vi, vj, beta[itype][jtype]);
        else if (!(mask[i] & groupbit) && (mask[j] & groupbit))
          vj = noslip(xcm[groupid], rsphere[itype][jtype], {x[j][0], x[j][1], x[j][2]}, {xtmp, ytmp, ztmp},
                      vcm[groupid], omega[groupid], vj, vi, beta[itype][jtype]);
        const vector<double, 3> vel = vi - vj;

        // kernel

        const double r = sqrt(rsq);
        wfd = kernel[itype][jtype].kernel_gradient(r) / r;

        // compute pressure  of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        fj = tmp * tmp * tmp;
        fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);

        velx = vel[0];
        vely = vel[1];
        velz = vel[2];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        // Morris Viscosity (Morris, 1996)

        fvisc = 10 * viscosity[itype][jtype] / (rho[i] * rho[j]);

        fvisc *= imass * jmass * wfd * delVdotDelR / rsq;

        // total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj) * wfd;
        fpair += fvisc;
        deltaE = -0.5 * fpair * delVdotDelR;

       // printf("testvar= %f, %f \n", delx, dely);

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        de[i] += deltaE;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
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

void PairSPHTaitwaterAngularNoslipSphere::allocate() {
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
  memory->create(rsphere, n + 1, n + 1, "pair:rsphere");
  memory->create(beta, n + 1, n + 1, "pair:beta");

  kernel.resize(n + 1);
  for (int i = 1; i <= n; ++i)
    kernel[i].resize(n + 1);
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHTaitwaterAngularNoslipSphere::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/taitwater/morris");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHTaitwaterAngularNoslipSphere::coeff(int narg, char **arg) {
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
  double rsphere_one = force->numeric(FLERR,arg[8]);
  double beta_one = force->numeric(FLERR,arg[9]);
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
      kernel[i][j].init(domain->dimension, kernel_one, cut_one);
      rigid_group[i][j] = rigid_group_one;
      rsphere[i][j] = rsphere_one;
      beta[i][j] = beta_one;

      setflag[i][j] = 1;
      count++;
    }
  }
  if (masstotal.find(rigid_group_one) == masstotal.end())
    masstotal.insert(std::pair<int, double>(rigid_group_one, group->mass(rigid_group_one)));

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHTaitwaterAngularNoslipSphere::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/taitwater/morris coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];
  kernel[j][i].init(domain->dimension, kernel[i][j].get_kernel_type(), cut[j][i]);
  rigid_group[j][i] = rigid_group[i][j];
  rsphere[j][i] = rsphere[i][j];
  beta[j][i] = beta[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHTaitwaterAngularNoslipSphere::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
