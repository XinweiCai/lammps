 /* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 Reference:Adami2012
 Date:2022-04-01
 by xinwei
 ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "pair_sph_freesurface.h"
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
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHFreeSurface::PairSPHFreeSurface(LAMMPS *lmp) : Pair(lmp)
{
	restartinfo = 0;
	first = 1;
}



PairSPHFreeSurface::~PairSPHFreeSurface() {
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutsq);
		memory->destroy(cut);
		memory->destroy(soundspeed);
		memory->destroy(B);
		memory->destroy(viscosity);
		memory->destroy(rho0);
	}
}

/* ---------------------------------------------------------------------- */
// get the quintic spline gradient
inline double getQuinticGradient(double r, double cutoff, int dimension)
{
	double h, ih, ihsq, s, w_g;
	double s1, s2, s3, s1sq, s2sq, s3sq;
	double s1_4, s2_4, s3_4;

	h = cutoff / 3;
	ih = 1.0 / h;
	ihsq = ih * ih;

	s  = r/h;
	s1 = 1.0 - s;
	s1sq = s1 * s1;
	s2 = 2.0 - s;
	s2sq = s2 * s2;
	s3 = 3.0 - s;
	s3sq = s3 * s3;
	s1_4 = -75.0  * s1sq * s1sq;
	s2_4 = 30.0 * s2sq * s2sq;
	s3_4 = -5.0 *   s3sq * s3sq;

	if ( s < 1.0 ) w_g = s3_4 + s2_4 + s1_4;
	else if ( s < 2.0 ) w_g = s3_4 + s2_4;
	else if ( s < 3.0 ) w_g = s3_4;
	else w_g = 0.0;

	if (dimension == 3) {
		w_g *= 0.002652582384864922 * ihsq * ihsq;	//w_g *= 1.0/120.0/PI* ih * ih * ih * ih
	}
	else if (dimension == 2) {
		w_g *= 0.00466144184787978 *  ihsq * ih;	//w_g *= 7.0/478.0/PI *  ih * ih * ih
	}
	else {
		w_g *= 0.008333333333333333 * ihsq;	//w_g *= 1.0/120.0 *  ih * ih
	}

	return w_g;
}

/* ----------------------------compute momentum------------------------------------------ */

void PairSPHFreeSurface::compute(int eflag, int vflag) {
	int i, j, ii, jj, inum, jnum, numf, itype, jtype;
	double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

	int *ilist, *jlist, *numneigh, **firstneigh;
	double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, pi, pj, pij, fvisc, velx, vely, velz;
	double h, hsq, rsq, tmp, wfd, delVdotDelR, deltaE;

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
		// compute pressure of atom i with Tait EOS
		tmp = rho[i] / rho0[itype];
		pi = tmp * tmp * tmp;
		pi = B[itype] * (pi * pi * tmp - 1.0);
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
				vector<double, 3> vi = { vxtmp, vytmp, vztmp };
				vector<double, 3> vj = { v[j][0], v[j][1], v[j][2] };
				const vector<double, 3> vel = vi - vj;
				velx = vel[0];
				vely = vel[1];
				velz = vel[2];
				// dot product of velocity delta and distance vector
				delVdotDelR = delx * velx + dely * vely + delz * velz;

				// compute pressure of atom j with Tait EOS
				tmp = rho[j] / rho0[jtype];
				pj = tmp * tmp * tmp;
				pj = B[jtype] * (pj * pj * tmp - 1.0);

				const double r = sqrt(rsq);
				wfd = getQuinticGradient(r, cut[itype][jtype], domain->dimension) / r;

				// using the density-weighted inter-particle averaged pressure(X.Y. Hu,J.C.P,2007)
				pij = (rho[j] * pi + rho[i] * pj) / (rho[i] + rho[j]);

				// the acceleration of particle a caused by shear forces(adami2012)
				//artifical viscosity 
				h = cut[itype][jtype] / 3;
				hsq = h * h;
				fvisc = imass * jmass * viscosity[itype][jtype] * h * (soundspeed[itype] + soundspeed[jtype]) * delVdotDelR * wfd;
				fvisc /= (rho[i] + rho[j]) * (rsq + 0.01 * hsq);
				
				// total pair force 
				fpair = -((imass / rho[i])*(imass / rho[i]) + (jmass / rho[j])*(jmass / rho[j])) * pij * wfd;
				fpair += fvisc;
				// thermal energy increment
				deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely * vely + velz * velz));

				f[i][0] += delx * fpair;
				f[i][1] += dely * fpair;
				f[i][2] += delz * fpair;
			
				// and change in density
				drho[i] += rho[i] * (jmass / rho[j]) * delVdotDelR * wfd;
				// change in thermal energy
				de[i] += deltaE;

				if (newton_pair || j < nlocal) {
					f[j][0] -= delx * fpair;
					f[j][1] -= dely * fpair;
					f[j][2] -= delz * fpair;
					de[j] += deltaE;
					drho[j] += rho[j] * (imass / rho[i]) * delVdotDelR * wfd;
					
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

void PairSPHFreeSurface::allocate() {
	allocated = 1;
	int n = atom->ntypes;

	memory->create(setflag, n + 1, n + 1, "pair:setflag");
	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			setflag[i][j] = 0;

	memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
	memory->create(soundspeed, n + 1, "pair:soundspeed");
	memory->create(B, n + 1, "pair:B");
	memory->create(cut, n + 1, n + 1, "pair:cut");
	memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
	memory->create(rho0, n + 1, "pair:rho");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHFreeSurface::settings(int narg, char **arg) {
	if (narg != 0)
		error->all(FLERR,
			"Illegal number of setting arguments for pair_style sph/adami/noslip/freesurface");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHFreeSurface::coeff(int narg, char **arg) {
	if (narg != 6)
		error->all(FLERR,
			"Incorrect args for pair_style sph/adami/noslip/freesurface coefficients");
	if (!allocated)
		allocate();

	int ilo, ihi, jlo, jhi;
	force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
	force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

	double rho0_one = force->numeric(FLERR, arg[2]);
	double soundspeed_one = force->numeric(FLERR, arg[3]);
	double viscosity_one = force->numeric(FLERR, arg[4]);
	double cut_one = force->numeric(FLERR, arg[5]);
	double B_one = soundspeed_one * soundspeed_one * rho0_one / 7;

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		rho0[i] = rho0_one;
		soundspeed[i] = soundspeed_one;
		B[i] = B_one;
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			rho0[j] = rho0_one;
			soundspeed[j] = soundspeed_one;
			viscosity[i][j] = viscosity_one;
			//printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
			cut[i][j] = cut_one;

			setflag[i][j] = 1;
			setflag[i][i] = 1;
			setflag[j][j] = 1;
			count++;
		}
	}

	if (count == 0)
		error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHFreeSurface::init_one(int i, int j) {
	if (setflag[i][j] == 0) {
		error->all(FLERR, "Not all pair coeffs are not set");
	}
	cut[j][i] = cut[i][j];
	viscosity[j][i] = viscosity[i][j];

	return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHFreeSurface::single(int i, int j, int itype, int jtype,
	double rsq, double factor_coul, double factor_lj, double &fforce) {
	fforce = 0.0;

	return 0.0;
}

