 /* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 Reference:Adami et al. A transport-velocity formulation for smoothed particle hydrodynamics.JCP2013
 Date:2021-06
 by xinwei
 ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "pair_sph_adami2013.h"
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

PairSPHAdami2013::PairSPHAdami2013(LAMMPS *lmp) : Pair(lmp)
{
	restartinfo = 0;
	first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHAdami2013::~PairSPHAdami2013() {
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutsq);
		memory->destroy(cut);
		memory->destroy(soundspeed);
		memory->destroy(B);
		memory->destroy(viscosity);
		memory->destroy(rho0);
		//memory->destroy(dp);

	}
}

/* ---------------------------------------------------------------------- */

void PairSPHAdami2013::compute(int eflag, int vflag) {
	int i, j, ii, jj, inum, jnum, numf, itype, jtype;
	double xtmp, ytmp, ztmp, delx, dely, delz, fpair, fpairpb;

	int *ilist, *jlist, *numneigh, **firstneigh;
	double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, pi, pj, pij, fvisc, velx, vely, velz,vel1x,vel1y,vel1z, vtrxtmp, vtrytmp, vtrztmp, ax,ay,az;
	double rsq, tmp, wf, wfd, sumwf, delVdotDelR, deltaE;
	double fpairx, fpairy, fpairz, fpairpbx, fpairpby, fpairpbz;

	if (eflag || vflag)
		ev_setup(eflag, vflag);
	else
		evflag = vflag_fdotr = 0;

	double **vtr = atom->v; 
	double **vv = atom->vv;  // extrapolated velocity of the wall particle
	double **v = atom->vest; 
	double **x = atom->x;
	double **f = atom->f;
	double **fpb = atom->fpb;
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
			vtrxtmp = vtr[i][0];
			vtrytmp = vtr[i][1];
			vtrztmp = vtr[i][2];
			// compute pressure of atom i 
			//f[i][0] += imass * 0.00000015;
			pi = B[itype] * (rho[i] / rho0[itype] - 1.0);


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
					//const vector<double, 3> vel1 = vi - vj;
					vector<double, 3> vi1 = { vtrxtmp, vtrytmp, vtrztmp };
					vector<double, 3> vj1 = { vtr[j][0], vtr[j][1], vtr[j][2] };
					const vector<double, 3> vel1 = vi1 - vj1;

					vel1x = vel1[0];
					vel1y = vel1[1];
					vel1z = vel1[2];

					//dot product of velocity delta and distance vector
					delVdotDelR = delx * vel1x + dely * vel1y + delz * vel1z;


					// no slip
					
					if ((mask[j] & groupbit) && !(mask[i] & groupbit)) {
						vector<double, 3> vvw = { vv[j][0], vv[j][1], vv[j][2] };
						vj = 2 * vj - vvw;
					}
					else if ((mask[i] & groupbit) && !(mask[j] & groupbit))
					{
						vector<double, 3> vvw = { vv[i][0], vv[i][1], vv[i][2] };
						vi = 2 * vi - vvw;
					}

					//adami2013 
					vector<double, 3> vitransport = { vtrxtmp, vtrytmp, vtrztmp };
					vector<double, 3> vjtransport = { vtr[j][0], vtr[j][1], vtr[j][2] };
					vector<double, 3> delvi = vitransport - vi;
					vector<double, 3> delvj = vjtransport - vj;
					
					ax = 0.5 * ((vi[0] * delvi[0] * rho[i] + vj[0] * delvj[0] * rho[j]) * delx + (vi[0] * delvi[1] * rho[i] + vj[0] * delvj[1] * rho[j]) * dely + (vi[0] * delvi[2] *rho[i] + vj[0] * delvj[2] * rho[j]) * delz );
					ay = 0.5 * ((vi[1] * delvi[0] * rho[i] + vj[1] * delvj[0] * rho[j]) * delx + (vi[1] * delvi[1] * rho[i] + vj[1] * delvj[1] * rho[j]) * dely + (vi[1] * delvi[2] *rho[i] + vj[1] * delvj[2] * rho[j]) * delz );
                    az = 0.5 * ((vi[2] * delvi[0] * rho[i] + vj[2] * delvj[0] * rho[j]) * delx + (vi[2] * delvi[1] * rho[i] + vj[2] * delvj[1] * rho[j]) * dely + (vi[2] * delvi[2] *rho[i] + vj[2] * delvj[2] * rho[j]) * delz );
					//ax = ay = az = 0;
					//const double aveaij = 0.5*(rho[i] * dot(vi, delvi) + rho[j] * dot(vj, delvj));



					const vector<double, 3> vel = vi - vj;
					velx = vel[0];
					vely = vel[1];
					velz = vel[2];

				// compute pressure of atom j
					pj = B[jtype] * (rho[j] / rho0[jtype] - 1.0);

					const double r = sqrt(rsq);
					wfd = kernel[itype][jtype].kernel_gradient(r) / r;

					// using the density-weighted inter-particle averaged pressure(X.Y. Hu,J.C.P,2007)
					pij = (rho[j] * pi + rho[i] * pj) / (rho[i] + rho[j]);



					// the acceleration of particle a caused by shear forces(adami2012)
					fvisc = viscosity[itype][jtype] * ((imass / rho[i]) * (imass / rho[i]) + (jmass / rho[j]) * (jmass / rho[j])) * wfd;					
				    
					fpair =  - ((imass / rho[i])*(imass / rho[i]) + (jmass / rho[j])*(jmass / rho[j])) * pij * wfd;

					fpairx = ((imass / rho[i])*(imass / rho[i]) + (jmass / rho[j])*(jmass / rho[j])) * ax * wfd;
					fpairy = ((imass / rho[i])*(imass / rho[i]) + (jmass / rho[j])*(jmass / rho[j])) * ay * wfd;
					fpairz = ((imass / rho[i])*(imass / rho[i]) + (jmass / rho[j])*(jmass / rho[j])) *  az * wfd;

					fpairpb = - ((imass / rho[i])*(imass / rho[i]) + (jmass / rho[j])*(jmass / rho[j])) * pb  * wfd;

					deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely * vely + velz * velz));

					f[i][0] += delx * fpair + velx  *  fvisc + fpairx;
					f[i][1] += dely * fpair + vely  *  fvisc+ fpairy;
					f[i][2] += delz * fpair + velz   *  fvisc+ fpairz;

					//fpb
					fpb[i][0] += delx * fpairpb;
				    fpb[i][1] += dely * fpairpb ;
				    fpb[i][2] += delz * fpairpb ;

					// change in thermal energy
					de[i] += deltaE;
					drho[i] += rho[i] * (jmass / rho[j]) * delVdotDelR * wfd;


					if (newton_pair || j < nlocal) {
						f[j][0] -= delx * fpair + velx * fvisc + fpairx;
						f[j][1] -= dely * fpair + vely * fvisc + fpairy;
						f[j][2] -= delz * fpair + velz * fvisc + fpairz;
						fpb[j][0] -= delx * fpairpb;
						fpb[j][1] -= dely * fpairpb;
						fpb[j][2] -= delz * fpairpb ;
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

void PairSPHAdami2013::allocate() {
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
	//memory->create(dp, n + 1, "pair:dp");

	kernel.resize(n + 1);
	for (int i = 1; i <= n; ++i)
		kernel[i].resize(n + 1);
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHAdami2013::settings(int narg, char **arg) {
	int rigid_group;
	if (narg != 2)
		error->all(FLERR,
			"Illegal number of setting arguments for pair_style sph/taitwater/morris");
	rigid_group = group->find(arg[0]);
	groupbit = group->bitmask[rigid_group];
	pb = force->numeric(FLERR, arg[1]);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHAdami2013::coeff(int narg, char **arg) {
	if (narg != 7)
		error->all(FLERR,
			"Incorrect args for pair_style sph/taitwater/morris coefficients");
	if (!allocated)
		allocate();

	int ilo, ihi, jlo, jhi;
	force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
	force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

	double rho0_one = force->numeric(FLERR, arg[2]);
	double soundspeed_one = force->numeric(FLERR, arg[3]);
	//double dp_one = force->numeric(FLERR, arg[4]);
	double viscosity_one = force->numeric(FLERR, arg[4]);
	double cut_one = force->numeric(FLERR, arg[5]);
	int kernel_one = force->inumeric(FLERR, arg[6]);
	double B_one = soundspeed_one * soundspeed_one * rho0_one;

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		rho0[i] = rho0_one;
		soundspeed[i] = soundspeed_one;
		//dp[i] = dp_one;
		B[i] = B_one;
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			rho0[j] = rho0_one;
			soundspeed[j] = soundspeed_one;
			B[j] = B_one;
			//dp[j] = dp_one;
			viscosity[i][j] = viscosity_one;
			//printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
			cut[i][j] = cut_one;
			kernel[i][j].init(domain->dimension, kernel_one, cut_one);

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

double PairSPHAdami2013::init_one(int i, int j) {

	if (setflag[i][j] == 0) {
		error->all(FLERR, "Not all pair sph/taitwater/morris coeffs are not set");
	}

	cut[j][i] = cut[i][j];
	viscosity[j][i] = viscosity[i][j];
	kernel[j][i].init(domain->dimension, kernel[i][j].get_kernel_type(), cut[j][i]);


	return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHAdami2013::single(int i, int j, int itype, int jtype,
	double rsq, double factor_coul, double factor_lj, double &fforce) {
	fforce = 0.0;

	return 0.0;
}

