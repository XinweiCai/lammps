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
#include "pair_sph_adami_boundary_stationary.h"
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
#include "sdpd_kernel.h"
#include "math_vector2.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHABStationary::PairSPHABStationary(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair

  comm_forward = 4;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHABStationary::~PairSPHABStationary() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairSPHABStationary::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHABStationary::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype, numf;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, h, ih, ihsq;
  int *jlist;
  double wf, sumwf, psum, pi, pj, tmp;
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->vest; 
  double **vv = atom->vv; //
  double *rho = atom->rho;
  double **f = atom->f;
  int *type = atom->type;
  double *mass = atom->mass;
  int *mask = atom->mask;

  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 0.0) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact, but not all of their single particle properties are set.\n",
                  i, j);
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


  for (ii = 0; ii < inum; ii++) {
	  i = ilist[ii];
	  if (mask[i] & groupbit) { // 若粒子i为边界粒子
		  xtmp = x[i][0];
		  ytmp = x[i][1];
		  ztmp = x[i][2];
		  itype = type[i];
		  jlist = firstneigh[i];
		  jnum = numneigh[i];
		  psum = 0;
		  sumwf = 0;
		  numf = 0;
		  vector<double, 3> vtemp = { 0,0,0 };
		  for (jj = 0; jj < jnum; jj++) {
			  j = jlist[jj];
			  j &= NEIGHMASK;
			  jtype = type[j];
			  if (!(mask[j] & groupbit)) { //若粒子j为流体粒子
				  vector<double, 3> xw = { xtmp, ytmp, ztmp }; // wall particle position
				  vector<double, 3> xf = { x[j][0],x[j][1], x[j][2] }; // fulid particle position 
				  vector<double, 3> vf = { v[j][0],v[j][1], v[j][2] }; //fluid particle velocity
				  vector<double, 3> rwf = xw - xf; // wall particle -> fluid particle 位置矢量
				  vector<double, 3> g = { gx, gy, gz }; //body force
				  double dwf = norm(rwf); //fluid and wall particle distance

				  if (dwf < cut[jtype][itype]) {
            pj = rp * (pow(rho[j]/rho0, ga)-1) + dp; // get the fluid particle  pressure from EOC
					  wf = kernel[itype][jtype].kernel(dwf); // get the kernal function
					  psum += pj * wf + dot(g, rwf) * rho[j] * wf; 
					  sumwf += wf;
					  vtemp += wf * vf;
					  numf++;
				  }
			  }
		  }
		  if (numf != 0) { // numf 是为了保证边界粒子的cutoff中含有fluid particle
			  vector<double, 3> vvw = vtemp / sumwf;
			  pi = psum / sumwf;
			  //rho[i] = rho0 * pow((pi - dp) / rp + 1, 0.14285714285714285);
        rho[i] = rho0 * pow((pi - dp) / rp + 1, 1.0 / ga); // 得到wall particle 的 密度
        //边界粒子的虚拟速度
			  vv[i][0] = vvw[0];
			  vv[i][1] = vvw[1];
			  vv[i][2] = vvw[2];
		  }
		  else 
		  { // 若wall particle 附近没有fluid particle ,则密度设为初始密度，虚拟速度设为0
			  rho[i] = rho0;
			  vv[i][0] = 0;
			  vv[i][1] = 0;
			  vv[i][2] = 0;
		  }
	  }
  }


  

  // communicate densities
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHABStationary::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");

  kernel.resize(n + 1);
  for (int i = 1; i <= n; ++i)
    kernel[i].resize(n + 1);
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHABStationary::settings(int narg, char **arg) {
  int rigid_group;
  if (narg != 8)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/rhosum/rigid");
  rigid_group = group->find(arg[0]); // 将所有的边界粒子变成一个group
  groupbit = group->bitmask[rigid_group];
  // g : body force
  gx = force->numeric(FLERR, arg[1]);
  gy = force->numeric(FLERR, arg[2]);
  gz = force->numeric(FLERR, arg[3]);
  dp = force->numeric(FLERR, arg[4]);
  rho0 = force->numeric(FLERR, arg[5]);
  soundspeed = force->numeric(FLERR, arg[6]);
  // ga 为参数，一般在Tait's EOS 中取7，而Adami2013取了1
  ga =  force->numeric(FLERR, arg[7]);
  //rp 背景pressure
  rp = soundspeed * soundspeed * rho0 / ga;

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHABStationary::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for sph/rhosum coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR,arg[2]);
  int kernel_one = force->inumeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      kernel[i][j].init(domain->dimension, kernel_one, cut_one);
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

double PairSPHABStationary::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/rhosum coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  kernel[j][i].init(domain->dimension, kernel[i][j].get_kernel_type(), cut[j][i]);

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHABStationary::single(int i, int j, int itype, int jtype, double rsq,
    double factor_coul, double factor_lj, double &fforce) {
   fforce = 0.0;

   return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHABStationary::pack_forward_comm(int n, int *list, double *buf,
                                     int pbc_flag, int *pbc) {
  int i, j, m;
  double *rho = atom->rho;
  double **vv = atom->vv;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
	buf[m++] = vv[j][0];
	buf[m++] = vv[j][1];
	buf[m++] = vv[j][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHABStationary::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *rho = atom->rho;
  double **vv = atom->vv;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
      rho[i] = buf[m++];
	  vv[i][0] = buf[m++];
	  vv[i][1] = buf[m++];
	  vv[i][2] = buf[m++];  
  }
}
