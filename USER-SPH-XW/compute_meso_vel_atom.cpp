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

#include "compute_meso_vel_atom.h"
#include <math.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoVelAtom::ComputeMesoVelAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute meso/vel/atom command");
  if (atom->rho_flag != 1) error->all(FLERR,"compute meso/rho/atom command requires atom_style with density (e.g. meso)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  velVector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMesoVelAtom::~ComputeMesoVelAtom()
{
  memory->sfree(velVector);
}

/* ---------------------------------------------------------------------- */

void ComputeMesoVelAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"velVector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute velVector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMesoVelAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow velVector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(velVector);
    nmax = atom->nmax;
    velVector = (double *) memory->smalloc(nmax*sizeof(double),"atom:velVector");
    vector_atom = velVector;
  }

  // compute velocities for each atom in group

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              velVector[i] = pow(v[i][0]* v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2], 0.5);
      }
      else {
              velVector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoVelAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
