/*
  (C) 2010 Iñaki García Etxebarria

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "cohomology.h"
#include "chomp.h"

/* Computes n! */
static long long unsigned int factorial(int n)
{
  if (n<=1)
    return 1;

  return n*factorial(n-1);
}

/* Compute n choose m = n!/(m!(n-m)!)*/
static long long unsigned int choose_m_in_n(int m, int n)
{
  int i;
  long long unsigned int result = 1;

  /* We define (n 0) to be 1. */
  if (m == 0)
    return 1;

  for (i=0; i<m; i++)
    result *= n-i;

  return (result / factorial(m));
}

/* Returns nonzero iff the given vector contains the given value */
static int contains(int *v, int len, int t)
{
  int i;

  for (i=0; i<len; i++) {
    if (v[i] == t)
      return 1;
  }

  return 0;
}

/* Returns the cone given by the rays contained in all the given
   cones. This allocates a new result. */
static cone_t *intersect_cones(cone_t **cones, int ncones, int id)
{
  int i;
  cone_t *result = malloc(sizeof(cone_t));

  result->id = id;
  result->nintersections = ncones;
  result->intersections = malloc(ncones*sizeof(int));

  /* This is slightly memory inefficient, but it saves reallocs. */
  result->rays = malloc(cones[0]->nrays*sizeof(int));
  result->nrays = 0;

  for (i=0; i<ncones; i++) {
    result->intersections[i] = cones[i]->id;
  }

  for (i=0; i<cones[0]->nrays; i++) {
    int j;

    for (j=1; j<ncones; j++) {
      if (!contains(cones[j]->rays, cones[j]->nrays, cones[0]->rays[i]))
	break;
    }

    /* If j == cones the ray is contained in all cones, otherwise the
       previous loop would have stopped early. */
    if (j == ncones) {
      result->rays[result->nrays++] = cones[0]->rays[i];
    }
  }

  return result;
}

/* We need to pass quite a few parameters to choose_k_cones, so we
   pack them neatley here. */
typedef struct {
  cone_t **Cech;
  /* Number of elements added to the Cech complex */
  int nitems;

  int *sign_pattern;

  cone_t **result;
  int nresult; /* Number of cones in the result */
} choose_params_t;

/*
  Chooses k different cones from the set of cones, storing them in
  result.
*/
static void choose_k_cones(cone_t **cones, int ncones, int k,
			   choose_params_t *par)
{
  int i;

  for (i=0; i<=(ncones-k); i++) {
    if (k==1) {
      cone_t *intersection;
      int j;

      /* Done choosing, we have all desired cones. Process the
	 intersection, and add it to the Cech complex if no rays in
	 the intersection are negative. */
      par->result[par->nresult] = cones[i];

      intersection = intersect_cones(par->result, par->nresult+1, par->nitems);

      for (j=0; j<intersection->nrays; j++) {
	if (par->sign_pattern[intersection->rays[j]] < 0) {
	  /* The monomial is not well defined in the intersection, do
	     not bother with it. */
	  break;
	}
      }

      if (j == intersection->nrays) {
	/* The monomial is well defined in this patch, add it to the
	   complex. */
	par->Cech[par->nitems++] = intersection;
      } else
	free_cone(intersection);
    } else {
      /* Not done choosing. Store this cone and select the next
	 element. */
      par->result[par->nresult++] = cones[i];
      choose_k_cones(&cones[i+1], ncones-(i+1), k-1, par);
      par->nresult--;
    }
  }
}

/* Populate one entry of the Chech cochain. The elements of the k-th
   entry are intersections of k+1 maximal dimensional cones, such that
   the monomial is well defined for all rays in the cone.  Returns how
   many elements are in the entry.
*/
static int populate_cech(cone_t **Cech, int degree,
			 int *sign_pattern, cone_t **cones, int ncones)
{
  cone_t *cones_to_intersect[degree+1];
  choose_params_t par = {
    .Cech = Cech,
    .nitems = 0,
    .sign_pattern = sign_pattern,
    .result = cones_to_intersect,
    .nresult = 0
  };

  /* C^{-1} is slightly special, it just denotes the empty set,
     with zero differential. */
  if (degree == -1) {
    par.nitems = 0;
  } else
    choose_k_cones(cones, ncones, degree+1, &par);
  
  return par.nitems;
}

/* Return +1 if a>b, -1 if a<b, and 0 if they are equal. The
   comparison proceeds term by term. */
static int compare(int *a, int *b, int dim)
{
  int j;

  for (j=0; j<dim; j++) {
    if (a[j]<b[j])
      return -1;
    else if (a[j]>b[j])
      return +1;
  }

  return 0;
}

/* Find the given intersection in the given list of cones, starting
   from start (which is excluded from the search). This assumes that
   the list of cones is ordered by intersection, which is true by
   construction in our case. */
static int find_intersection(int *intersection, int dim,
			     cone_t **cones, int ncones, int origin)
{
  int start = origin+1, end = ncones-1;

  /* The algorithm below simplifies slightly if we can tell for sure
     that the boundaries of the interval do not match. */
  if (compare(intersection, cones[start]->intersections, dim) == 0)
    return start;
  
  if (compare(intersection, cones[end]->intersections, dim) == 0)
    return end;

  while (1) {
    int midpoint;

    assert(end-start > 1);

    /* Bisect and compare */
    midpoint = (start+end)/2;

    switch (compare(intersection, cones[midpoint]->intersections, dim)) {
    case 1:
      start = midpoint;
      break;
    case -1:
      end = midpoint;
      break;
    default:
      return midpoint;
      break;
    }
  }
}

/* Build the appropriate differentials in the Cech complex, and
   output them to CHomP. */
static void build_differentials(cone_t **from, int nfrom, cone_t **to, int nto,
				int ncones, FILE *chomp)
{
  int i;

  for (i=0; i<nfrom; i++) {
    int j;
    int target[from[i]->nintersections+1];
    int signature = +1;
    int last_find = -1;

    fprintf(chomp, "   boundary %d =", i+1);

    for (j=0; j<=from[i]->nintersections; j++) {
      int k;
      int min, max;

      if (j>0)
	min = from[i]->intersections[j-1]+1;
      else
	min = 0;

      if (j<from[i]->nintersections)
	max = from[i]->intersections[j]-1;
      else
	max = ncones-1;

      for (k=min; k<=max; k++) {
	int t;

	for (t=0; t<j; t++)
	  target[t] = from[i]->intersections[t];
	target[t] = k;
	for (t=j; t<from[i]->nintersections; t++)
	  target[t+1] = from[i]->intersections[t];

	/* We now have built the target differential in
	   target. Search for it in the destination module. */
	last_find = find_intersection(target, from[i]->nintersections+1,
				      to, nto, last_find);

	fprintf(chomp, " %s %d", (signature>0)?"+":"-", last_find+1);
      }

      signature = -signature;
    }

    fprintf(chomp, "\n");
  }
}

int compute_kth_cohomology(int k, int *sign_pattern,
			   cone_t **cones, int ncones)
{
  FILE *chomp;
  int i;
  cone_t **Cech[3]; /* C^{k-1}, C^{k}, C^{k+1} */
  int nCech[3]; /* Elements in the k-th Cech complex */
  int result;
  char *fname = "__chomp_input__";

  /* CHomP computes things in terms of chain complexes, instead of
     cochain complexes. So we have to invert the sense of the arrows,
     but otherwise things map easily. In particular, we start from
     the highest node in the Cech complex, C^{k+1}, and call that the
     zero dimensional space. */
  chomp = fopen(fname, "wb");

  fprintf(chomp, "chain complex\n\n");

  fprintf(chomp, "max dimension = 2\n\n");

  /* Populate the Cech patches */
  for (i=0; i<3;i++) {
    /* Maximum possible number of patches. In general not all will be
       properly defined, this is just the theoretical maximum so we
       just need to allocate things once. */
    int nmax = choose_m_in_n(k+i, ncones);

    Cech[i] = malloc(sizeof(cone_t*)*nmax);

    /* Generate all the possible combinations (in a well defined order
       so we can do binary searches later on). */
    nCech[i] = populate_cech(Cech[i], k-1+i, sign_pattern,
			     cones, ncones);
  }

  /* We now have the elements of the Cech complex in place, let us
     build the differentials. Recall that we are building a chain
     complex, so the ordering is different to one in which the Cech
     complex is usually presented. */
  for (i=0; i<3; i++) {
    fprintf(chomp, "dimension %d: %d\n", i, nCech[2-i]);

    /* The lowest complex has always zero differential */
    if (i==0) {
      int j;

      for (j=0; j<nCech[2-i]; j++) {
	fprintf(chomp, "   boundary %d = 0\n", j+1);
      }
    } else {
      build_differentials(Cech[2-i], nCech[2-i], Cech[3-i], nCech[3-i],
			  ncones, chomp);
    }

    fprintf(chomp, "\n");
  }

  fclose(chomp);

  /* We have all the data in fname, pass this on to CHomP */
  result = run_chomp(fname);

  /* Done with the file, clean up. */
  unlink(fname);

  for (i=0; i<3; i++) {
    int k;

    for (k=0; k<nCech[i]; k++) {
      free_cone(Cech[i][k]);
    }

    free(Cech[i]);
  }

  return result;
}
