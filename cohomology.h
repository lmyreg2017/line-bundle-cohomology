#ifndef __COHOMOLOGY_H__
#define __COHOMOLOGY_H__

/* A cone. We could assume that the triangulation is simplicial, but
   it is not hard to keep it general. */
typedef struct {
  int id; /* Numerical identifier for the cone in its Cech element. */

  int *rays;
  int nrays;

  /* We generally construct the cones by intersection of basic
     cones. We store in here which intersections defined the current
     cone (this gives a well defined ordering of cones). */
  int *intersections;
  /* Number of cones intersected to obtain this cone. Tautologically,
     top dimensional cones have 1 here. */
  int nintersections;
} cone_t;

int compute_kth_cohomology(int k, int *sign_pattern,
			   cone_t **cones, int ncones);

/* Frees the memory associated with the cone, including the pointer
   to the structure itself. */
void free_cone(cone_t *cone);

#endif
