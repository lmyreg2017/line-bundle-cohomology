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
#include "cohomology.h"

#define wrong_input(buf) do {\
  fprintf(stderr, "[%s:%d] Wrong input!!\n", __FILE__, __LINE__);\
  fprintf(stderr, "<<%s>>\n", buf);				 \
  abort();} while(0)

static int point_count = 0;

/* List of sign patterns found */
struct {
  int *pattern;

  /* Number of points in the region. */
  int npoints;

  /* Whether it's at the boundary of the cone. */
  int boundary;
} *patterns = NULL;

/* Total number of patterns */
static int npatterns = 0;

static int dot(int *a, int *b, int dim)
{
  int i, result=0;

  for (i=0; i<dim; i++) {
    result += a[i]*b[i];
  }

  return result;
}

/* Returns 0 nonzero iff equal. */
static int equal(int *a, int *b, int dim)
{
  int j;

  for (j=0; j<dim; j++) {
    if (a[j] != b[j])
      return 0;
  }

  return 1;
}

/* Find a sign pattern in the list of sign patterns. Returns -1 if
   not found in the list. */
static int find_pattern(int *pattern, int nrays)
{
  int i;

  /* We are doing a linear search. If we kept the array sorted we
     could do a binary search and save some time. */
  for (i=0; i<npatterns; i++) {
    if (equal(pattern, patterns[i].pattern, nrays))
      return i; /* Found a match */
  }

  /* No match */
  return -1;
}

/* Adds a pattern to the list, returning its position. */
static int add_pattern(int *pattern, int nrays)
{
  /* Found no match, add the pattern to the list of patterns. */
  patterns = realloc(patterns, (npatterns+1)*sizeof(patterns[0]));
  patterns[npatterns].pattern = malloc(sizeof(int)*nrays);
  memcpy(patterns[npatterns].pattern, pattern, sizeof(int)*nrays);

  patterns[npatterns].npoints = 0;

  patterns[npatterns].boundary = 0;

  npatterns++;

  return (npatterns-1);
}

/* Nonzero iff m is on the boundary of the box */
static int m_in_boundary(int **box, int *m, int dim)
{
  int i;

  for (i=0; i<dim; i++) {
    if ((m[i] == box[i][0]) || (m[i] == box[i][1]))
      return 1;
  }

  return 0;
}

static void traverse_box(int **box, int dim, int **rays, int nrays,
			 int *divisor, int k, int *m)
{
  int i;
  int sign_pattern[nrays];
  /* Pattern for the previous point analyzed. Since neighboring points
     will generally have the same pattern this saves some searches in
     the pattern list. */
  int last_pattern[nrays];
  int last_pattern_index = 0xdeadbeef;

  memset(last_pattern, 0, nrays*sizeof(int));

  for (i=box[k][0]; i<=box[k][1]; i++) {
    m[k] = i;

    if (k<(dim-1)) {
      traverse_box(box, dim, rays, nrays, divisor, k+1, m);
    } else {
      int j;

      /* We have a point in m, analyze it */
      point_count++;
      
      for (j=0; j<nrays; j++) {
	if (dot(m, rays[j], dim) >= -divisor[j]) {
	  sign_pattern[j] = +1;
	} else {
	  sign_pattern[j] = -1;
	}
      }

      /* Find this pattern in the list of patterns. */
      if (!equal(sign_pattern, last_pattern, nrays)) {
	last_pattern_index = find_pattern(sign_pattern, nrays);

	/* We hadn't encountered this sign pattern before, add it to
	   the list. */
	if (last_pattern_index < 0)
	  last_pattern_index = add_pattern(sign_pattern, nrays);

	memcpy(last_pattern, sign_pattern, sizeof(sign_pattern));
      }

      patterns[last_pattern_index].npoints++;
      if (m_in_boundary(box, m, dim))
	patterns[last_pattern_index].boundary = 1;
    }
  }
}

void free_cone(cone_t *cone)
{
  free(cone->intersections);
  free(cone->rays);
  free(cone);
}

/* Read the info for the cohomology to compute from the input
   file. dim is the dimension of the M lattice, and k the cohomology
   we are interested in. */
static void scan_box_info(FILE *fd, int dim, int k)
{
  int **box;
  char *line = NULL;
  size_t nline = 0;
  int i;
  int nrays;
  int **rays;
  int m[dim];
  int *divisor;
  char *ptr;
  cone_t **cones;
  int ncones;
  int result = 0;

  box = malloc(dim*sizeof(int*));

  for (i=0; i<dim; i++) {
    box[i] = malloc(2*sizeof(int));

    if (getline(&line, &nline, fd) < 0)
      wrong_input(line);

    if (sscanf(line, "%d %d", &box[i][0], &box[i][1]) != 2)
      wrong_input(line);
  }

  /* Got the dimension of the box, read the rays */
  if (getline(&line, &nline, fd) < 0)
    wrong_input(line);

  if (sscanf(line, "%d\n", &nrays) != 1)
    wrong_input(line);

  rays = malloc(nrays*sizeof(int*));

  for (i=0; i<nrays; i++) {
    int j;

    rays[i] = malloc(dim*sizeof(int));

    if (getline(&line, &nline, fd) < 0)
      wrong_input(line);

    ptr = line;

    for (j=0; j<dim; j++) {
      char *p;

      rays[i][j] = strtol(ptr, &p, 10);

      if (ptr == p)
	wrong_input(line);

      ptr = p;
    }
  }

  /* Information for the divisor */
  divisor = malloc(nrays*sizeof(int));

  if (getline(&line, &nline, fd) < 0)
    wrong_input(line);

  ptr = line;

  for (i=0; i<nrays; i++) {
    char *p;

    divisor[i] = strtol(ptr, &p, 10);
    
    if (ptr == p)
      wrong_input(line);
    
    ptr = p;
  }

  /* The cones. */
  if (getline(&line, &nline, fd) < 0)
    wrong_input(line);

  if (sscanf(line, "%d\n", &ncones) != 1)
    wrong_input(line);

  cones = malloc(ncones*sizeof(cone_t*));

  for (i=0; i<ncones; i++) {
    int j;

    cones[i] = malloc(sizeof(cone_t));

    if (getline(&line, &nline, fd) < 0)
      wrong_input(line);

    if (sscanf(line, "%d\n", &cones[i]->nrays) != 1)
      wrong_input(line);

    /* Top dimensional cones are the intersection with themselves. */
    cones[i]->id = i;
    cones[i]->nintersections = 1;
    cones[i]->intersections = malloc(sizeof(int));
    cones[i]->intersections[0] = i;

    cones[i]->rays = malloc(cones[i]->nrays*sizeof(int));

    if (getline(&line, &nline, fd) < 0)
      wrong_input(line);

    ptr = line;

    for (j=0; j<cones[i]->nrays; j++) {
      char *p;

      cones[i]->rays[j] = strtol(ptr, &p, 10);

      if (ptr == p)
	wrong_input(line);

      ptr = p;
    }
  }

  /* We read all the information successfully, traverse the box */
  traverse_box(box, dim, rays, nrays, divisor, 0, m);

  /* Compute the cohomology for each compact region. */
  for (i=0; i<npatterns; i++) {
    if (!patterns[i].boundary) {
      result += compute_kth_cohomology(k, patterns[i].pattern, cones, ncones) * patterns[i].npoints;
    }
  }

  for (i=0; i<ncones; i++) {
    free_cone(cones[i]);
  }
  free(cones);

  free(divisor);

  for (i=0; i<nrays; i++)
    free(rays[i]);
  free(rays);

  for (i=0; i<dim; i++)
    free(box[i]);
  free(box);

  for (i=0; i<npatterns; i++) {
    free(patterns[i].pattern);
  }
  free(patterns);

  free(line);

  /* Done, print the result. */
  printf("%d\n", result);
}

int main(int argc, char *argv[])
{
  FILE *fd;
  char *line = NULL, *p;
  size_t nline = 0;
  int dim; /* Dimension of the M lattice */
  int k;

  if (argc != 3) {
    printf("Usage: %s box_info k\n", argv[0]);
    printf("\twhere box_info is the path to a file holding the information\n");
    printf("\tabout the box and the divisors, and k tells the program to\n");
    printf("\tcompute the chain complex relevant for H^k.\n");
    return -1;
  }

  fd = fopen(argv[1], "r");
  if (fd == NULL) {
    perror("fopen");
    printf("ERROR: could not open input file '%s'.\n", argv[1]);
    return -1;
  }

  k = strtol(argv[2], &p, 10);

  if (p == argv[2])
    wrong_input(argv[2]);

  if (getline(&line, &nline, fd) < 0)
    wrong_input(line);

  if (sscanf(line, "%d", &dim) != 1)
    wrong_input(line);

  scan_box_info(fd, dim, k);

  free(line);

  fclose(fd);

  return 0;
}
