#ifndef __DEFINITIONS_H__
#define __DEFINITIONS_H__

typedef struct {
  int *pattern;

  /* Number of points in the region. */
  int npoints;

  /* Whether it's at the boundary of the cone. */
  int boundary;
} sign_pattern_t;

#endif
