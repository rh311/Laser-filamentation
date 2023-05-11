#ifndef FILAMENTATION_REFACTORING_GRID_H
#define FILAMENTATION_REFACTORING_GRID_H

struct grid_t_r_z {
  int Nr;
  int Nt;
  int Nz;
  double *r;
  double *t;
  double *z;
};

#endif
