#ifndef ASA047_H 
#define ASA047_H 1

void nelmin ( double fn ( double x[] ), int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault );
void timestamp ( );

#endif
