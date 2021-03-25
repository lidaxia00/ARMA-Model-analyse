/*  mattr.c    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
void mattr(double *a, double *b, int m, int n)
{
  double *p;
  int i, j;
  for (i = 0; i < n; ++i, ++b)
    for (j = 0, p = b; j < m; ++j, p += n)
      *a++ = *p;
}
/*
mattr

     Transpose an m by n matrix A = B~.

     void mattr(double *a,double *b,int m,int n)
       a = pointer to array containing output n by m matrix 
       b = pointer to array containing input m by n matrix
            (matrices stored in row order)
       m,n = dimension parameters (dim(a)=dim(b)=n*m)

     ------------------------------------------------------------

*/