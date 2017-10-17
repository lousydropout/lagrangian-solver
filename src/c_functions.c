#include <cs.h>

cs_dl *from_arrays (long m, long n, long nz, long *i, long *j, double *x)
{
  cs_dl *T, *A;
  long k;

  T = cs_dl_spalloc (m, n, nz, 1, 1);
  T->m = m; T->n = n; T->nz = nz; T->nzmax = nz;
  for (k = 0; k < nz; k++)
    {
      T->i [k] = i [k] - 1;  // C is zero-based
      T->p [k] = j [k] - 1;
      T->x [k] = x [k];
    }

  A = cs_dl_compress (T);  // create csc A
  cs_dl_spfree (T);        // free triplet T
  cs_dl_dupl (A);          // remove duplicate entries from A
  cs_dl_dropzeros (A);     // drop zeros from A

  return A;
}

cs_dl *to_cs (long m, long n, long nzmax, long *ivec, long *pvec, double *xvec)
{
  cs_dl *A;
  long k;
  A = cs_dl_spalloc (m, n, nzmax, 1, 0);
  for (k = 0; k <= n; k++)    A->p [k] = pvec [k] - 1;
  for (k = 0; k < nzmax; k++) A->i [k] = ivec [k] - 1;
  for (k = 0; k < nzmax; k++) A->x [k] = xvec [k];

  return A;
}


double *solve_cs(long n, cs_dls *S, cs_dln *N, double *b, int err)
{
  double *x, *y;
  long ok;
  err = 0;
  y = cs_dl_malloc (n, sizeof (double)); //work space
  
  ok = (S && N && y);
  if (!S) {err = 1; printf ("Error in Symbolic type");}
  if (!N) {err = 1; printf ("Error in Numeric type");}
  if (!y) {err = 1; printf ("Error: cannot allocate workspace");}
  if (ok)
    {
      x = b; 
      cs_dl_ipvec (N->pinv, x, y, n) ;       /*    y = r(p) */
      cs_dl_lsolve (N->L, y) ;               /*    y = L\y  */
      cs_dl_usolve (N->U, y) ;               /*    y = U\y  */
      cs_dl_ipvec (S->q, y, x, n) ;          /* r(q) = y    */
    }
  else
    {
      printf ("\n\nError in Solve_CS\n\n");
    }
  cs_free (y); // free workspace
  return x;
}
