#include <cs.h>
#include <stdio.h>


double *dvec (long n) {
  long i;
  double *result;
  result = (double *) malloc (n * sizeof (double));
  for (i = 0; i < n; i++) result [i] = 0.314159 + (double) i;
  return result;
}

cs_dl *from_arrays (long m, long n, long nz, long *i, long *j, double *x)
{
  cs_dl *T, *A;
  long k;

  T = cs_dl_spalloc (m, n, nz, 1, 1);
  T->m = m; T->n = n; T->nz = nz; T->nzmax = nz;
  for (k = 0; k < nz; k++) {
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

cs_dl *to_cs (long m, long n, long nzmax, long *ivec, long *pvec, double *xvec) {
  cs_dl *A;
  long k;
  A = cs_dl_spalloc (m, n, nzmax, 1, 0);
  for (k = 0; k <= n; k++)    A->p [k] = pvec [k] - 1;
  for (k = 0; k < nzmax; k++) A->i [k] = ivec [k] - 1;
  for (k = 0; k < nzmax; k++) A->x [k] = xvec [k];

  /** Not safe when inputs are local!
      A = cs_dl_spalloc (m, n, 0, 1, 0);
      A->nzmax = nzmax;
      A->m = m;
      A->n = n;
      A->i = ivec;
      A->p = pvec;
      A->x = xvec; 
  **/
  
  /* for (k = 0; k < nzmax; k++) printf ("%g   %g\n", xvec [k], A->x [k]); */
  return A;
}


cs_dl *load (FILE *f)
{
  long i, j ;   /* use double for integers to avoid csi conflicts */
  double x ;
  cs_dl *T ;
  if (!f) return (NULL) ;                             /* check inputs */
  T = cs_dl_spalloc (0, 0, 1, 1, 1) ;                    /* allocate result */
  while (fscanf (f, "%ld %ld %lg\n", &i, &j, &x) == 3)
    {
      /* printf("%ld %ld %lg\n", i, j, x); */
      if (!cs_dl_entry (T, i, j, x)) return (cs_dl_spfree (T)) ;
    }
  return (T) ;
}


/* cs_dl_demo2: read a matrix and solve a linear system */
cs_dl *get_cs_dl_prob (FILE *f)
{
  cs_dl *T, *A;
  /* long m, n, mn, nz1, nz2; */
  T = load (f);
  A = cs_dl_compress (T);
  
  cs_dl_spfree (T);
  return A;
}

void print_cs (cs_dl *A)
{
  long j;
  long m, n, nz, nzmax;
  m=A->m; n=A->n; nz=A->nz; nzmax=A->nzmax;
  printf ("A.nz= %ld\n",nz);
  printf ("A.m= %ld\n",m);
  printf ("A.n= %ld\n",n);
  printf ("A.nzmax= %ld\n",nzmax);

  printf ("Printing the values of matrix A\n");
  printf ("---------- jcol ------------\n");
  for (j = 0; j < A->n + 1; j++)
    printf ("%ld, ",A->p [j]);
  printf ("\b\b\n---------- i ------------\n");
  for (j = 0; j < nzmax; j++)
    printf ("%ld, ",A->i [j]);
  printf ("\b\b\n---------- x ------------\n");
  for (j = 0; j < A->nzmax; j++)
    printf ("%g, ",A->x [j]);
  printf ("\b\b\n----------------------------\n");
}

void print_sol (cs_dl *Prob, double *x) {
  long m = Prob->m;
  long i;
  printf ("Solution \n");
  for (i = 0; i < m; i++)
    printf ("%g\n",x [i]);
}


double *solve_cs(long n, cs_dls *S, cs_dln *N, double *b)
{
  double *x, *y;
  long ok;
  y = cs_dl_malloc (n, sizeof (double)); //work space
  
  ok = (S && N && y);
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
  return x;
}



cs_dl *get_prob()
{
  FILE *f = fopen ("t1","r");
  cs_dl *Prob = get_cs_dl_prob (f);
  fclose (f);
  return Prob;
}

