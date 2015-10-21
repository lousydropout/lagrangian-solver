#include <cs.h>
#include <stdio.h>


cs *from_arrays(int *ivec, int *jvec, double *xvec, int nz)
{
  cs *T, *A;
  int k;
  T = cs_spalloc(0, 0, 1, 1, 1);
  for ( k = 0; k < nz; k++)
    {
      if (!cs_entry (T, ivec[k], jvec[k], xvec[k]))
	return (cs_spfree (T)) ;
    }
  A = cs_compress(T);
  cs_dupl(A);
  cs_dropzeros(A);
  return A;
}

cs *to_cs (int m, int n, int nzmax, int *ivec, int *pvec, double *xvec) {
  cs *A;

  A = cs_spalloc (m, n, 0, 1, 0);
  A->nzmax = nzmax;
  A->m = m;
  A->n = n;
  A->i = ivec;
  A->p = pvec;
  A->x = xvec; 

  /* for (k = 0; k < nzmax; k++) printf ("%g   %g\n", xvec [k], A->x [k]); */
  return A;
}


cs *load (FILE *f)
{
  int i, j ;   /* use double for integers to avoid csi conflicts */
  double x ;
  cs *T ;
  if (!f) return (NULL) ;                             /* check inputs */
  T = cs_spalloc (0, 0, 1, 1, 1) ;                    /* allocate result */
  while (fscanf (f, "%d %d %lg\n", &i, &j, &x) == 3)
    {
      /* printf("%d %d %lg\n", i, j, x); */
      if (!cs_entry (T, i, j, x)) return (cs_spfree (T)) ;
    }
  return (T) ;
}


/* cs_demo2: read a matrix and solve a linear system */
cs *get_cs_prob(FILE *f)
{
  cs *T, *A;
  /* int m, n, mn, nz1, nz2; */
  T = load(f);
  A = cs_compress(T);
  
  cs_spfree(T);
  return A;
}

void print_cs(cs *A)
{
  int j;
  int m, n, nz, nzmax;
  m=A->m; n=A->n; nz=A->nz; nzmax=A->nzmax;
  printf("A.nz= %d\n",nz);
  printf("A.m= %d\n",m);
  printf("A.n= %d\n",n);
  printf("A.nzmax= %d\n",nzmax);

  printf("Printing the values of matrix A\n");
  printf("---------- jcol ------------\n");
  for(j=0; j<A->n+1; j++)
    printf("%d, ",A->p[j]);
  printf("\b\b\n---------- i ------------\n");
  for(j=0; j<nzmax; j++)
    printf("%d, ",A->i[j]);
  printf("\b\b\n---------- x ------------\n");
  for(j=0; j<A->nzmax; j++)
    printf("%g, ",A->x[j]);
  printf("\b\b\n----------------------------\n");
}

void print_sol(cs *Prob, double *x) {
  int m = Prob->m;
  int i;
  printf("Solution \n");
  for(i=0; i<m; i++)
    printf("%g\n",x[i]);
}


double *solve_cs(int n, css *S, csn *N, double *b)
{
  double *x, *y;
  int ok;
  y = cs_malloc(n, sizeof(double)); //work space
  
  ok = (S && N && y);
  if (ok)
    {
      x = b; 
      cs_ipvec (N->pinv, x, y, n) ;       /*    y = r(p) */
      cs_lsolve (N->L, y) ;               /*    y = L\y  */
      cs_usolve (N->U, y) ;               /*    y = U\y  */
      cs_ipvec (S->q, y, x, n) ;          /* r(q) = y    */
    }
  else
    {
      printf("\n\nError in Solve_CS\n\n");
    }
  return x;
}



cs *get_prob()
{
  FILE *f = fopen("t1","r");
  //FILE *f = fopen(name,"r");
  cs *Prob = get_cs_prob(f);
  fclose(f);
  return Prob;
}

