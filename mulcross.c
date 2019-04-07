/* MULCROSS.C */  /* 5/3/93 */                       
/* generate test data sets for outlier and cluster detection */
/* copyright 1993,94,96 by David L. Woodruff and David M. Rocke */
/* 3/9/93 contam some number of sqrt(chi^2(p;0.001)/p) on each dim */
/* converted from gencross to create multiple clusters */
/* independant N(0,1) entries */
/* pollute with with N[PollMean, y] data entries */
/* where y is set to be at the "crossover" point */
/* This software may be freely distributed for non-commercial use */

/* Compilation Notes: */
/* 1. This is a single file C program using only standard libraries */
/* 2. you will need the -lm linker switch on most Unix and related systems */

/* Other Notes: */
/* 1. see the function Give_Info() for usage instructions */
/* 2. see the appendix of the paper */
/*    "Idenfication of Outliers in Multivariate Data" */
/*    by Rocke and Woodruff for detailed explanations */
/* 3. the program relpart.c can analyze the test data generated */
/* 4. send questions to dlwoodruff@ucdavis.edu; (916)752-0515 */

#include <stdio.h>       /* magical incantations */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <memory.h>

#define BANNER "mulcross version 1.00\nCopyright 1993,94,96 by David L. Woodruff and David M. Rocke\n"

#define DATAFILE "MULCROSS.DAT"   /* data output file */
#define STATSFILE "MULCROSS.STT"    /* statistics */
#define PARMSFILE "GEN.PRM"   /* problem parameters */
#define SEEDFILE "SEED.DAT"   /* random number seed */
/*----------------------------------------------------------------------------*/
void Give_Info()
{
  printf("This program reads parameters from a file named exactly %s\n",PARMSFILE);
  printf("It reads a integer random number seed with 1 to 5 digits from %s\n",SEEDFILE);
  printf("This seed file is overwritten with a psuedo-random seed by the program\n");
  printf("Data is output to %s and statistics to %s\n",DATAFILE, STATSFILE);
  printf("The parameters file (%s) should contain P N A D\n",PARMSFILE);
  printf("Where P is the dimension, N is the number of points\n");
  printf("A is the fraction of bad data and D controls its distance\n");
  printf("See Rocke and Woodruff: Identification of Outliers in Multivariate Data\n");
  printf("Sample parameters file: 10 100 0.2 2\n");
}

#define ALLCHK(x) if (x == NULL) {printf ("allocation error\n"); assert(x);}

/* global data */
int VectLen;        /* number of attributes */
int XCnt;       /* number of observations */
int PollCnt;        /* number to pollute */
double PollFrac;      /* fraction to pollute */
double NumUnits;      /* location of polluted means */
int NumClusters;      /* number of outliers clusters */
double *X;        /* The observations */
#define Xof(i,j) *(X+(i-1)*VectLen+j-1) /* X[i,j] */
int *JBits;                       /* indicators for J set (squander bits) */
double *XBarJ;                    /* x bar values for the J set */
#define XBarJof(i) *(XBarJ+i-1)   /* XBarJ[i] */
double *XJ;                       /* data matrix corresp. to J set */
#define XJof(i,j) *(XJ+(i-1)*VectLen+j-1) /* XJ[i,j] */

/* covariance matrix, rectangular space to facilitate GJ inversion */
/* the inverse is going to be in the right half of the C matrix, so */
/* space is not allocated, just a macro to get at the data */
double *C;                        /* COLUMN MAJOR covariance matrix */
#define Cof(i,j) (*(C+(i-1)+(j-1)*VectLen))  /* C[i,j] */
#define C_1of(i,j) *(C+VectLen*VectLen+(i-1)+(j-1)*VectLen)     /* C-1[i,j] */
double *Chol;                     /* COLUMN MAJOR "left" cholesky factor */
#define Cholof(i,j) (*(Chol+(i-1)+(j-1)*VectLen))  /* Chol[i,j] */
double *SqResiduals;   /* squared distances, zero based */
double *XSave;         /* scratch space */
#   define XSaveof(i,j) *(XSave+(i-1)*VectLen+j-1) /* col major...*/

long seed;      /* URan seed */
double Lambda;              /* covariance matrix multiplier */

/*----------------------------------------------------------------------------*/
double URan( long *seed)
/*return a uniform 0,1 rv using and returning the seed*/
{
#define c (long)2147483647
   do *seed = ((long long)16807 * *seed) % c; while (*seed == 0);
   return((double)*seed / (double)c);
#undef c
}

/*--------------------------------------------------------------------------*/
double Norm(double mu, double sd, long *z)
/* return a normal deviate with mean mu and std dev sd,
   use seed z */
/* just waste a deviate... */

{
double v1,v2,r;

   do {
     v1 = 2.0 * URan(z) - 1.;
     v2 = 2.0 * URan(z) - 1.;
     r = v1*v1 + v2*v2;
   } while (r >= 1.0);
   return(mu + sd * (v1 * sqrt(-2.0 * log(r)/r)));
}

/*---------------------------------------------------------------------------*/
double Chi2_At_Pt001(int v)
/* return the v deg. of freed. 0.001 point of a chi square*/
/* note that  in this case 0.001 means just that in Bowker and Liberman
   or Bickel and Doksum */
{
    switch (v) {
  case 1: return(10.8);
  case 2: return(13.8);
  case 3: return(16.3);
  case 4: return(18.5);
  case 5: return(20.5);
  case 6: return(22.5);
  case 7: return(24.3);
  case 8: return(26.1);
  case 9: return(27.9);
  case 10: return(29.6);
  case 11: return(31.3);
  case 12: return(32.9);
  case 13: return(34.5);
  case 14: return(36.1);
  case 15: return(37.7);
  case 16: return(39.3);
  case 17: return(40.8);
  case 18: return(43.3);
  case 19: return(43.8);
  case 20: return(45.3);
  case 25: return(52.6);
  case 50: return(86.7);
  case 100: return(149.4);
  default: return(v *
       pow(1-(double)2/(double)(9*v)+0.6745*sqrt(2./(9.*(double)v)),3));
       /* Bowker and Liberman pp 601-602 */
   }
}

/*-------------------------------------------------------------------------*/
void Load_Parms()
/* read in the parameters */
{
    FILE *f;      /* stream record pointer */

  if ((f = fopen(PARMSFILE,"r")) == NULL) {
    Give_Info();
    printf("mulcross could not open paramaters file %s for read\n",PARMSFILE);
    exit(1);
  }
  NumClusters = 0;  // to facilitate gencross input files
  fscanf(f,"%d %d %lf %lf %d",&VectLen, &XCnt, &PollFrac, &NumUnits, &NumClusters);
  fclose(f);
  if (!NumClusters) NumClusters = 1;
  PollCnt = (int)(PollFrac * (double)XCnt + 0.5);
  if (XCnt <= VectLen) {
    printf("For vectors of length %d, there must be at least %d points\n",
      VectLen, VectLen+1);
    exit(1);
  }
  if (PollCnt > XCnt / 2) {
    printf("You have %d data points, so you cannot pollute %d\n",
      XCnt, PollCnt);
    exit(1);
  }
  if (PollCnt <=0 ) {
    printf("You must have at least one bad point (see PollFrac in %s)\n",
           PARMSFILE);
    exit(1);
  }
  if ((f = fopen(SEEDFILE,"r")) == NULL) {
    Give_Info();
    printf("Could not open seed file %s for read\n",SEEDFILE);
    exit(1);
  }
  Lambda = (1. - PollFrac) * (PollFrac * (double)VectLen - (1. - PollFrac)) /
           (PollFrac * ((1.-PollFrac)*(double)VectLen - PollFrac));
  if (Lambda < 0.01) {
     Lambda = 0.01;
     printf("warning: setting lambda to %lf\n",Lambda);
   }

  fscanf(f,"%ld",&seed);
  fclose(f);
}

/*-------------------------------------------------------------------------*/
void Generate_Data()
/* generate the data (see comments at top of file) */
{
    int row, col;       /* to loop */

  for (row = 1; row <= XCnt; row++)
   for (col = 1; col <= VectLen; col++)
    Xof(row,col) =  Norm((double)0.,(double)1.,&seed);
}

/*-------------------------------------------------------------------------*/
void Pollute()
/* find the sample number of the computed sample closest to 0 */
/* no longer write the sample number to the file NEARFILE */
{
  int Samp;     /* sample number, not zero based */
  int i;      /* loop index for distance computation */
  double PollDist;            /* how far out on each dimension */
  double *DirVect;            /* direction vector */
  int SampStart;              /* first polluted sample index */
                 
  DirVect = malloc((VectLen+1)*sizeof(double)); ALLCHK(DirVect)               
  PollDist = NumUnits * sqrt(Chi2_At_Pt001(VectLen) / VectLen);
  SampStart = XCnt - PollCnt + 1;
  for (Samp = SampStart; Samp <=XCnt; Samp++) {
    /* make sure you get enough bad points, even if clusters have differing sizes */
    if (!((Samp-SampStart) % (PollCnt / NumClusters))) {
      for (i=1; i<=VectLen; i++) 
        if (URan(&seed) > 0.5) DirVect[i] = -1.; else DirVect[i] = 1.;
    }
    for (i=1; i<=VectLen; i++) {
        Xof(Samp, i) = Norm(PollDist*DirVect[i], sqrt(Lambda), &seed);
    }
  }
}

/*-------------------------------------------------------------------------*/
void Write_Data()
/* output the X Array */
{
  FILE *f;        /* input file stream record */
  int row,col;      /* to loop */

  if ((f = fopen(DATAFILE,"w")) == NULL) {
    Give_Info();
    printf("Could not open %s for write\n",DATAFILE);
    exit(1);
  }
  fprintf(f,"%d %d\n",VectLen, XCnt);

  for (row = 1; row <= XCnt; row++) {
    for (col = 1; col <= VectLen; col++) fprintf(f," %14.11lf",Xof(row,col));
    fprintf(f,"\n");
  }
  fclose(f);
}

/*-------------------------------------------------------------------------*/
void Dump_Data(char *msg)
/* dump the X array */
{
  int row, col;

  printf("Dump of X data array %s\n",msg);
  for (row=1; row<=XCnt; row++) {
    for (col=1; col<=VectLen; col++) printf("%7.5lf ",Xof(row,col));
    printf("\n");
  }
}

/*-------------------------------------------------------------------------*/
void Make_Room()
/* allocate space for global data structures */
{
  X = malloc(VectLen*XCnt*sizeof(double)); ALLCHK(X)
  JBits = malloc((XCnt)*sizeof(int)); ALLCHK(JBits)
  XBarJ = malloc(VectLen*sizeof(double)); ALLCHK(XBarJ)
  XJ = malloc((XCnt)*VectLen*sizeof(double)); ALLCHK(XJ)
  C = malloc(VectLen*VectLen*2*sizeof(double)); ALLCHK(C)
  Chol = malloc(VectLen*VectLen*sizeof(double)); ALLCHK(Chol)
  XSave = malloc(VectLen*XCnt*sizeof(double)); ALLCHK(XSave)
  SqResiduals = malloc((XCnt)*sizeof(double)); ALLCHK(SqResiduals)
}
/*-------------------------------------------------------------------------*/
void InvertC(double *C, int VectLen, double *Determinant)
/* things are getting a bit hacked up here..... */
/* assume C is VectLen by VectLen */
/* C must have room for another square on the right */
/* do Gauss-Jordon Elimination on C*x(j) = e(j); Strang 29 */
/* see also page 38 (no use exchanging rows since PDS) */
/* keep a running product of the pivots as the determinant of C */
/* (later, try to do something more efficient (e.g. use Symm)) */
/* COLUMN MAJOR arrays, look out... */
/* pivots are left in original C space, result is on ``right'' side */

{
  /* "locals" */
  int pivotrow, row, col;     /* to loop */
  double Pivot;               /* C[pivotrow, pivotrow]        */
  double m;                   /* multiplier */

  /* put I in the right half */
  for (row=1; row <= VectLen; row++)
    for (col = VectLen+1; col <= 2*VectLen; col++)
      if (col == row+VectLen) Cof(row,col) = 1.; else Cof(row,col) = 0.;

  /* forward elimination */
  *Determinant = 1;   /* running product */
  for (pivotrow = 1; pivotrow < VectLen; pivotrow++) {
    if (!Cof(pivotrow,pivotrow)) {
      /* printf("zero pivot at row %d, you lose\n",pivotrow); */
      *Determinant = 0.;
      return;
    } else {
      /* non-zero pivot in C[row,row] */
      Pivot = Cof(pivotrow, pivotrow);   /* two uses for pivot */
      for (row = pivotrow+1; row <= VectLen; row++) {
        m = Cof(row, pivotrow) / Pivot;
        for (col = pivotrow; col <= VectLen+pivotrow; col++) { /* assumes no zero pivots shuffles!!!*/
          Cof(row,col) = Cof(row,col) - m * Cof(pivotrow,col);
        }
      }
    }
    *Determinant = *Determinant * Pivot;

  }
  *Determinant = *Determinant * Cof(VectLen, VectLen);
  if (!*Determinant) {return;}
    
  /* back substitution, pivots are non-zero, so use them */
  for (pivotrow = 2; pivotrow <=VectLen; pivotrow++){
    Pivot = Cof(pivotrow, pivotrow);
    for (row = 1; row < pivotrow; row++) {
      if (Cof(row,pivotrow)) {
        m = Cof(row,pivotrow) / Pivot;
          for (col = pivotrow; col <= VectLen+pivotrow; col++) {
            Cof(row,col) = Cof(row,col) - m*Cof(pivotrow, col);
          }
      }
    }
  }

  /* finally, divide the rows in the inverse by the pivots */
  for (row = 1; row <= VectLen; row++)
    for (col = VectLen+1; col <= 2*VectLen; col++)
      Cof(row, col) = Cof(row,col) / Cof(row,row);
}

/*-------------------------------------------------------------------------*/
void SqrtC(double *C, double *Chol)
/* cholesky factorization */
/* assume C is VectLen by VectLen */
/* do Gauss-Jordon; Strang 241, 29 */
/* see also page 38 (no use exchanging rows since PDS) */
/* we need L root D (not U) for random deviates so transpose */
/* COLUMN MAJOR arrays, look out... */
/* pivots are left in original C space, result is on ``right'' side */

{
  /* "locals" */
  int pivotrow, row, col;     /* to loop */
  double Pivot;               /* C[pivotrow, pivotrow]        */
  double m;                   /* multiplier */

  for (row=1; row <= VectLen; row++) for (col=1; col <= VectLen; col++)
    Cholof(row,col) = Cof(row,col);

  /* forward elimination to produce U and also mult by pivot matrix */
  for (pivotrow = 1; pivotrow < VectLen; pivotrow++) {
    if (!Cholof(pivotrow,pivotrow)) {
      printf("zero pivot at row %d, you lose\n",pivotrow);
      exit(1);
    } else {
      /* non-zero pivot in C[row,row] */
      Pivot = Cholof(pivotrow, pivotrow);   /* two uses for pivot */
      for (row = pivotrow+1; row <= VectLen; row++) {
        m = Cholof(row, pivotrow) / Pivot;
        for (col = pivotrow+1; col <= VectLen; col++) { 
          Cholof(row,col) = Cholof(row,col) - m*Cholof(pivotrow,col);
        }
      }
      Cholof(pivotrow, pivotrow) = sqrt(Pivot);
      for (col = pivotrow+1; col <=VectLen; col++)
       Cholof(pivotrow,col) = Cholof(pivotrow,col) * sqrt(Pivot) / Pivot;
    }
  }
  Cholof(VectLen, VectLen) = sqrt(Cholof(VectLen, VectLen));
  /* do the transpose */
  /* (have upper, want lower) */
  for (row=1; row<=VectLen;row++) {
    for (col=row+1; col<=VectLen;col++) {
      Cholof(col,row) = Cholof(row,col);
      Cholof(row,col) = 0.;
    }
  }
}

/*---------------------------------------------------------------------------*/
void Standardize_X()
/* transform X so that it has unit mean and covariance I */
/* (ab)use various and sundry global data structures */
/* in other words, call this only early in the program */
/* for now, just throw away the transform vector and matrix */
{
  int i,j,k;        /* to loop */
  double Det;       /* to throw away determinant */

  for (j = 1; j <= VectLen; j++) {
    XBarJof(j) = 0;
    for (i = 1; i <= XCnt; i++) XBarJof(j) += Xof(i,j);
    XBarJof(j) /= XCnt;
  }
  for (i = 1; i <= XCnt; i++)
    for (j = 1; j <= VectLen; j++) Xof(i,j) -= XBarJof(j);

  /* now find the appropriate rotation (sqrt(XXt)^-1) */
  /* put XXt into C (but this is XtX in this program... )*/
  for (i=1; i<=VectLen; i++) {
    for (j=1; j<=VectLen; j++) {
      Cof(i,j) = 0;
      for (k=1; k<=XCnt; k++) {
        Cof(i,j) += Xof(k,i) * Xof(k,j);
      }
    }
  }
  InvertC(C, VectLen, &Det);
  if (Det <= 0.0) {printf("Unexpected singularity in mulcross\n"); exit(1);}
  for (i=1; i<=VectLen;i++) for (j=1; j<=VectLen;j++) Cof(i,j) = C_1of(i,j);
  SqrtC(C, Chol);
  /* put sqrt times X (or X times sqrtt here) into XSave then move it to X */
  /* becuase X here is really Xt */
  /* but the sqrt is transposed above .... */
  for (i=1; i<=XCnt; i++) {
    for (j=1; j<=VectLen; j++) {
      XSaveof(i,j) = 0;
      for (k=1; k<=VectLen; k++) {
        XSaveof(i,j) += Xof(i,k) * Cholof(k,j);
      }
    }
  }
  for (i=1; i<=XCnt; i++) for (j=1; j<=VectLen; j++) 
  Xof(i,j) = XSaveof(i,j) * sqrt((double)(XCnt-1)); /* make C near I (est = I)*/
}

/*-------------------------------------------------------------------------*/
void Form_XJ()
/* put rows in XJ corresponding to the J set indexes */
/* note that the indicator set is zero based; also note row data vectors*/
{
  int Xrow, XJrow, col;                /* to loop */

  XJrow = 1;
  for (Xrow = 1; Xrow <= XCnt; Xrow++) {
    if (*(JBits+Xrow-1)) {
      for (col = 1; col <= VectLen; col++)
        XJof(XJrow, col) = Xof(Xrow, col);
      /* the row index is in the J set */
      XJrow++;
    }
  }
}

/*-------------------------------------------------------------------------*/
void Compute_XBarJ(const int JCnt)
/* compute XBarJ values (observations are "row" vectors) */
/* JCnt is the size of the sub-sample */
/* very simple */
{
  int row, col;                   /* to loop */
# define XBCOL *(XBarJ+col-1)   /* typing aid */

  for (col=1; col<=VectLen; col++) {
    XBCOL = XJof(1,col);
    for (row=2; row<=JCnt; row++) XBCOL += XJof(row,col);
    XBCOL = XBCOL / (JCnt);
  }
#undef XBCOL
}

/*-------------------------------------------------------------------------*/
void Dump_XBarJ(char *msg)
/* display the estimate of the center */
{
    int col;                /* to loop */

    printf("Dump of Mean for current J set: %s\n",msg);
    for (col=1; col <= VectLen; col++) printf("  %E\n", XBarJof(col));
}

/*-------------------------------------------------------------------------*/
void Form_C(const int JCnt)
/* form the covariance matrix for XJ */
/* JCnt is the size of the sub-sample */

{
  int i, j;                  /* current cell in C */
  int k;                     /* to loop through samples */

  for (i=1; i<=VectLen; i++) {
    for (j=1; j<=i; j++) {
      Cof(i,j) = 0;
      for (k=1; k<=JCnt; k++) {
        Cof(i,j) += (XJof(k,i) - XBarJof(i))
        * (XJof(k,j) - XBarJof(j));
      }
      Cof(i,j) = Cof(i,j) / (double)(JCnt - 1);
    }
  }
  for (i=1; i<=VectLen;i++) for (j=i+1; j<=VectLen;j++) Cof(i,j) = Cof(j,i);
}


/*--------------------------------------------------------------------------*/
void Compute_Distance_Vector()
/* compute a squared distance vector (called SqResiduals) for the current C_1
 and sub-sample
*/
{
  double rowsum;            /* to accumulate rightmost mult first */
  int i;                    /* index into vector being formed */
  int row, col;             /* indexes for vector matrix mult. */
            /* (row and col refer to c-1) */

  for (i=1; i<=XCnt; i++) {
    *(SqResiduals+i-1) = 0.;
    for (row = 1; row <= VectLen; row++) {
     rowsum = 0;
     for (col = 1; col <= VectLen; col++) {
       rowsum += (Xof(i,col) - XBarJof(col)) * C_1of(row,col);
     }
     *(SqResiduals+i-1) += rowsum * (Xof(i,row) - XBarJof(row));
    }
  }
}
                     
/*------------------------------------------------------------------*/
void Generate_Permutation(int N, int *p)
/* generate permutation of length N in p */
/* p is zero based, but the permutation is not */
{
  int i,j;         /* to loop */
  int lspot;         /* offset in l */
  int *FullList;       /* unpermuted */
  int *l;          /* left to be used */

  FullList = malloc(N * sizeof(int));
  ALLCHK(FullList)
  for (i=0; i<N; i++) *(FullList+i) = i+1;
  l = malloc(N * sizeof(int));
  ALLCHK(l)

  memcpy(l, FullList, sizeof(int)*N);
  for (i=0; i < N; i++) {
    lspot = (int)(URan(&seed) * (N - i));
    *(p+i) = *(l+lspot);
    for (j=lspot; j<N-i; j++) *(l+j) = *(l+j+1);
  }
  free(l); free(FullList);
}
                     
/**************************************************************************/
int main()
{
  FILE *f;           /* to update seed file */
  int i,j,k, row, col; /* loop indexes (cannablized code)*/
  int *Permutation;
  int JCnt;   /* number in basis */
  double Det; /* to ignore the determinant...*/

  printf(BANNER);
  Load_Parms();
  Make_Room();
  Permutation = malloc(XCnt * sizeof(int)); ALLCHK(Permutation)
  Generate_Data();
  Pollute();      
  Standardize_X();
  JCnt = XCnt-PollCnt;
  for (i=0; i<JCnt; i++) *(JBits+i) = 1;
  for (i=JCnt; i<XCnt; i++) *(JBits+i) = 0;
  Form_XJ();
  Compute_XBarJ(JCnt);
  Form_C(JCnt);
  InvertC(C, VectLen, &Det);
    
  memcpy(XSave, X, XCnt * VectLen * sizeof(double));
  /* randomize the rows of X to make X working */
  /* column major is a pain....*/
  Generate_Permutation(XCnt, Permutation);
  for (j=1; j<=XCnt; j++) for (k=1; k<=VectLen; k++) 
      Xof(j,k) = XSaveof((*(Permutation+j-1)),k);
  Compute_Distance_Vector();
  Write_Data();
  if (!(f=fopen(STATSFILE,"w"))) {
    Give_Info();
    printf("Could not open %s for write\n",STATSFILE);
    exit(1);
  }
  fprintf(f,BANNER);
  fprintf(f,"Robust Mean:\n");                                         
  for (col=1; col <= VectLen; col++) fprintf(f,"  %E\n", XBarJof(col));
  fprintf(f,"Robust Mahalonobis distances\n");
  for (row=0; row<XCnt; row++)
     fprintf(f,"point %4d: %.4E\n", row+1, *(SqResiduals+row));
  fclose(f);
  if ((f = fopen(SEEDFILE,"w")) == NULL) {
      printf("Could not open %s for write\n",SEEDFILE);
      exit(1);
  }
  fprintf(f,"%ld\n",seed);
  fclose(f);
  printf("done.\n");
  printf("Parameters read from %s; seed file, %s, read and updated\n",
           PARMSFILE, SEEDFILE);
  printf("Data is in %s, statistics are in %s\n",DATAFILE, STATSFILE);
}
