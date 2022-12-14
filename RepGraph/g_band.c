/**
>HEADER
    Copyright (c) 2004 Haixu Tang hatang@indiana.edu

    This file is part of the RepGraph package.

    RepGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RepGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RepGraph.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

/* A PACKAGE FOR GLOBALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call ALIGN0(A,B,M,N,L,U,W,G,H,S,MW,MX).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	L : lower bound of the band
	U : upper bound of the band
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	S : script for DISPLAY routine
	MW : maximum window size
	MX : maximum length sequence M to be aligned
*/

#include <stdio.h>
#include <perdef.h>
/*#include <extfunc.h> */

extern char na_name[16];

int ALIGN0(char *A, char *B, int M, int N, int low, int up, int W[][15],
	  int G, int H, int *S, int MW, int MX);

#define MININT -999999999

#define DIGIT 10.0

int *CC, *DD;

/* pointer to the previous crossing point */
static int *CP, *DP;

static int IP;

static int *MP[3];		/* save crossing points */
static char *MT[3];		/* 0: rep, 1: del, 2: ins */
static int *FP;			/* forward dividing points */
static char *FT;

static int (*w)[15];				/* w = W */
static int g, h, m;				/* g = G, h = H, m = g+h */

#define gap(k)  ((k) <= 0 ? 0 : (g+h*(k)))	/* k-symbol indel cost */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

						/* Append "Delete k" op */
#define DEL(k)				\
{ 					\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ 					\
  if (last > 0)				\
    last = sapp[-1] += (k);		\
  else					\
    last = *sapp++ = (k);		\
}

						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
}

/* align(A,B,M,N,up,low,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] and appends such a conversion to the current script.
   tb(te)= 1  no gap-open penalty if the conversion begins(ends) with a delete.
   tb(te)= 2  no gap-open penalty if the conversion begins(ends) with an insert.
*/
static int align(char *A, char *B, int M, int N, int low, int up, char tb,
	   char te)
{
	int rmid, k, l, r, v, kt;
	int t1, t2, t3;

   {	int band, midd;
	int leftd, rightd;	/* for CC, DD, CP and DP */
	register int midc;
	register int curd;	/* current index for CC, DD CP and DP */
	register int i, j;
	register int c, d, e;
	int t, fr, *wa, ib;

	/* Boundary cases: M <= 0 , N <= 0, or up-low <= 0 */
	if (N <= 0) { 
		if (M > 0) DEL(M)
		return;
	}
	if (M <= 0) {
		INS(N)
		return;
	}
	if ((band = up-low+1) <= 1) {
		for (i = 1; i <= M; i++) REP
		return;
	}

	/* Divide: Find all crossing points */

	/* Initialization */
	midd = band/2 + 1;
	rmid = low + midd - 1;
	leftd = 1-low;
	rightd = up-low+1;
	if (leftd < midd) {
		fr = -1;
		for (j = 0; j < midd; j++) 
		    CP[j] = DP[j] = -1;
		for (j = midd; j <= rightd; j++) {
		    CP[j] = DP[j] = 0;
		}
		MP[0][0] = -1;
		MP[1][0] = -1;
		MP[2][0] = -1;
	} else if (leftd > midd) {
		fr = leftd-midd;
		for (j = 0; j <= midd; j++) {
		    CP[j] = DP[j] = fr;
		}
		for (j = midd+1; j <= rightd; j++) 
		    CP[j] = DP[j] = -1;
		MP[0][fr] = -1;
		MP[1][fr] = -1;
		MP[2][fr] = -1;
	} else {
		fr = 0;
		for (j = 0; j < midd; j++) {
		    CP[j] = DP[j] = 0;
		}
		for (j = midd; j <= rightd; j++) {
		    CP[j] = DP[j] = 0;
		}
		MP[0][0] = -1;
		MP[1][0] = -1;
		MP[2][0] = -1;
	}

	CC[leftd] = 0;
	if (tb == 2) t = 0;
	else t = -g;
	for (j = leftd+1; j <= rightd; j++) {
		CC[j] = t = t-h;
		DD[j] = t-g;
	}
	CC[rightd+1] = MININT;
	DD[rightd+1] = MININT;
	if (tb == 1) DD[leftd] = 0;
	else DD[leftd] = -g;
	CC[leftd-1] = MININT;
	for (i = 1; i <= M; i++) {
	    if (i > N-up) rightd--;
	    if (leftd > 1) leftd--;
	    wa = w[A[i]];
	    if ((c = CC[leftd+1]-m) > (d = DD[leftd+1]-h)) {
		d = c;
		DP[leftd] = CP[leftd+1];
	    } else DP[leftd] = DP[leftd+1];
	    if ((ib = leftd+low-1+i) > 0) c = CC[leftd]+wa[B[ib]];
	    if (d > c || ib <= 0) {
		c = d;
		CP[leftd] = DP[leftd];
	    }
	    e = c-g;
	    DD[leftd] = d;
	    CC[leftd] = c;
	    IP = CP[leftd];
	    if (leftd == midd) CP[leftd] = DP[leftd] = IP = i;
	    for (curd=leftd+1; curd <= rightd; curd++) {
	       if (curd != midd) {
		   if ((c = c-m) > (e = e-h)) {
		      e = c;
		      IP = CP[curd-1];
		   }  /* otherwise, IP is unchanged */
		   if ((c = CC[curd+1]-m) > (d = DD[curd+1]-h)) {
		      d = c;
		      DP[curd] = CP[curd+1];
		   } else {
		      DP[curd] = DP[curd+1];
		   }
		   c = CC[curd] + wa[B[curd+low-1+i]];
		   if (c < d || c < e) {
		      if (e > d) {
		         c = e;
		         CP[curd] = IP;
		      } else {
		         c = d;
		         CP[curd] = DP[curd];
		      }
		   } /* otherwise, CP is unchanged */
		   CC[curd] = c;
		   DD[curd] = d;
		} else { /* j == midc */
		   if ((c = c-m) > (e = e-h)) {
		      e = c;
		      MP[1][i] = CP[curd-1];
		      MT[1][i] = 2;
		   } else {
		      MP[1][i] = IP;
		      MT[1][i] = 2;
		   }
		   if ((c = CC[curd+1]-m) > (d = DD[curd+1]-h)) {
		      d = c;
		      MP[2][i] = CP[curd+1];
		      MT[2][i] = 1;
		   } else {
		      MP[2][i] = DP[curd+1];
		      MT[2][i] = 1;
		   }
		   c = CC[curd] + wa[B[curd+low-1+i]];
		   if (c < d || c < e) {
		      if (e > d) {
		         c = e;
		         MP[0][i] = MP[1][i];
		         MT[0][i] = 2;
		      } else {
		         c = d;
		         MP[0][i] = MP[2][i];
		         MT[0][i] = 1;
		      }
		   } else {
			 MP[0][i] = i-1;
			 MT[0][i] = 0;
		   }
		   if (c-g > e) {
			MP[1][i] = MP[0][i];
			MT[1][i] = MT[0][i];
		   }
		   if (c-g > d) {
			MP[2][i] = MP[0][i];
			MT[2][i] = MT[0][i];
		   }
		   CP[curd] = DP[curd] = IP = i;
		   CC[curd] = c;
		   DD[curd] = d;
		}
	    }
	}

	/* decide which path to be traced back */
	if (te == 1 && d+g > c) {
		k = DP[rightd];
		l = 2;
	} else if (te == 2 && e+g > c) {
		k = IP;
		l = 1;
	} else {
		k = CP[rightd];
		l = 0;
	}
	if (rmid > N-M) l = 2;
	else if (rmid < N-M) l = 1;
	v = c;
   }
	/* Conquer: Solve subproblems recursively */

	/* trace back */
	r = -1;	
	for (; k > -1; r=k, k=MP[l][r], l=MT[l][r]){
		FP[k] = r;
		FT[k] = l;
	}
	/* forward dividing */
	if (r == -1) { /* optimal alignment did not cross the middle diagonal */
	   if (rmid < 0) align(A,B,M,N,rmid+1,up,tb,te);
	   else align(A,B,M,N,low,rmid-1,tb,te);
	} else {
	   k = r;
	   l = FP[k];
	   kt = FT[k];

	   /* first block */
	   if (rmid < 0) {
		align(A,B,r-1,r+rmid,rmid+1,min(up,r+rmid),tb,1);
		DEL(1)
	   } else if (rmid > 0) {
		align(A,B,r,r+rmid-1,max(-r,low),rmid-1,tb,2);
		INS(1)
	   }

	   /* intermediate blocks */
	   t2 = up-rmid-1;
	   t3 = low-rmid+1;
	   for (; l > -1; k = l, l = FP[k], kt = FT[k]) {
		if (kt == 0) REP
		else if (kt == 1) { /* right-hand side triangle */
		    INS(1)
		    t1 = l-k-1;
		    align(A+k,B+k+rmid+1,t1,t1,0,min(t1,t2),2,1);
		    DEL(1)
		} else { /* kt == 2, left-hand side triangle */
		    DEL(1)
		    t1 = l-k-1;
		    align(A+k+1,B+k+rmid,t1,t1,max(-t1,t3),0,1,2);
		    INS(1)
		}
	   }

	   /* last block */
	   if (N-M > rmid) {
		INS(1)
		t1 = k+rmid+1;
		align(A+k,B+t1,M-k,N-t1,0,min(N-t1,t2),2,te);
	   } else if (N-M < rmid) {
		DEL(1)
		t1 = M-(k+1);
		align(A+k+1,B+k+rmid,t1,N-(k+rmid),max(-t1,t3),0,1,te);
	   }
	}
	return(v);
}

#ifndef FAR_PTR
#define FCKALLOC ckalloc
#else
#define FCKALLOC fckalloc
#endif

static int CHECK_SCORE0();

int ALIGN0(char *A, char *B, int M, int N, int low, int up, int W[][15],
	  int G, int H, int *S, int MW, int MX)

{ 
	int c, i, j;
	int band;
	int check_score;

	w = W;			/* Setup global parameters */
	g = G;
	h = H;
	m = g+h;
	sapp = S;
	last = 0;
	low = min(max(-M, low),min(N-M,0));
	up = max(min(N, up),max(N-M,0));

	if (N <= 0) { 
		if (M > 0) DEL(M);
		return -gap(M);
	}
	if (M <= 0) {
		INS(N);
		return -gap(N);
	}
	if ((band = up-low+1) <= 1) {
		c = 0;
		for (i = 1; i <= M; i++) {
			REP;
			c += w[A[i]][B[i]];
		}
		return c;
	}

	j = (MW+1+2) * sizeof(int);
	CC = (int *) ckalloc(j);
	DD = (int *) ckalloc(j);
	CP = (int *) ckalloc(j);
	DP = (int *) ckalloc(j);
	
	j = MX+1;
	MT[0] = (char *) ckalloc(j);
	MT[1] = (char *) ckalloc(j);
	MT[2] = (char *) ckalloc(j);
	FT = (char *) ckalloc(j);

	j *= sizeof(int);
	MP[0] = (int *) ckalloc(j);
	MP[1] = (int *) ckalloc(j);
	MP[2] = (int *) ckalloc(j);
	FP = (int *) ckalloc(j);

  	c = align(A,B,M,N,low,up,0,0);

	check_score = CHECK_SCORE0(A,B,M,N,S);
	if (check_score != c) {
		printf("M %d N %d\n", M, N);
		printf("\nCheck_score=%d c = %d\n", check_score, c);
	}

	free((void *) CC);
	free((void *) DD);
	free((void *) CP);
	free((void *) DP);
	for(i = 0; i < 3; i ++)	{
		free((void *) MT[i]);
		free((void *) MP[i]);
	}
	free((void *) FT);
	free((void *) FP);

	return c;
}


/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];

int DISPLAY0(char *A, char *B, int M, int N, int *S, int AP, int BP)
{ register char *a, *b, *c;
  register int   i,  j, op;
           int   lines, ap, bp;

  i = j = op = lines = 0;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  while (i < M || j < N)
    { if (op == 0 && *S == 0)
        { op = *S++;
          *a = na_name[A[++i]];
          *b = na_name[B[++j]];
          *c++ = (*a++ == *b++) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = na_name[B[++j]];
              op--;
            }
          else
            { *a++ = na_name[A[++i]];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N)
        { *a = *b = *c = '\0';
          printf("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
            printf("    .    :");
          if (b <= a+5)
            printf("    .");
          printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}

/* CHECK_SCORE - return the score of the alignment stored in S */

static int CHECK_SCORE0(char *A, char *B, int M, int N, int *S)
{ 
  register int   i,  j, op;
  int score;

  score = i = j = op = 0;
  while (i < M || j < N) {
	op = *S++;
	if (op == 0) {
		score = w[A[++i]][B[++j]] + score;
	} else if (op > 0) {
		score = score - (g+op*h);
		j = j+op;
	} else {
		score = score - (g-op*h);
		i = i-op;
	}
  }
  return(score);
}
