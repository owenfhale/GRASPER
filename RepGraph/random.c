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

#include <stdinc.h>

int  ran_number(int n, int *idum);
double random1(int *idum);

#define MBIG 1000000000
#define MSEED 16183398
#define MZ 0
#define FAC (1.0/MBIG)

double random1(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj, mk;
	int i, ii,k;

	if(*idum < 0 || iff == 0)	{ /* initialization */
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for(i = 1; i < 54; i ++)	{
			ii = (21 + i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if(mk < MZ)	mk += MBIG;
			mj = ma[ii];
		}
		for(k = 1; k <= 4; k ++)	{
			for(i = 1; i <= 55; i ++)	{
				ma[i] -= ma[1 + (i + 30) % 55];
				if(ma[i] < MZ)	ma[i] += MBIG;
			}
		}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if( ++inext == 56) inext = 1;
	if( ++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if(mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return(mj * FAC);
}

int  ran_number(int n, int *idum)
{
	double	     t;
	int          p;

	t = random1(idum);
	if(t == 1.0)	{
		t = t - 1.0e-10;
	}
	p = ( int ) (n * t);
	return( p );
}
