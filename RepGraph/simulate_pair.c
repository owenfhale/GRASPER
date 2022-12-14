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
#include <extvab.h>
#include <extfunc.h>

double simulate_pair(int length1, int length2, int copynum1, int copynum2, int lengththresh,
	double *distchr, int *num_chrseg, int num_chr, int num_seg, int seglength);
double simulate_pair_all(int length1, int length2, int copynum1, int copynum2, int lengththresh,
	double *distchr, int *num_chrseg, int num_chr, int num_seg, int seglength, int simtimes);

double simulate_pair_all(int length1, int length2, int copynum1, int copynum2, int lengththresh,
	double *distchr, int *num_chrseg, int num_chr, int num_seg, int seglength, int simtimes)
{
	int	i, j, k, l, m, n;
	double	d;

	d = 0;
	for(i = 0; i < simtimes; i ++)	{
		d += simulate_pair(length1, length2, copynum1, copynum2, lengththresh, distchr, num_chrseg, num_chr, num_seg, seglength);
	}
	d /= simtimes;
	return(d);
}

double simulate_pair(int length1, int length2, int copynum1, int copynum2, int lengththresh,
	double *distchr, int *num_chrseg, int num_chr, int num_seg, int seglength)
{
	int i, j, k, l, m, n;
	int k1, k2, *d1, *d2;
	double s1, s2, s0, r;
	int nclosepair;

	if(lengththresh > 2 * seglength)	{
		return(1.0);
	}

	d1 = (int *) ckalloc(copynum1 * sizeof(int));
	d2 = (int *) ckalloc(copynum2 * sizeof(int));

	for(j = 0; j < copynum1; j ++)	{
		s1 = random1(&idum);
		s0 = 0;
		for(i = 0; i < num_seg; i ++)	{
			if(s1 > s0)	{
				d1[j] = i;
				break;
			}
			s0 += distchr[i];
		}
	}
	for(j = 0; j < copynum2; j ++)	{
		s2 = random1(&idum);
		s0 = 0;
		for(i = 0; i < num_seg; i ++)	{
			if(s2 > s0)	{
				d2[j] = i;
				break;
			}
			s0 += distchr[i];
		}
	}
	r = 0;
	for(i = 0; i < copynum1; i ++)	{
		n = k1 = 0;
		for(k = 0; k < num_chr; k ++)	{
			if(d1[i] == n)	{
				k1 = 1;
			}
			n = num_chrseg[k];
		}
		for(j = 0; j < copynum2; j ++)	{
			n = k2 = 0;
			for(k = 0; k < num_chr; k ++)	{
				if(d2[j] == n)	{
					k2 = 1;
				}
				n = num_chrseg[k];
			}
			if(d1[i] == d2[j])	{
				r += ((double) lengththresh) / (seglength - length1 - length2);
			} else if(d1[i] - d2[j] == 1 && k1 == 0 || d2[j] - d1[i] == 1 && k2 == 0)	{
				r += ((double) lengththresh) / (seglength - length1) * ((double) lengththresh) / (seglength - length2);
			}
		}
	}
	free((void *) d1);
	free((void *) d2);
	return(r);
}
