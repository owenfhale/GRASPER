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

double Kimura2D(double TransitionRate, double TranversionRate);
double Align2Dist(char *seq1, char *seq2, int len1, int len2);

double Kimura2D(double TransitionRate, double TranversionRate)
{
	double EvoDist, Ratio;
	double di, dj;

	di = 1 - 2 * TranversionRate;
	dj = 1 - 2 * TransitionRate - TranversionRate;
	EvoDist = -0.25 * log(di * dj * dj);
	if(TranversionRate > 0.0)	{
		Ratio = -0.5 + log(dj) / log(di);
	}
	return(EvoDist);
}

double Align2Dist(char *seq1, char *seq2, int len1, int len2)
{
	int	i, j, k, l, n1, n2;
	int	tot_site, gscore;
	int	windows;
	double	EvoDist, Dist1, Dist2;
	int	*sapp;

	sapp = (int *) ckalloc((len1 + len2 + 1) * sizeof(int));
	windows = abs(len1 - len2) + 1000;
	gscore = ALIGN0(&seq1[-1], &seq2[-1], len1, len2, -windows, windows, W, 16, 4, sapp,
		 len1 + len2 + 1, len1 + len2 + 1);
	i = j = k = 0;
	tot_site = n1 = n2 = 0;
	while(i < len1 && j < len2)	{
		if(sapp[k] == 0)	{
			if(seq1[i] != seq2[j])	{
				if(seq1[i] < 2 && seq2[j] < 2 || seq1[i] >= 2 && seq2[j] >= 2)	{
					n1 ++;
				} else	{
					n2 ++;
				}
			}
			tot_site ++;
			i ++;
			j ++;
		} else if(sapp[k] > 0)	{
			j += sapp[k];
		} else	{
			i -= sapp[k];
		}
		k ++;
	}
	Dist1 = ((double) n1) / tot_site;
	Dist2 = ((double) n2) / tot_site;
	EvoDist = Kimura2D(Dist1, Dist2);
	return(EvoDist);
}
