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

int count_comp(double *prob, int *visit, double probthresh, int num_subrep, int index);

int count_comp(double *prob, int *visit, double probthresh, int num_subrep, int index)
{
	int	i, j, k, l, m, n;

	k = 0;
	for(i = 0; i < num_subrep; i ++)	{
		if(visit[i] == 1)	continue;
		visit[i] = 1;
		if(index < i)	{
			n = numc(index, i);
		} else {
			n = numc(i, index);
		}
		if(prob[n] > probthresh)	{
			k += count_comp(prob, visit, probthresh, num_subrep, i);
			k ++;
		}
	}
	return(k);
}
