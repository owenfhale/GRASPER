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

int comput_bp(int *rep1, int length1, int *rep2, int length2);

int comput_bp(int *rep1, int length1, int *rep2, int length2)
{
	int	i, j, k, l, m, n;
	int	*pointlist1, *pointlist2;

	pointlist1 = (int *) ckalloc(2 * length1 * sizeof(int));
	pointlist2 = (int *) ckalloc(2 * length2 * sizeof(int));
	for(i = 0; i < length1 - 1; i ++)	{
		pointlist1[2 * i] = rep1[i];
		pointlist1[2 * i + 1] = rep1[i + 1];
	}
	for(i = 0; i < length2 - 1; i ++)	{
		pointlist2[2 * i] = rep2[i];
		pointlist2[2 * i + 1] = rep2[i + 1];
	}
	n = 0;
	for(i = 0; i < length1 - 1; i ++)	{
		for(j = 0; j < length2 - 1; j ++)	{
			if((pointlist1[i * 2] == pointlist2[j * 2] && 
			   pointlist1[i * 2 + 1] == pointlist2[j * 2 + 1]) || 
			   (pointlist1[i * 2] == -pointlist2[j * 2 + 1] && 
			   pointlist1[i * 2 + 1] == -pointlist2[j * 2]))	{
				n ++;
				break;
			}
		}
	}
	free((void *) pointlist1);
	free((void *) pointlist2);
	return(n);
}
