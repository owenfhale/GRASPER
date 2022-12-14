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

int trans_loc(int loc, int *range, int num);

int trans_loc(int loc, int *range, int num)
{
        int     i, j, k, l;

        k = 0;
        if(loc <= range[0])      return(loc);
        k += SEGLEN;
        for(i = 0; i < num - 1; i ++)   {
                if(loc >= range[2 * i + 1] && loc <= range[2 * i + 2])    {
                        return(k + loc - range[2 * i + 1] - 1);
                }
                k += range[2 * i + 2] - range[2 * i + 1];
                k += SEGLEN;
        }
        if(loc >= range[2 * i + 1])      return(k + loc - range[2 * i + 1] - 1);
        printf("Location not found: %d.\n", loc);
        for(i = 0; i < num; i ++)   {
		printf("range %d %d\n", range[2 * i], range[2 * i + 1]);
	}
        exit(-1);
}
