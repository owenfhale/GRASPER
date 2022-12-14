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

int chk_consist(PATH *startpath, PATH *midpath, int num_midpath, int *n);

int chk_consist(PATH *startpath, PATH *midpath, int num_midpath, int *n)
{
	int	i, j, k, l, c;

	c = -1;
	for(i = 0; i < num_midpath; i ++)	{
		if(startpath -> len_path <= midpath[i].len_path)	{
			for(j = 0; j < startpath -> len_path; j ++)	{
				if(startpath -> edge[j] != midpath[i].edge[j])	{
					break;
				}
			}
			if(j == startpath -> len_path)	{
				c = 0;
				*n = i;
				break;
			}
		} else	{
			for(j = 0; j < midpath[i].len_path; j ++)	{
				if(startpath -> edge[j] != midpath[i].edge[j])	{
					break;
				}
			}
			if(j == midpath[i].len_path)	{
				c = i + 1;
				break;
			}
		}
	}

	return(c);
}
