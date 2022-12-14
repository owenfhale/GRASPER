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


void reducepath(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2, EDGE *newedge, 
		int *edgematch1, int *edgematch2);

void reducepath(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2, EDGE *newedge, 
		int *edgematch1, int *edgematch2)
{
	int	i, j, k, l, num, m, c, n1, n2, n;
	int	*del_path;
	NODES	*begin, *end;

	del_path = (int *) ckalloc((MAX_BRA + vertex -> num_path) * sizeof(int));

	n1 = n2 = 0;
	num = vertex -> num_path;
	for(i = 0; i < num; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path > 0 && j == path[k].len_path)	{
			if(path[k].edge[j - 1] == edge1)	{
				if(edgematch1[n1])	{
					path[k].edge[j - 1] = newedge;
					if(vertex != newedge -> end)	{
						del_path[i] = 1;
					}
					add_path(newedge -> end, k, j);
				}
				n1 ++;
			}
		} else if(path[k].len_path > 0 && j == 0)	{
			if(path[k].edge[j] == edge2)	{
				if(edgematch2[n2])	{
					path[k].edge[j] = newedge;
					if(vertex != newedge -> begin)	{
						del_path[i] = 1;
					}
					add_path(newedge -> begin, k, j);
				}
				n2 ++;
			}
		} else if(path[k].len_path > 0 && path[k].edge[j - 1] == edge1 && path[k].edge[j] == edge2)	{
			del_path[i] = 1;
			path[k].edge[j - 1] = newedge;
			remove_edge(&path[k], k, j);
		}
	}

	n = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		if(del_path[i] == 0)	{
			vertex -> path_index[n] = vertex -> path_index[i];
			vertex -> path_pos[n] = vertex -> path_pos[i];
			n ++;
		}
	}
	vertex -> num_path = n;

	free((void *) del_path);
}
