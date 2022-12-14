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

int readpath(NODES **start_node, PATH *path, int num_seq);
int singlepath(NODES *start_node, PATH *path, int reads, int begin);

int readpath(NODES **start_node, PATH *path, int num_seq)
{
	int	i, j, k, l, m, n;
	int	num_path;
	int	nch;
	int	reads;
	EDGE	*edge, *bal_edge;
	NODES	*v;

/*	Define read paths	*/

	for(i = 0; i < num_seq; i ++)	{
		if(!start_node[i])	{
			continue;
		}
		m = singlepath(start_node[i], path, i, 0);
		l = path[i + num_seq].len_path = path[i].len_path;
		for(j = 0; j < path[i + num_seq].len_path; j ++)	{
			path[i + num_seq].edge[j] = path[i].edge[l - 1 - j] -> bal_edge;
		}
	}
	num_path = 2 * num_seq;
	return(num_path);
}

int singlepath(NODES *start_node, PATH *path, int reads, int begin)
{
	int	i, j, k, l;
	EDGE	*edge;

	for(i = 0; i < start_node -> num_nextedge; i ++)	{
		edge = start_node -> nextedge[i];
		for(j = 0; j < edge -> multip; j ++)	{
			if(edge -> readinterval[j].eq_read == reads &&
			   edge -> readinterval[j].begin == begin)	{
				path[reads].edge[path[reads].len_path ++] = edge;
				k = singlepath(edge -> end, path, reads, begin + edge -> readinterval[j].length - 1);
				return(begin + edge -> readinterval[j].length - 1);
			}
		}
	}
/*
	printf("reads %d begin %d \n", reads, begin);
	for(i = 0; i < start_node -> num_nextedge; i ++)	{
		edge = start_node -> nextedge[i];
		printf("edge %d length %d\n", edge, edge -> length);
		for(j = 0; j < edge -> multip; j ++)	{
			printf("j %d reads %d begin %d length %d\n", j, edge -> readinterval[j].eq_read,
				edge -> readinterval[j].begin, edge -> readinterval[j].length);
		}
	}
*/
	return(k);
}
