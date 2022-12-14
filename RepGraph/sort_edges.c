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

static int bincompar(BINDEX *a, BINDEX *b);
int sort_edges(EDGE **impedges, int num_imp, BINDEX *index, int pos);
int getedges(NODES **vertex, int num_vertex, EDGE **edge);

int getedges(NODES **vertex, int num_vertex, EDGE **edge)
{
	int	i, j, k, l;

	k = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge[k ++] = vertex[i] -> nextedge[j];
		}
	}
	return(k);
}

int sort_edges(EDGE **impedges, int num_imp, BINDEX *index, int pos)
{
	int	i, j, k, l;
	int	n;
	EDGE	*edge;

	n = 0;
	for(i = 0; i < num_imp; i ++)	{
		for(j = 0; j < impedges[i] -> multip; j ++)	{
			if(impedges[i] -> readinterval[j].eq_read == pos)	{
				index[n].index = i;
				index[n].index_mul = j;
				index[n ++].begin = impedges[i] -> readinterval[j].begin;
			}
		}
	}

	qsort(index, n, sizeof(BINDEX), (void *) bincompar);
	return(n);
}

static int bincompar(BINDEX *a, BINDEX *b)
{
	if(a -> begin > b -> begin)	return(1);
	else				return(-1);
}
