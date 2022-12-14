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

int create_impedges(EDGE **impedges, EDGE *edge, NODES **vertex_new, int num_vertex_new, NODES **vertex, int *num_imp);
NODES *initiate_node(NODES *node1, NODES *node);
EDGE *initiate_edge(EDGE *edge1, EDGE *edge);

int create_impedges(EDGE **impedges, EDGE *edge, NODES **vertex_new, int num_vertex_new, NODES **vertex, int *num_imp)
{
	int	i, j, k, l, m, n;
	NODES	*begin, *end;
	EDGE	*edge1, *edge2, *bal_edge;

	bal_edge = edge -> bal_edge;
/*
printf("edge %d bal_edge %d\n", edge -> multip, bal_edge -> multip);
*/

	edge1 = initiate_edge(edge1, edge);
	if(edge -> begin -> visit > 0)	{
		begin = vertex_new[edge -> begin -> visit - 1];
	} else	{
		begin = initiate_node(begin, edge -> begin);
		vertex_new[num_vertex_new ++] = begin;
		edge -> begin -> visit = num_vertex_new;
	}
	begin -> nextedge[begin -> num_nextedge ++] = edge1;
	edge1 -> begin = begin;
	if(edge -> end -> visit > 0)	{
		end = vertex_new[edge -> end -> visit - 1];
	} else	{
		end = initiate_node(end, edge -> end);
		vertex_new[num_vertex_new ++] = end;
		edge -> end -> visit = num_vertex_new;
	}
	end -> lastedge[end -> num_lastedge ++] = edge1;
	edge1 -> end = end;
	impedges[*num_imp] = edge1;
	(*num_imp) ++;
	edge -> visit = 1;

	if(edge == bal_edge)	{
		edge1 -> bal_edge = edge1;
		edge1 -> begin -> bal_node = edge1 -> end;
		edge1 -> end -> bal_node = edge1 -> begin;
	} else	{
		edge2 = initiate_edge(edge2, bal_edge);
		if(bal_edge -> begin -> visit > 0)	{
			begin = vertex_new[bal_edge -> begin -> visit - 1];
		} else	{
			begin = initiate_node(begin, bal_edge -> begin);
			vertex_new[num_vertex_new ++] = begin;
			bal_edge -> begin -> visit = num_vertex_new;
		}
		begin -> nextedge[begin -> num_nextedge ++] = edge2;
		edge2 -> begin = begin;
		if(bal_edge -> end -> visit > 0)	{
			end = vertex_new[bal_edge -> end -> visit - 1];
		} else	{
			end = initiate_node(end, bal_edge -> end);
			vertex_new[num_vertex_new ++] = end;
			bal_edge -> end -> visit = num_vertex_new;
		}
		end -> lastedge[end -> num_lastedge ++] = edge2;
		edge2 -> end = end;
		impedges[*num_imp] = edge2;
		(*num_imp) ++;
		bal_edge -> visit = 1;
		edge1 -> bal_edge = edge2;
		edge2 -> bal_edge = edge1;
		edge1 -> begin -> bal_node = edge2 -> end;
		edge1 -> end -> bal_node = edge2 -> begin;
		edge2 -> begin -> bal_node = edge1 -> end;
		edge2 -> end -> bal_node = edge1 -> begin;
	}
	return(num_vertex_new);
}

EDGE *initiate_edge(EDGE *edge1, EDGE *edge)
{
	int	i;

	edge1 = (EDGE *) ckalloc(1 * sizeof(EDGE));
	edge1 -> length = edge -> length;
	edge1 -> multip = edge -> multip;
	edge1 -> readinterval = (READINTERVAL *) ckalloc(edge1 -> multip * sizeof(READINTERVAL));
	for(i = 0; i < edge1 -> multip; i ++)	{
		edge1 -> readinterval[i] = edge -> readinterval[i];
	}
	return(edge1);
}

NODES *initiate_node(NODES *node1, NODES *node)
{
	node1 = (NODES *) ckalloc(1 * sizeof(NODES));
	node1 -> lastedge = (EDGE **) ckalloc((MAX_BRA + node -> num_lastedge) * sizeof(EDGE *));
	node1 -> nextedge = (EDGE **) ckalloc((MAX_BRA + node -> num_nextedge) * sizeof(EDGE *));
	return(node1);
}
