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

int n1 = 0, n2 = 0;

int add_int_edge(EDGE **impedges, BINDEX *index, int num, NODES **vertex, int num_vertex, int pos, int num_seq, int len_seq, NODES **start_node);
void buildedge(NODES *node1, NODES *node2, int index, int sp, int ep, int len_seq, int num_seq);
EDGE *single_edge(NODES *node1, NODES *node2, int index, int sp, int ep);

int add_int_edge(EDGE **impedges, BINDEX *index, int num, NODES **vertex, int num_vertex, int pos, int num_seq, int len_seq, NODES **start_node)
{
	int	i, j, k, l, n, m;
	NODES	*node, *bal_node;
	EDGE	*edge1, *edge, *edge2, *bal_edge;
	READINTERVAL *readinterval, *readinterval1;

	if(num == 0)	return(num_vertex);

	if(index[0].begin > 0)	{
		edge2 = impedges[index[0].index];
		node = (NODES *) ckalloc(1 * sizeof(NODES));
		node -> lastedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		node -> nextedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		bal_node = (NODES *) ckalloc(1 * sizeof(NODES));
		bal_node -> lastedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		bal_node -> nextedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		node -> bal_node = bal_node;
		bal_node -> bal_node = node;
		vertex[num_vertex ++] = node;
		vertex[num_vertex ++] = bal_node;
		start_node[pos] = node;
		buildedge(node, edge2 -> begin, pos, 0, index[0].begin, len_seq, num_seq);
	} else	{
		start_node[pos] = impedges[index[0].index] -> begin;
	}
	for(i = 0; i < num - 1; i ++)	{
		readinterval = &(impedges[index[i].index] -> readinterval[index[i].index_mul]);
		readinterval1 = &(impedges[index[i + 1].index] -> readinterval[index[i + 1].index_mul]);
		if(readinterval -> begin + readinterval -> length - 1 < readinterval1 -> begin)	{
			edge1 = impedges[index[i].index];
			edge2 = impedges[index[i + 1].index];
/*
printf("edge1 %d %d edge2 %d %d\n", edge1 -> end -> num_lastedge, edge2 -> end -> num_nextedge,
 edge2 -> begin -> num_lastedge, edge2 -> begin -> num_nextedge);
*/
			buildedge(edge1 -> end, edge2 -> begin, pos, readinterval -> begin + readinterval -> length - 1,
				readinterval1 -> begin, len_seq, num_seq);
		}
	}
	readinterval = &(impedges[index[num - 1].index] -> readinterval[index[num - 1].index_mul]);
	if(readinterval -> begin + readinterval -> length < len_seq)	{
		edge1 = impedges[index[num - 1].index];
		node = (NODES *) ckalloc(1 * sizeof(NODES));
		node -> lastedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		node -> nextedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		bal_node = (NODES *) ckalloc(1 * sizeof(NODES));
		bal_node -> lastedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		bal_node -> nextedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		node -> bal_node = bal_node;
		bal_node -> bal_node = node;
		vertex[num_vertex ++] = node;
		vertex[num_vertex ++] = bal_node;
		buildedge(edge1 -> end, node, pos, readinterval -> begin + readinterval -> length - 1, len_seq - 1, len_seq, num_seq);
	}
	return(num_vertex);
}

void buildedge(NODES *node1, NODES *node2, int index, int sp, int ep, int len_seq, int num_seq)
{
	NODES	*bal_node1, *bal_node2;
	EDGE	*edge, *bal_edge;

	edge = single_edge(node1, node2, index, sp, ep);
	node2 -> lastedge[node2 -> num_lastedge ++] = edge;
	node1 -> nextedge[node1 -> num_nextedge ++] = edge;
	bal_node1 = node2 -> bal_node;
	bal_node2 = node1 -> bal_node;
	bal_edge = single_edge(bal_node1, bal_node2, index + num_seq, len_seq - ep - 1, len_seq - sp - 1);
	bal_node1 -> nextedge[bal_node1 -> num_nextedge ++] = bal_edge;
	bal_node2 -> lastedge[bal_node2 -> num_lastedge ++] = bal_edge;
	edge -> bal_edge = bal_edge;
	bal_edge -> bal_edge = edge;
}

EDGE *single_edge(NODES *node1, NODES *node2, int index, int sp, int ep)
{
	EDGE	*edge;
	READINTERVAL *readinterval;

	edge = (EDGE *) ckalloc(1 * sizeof(EDGE));
	edge -> length = ep - sp + 1;
	edge -> multip = 1;
	readinterval = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
	readinterval -> eq_read = index;
	readinterval -> offset = 0;
	readinterval -> begin = sp;
	readinterval -> length = edge -> length;
	edge -> readinterval = readinterval;
	edge -> begin = node1;
	edge -> end = node2;
	return(edge);
}
