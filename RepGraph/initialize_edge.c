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

void initial_edge(NODES **vertex, int num_vertex, char **src_seq, int num_seq);
void initedge(EDGE *edge, char **src_seq);

void initial_edge(NODES **vertex, int num_vertex, char **src_seq, int num_seq)
{
	int	i, j, k, l;
	EDGE	*edge, *bal_edge;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			if(!(edge -> seq))	{
				bal_edge = edge -> bal_edge;
				initedge(edge, src_seq);
				if(bal_edge != edge)	{
				   bal_edge -> seq = (char *) ckalloc(bal_edge -> length * sizeof(char));
				   for(k = 0; k < edge -> length; k ++)	{
					bal_edge -> seq[edge -> length - k - 1] = rev(edge -> seq[k]);
				   }
				}
			}
		}
	}
}

void initedge(EDGE *edge, char **src_seq)
{
	int	i, j, k, l, len;
	int	reads, pos;
	char	*seq;

	edge -> seq = (char *) ckalloc(edge -> length * sizeof(char));
	seq = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));

	sortreadinterval(edge -> readinterval, edge -> multip);
	for(i = 0; i < edge -> multip; i ++)	{
		reads = edge -> readinterval[i].eq_read;
		pos = edge -> readinterval[i].offset;
		l = edge -> readinterval[i].length;
		for(j = 0; j < l; j ++)	{
			seq[pos + j] = src_seq[reads][edge -> readinterval[i].begin + j];
		}
	}
	for(i = 0; i < edge -> length; i ++)	{
		edge -> seq[i] = seq[i];
	}
	free((void *) seq);
}
