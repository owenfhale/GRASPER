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

void write_graph(NODES **vertex, int num_vertex, FILE *fp, FILE *fp1);

void write_thread(NODES **vertex, int num_vertex, FILE *fp)
{
  int i, j, k , l, n;
  char temp[100];
  EDGE *edge, *bal_edge;
    
  fprintf(fp, "#EdgeNum\tParrallelEdgeNum\tMultiplicity\tLength\n");
}

void write_graph(NODES **vertex, int num_vertex, FILE *fp, FILE *fp1)
{
	int	i, j, k, l, n;
	char	temp[100];
	EDGE	*edge, *bal_edge;

	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			sprintf(temp, "edge%d %d", n + 1, edge -> multip);
			if(edge -> seq)	{
				writeseq(fp, edge -> seq, temp, edge -> length);
			}
			edge -> start_cover = n;
			n ++;
		}
		fflush(fp);
	}
	fprintf(fp1, "Number_of_Vertex %d\n", num_vertex);
	for(i = 0; i < num_vertex; i ++)	{
		fprintf(fp1, "Vertex %d %d %d\n", i, vertex[i] -> num_nextedge, vertex[i] -> num_lastedge);
		fprintf(fp1, "Last_edge ");
		for(j = 0; j < vertex[i] -> num_lastedge; j ++)	{
			bal_edge = vertex[i] -> lastedge[j] -> bal_edge;
			if(bal_edge)	{
				fprintf(fp1, "%d %d ", vertex[i] -> lastedge[j] -> start_cover,
					 vertex[i] -> lastedge[j] -> bal_edge -> start_cover);
			} else	{
				fprintf(fp1, "%d -1 ", vertex[i] -> lastedge[j] -> start_cover);
			}
		}
		fprintf(fp1, "\n");
		fprintf(fp1, "Next_edge ");
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			bal_edge = vertex[i] -> nextedge[j] -> bal_edge;
			if(bal_edge)	{
				fprintf(fp1, "%d %d ", vertex[i] -> nextedge[j] -> start_cover,
					 vertex[i] -> nextedge[j] -> bal_edge -> start_cover);
			} else	{
				fprintf(fp1, "%d -1 ", vertex[i] -> nextedge[j] -> start_cover);
			}
		}
		fprintf(fp1, "\n");
	}
	fflush(fp1);
}


