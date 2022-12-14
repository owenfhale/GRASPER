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

void write_interval(NODES **vertex, int num_vertex, FILE *fp);
void read_interval(NODES **vertex, int num_vertex, FILE *fp);
void write_interval_single(NODES **vertex, int num_vertex, FILE *fp);
void write_edge_inteval(FILE *fp, EDGE *edge);

void write_interval_single(NODES **vertex, int num_vertex, FILE *fp)
{
	int	i, j, k, l, m, n;
	EDGE	*edge;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			if(edge -> multip == 1)	{
				edge -> visit = 1;
			} else	{
				edge -> visit = 0;
			}
		}
	}
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			if(edge -> visit == 1)	{
				continue;
			}
			write_edge_inteval(fp, edge);
		}
	}
}

void write_edge_inteval(FILE *fp, EDGE *edge)
{
	int	i, k;
	EDGE	*edge0;

	fprintf(fp, "EDGE %d Length %d Multiplicity %d.\n", edge -> start_cover + 1, edge -> length, edge -> multip);
	for(k = 0; k < edge -> multip; k ++)	{
		fprintf(fp, "INTV %d %d %d %d\n", edge -> readinterval[k].eq_read, edge -> readinterval[k].begin,
			edge -> readinterval[k].length, edge -> readinterval[k].offset);
	}
	edge -> visit = 1;
	edge -> bal_edge -> visit = 1;
	for(i = 0; i < edge -> begin -> num_lastedge; i ++)	{
		edge0 = edge -> begin -> lastedge[i];
		if(edge0 -> visit == 0)	{
			write_edge_inteval(fp, edge0);
		}
	}
	for(i = 0; i < edge -> end -> num_nextedge; i ++)	{
		edge0 = edge -> end -> nextedge[i];
		if(edge0 -> visit == 0)	{
			write_edge_inteval(fp, edge0);
		}
	}
}

void write_interval(NODES **vertex, int num_vertex, FILE *fp)
{
	int	i, j, k, l, m, n;
	EDGE	*edge;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			fprintf(fp, "EDGE %d Length %d Multiplicity %d.\n", edge -> start_cover + 1, edge -> length, edge -> multip);
			for(k = 0; k < edge -> multip; k ++)	{
				fprintf(fp, "INTV %d %d %d %d\n", edge -> readinterval[k].eq_read, edge -> readinterval[k].begin,
					edge -> readinterval[k].length, edge -> readinterval[k].offset);
			}
		}
	}
}

void read_interval(NODES **vertex, int num_vertex, FILE *fp)
{
	int	i, j, k, l, m, n, multip;
	char	str[500];
	EDGE	*edge;

	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			fgets(str, 400, fp);
			sscanf(str, "%*s%*d%*s%*d%*s%d\n", &multip);
			if(multip != edge -> multip)	{
				if(edge -> multip == 1)	{
					edge -> readinterval = (READINTERVAL *) ckalloc(edge -> multip * sizeof(READINTERVAL));
					continue;
				} else	{
printf("edge -> multip %d multip %d\n", edge -> multip, multip);
					printf("Interval and edge files are not consistent!\n");
					exit(0);
				}
			}
			edge -> readinterval = (READINTERVAL *) ckalloc(edge -> multip * sizeof(READINTERVAL));
			for(k = 0; k < edge -> multip; k ++)	{
				fgets(str, 400, fp);
				sscanf(str, "%*s%d%d%d%d\n", &(edge -> readinterval[k].eq_read), &(edge -> readinterval[k].begin),
					&(edge -> readinterval[k].length), &(edge -> readinterval[k].offset));
			}
			n ++;
		}
	}
}
