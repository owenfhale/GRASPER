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


int rem_cycle(NODES **vertex, int num_vertex);

int rem_cycle(NODES **vertex, int num_vertex)
{
	int	i, j, k, l, n, m, n1, n2, n3, nl, nall;
	int	*distv;
	int	tot_edge;
	int	nbul;
	NODES	*begin, *end, **vertexl;
	EDGE	**edge, *edge1;

	l = num_vertex * (num_vertex - 1) / 2;
	distv = (int *) ckalloc(l * sizeof(int));

/*	Get the list of edges	*/
	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}
	edge = (EDGE **) ckalloc(tot_edge * sizeof(EDGE *));
	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		vertex[i] -> visit = i;
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge[tot_edge ++] = vertex[i] -> nextedge[j];
		}
	}

/*	Sort the edges with descent multiplicity	*/

	for(i = 0; i < tot_edge; i ++)	{
		for(j = i + 1; j < tot_edge; j ++)	{
			if(edge[i] -> multip < edge[j] -> multip)	{
				edge1 = edge[j];
				edge[j] = edge[i];
				edge[i] = edge1;
			}
		}
	}

/*	Initialize pairwise distance between vertics with SHORTCYC + 1	*/

	for(i = 0; i < l; i ++)		distv[i] = SHORTCYC + 1;

/*	Add the edge the graph with the descent order of
	multiplicy; if an edge form a cycle shorter
	than the threshold, mark this edge and don't add it back
	edge -> visit: 0 	not decided
		       1	added
		       2	removed
*/

	for(i = 0; i < tot_edge; i ++)	{
		begin = edge[i] -> begin;
		end = edge[i] -> end;
/*	a loop in the graph	*/
		if(begin == end)	{
			if(edge[i] -> length <= SHORTCYC)	{
				edge[i] -> visit = 2;
			} else	{
				edge[i] -> visit = 1;
			}
			continue;
		}
		n = numc(begin -> visit, end -> visit);
/*	Should edge[i] be added?	*/
		if(edge[i] -> length + distv[n] - 1 > SHORTCYC)	{
/*	Yes. But lets check its reverse complent first	*/
			if(!(edge[i] -> bal_edge) || edge[i] -> bal_edge -> visit != 2)	{
/*	Add it, if it is not removed	*/
				edge[i] -> visit = 1;
/*	Update its pairwise distances	*/
				distv[n] = min(edge[i] -> length, distv[n]);
				for(j = 0; j < num_vertex; j ++)	{
					if(vertex[j] == begin || vertex[j] == end)	continue;
					n = numc(end -> visit, j);
					m = numc(begin -> visit, j);
					distv[m] = min(distv[m], distv[n] + edge[i] -> length - 1);
					distv[n] = min(distv[n], distv[m] + edge[i] -> length - 1);
					for(k = j + 1; k < num_vertex; k ++)	{
						if(vertex[k] == begin || vertex[k] == end)	continue;
						n1 = numc(end -> visit, k);
						n2 = numc(begin -> visit, k);
						n3 = numc(j, k);
						nall = min(distv[n] + distv[n2] + edge[i] -> length - 1,
							distv[m] + distv[n1] + edge[i] -> length - 1);
						distv[n3] = min(distv[n3], nall);
					}
				}
			} else	{
/*	Otherwise remove this edge too.	*/
				edge[i] -> visit = 2;
			}
		} else	{
/*	No. Remove its reverse complement too	*/
			edge[i] -> visit = 2;
			if(edge[i] -> bal_edge)
				edge[i] -> bal_edge -> visit = 2;
		}
	}

/*	Remove the edges that are not added to the graph	*/

	vertexl = (NODES **) ckalloc(num_vertex * sizeof(NODES *));
	nbul = nl = 0;
	for(i = 0; i < tot_edge; i ++)	{
		if(edge[i] -> visit == 2)	{
			if(edge[i] -> begin -> num_nextedge == 1)	{
				vertexl[nl ++] = edge[i] -> begin;
			}
			if(edge[i] -> end -> num_lastedge == 1)	{
				vertexl[nl ++] = edge[i] -> end;
			}
			erasedge(edge[i]);
			nbul ++;
		}
	}

	for(i = 0; i < nl; i ++)	{
		if(vertexl[i] -> num_lastedge == 0 && vertexl[i] -> num_nextedge == 0)	{
			vertexl[i] = vertexl[nl - 1];
			nl --;
		}
	}
	num_vertex = merge_graph(vertex, num_vertex);
	printf("# of new created shortends: %d\n", nl);

/*	Remove the shortend created by removing the cycles	*/

	for(i = 0; i < nl ; i ++)	{
		if(vertexl[i] -> num_nextedge == 0)	{
			while(vertexl[i] -> num_lastedge > 0)	{
				erasedge(vertexl[i] -> lastedge[0]);
			}
		}
		if(vertexl[i] -> num_lastedge == 0)	{
			while(vertexl[i] -> num_nextedge > 0)	{
				erasedge(vertexl[i] -> nextedge[0]);
			}
		}
	}
	num_vertex = merge_graph(vertex, num_vertex);

	free((void **) vertexl);

	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}
	printf("%d bulges removed, %d vertics %d edges left.\n",
		nbul, num_vertex, tot_edge);

	free((void **) edge);
	free((void *) distv);
	return(num_vertex);
}
