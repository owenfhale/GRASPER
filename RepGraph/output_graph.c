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

#define components 0

int max_multip;
int diff = 0;
int min_length;

int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength, *lmultip;
char partmark;

void output_graph(NODES **vertex, int num_vertex, FILE *fp);
int output_edge(FILE *fp, EDGE *edge0, int num_source);

void output_graph(NODES **vertex, int num_vertex, FILE *fp)
{
	int	i, j, k, l, m, n, len, num_source;
	int	tot_edge;
	int	**num_pa, disttangle[8][7], num_tangle;
	int	*label;
	EDGE	*edge0;
	NODES	*v;

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}

	max_multip = 0;
	min_length = 1000;
	for(i = 0; i < num_vertex; i ++)	{
		vertex[i] -> visit = i;
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge0 = vertex[i] -> nextedge[j];
			edge0 -> visit = 0;
			if(edge0 -> multip > max_multip && edge0 -> length > min_length)	{
				max_multip = edge0 -> multip;
			}
		}
	}

	num_source = 1;
	fprintf(fp, "digraph G {\n");
	fprintf(fp, "\tsize=\"8,8\";\n");
	for(m = 0; m < num_vertex; m ++)	{
		for(j = 0; j < vertex[m] -> num_nextedge; j ++)	{
			edge0 = vertex[m] -> nextedge[j];
			if(edge0 -> visit == 0)	{
				num_source = output_edge(fp, edge0, num_source);
			}
		}
	}
	fprintf(fp, "}\n");

	for(m = 0; m < num_vertex; m ++)	{
		for(j = 0; j < vertex[m] -> num_nextedge; j ++)	{
			edge0 = vertex[m] -> nextedge[j];
			edge0 -> visit = 0;
		}
	}

	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}
	numtangle = (int *) ckalloc(tot_edge * sizeof(int));;
	maxmultip = (int *) ckalloc(tot_edge * sizeof(int));;
	mlength = (int *) ckalloc(tot_edge * sizeof(int));;
	maxlength = (int *) ckalloc(tot_edge * sizeof(int));;
	avelength = (int *) ckalloc(tot_edge * sizeof(int));;
	lmultip = (int *) ckalloc(tot_edge * sizeof(int));;
	printf("---------------------------------------------------------------------------------\n");
	printf("Vertices Edge Source Sink Tangles Super-tangles Overall-length\n");
	printf("---------------------------------------------------------------------------------\n");
	nsuper = 0;
	len = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> start_cover = vertex[i] -> nextedge[j] -> multip;
			len += vertex[i] -> nextedge[j] -> length - 1;
		}
	}
	num_tangle = count_tangle(vertex, num_vertex, disttangle);
	tot_edge = count_edge(vertex, num_vertex, num_pa);
	printf("%-8d %-4d %-6d %-4d %-7d %-13d %-14d\n",
		num_vertex, tot_edge, num_pa[0][1], num_pa[1][0], num_tangle, nsuper, len);
	printf("---------------------------------------------------------------------------------\n");

	if(nsuper > 0)	{
		printf("\nStatistics of the super-tangles:\n");
		printf("---------------------------------------------------------------------------------\n");
		printf("Supertangle #tangles(total len) length of max-multip(multip)    max-length(multip) \n");
		printf("---------------------------------------------------------------------------------\n");
	}
	for(j = 0; j < nsuper; j ++)	{
		printf("%-11d %-8d(%-9d) %-15d(%-8d) %10d(%11d)\n", j + 1, numtangle[j], maxlength[j], mlength[j], maxmultip[j],
			  avelength[j], lmultip[j]);
		numtangle[j] = maxmultip[j] = mlength[j] = 0;
	}
	printf("---------------------------------------------------------------------------------\n\n");
	printf("Distribution of vertex degrees:\n");
	printf("---------------------------------------------------------------------------------\n");
	printf("         \\ Indegree     0    1    2    3    4    5    >6\n");
	printf("Outdegree \\       \n");
	printf("---------------------------------------------------------------------------------\n");
	printf("   0                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[0][0],
		num_pa[0][1], num_pa[0][2], num_pa[0][3], num_pa[0][4], num_pa[0][5], num_pa[0][6]);
	printf("   1                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[1][0],
		num_pa[1][1], num_pa[1][2], num_pa[1][3], num_pa[1][4], num_pa[1][5], num_pa[1][6]);
	printf("   2                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[2][0],
		num_pa[2][1], num_pa[2][2], num_pa[2][3], num_pa[2][4], num_pa[2][5], num_pa[2][6]);
	printf("   3                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[3][0],
		num_pa[3][1], num_pa[3][2], num_pa[3][3], num_pa[3][4], num_pa[3][5], num_pa[3][6]);
	printf("   4                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[4][0],
		num_pa[4][1], num_pa[4][2], num_pa[4][3], num_pa[4][4], num_pa[4][5], num_pa[4][6]);
	printf("   5                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[5][0],
		num_pa[5][1], num_pa[5][2], num_pa[5][3], num_pa[5][4], num_pa[5][5], num_pa[5][6]);
	printf("  >6                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[6][0],
		num_pa[6][1], num_pa[6][2], num_pa[6][3], num_pa[6][4], num_pa[6][5], num_pa[6][6]);
	printf("------------------------------------------------------------------------\n");
	printf("\nNumber of tangles (repeat edges): %d.\n", num_tangle);
	printf("Distribution of tangle multiplicities:\n");
	printf("---------------------------------------------------------------------------------\n");
               printf("            \\ Length     <500     <1000    <5000   <10000   >10000   Total\n");
	printf("Multiplicity \\   \n");
	printf("---------------------------------------------------------------------------------\n");
	printf("    1 (non-repeated) %8d %8d %8d %8d %8d %8d\n", disttangle[0][0], disttangle[0][1],
		disttangle[0][2], disttangle[0][3], disttangle[0][4], disttangle[0][5]);
	printf("    2                %8d %8d %8d %8d %8d %8d\n", disttangle[1][0], disttangle[1][1],
		disttangle[1][2], disttangle[1][3], disttangle[1][4], disttangle[1][5]);
	printf("    3                %8d %8d %8d %8d %8d %8d\n", disttangle[2][0], disttangle[2][1],
		disttangle[2][2], disttangle[2][3], disttangle[2][4], disttangle[2][5]);
	printf("    4                %8d %8d %8d %8d %8d %8d\n", disttangle[3][0], disttangle[3][1],
		disttangle[3][2], disttangle[3][3], disttangle[3][4], disttangle[3][5]);
	printf("    5                %8d %8d %8d %8d %8d %8d\n", disttangle[4][0], disttangle[4][1],
		disttangle[4][2], disttangle[4][3], disttangle[4][4], disttangle[4][5]);
	printf("    6                %8d %8d %8d %8d %8d %8d\n", disttangle[5][0], disttangle[5][1],
		disttangle[5][2], disttangle[5][3], disttangle[5][4], disttangle[5][5]);
	printf("   >6                %8d %8d %8d %8d %8d %8d\n", disttangle[6][0], disttangle[6][1],
		disttangle[6][2], disttangle[6][3], disttangle[6][4], disttangle[6][5]);
	printf("Total                %8d %8d %8d %8d %8d %8d\n", disttangle[7][0], disttangle[7][1],
		disttangle[7][2], disttangle[7][3], disttangle[7][4], disttangle[7][5]);
	printf("---------------------------------------------------------------------------------\n");
	printf("The overall length of the edges is: %d.\n", len);

	free((void *) numtangle);
	free((void *) maxmultip);
	free((void *) mlength);
	free((void *) maxlength);
	free((void *) avelength);
	free((void *) lmultip);
	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
}

int output_edge(FILE *fp, EDGE *edge0, int num_source)
{
	int	i;
	EDGE	*edge;

/*
	if((edge0 -> end -> num_nextedge == 0 || edge0 -> begin -> num_lastedge == 0) &&
	   edge0 -> multip == 1)	{
		fprintf(fp, "\t%d -> %d [", edge0 -> begin -> visit, edge0 -> end -> visit);
		fprintf(fp, "label = \"");
		fprintf(fp, "\"];\n");
	} else if(edge0 -> length == 503)	{
		fprintf(fp, "\t%d -> e%d [", edge0 -> begin -> visit, num_source);
		fprintf(fp, "label = \"");
		fprintf(fp, "\"];\n");
		fprintf(fp, "\tb%d -> %d [", num_source, edge0 -> end -> visit);
		fprintf(fp, "label = \"");
		fprintf(fp, "\"];\n");
		num_source ++;
	} else	{
*/
	edge0 -> visit = 1;
	if(edge0 -> multip > 1)	{
		if(edge0 -> length > min_length && abs(edge0 -> multip - max_multip) <= diff)	{
			fprintf(fp, "\t%d -> %d [style=bold,color=red,", edge0 -> begin -> visit, edge0 -> end -> visit);
		} else	{
			fprintf(fp, "\t%d -> %d [", edge0 -> begin -> visit, edge0 -> end -> visit);
		}
		fprintf(fp, "label = \"s%d(%d,%d", edge0 -> start_cover + 1, edge0 -> length, edge0 -> multip);
		fprintf(fp, ")\"];\n");
	} else	{
		return(num_source);
	}
	if(partmark == 1)	{
		edge0 -> bal_edge -> visit = 1;
	}
	for(i = 0; i < edge0 -> begin -> num_lastedge; i ++)	{
		edge = edge0 -> begin -> lastedge[i];
		if(edge -> visit == 0)	{
			num_source = output_edge(fp, edge, num_source);
		}
	}
	for(i = 0; i < edge0 -> end -> num_nextedge; i ++)	{
		edge = edge0 -> end -> nextedge[i];
		if(edge -> visit == 0)	{
			num_source = output_edge(fp, edge, num_source);
		}
	}
	return(num_source);
}
