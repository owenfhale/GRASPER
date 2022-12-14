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

#define MAX_NODES 40000
#define MAX_EDGE  80000

int read_graph(NODES **vertex, EDGE **edge, int *num_edge, FILE *fp, FILE *fp1);

int read_graph(NODES **vertex, EDGE **edge, int *num_edge, FILE *fp, FILE *fp1)
{
	int	i, j, k, l, n, m, d, num_vertex, p1, p2, multip;
	char	targ_name[100], *targ_seq, str[5000];

	targ_seq = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));

	k = -1;
	while(fgets(str, 4990, fp))	{
		if(str[0] == '>')	{
			sscanf(&str[2], "%s%d", targ_name, &multip);
			if(k >= 0)	{
				edge[k] -> seq = (char *) ckalloc(n * sizeof(char));
				for(i = 0; i < n; i ++)	{
					edge[k] -> seq[i] = targ_seq[i];
				}
				edge[k] -> length = n;
			}
			k ++;
			edge[k] = (EDGE *) ckalloc(1 * sizeof(EDGE));
			edge[k] -> multip = multip;
			n = 0;
		} else {
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					targ_seq[n ++] = char2int(str[i]);
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					targ_seq[n ++] = char2int(str[i] - 'A' + 'a');
				}
			}
		}
	}
	edge[k] -> seq = (char *) ckalloc(n * sizeof(char));
	for(i = 0; i < n; i ++)	{
		edge[k] -> seq[i] = targ_seq[i];
	}
	edge[k] -> length = n;
	k ++;

	*num_edge = k;

	p1 = p2 = 0;

	n = -1;
	while(fgets(str, 4990, fp1))	{
/*
printf("str %s\n", str);
*/
		if(!strncmp(str, "Number_of_Vertex", 16))	{
			sscanf(str, "%*s%d", &num_vertex);
		} else if(!strncmp(str, "Vertex", 6))	{
			n ++;
			vertex[n] = (NODES *) ckalloc(1 * sizeof(NODES));
			sscanf(str, "%*s%*d%d%d", &(vertex[n] -> num_nextedge), &(vertex[n] -> num_lastedge));
			vertex[n] -> lastedge = (EDGE **) ckalloc((vertex[n] -> num_lastedge) * sizeof(EDGE *));
			vertex[n] -> nextedge = (EDGE **) ckalloc((vertex[n] -> num_nextedge) * sizeof(EDGE *));
			p1 += vertex[n] -> num_lastedge;
			p2 += vertex[n] -> num_nextedge;
		} else if(!strncmp(str, "Last_edge", 9))	{
			l = 9;
			for(j = 0; j < vertex[n] -> num_lastedge; j ++)	{
				sscanf(&str[l + 1], "%d %d", &m, &d);
				vertex[n] -> lastedge[j] = edge[m];
				if(d >= 0) 	{
					vertex[n] -> lastedge[j] -> bal_edge = edge[d];
					if(vertex[n] -> bal_node == NULL)
						vertex[n] -> bal_node = edge[d] -> begin;
				}
				edge[m] -> end = vertex[n];
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
			}
		} else if(!strncmp(str, "Next_edge", 9))	{
			l = 9;
			for(j = 0; j < vertex[n] -> num_nextedge; j ++)	{
				sscanf(&str[l + 1], "%d %d", &m, &d);
				vertex[n] -> nextedge[j] = edge[m];
				if(d >= 0)
					vertex[n] -> nextedge[j] -> bal_edge = edge[d];
				edge[m] -> begin = vertex[n];
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
			}
		}
	}

	num_vertex = n + 1;
        for(i = 0; i < num_vertex; i ++)        {
                if(vertex[i] -> num_lastedge > 0)       {
                        vertex[i] -> bal_node = vertex[i] -> lastedge[0] -> bal_edge -> begin;
                } else  {
                        vertex[i] -> bal_node = vertex[i] -> nextedge[0] -> bal_edge -> end;
                }
        }
	free((void *) targ_seq);
	return(num_vertex);
}


/**********************************************************************
 * Read the graph files 
 *
 * TODO:
 * 1. adjust allocation sizes to input, instead of assuming fixed sizes
 * 2. Return a graph structure instead of individual parameters
 **********************************************************************/

void read_graph_file(char *edgefile,
		     char *graphfile,
		     int *num_vertex,
		     NODES ***vertex,
		     int *num_edge,
		     EDGE ***edge)
{
  FILE *fp;
  FILE *fp1;

  fp = ckopen(edgefile, "r");
  fp1 = ckopen(graphfile, "r");
  *vertex = (NODES **) ckalloc(MAX_NODES * sizeof(NODES *));
  *edge = (EDGE **) ckalloc(MAX_EDGE * sizeof(EDGE *));
  *num_vertex = read_graph(*vertex, *edge, num_edge, fp, fp1);
  printf("Input graph ... done. %d vertices and %d edges.\n",
	 *num_vertex, *num_edge);
  fclose(fp1);
  fclose(fp);
}
