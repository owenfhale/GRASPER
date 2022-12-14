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

int graph(int num_seq, int *len_seq, LIST **list, EDGE **edge);
void movereadinterval(EDGE *edge1, EDGE *edge);

int graph(int num_seq, int *len_seq, LIST **list, EDGE **edge)
{
	int	i, j, k, l, m, n;
	char	c, c1;
	char	*label;
	int	num_edge;
	int	pos1, pos2, *dist;
	POSITION	*position;
	EDGE	*edge1, *edge2;
	NODES	*node, *node_next, *node1, *node2, *node0;

	dist = (int *) ckalloc(10 * sizeof(int));
	for(i = 0; i < 2 * num_seq; i ++)	{
		n = 0;
		for(j = 0; j < len_seq[i]; j ++)	{
			if(list[i][j].node -> num_path <= 1)	{
				break;
			}
		}
		if(j < len_seq[i])	{
			node = list[i][j].node;
			pos1 = j;
		}
		while(j < len_seq[i])	{
			for(j ++; j < len_seq[i]; j ++)	{
				if(list[i][j].node -> num_path <= 1)	{
					break;
				}
			}
			if(j == len_seq[i])	continue;
			pos2 = j;
			node_next = list[i][j].node;
			edge1 = insert_edge(node, node_next, i, pos1, pos2);
			if(pos2 - pos1 > 1)	{
				n ++;
			}
			node = node_next;
			pos1 = pos2;
		}
		if(n >= 7)	{
			dist[7] ++;
		} else	{
			dist[n] ++;
		}
		dist[8] += n;
	}
	printf("Edges assigned.\n");
	printf("Distribution of # assigned long edges in one read.\n");
	printf("-------------------------------------------------------------------------------------\n");
	printf("        0        1        2        3        4        5        6       7+      All\n");
	printf("-------------------------------------------------------------------------------------\n");
	for(i = 0; i < 9; i ++)	{
		printf("%9d", dist[i]);
	}
	printf("\n");
	printf("-------------------------------------------------------------------------------------\n");
	free((void *) dist);

/*	Free nodes with width >= 2	*/

	n = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			node = list[i][j].node;
			if(node && node -> num_lastedge == 0 && node -> num_nextedge == 0)	{
				position = node -> position;
				while(position)	{
					if(position -> position >= 0)	{
						list[position -> readindex][position -> position].node = (NODES *) NULL;
					}
					position = position -> next;
				}
				if(node -> num_path <= 1)	{
					printf("node %d num_path %d edges %d %d\n", node, node -> num_path,
						node -> num_lastedge, node -> num_nextedge);
				}
				free_nodes(node);
				n ++;
			}
		}
	}
	printf("%d nodes are removed.\n", n);

/*	Determine long edges and remove 1-in-1-out nodes	*/

	num_edge = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		node = list[i][0].node;
		if(node && node -> visit == 0)	{
			num_edge = makedge(node, edge, num_edge, list);
		}
	}
	printf("Edge made.\n");

	for(i = 0; i < 2 * num_seq; i ++)	{
		node = list[i][0].node;
		if(node)	{
			node_next = node -> bal_node;
			if(node -> num_lastedge != node_next -> num_nextedge ||
			   node_next -> num_lastedge != node -> num_nextedge)	{
				printf("%d(%d-%d) %d(%d-%d)\n", node, node -> num_lastedge,
					node -> num_nextedge, node_next, node_next -> num_lastedge,
					node_next -> num_nextedge);
				getchar();
			}
		}
	}
	return(num_edge);
}

void movereadinterval(EDGE *edge1, EDGE *edge)
{
	int	i, j, k, l, m;
	READINTERVAL	*tmpreadinterval;

	if(!(edge1 -> readinterval) || !(edge -> readinterval))	{
		edge1 -> multip += edge -> multip;
		return;
	}
	tmpreadinterval = (READINTERVAL *) ckalloc(edge1 -> multip * sizeof(READINTERVAL));
	for(m = 0; m < edge1 -> multip; m ++)	{
		tmpreadinterval[m] = edge1 -> readinterval[m];
	}
	free((void *) edge1 -> readinterval);
	l = edge1 -> multip + edge -> multip;
	edge1 -> readinterval = (READINTERVAL *) ckalloc(l * sizeof(READINTERVAL));
	for(m = 0; m < edge1 -> multip; m ++)	{
		edge1 -> readinterval[m] = tmpreadinterval[m];
	}
	for(m = 0; m < edge -> multip; m ++)	{
		edge1 -> readinterval[edge1 -> multip ++] = edge -> readinterval[m];
	}
	free((void *) tmpreadinterval);
}
