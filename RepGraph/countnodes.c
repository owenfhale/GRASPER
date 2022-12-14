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

int singstatnode(NODES *node, int **nstat);
void statnode(LIST **list, int *len_seq, char **src_seq, int num_seq);
int countnode(LIST **list, int *len_seq, int num_seq);
int cleannode(LIST **list, int *len_seq, int num_seq);

int countnode(LIST **list, int *len_seq, int num_seq)
{
	int	i, j, k, l, n;
	POSITION *position;

	n = 0;
	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			if(list[i][j].node && list[i][j].node -> visit == 0)	{
				n ++;
				list[i][j].node -> visit = 1;
			}
		}
	}
	cleannode(list, len_seq, num_seq);
	return(n);
}

void statnode(LIST **list, int *len_seq, char **src_seq, int num_seq)
{
	int	i, j, k, l, n, max_l;
	int	**stat, **nstat, *lax;
	int	width;
	POSITION *position;

	stat = (int **) ckalloc(500 * sizeof(int *));
	nstat = (int **) ckalloc(500 * sizeof(int *));
	lax = (int *) ckalloc(500 * sizeof(int));
	for(i = 0; i < 500; i ++)	{
		stat[i] = (int *) ckalloc(500 * sizeof(int));
		nstat[i] = (int *) ckalloc(500 * sizeof(int));
	}

	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			for(k = 0; k < 5; k ++)	{
				lax[k] = 0;
			}
			if(list[i][j].node && list[i][j].node -> visit == 0)	{
				list[i][j].node -> visit = 1;
				width = singstatnode(list[i][j].node, nstat);
				if(width > 1)	continue;
				position = list[i][j].node -> position;
				while(position)	{
					if(position -> position >= 0)	{
						lax[src_seq[position -> readindex][position -> position]] ++;
					} else	{
						lax[4] ++;
					}
					position = position -> next;
				}
				max_l = lax[0];
				for(k = 1; k < 5; k ++)	{
					if(lax[k] > max_l)	{
						max_l = lax[k];
					}
				}
				if(max_l < 10 && list[i][j].node -> npos < 30)	{
					stat[list[i][j].node -> npos][max_l] ++;
				} else if(max_l < 10)	{
					stat[30][max_l] ++;
				} else if(list[i][j].node -> npos < 30)	{
					stat[list[i][j].node -> npos][10] ++;
				} else	{
					stat[30][10] ++;
				}
			}
		}
	}
	printf("----------------------------------------------------------------------------------------\n");
	printf("# of supernodes that contains n nodes and the major node has multiplicity k\n");
	printf("----------------------------------------------------------------------------------------\n");
	printf("major   1      2      3      4      5      6      7      8      9      10+     Total\n");
	printf("----------------------------------------------------------------------------------------\n");
	for(i = 1; i <= 30; i ++)	{
		if(i < 30)	{
			printf("%-7d ", i);
		} else	{
			printf("30+     ", i);
		}
		k = 0;
		for(j = 1; j <= 10; j ++)	{
			printf("%-6d ", stat[i][j]);
			k += stat[i][j];
		}
		printf("%-6d", k);
		printf("\n");
	}
	printf("----------------------------------------------------------------------------------------\n");
	printf("# of supernodes with width k and thickness n\n");
	printf("----------------------------------------------------------------------------------------\n");
	printf("          2         3         4         5         6        7+\n");
	printf("----------------------------------------------------------------------------------------\n");
	for(i = 2; i <= 50; i ++)	{
		if(i < 50)	{
			printf("%-9d ", i);
		} else	{
			printf("50+       ");
		}
		for(j = 2; j <= 7; j ++)	{
			printf("%-9d ", nstat[i][j]);
		}
		printf("\n");
	}
	printf("----------------------------------------------------------------------------------------\n");
	cleannode(list, len_seq, num_seq);
	for(i = 0; i < 500; i ++)	{
		free((void *) stat[i]);
		free((void *) nstat[i]);
	}
	free((void **) stat);
	free((void **) nstat);
	free((void *) lax);
}

int singstatnode(NODES *node, int **nstat)
{
	int	i, j, k, l, n, width, thickness;

	width = countwidth(node);
	thickness = countthickness(node);
	if(width < 7)	{
		if(thickness <= 50)	{
			nstat[thickness][width] ++;
		} else	{
			nstat[50][width] ++;
		}
	} else	{
		if(thickness <= 50)	{
			nstat[thickness][7] ++;
		} else	{
			nstat[50][7] ++;
		}
	}
	return(width);
}

int cleannode(LIST **list, int *len_seq, int num_seq)
{
	int	i, j, k, l, n;

	n = 0;
	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			if(list[i][j].node)	{
				list[i][j].node -> visit = 0;
			}
		}
	}
	return(n);
}
