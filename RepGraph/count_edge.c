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

#define MIN_LEG 100

int srindex = 18;
extern int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength,
	   *lmultip;;

/* for test only	*/
int thresh1 = 2;
FILE *fp, *fp1;

int count_edge(NODES **vertex, int num_vertex, int **num_pa);
void count_multip(NODES **vertex, int num_vertex);
int count_edge_simp(NODES **vertex, int num_vertex, int **num_pa);
int forwardtanglelen(EDGE *edge, int *ave, int *multip);
int backtanglelen(EDGE *edge, int *ave, int *multip);
int count_tangle(NODES **vertex, int num_vertex, int disttangle[][7]);
int count_bal(NODES *vertex);
void getmaxedge(EDGE *edge, EDGE **maxedge);

int count_edge_simp(NODES **vertex, int num_vertex, int **num_pa)
{
	int	i, j, k, l, m;
	int	l1, l2;

	for(i = 0; i < MAX_BRA; i ++)	{
		for(j = 0; j < MAX_BRA; j ++)	{
			num_pa[i][j] = 0;
		}
	}

	l1 = l2 = 0;
	for(i = 0; i < num_vertex; i ++)	{
		num_pa[vertex[i] -> num_lastedge][vertex[i] -> num_nextedge] ++;
		l1 += vertex[i] -> num_nextedge;
		l2 += vertex[i] -> num_lastedge;
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> visit = 0;
		}
	}
	if(l1 != l2)	{
		printf("edge not balanced %d %d\n", l1, l2);
		exit(-1);
	}
	return(l1);
}

int count_edge(NODES **vertex, int num_vertex, int **num_pa)
{
	int	i, j, k, l, m;
	int	l1, l2;
	int	ave, multip;
	EDGE	**maxedge;
	char temp[100];

	nsuper = 0;

	for(i = 0; i < MAX_BRA; i ++)	{
		for(j = 0; j < MAX_BRA; j ++)	{
			num_pa[i][j] = 0;
		}
	}

	l1 = l2 = 0;
	for(i = 0; i < num_vertex; i ++)	{
		num_pa[vertex[i] -> num_lastedge][vertex[i] -> num_nextedge] ++;
		l1 += vertex[i] -> num_nextedge;
		l2 += vertex[i] -> num_lastedge;
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> visit = 0;
		}
	}
	if(l1 != l2)	{
		printf("edge not balanced %d %d\n", l1, l2);
		exit(-1);
	}

	maxedge = (EDGE **) ckalloc(l1 * sizeof(EDGE *));
	for(i = 0; i < num_vertex; i ++)	{
if(nsuper == srindex)	{
	sprintf(temp, "tmp%d.gvz", nsuper + 1);
	fp = ckopen(temp, "w");
	fprintf(fp, "digraph G {\n");
	fprintf(fp, "\tsize=\"8,8\";\n");
}
/*
*/
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			if(vertex[i] -> nextedge[j] -> visit == 0)	{
				getmaxedge(vertex[i] -> nextedge[j], maxedge);
			}
		}
		if(numtangle[nsuper] > 0)	{
if(nsuper == srindex)	{
	fprintf(fp, "}\n");
	fclose(fp);
}
			nsuper ++;
		}
/*
		if(numtangle[nsuper] > 0)	{
			fclose(fp1);
			fprintf(fp, "}\n");
			fclose(fp);
			srindex = 0;
			sprintf(temp, "tmp%d.intv", nsuper + 1);
			fp1 = ckopen(temp, "w");
			sprintf(temp, "tmp%d-%d.gvz", nsuper + 1, thresh1);
			fp = ckopen(temp, "w");
		}
*/
	}

/*	Output intevals into different components	*/

/*
	for(i = 0; i < num_vertex; i ++)	{
	}
	fclose(fp1);
	fprintf(fp, "}\n");
	fclose(fp);
*/

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> visit = 0;
		}
	}

	for(i = 0; i < nsuper; i ++)	{
		ave = maxedge[i] -> length * maxedge[i] -> start_cover;
		multip = maxedge[i] -> start_cover;
		m = backtanglelen(maxedge[i], &ave, &multip);
		m += maxedge[i] -> length;
		m += forwardtanglelen(maxedge[i], &ave, &multip);
		maxlength[i] = m;
/*
		avelength[i] = ave / multip;
*/
	}

	free((void **) maxedge);
	return(l1);
}

int backtanglelen(EDGE *edge, int *ave, int *multip)
{
	int	i, j, k, l;
	NODES	*vertex;

	vertex = edge -> begin;
	l = 0;
	for(j = 0; j < vertex -> num_lastedge; j ++)	{
		if(vertex -> lastedge[j] -> visit == 0)	{
			vertex -> lastedge[j] -> visit = 1;
		} else	{
			continue;
		}
		if(vertex -> lastedge[j] -> start_cover > 1)	{
			*ave += (vertex -> lastedge[j] -> length - 1) * vertex -> lastedge[j] -> start_cover;
			*multip += vertex -> lastedge[j] -> start_cover;
			k = vertex -> lastedge[j] -> length + backtanglelen(vertex -> lastedge[j], ave, multip) - 1;
			if(k > l)	 l = k;
		}
	}
	return(l);
}

int forwardtanglelen(EDGE *edge, int *ave, int *multip)
{
	int	i, j, k, l;
	NODES	*vertex;

	vertex = edge -> end;
	l = 0;
	for(j = 0; j < vertex -> num_nextedge; j ++)	{
		if(vertex -> nextedge[j] -> visit == 0)	{
			vertex -> nextedge[j] -> visit = 1;
		} else	{
			continue;
		}
		if(vertex -> nextedge[j] -> start_cover > 1)	{
			*ave += (vertex -> nextedge[j] -> length - 1) * vertex -> nextedge[j] -> start_cover;
			*multip += vertex -> nextedge[j] -> start_cover;
			k = vertex -> nextedge[j] -> length + forwardtanglelen(vertex -> nextedge[j], ave, multip) - 1;
			if(k > l)	 l = k;
		}
	}
	return(l);
}

void count_multip(NODES **vertex, int num_vertex)
{
	int	i, j, k, l, n1, n2;
	EDGE	*edge, *edge0;
	int	max_red, l_red;
	NODES	*v;
	int	num_unb, num_unb_copy;

	num_unb_copy = 2 * num_vertex + 1;
	num_unb = num_vertex * 2;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			edge -> start_cover = 1;
		}
	}
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i] -> num_lastedge == 0)	{
			edge = vertex[i] -> nextedge[0];
			v = edge -> end;
			if(v -> num_nextedge == 1)	{
				v -> nextedge[0] -> start_cover = v -> num_lastedge;
			}
		} else if(vertex[i] -> num_nextedge == 0)	{
			edge = vertex[i] -> lastedge[0];
			v = edge -> begin;
			if(v -> num_lastedge == 1)	{
				v -> lastedge[0] -> start_cover = v -> num_nextedge;
			}
		}
	}
	while(num_unb_copy > num_unb)	{
		num_unb_copy = num_unb;
		num_unb = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_nextedge != 0 && vertex[i] -> num_lastedge != 0)	{
				num_unb += abs(count_bal(vertex[i]));
			}
		}
		max_red = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_nextedge != 0 && vertex[i] -> num_lastedge != 0)	{
				for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
					edge = vertex[i] -> nextedge[j];
					n1 = count_bal(edge -> begin);
					n2 = count_bal(edge -> end);
					if(n1 > 0 && n2 < 0)	{
						if(max_red < 2)		{
							edge0 = edge;
							l_red = 1;
							max_red = 2;
						}
					} else if(n1 < 0 && n2 > 0 && edge -> start_cover > 1)	{
						if(max_red < 2)		{
							edge0 = edge;
							l_red = 0;
							max_red = 2;
						}
					}
				}
			}
			if(max_red == 2)	break;
		}
		if(max_red == 2)	{
			if(l_red == 1)	{
				edge0 -> start_cover ++;
			} else	{
				edge0 -> start_cover --;
			}
		}
	}
}

int count_tangle(NODES **vertex, int num_vertex, int disttangle[][7])
{
	int	i, j, k, l, n1, n2;
	EDGE	*edge, *edge0;
	int	num_tangle;

	for(i = 0; i < 8; i ++)	{
		for(j = 0; j < 7; j ++)	{
			disttangle[i][j] = 0;
		}
	}

	num_tangle = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			if(edge -> start_cover <= 6)	{
				n1 = edge -> start_cover - 1;
			} else	{
				n1 = 6;
			}
			if(edge -> length < 500)	{
				n2 = 0;
			} else if(edge -> length < 1000)	{
				n2 = 1;
			} else if(edge -> length < 5000)	{
				n2 = 2;
			} else if(edge -> length < 10000)	{
				n2 = 3;
			} else	{
				n2 = 4;
			}
			disttangle[n1][n2] ++;
			disttangle[7][n2] ++;
			disttangle[n1][5] ++;
			disttangle[7][5] ++;
			if(edge -> start_cover > 1)	{
				if(edge -> length >= MIN_LEG)	{
					num_tangle ++;
				}
			}
		}
	}
	return(num_tangle);
}

int count_bal(NODES *vertex)
{
	int	i, j, k, l, n1, n2;

	n1 = n2 = 0;
	for(j = 0; j < vertex -> num_lastedge; j ++)	{
		n1 += vertex -> lastedge[j] -> start_cover;
	}
	for(j = 0; j < vertex -> num_nextedge; j ++)	{
		n2 += vertex -> nextedge[j] -> start_cover;
	}

	return(n1 - n2);
}

void getmaxedge(EDGE *edge, EDGE **maxedge)
{
	int	i, j, k, l;
	char	c;

/*	Output intevals	*/
/*
	if(edge -> multip >= thresh1 && edge -> visit == 0 && edge -> bal_edge -> visit == 0)	{
		fprintf(fp1, "");
		fprintf(fp1, "EDGE Length %d Multiplicity %d.\n", edge -> length, edge -> multip);
		for(k = 0; k < edge -> multip; k ++)	{
			fprintf(fp1, "INTV %d %d %d %d\n", edge -> readinterval[k].eq_read, edge -> readinterval[k].begin,
				edge -> readinterval[k].length, edge -> readinterval[k].offset);
		}
		srindex ++;
		edge -> subrepeat = edge -> bal_edge -> subrepeat = srindex;
	}
*/

	if(nsuper == srindex)	{
		fprintf(fp, "\t%d -> %d [", edge -> begin -> visit, edge -> end -> visit);
		if(edge -> visit == 0 && edge -> bal_edge -> visit == 0)	{
			fprintf(fp, "label = \"(%d,%d", edge -> length, edge -> multip);
		} else if(edge -> visit == 0)	{
			fprintf(fp, "label = \"(%d,%d", edge -> length, edge -> multip);
		}
		fprintf(fp, ")\"];\n");
	}
/*
*/

	edge -> visit = 1;
//	edge -> bal_edge -> visit = 1;
	if(edge -> start_cover > 1)	{
		maxlength[nsuper] += edge -> length;
		if(edge -> length > avelength[nsuper])	{
			lmultip[nsuper] = edge -> multip;
			avelength[nsuper] = edge -> length;
		}
		if(edge -> length >= MIN_LEG)	{
			numtangle[nsuper] ++;
		}
	} else	{
		return;
	}
	if(edge -> start_cover > maxmultip[nsuper])	{
		maxmultip[nsuper] = edge -> start_cover;
		mlength[nsuper] = edge -> length;
		maxedge[nsuper] = edge;
	}
	for(j = 0; j < edge -> end -> num_nextedge; j ++)	{
		if(edge -> end -> nextedge[j] -> visit == 0)	{
			getmaxedge(edge -> end -> nextedge[j], maxedge);
		}
	}
	for(j = 0; j < edge -> begin -> num_lastedge; j ++)	{
		if(edge -> begin -> lastedge[j] -> visit == 0)	{
			getmaxedge(edge -> begin -> lastedge[j], maxedge);
		}
	}
}
