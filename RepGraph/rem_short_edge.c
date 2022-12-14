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

static int intvcompa(READINTERVAL *a, READINTERVAL *b);
int update_short_cycle(EDGE *edge, NODES **vertex, int num_vertex);
int update_short_cycle_back(EDGE *edge, NODES **vertex, int num_vertex);
void merge_inteval(NODES *vertex, READINTERVAL *readinterval);
void merge_inteval_back(NODES *vertex, READINTERVAL *readinterval);
int update_short(EDGE *edge, NODES **vertex, int num_vertex, int extra, int *len);
void split_merge_inteval(EDGE *edge, READINTERVAL *readinterval, int l1, int l2, int *len);
int rem_short_edge(NODES **vertex, int num_vertex, int *len);

int rem_short_edge(NODES **vertex, int num_vertex, int *len)
{
	int	i, j, k, l, m, n, l1, l2;
	NODES	*begin, *end, *vertex_new, *begin_bal, *end_bal;
	EDGE	*edge, *bal_edge;

	do	{
		n = num_vertex;
		for(i = 0; i < num_vertex; i ++)	{
		   for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			bal_edge = edge -> bal_edge;
			begin = edge -> begin;
			end = edge -> end;
			begin_bal = bal_edge -> begin;
			end_bal = bal_edge -> end;
			l = 0;
			for(m = 0; m < edge -> multip; m ++)	{
				if(l < edge -> readinterval[m].length)	{
					l = edge -> readinterval[m].length;
				}
			}
/*	If the longest read inteval is shorter than SHORTLEG,
	proceed with suppression	*/
			if(l < SHORTLEG)	{
				if(begin != end)	{
					num_vertex = update_short(edge, vertex, num_vertex, 1, len);
					if(bal_edge != edge)	{
						num_vertex = update_short(bal_edge, vertex, num_vertex, 0, len);
						vertex[num_vertex - 2] -> bal_node = vertex[num_vertex - 1];
						vertex[num_vertex - 1] -> bal_node = vertex[num_vertex - 2];
						free((void *) bal_edge -> readinterval);
						free((void *) bal_edge);
						begin_bal -> num_lastedge = begin_bal -> num_nextedge = 0;
						end_bal -> num_lastedge = end_bal -> num_nextedge = 0;
					} else	{
						vertex[num_vertex - 1] -> bal_node = vertex[num_vertex - 1];
					}
					free((void *) edge -> readinterval);
					free((void *) edge);
					begin -> num_lastedge = begin -> num_nextedge = 0;
					end -> num_lastedge = end -> num_nextedge = 0;
/*
				} else if(edge -> multip < 10)	{
*/
				} else	{
					num_vertex = update_short_cycle(edge, vertex, num_vertex);
					if(bal_edge != edge)	{
						num_vertex = update_short_cycle_back(bal_edge, vertex, num_vertex);
						erasedge(bal_edge);
					}
					erasedge(edge);
				}
			}
		  }
		}
		num_vertex = merge_graph(vertex, num_vertex);
	} while(num_vertex < n);

	return(num_vertex);
}

int update_short_cycle(EDGE *edge, NODES **vertex, int num_vertex)
{
	int	i, j, k, l, n;

	qsort(edge -> readinterval, edge -> multip, sizeof(READINTERVAL), (void *) intvcompa);
	for(i = edge -> multip - 1; i >= 0; i --)	{
		merge_inteval(edge -> begin, &(edge -> readinterval[i]));
	}
	return(num_vertex);
}

int update_short_cycle_back(EDGE *edge, NODES **vertex, int num_vertex)
{
	int	i, j, k, l, n;

	qsort(edge -> readinterval, edge -> multip, sizeof(READINTERVAL), (void *) intvcompa);
	for(i = 0; i < edge -> multip; i ++)	{
		merge_inteval_back(edge -> end, &(edge -> readinterval[i]));
	}
	return(num_vertex);
}

void merge_inteval(NODES *vertex, READINTERVAL *readinterval)
{
	int	i, j, k, l, n;
	EDGE	*edge, *bal_edge;

	for(i = 0; i < vertex -> num_lastedge; i ++)	{
		edge = vertex -> lastedge[i];
		for(j = 0; j < edge -> multip; j ++)	{
			if(edge -> readinterval[j].begin + edge -> readinterval[j].length - 1 ==
			   readinterval -> begin && edge -> readinterval[j].eq_read == readinterval -> eq_read)	{
/*
				edge -> length += readinterval -> length - 1;
*/
				edge -> readinterval[j].length += readinterval -> length - 1;
				return;
			}
		}
	}
	printf("Cycle edge read inteval not found: %d %d %d %d.\n",
		 readinterval -> begin, readinterval -> length, readinterval -> eq_read, readinterval -> offset);
	for(i = 0; i < vertex -> num_lastedge; i ++)	{
		edge = vertex -> lastedge[i];
		for(j = 0; j < edge -> multip; j ++)	{
			printf("i %d readinterval %d %d %d %d\n", i, edge -> readinterval[j].eq_read,
				edge -> readinterval[j].begin, edge -> readinterval[j].length,
				edge -> readinterval[j].offset);
		}
	}
	exit(-1);
}

void merge_inteval_back(NODES *vertex, READINTERVAL *readinterval)
{
	int	i, j, k, l, n;
	EDGE	*edge, *bal_edge;

	for(i = 0; i < vertex -> num_nextedge; i ++)	{
		edge = vertex -> nextedge[i];
		for(j = 0; j < edge -> multip; j ++)	{
			if(edge -> readinterval[j].begin == readinterval -> begin + readinterval -> length - 1
			 && edge -> readinterval[j].eq_read == readinterval -> eq_read)	{
/*
				edge -> length += (readinterval -> length - 1);
*/
				edge -> readinterval[j].begin -= (readinterval -> length - 1);
				edge -> readinterval[j].length += readinterval -> length - 1;
				return;
			}
		}
	}
	printf("Cycle edge back read inteval not found: %d %d %d %d.\n",
		 readinterval -> begin, readinterval -> length, readinterval -> eq_read, readinterval -> offset);
	exit(-1);
}

int update_short(EDGE *edge, NODES **vertex, int num_vertex, int extra, int *len)
{
	int	i, j, k, l, n, l1, l2;
	int	num_lastedge, num_nextedge;
	NODES	*begin, *end, *vertex_new;

	begin = edge -> begin;
	end = edge -> end;
	for(i = 0; i < edge -> multip; i ++)	{
		l = edge -> readinterval[i].length / 2;
		if(edge -> readinterval[i].length % 2 == 1)	{
			split_merge_inteval(edge, &(edge -> readinterval[i]), l, l, len);
		} else	{
			split_merge_inteval(edge, &(edge -> readinterval[i]), l - 1 + extra,
				 edge -> readinterval[i].length - l - extra, len);
		}
	}
	vertex_new = (NODES *) ckalloc(1 * sizeof(NODES));
	n = begin -> num_lastedge + end -> num_lastedge;
	vertex_new -> lastedge = (EDGE **) ckalloc(n * sizeof(EDGE *));
	n = begin -> num_nextedge + end -> num_nextedge;
	vertex_new -> nextedge = (EDGE **) ckalloc(n * sizeof(EDGE *));
	num_lastedge = num_nextedge = 0;
	for(i = 0; i < begin -> num_lastedge; i ++)	{
		vertex_new -> lastedge[num_lastedge ++] = begin -> lastedge[i];
		begin -> lastedge[i] -> end = vertex_new;
	}
	for(i = 0; i < end -> num_lastedge; i ++)	{
		if(end -> lastedge[i] != edge)	{
			vertex_new -> lastedge[num_lastedge ++] = end -> lastedge[i];
			end -> lastedge[i] -> end = vertex_new;
		}
	}
	for(i = 0; i < begin -> num_nextedge; i ++)	{
		if(begin -> nextedge[i] != edge)	{
			vertex_new -> nextedge[num_nextedge ++] = begin -> nextedge[i];
			begin -> nextedge[i] -> begin = vertex_new;
		}
	}
	for(i = 0; i < end -> num_nextedge; i ++)	{
		vertex_new -> nextedge[num_nextedge ++] = end -> nextedge[i];
		end -> nextedge[i] -> begin = vertex_new;
	}
	vertex_new -> num_lastedge = num_lastedge;
	vertex_new -> num_nextedge = num_nextedge;

	vertex[num_vertex] = vertex_new;
	return(num_vertex + 1);
}

void split_merge_inteval(EDGE *edge, READINTERVAL *readinterval, int l1, int l2, int *len)
{
	int	i, j, k, l;
	NODES	*begin, *end;
	EDGE	*tmpedge;

	begin = edge -> begin;
	end = edge -> end;
/*
printf("Read inteval: %d %d %d %d.\n",
	 readinterval -> begin, readinterval -> length, readinterval -> eq_read, readinterval -> offset);
*/
	if(readinterval -> begin > 0)	{
	  for(i = 0; i < begin -> num_lastedge; i ++)	{
		tmpedge = begin -> lastedge[i];
		for(j = 0; j < tmpedge -> multip; j ++)	{
/*
printf("1 Read inteval: %d %d %d %d.\n",
 tmpedge -> readinterval[j].begin, tmpedge -> readinterval[j].length,
 tmpedge -> readinterval[j].eq_read, tmpedge -> readinterval[j].offset);
*/
			if(tmpedge -> readinterval[j].begin + tmpedge -> readinterval[j].length - 1 ==
			   readinterval -> begin && tmpedge -> readinterval[j].eq_read == readinterval -> eq_read)	{
				if(readinterval -> begin + readinterval -> length - 1 == len[readinterval -> eq_read])	{
/*
					tmpedge -> length += readinterval -> length - 1;
*/
					tmpedge -> readinterval[j].length += readinterval -> length - 1;
					continue;
				} else	{
/*
					tmpedge -> length += l1;
*/
					tmpedge -> readinterval[j].length += l1;
					break;
				}
			}
		}
		if(j < tmpedge -> multip)		break;
	  }
	  if(i == begin -> num_lastedge)	{
		printf("1 Read inteval not found: %d %d %d %d.\n",
			 readinterval -> begin, readinterval -> length, readinterval -> eq_read, readinterval -> offset);
		exit(-1);
	  }
	}

	if(readinterval -> begin + readinterval -> length - 1 < len[readinterval -> eq_read])	{
	  for(i = 0; i < end -> num_nextedge; i ++)	{
		tmpedge = end -> nextedge[i];
		for(j = 0; j < tmpedge -> multip; j ++)	{
/*
printf("2 Read inteval: %d %d %d %d.\n",
 tmpedge -> readinterval[j].begin, tmpedge -> readinterval[j].length,
 tmpedge -> readinterval[j].eq_read, tmpedge -> readinterval[j].offset);
*/
			if(tmpedge -> readinterval[j].begin == readinterval -> begin + readinterval -> length - 1 &&
			   tmpedge -> readinterval[j].eq_read == readinterval -> eq_read)	{
				if(readinterval -> begin == 0)	{
/*
					tmpedge -> length += readinterval -> length - 1;
*/
					tmpedge -> readinterval[j].begin -= (readinterval -> length - 1);
					tmpedge -> readinterval[j].length += readinterval -> length - 1;
					continue;
				} else	{
/*
					tmpedge -> length += l2;
*/
					tmpedge -> readinterval[j].begin -= l2;
					tmpedge -> readinterval[j].length += l2;
					break;
				}
			}
		}
		if(j < tmpedge -> multip)		break;
	  }
	  if(i == end -> num_nextedge)	{
		printf("2 Read inteval not found: %d %d %d %d.\n",
			 readinterval -> begin, readinterval -> length, readinterval -> eq_read, readinterval -> offset);
		exit(-1);
	  }
	}
}

int intvcompa(READINTERVAL *a, READINTERVAL *b)
{
	if(a -> eq_read > b -> eq_read)	 return(1);
	if(a -> eq_read == b -> eq_read && a -> begin > b -> begin)	 return(1);

	return(-1);
}
