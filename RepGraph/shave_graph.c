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


int merge_graph(NODES **vertex, int num_vertex);
int shave_graph(NODES **vertex, int num_vertex);
EDGE *merge_vertex(NODES *vertex);
EDGE *new_edge(NODES *vertex, NODES *begin, NODES *end, EDGE *edge1, EDGE *edge2);
void combine_readinterval(EDGE *edge1, EDGE *edge2, EDGE *newedge);


int rem_short_end(NODES **vertex, int num_vertex)
{
	int	i, j, k, l, n, m, nsh;
	int	tot_edge;
	EDGE	*edge, *bal_edge;

/*	Remove end edges	*/

	nsh = 0;
	do	{
		n = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_lastedge != 0)	continue;
			j = 0;
			while(j < vertex[i] -> num_nextedge)	{
				edge = vertex[i] -> nextedge[j];
				if(edge -> end -> num_lastedge > 1 && edge -> length <= SHORTEND)	{
					bal_edge = edge -> bal_edge;
					erasedge(edge);
					if(bal_edge && edge != bal_edge)	{
						erasedge(bal_edge);
						n ++;
					}
					n ++;
				} else	{
					j ++;
				}
			}
		}
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_nextedge != 0)	continue;
			j = 0;
			while(j < vertex[i] -> num_lastedge)	{
				edge = vertex[i] -> lastedge[j];
				if(edge -> begin -> num_nextedge > 1 && edge -> length <= SHORTEND)	{
					bal_edge = edge -> bal_edge;
					erasedge(edge);
					if(bal_edge && edge != bal_edge)	{
						erasedge(bal_edge);
						n ++;
					}
					n ++;
				} else	{
					j ++;
				}
			}
		}
		nsh += n;
		num_vertex = merge_graph(vertex, num_vertex);
	} while(n > 0);

	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}
	printf("%d bulges removed, %d vertics %d edges left.\n",
		nsh, num_vertex, tot_edge);

	return(num_vertex);
}

int shave_graph(NODES **vertex, int num_vertex)
{
	int	i, j, k, l, n, m, n1, reads;
	int	tot_edge;
	int	nbul;
	int	maxmlt, maxl, maxk, multip;
	int	true_multip;
	NODES	*begin, *end, *bal_node;
	EDGE	*edge, *edge1, *edge2, *bal_edge, *bal_edge1, *bal_edge2;

/*	Remove bulges and shave edges linking sources & sinks	*/

	nbul = 0;
	do	{
		m = 0;
		for(i = 0; i < num_vertex; i ++)	{
			for(k = 0; k < vertex[i] -> num_nextedge; k ++)	{
				edge1 = vertex[i] -> nextedge[k];
				bal_edge1 = edge1 -> bal_edge;
				if(abs(edge1 -> length - 503) < 10)	{
					continue;
				}
				j = k + 1;
				while(j < vertex[i] -> num_nextedge)	{
					edge2 = vertex[i] -> nextedge[j];
					bal_edge2 = edge2 -> bal_edge;
					if(abs(edge1 -> length - 503) < 10)	{
						j ++;
						continue;
					}
					if(edge2 -> end == edge1 -> end && edge1 != bal_edge2 && edge2 -> length < MIN_INT &&
					   edge1 -> length < MIN_INT && abs(edge1 -> length - edge2 -> length) < MIN_INT) {
						if(edge1 -> readinterval && edge2 -> readinterval)
							movereadinterval(edge1, edge2);
						erasedge(edge2);
						if(bal_edge1 && bal_edge2 && bal_edge2 != edge2)	{
							movereadinterval(bal_edge1, bal_edge2);
							erasedge(bal_edge2);
							m ++;
						}
						m ++;
					} else if(edge2 -> end == edge1 -> end && edge1 != bal_edge2)	{
						if(edge1 -> length < 10 && edge2 -> length < MIN_OTH_LEN &&
						   edge2 -> length != 503)	{
							movereadinterval(edge1, edge2);
							erasedge(edge2);
							if(bal_edge1 && bal_edge2 && bal_edge2 != edge2)	{
								movereadinterval(bal_edge1, bal_edge2);
								erasedge(bal_edge2);
								m ++;
							}
							m ++;
						} else if(edge2 -> length < 10 && edge1 -> length < MIN_OTH_LEN &&
							  edge1 -> length != 503)	{
							movereadinterval(edge2, edge1);
							erasedge(edge1);
							if(bal_edge1 && bal_edge2 && bal_edge1 != edge1)	{
								movereadinterval(bal_edge2, bal_edge1);
								erasedge(bal_edge1);
								m ++;
							}
							m ++;
							break;
						} else {
							j ++;
						}
					} else	{
						j ++;
					}
				}
			}
		}
		num_vertex = merge_graph(vertex, num_vertex);
		nbul += m;
	} while(m > 0);

	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}
	printf("%d bulges removed, %d vertics %d edges left.\n",
		nbul, num_vertex, tot_edge);

	return(num_vertex);
}

int merge_graph(NODES **vertex, int num_vertex)
{
	int	i, j, k, l;
	int	num_vertex_old;
	NODES	*bal_node, *begin, *end;
	POSITION	*position;
	EDGE	*edge1, *edge2;

	do	{
		num_vertex_old = num_vertex;
		i = 0;
		while(i < num_vertex)	{
			if(vertex[i] -> num_nextedge == 0 && vertex[i] -> num_lastedge == 0)	{
				free_nodes(vertex[i]);
				for(j = i; j < num_vertex - 1; j ++)	{
					vertex[j] = vertex[j + 1];
				}
				num_vertex --;
			} else if(vertex[i] -> num_nextedge == 1 && vertex[i] -> num_lastedge == 1 &&
				  vertex[i] -> nextedge[0] != vertex[i] -> lastedge[0]) {
				bal_node = vertex[i] -> bal_node;
				edge1 = merge_vertex(vertex[i]);
				begin = edge1 -> begin;
				end = edge1 -> end;
				if(bal_node && bal_node != vertex[i])	{
					if(bal_node -> num_lastedge == 1 && bal_node -> num_nextedge == 1)	{
						edge2 = merge_vertex(bal_node);
						if((begin != bal_node && end != bal_node) || begin == end)	{
							edge1 -> bal_edge = edge2;
							edge2 -> bal_edge = edge1;
						} else 	{
							edge2 -> bal_edge = edge2;
						}
					} else	{
						printf("pos3 vertex %d(%d-%d) bal_node %d(%d-%d)\n",
							vertex[i], vertex[i] -> num_lastedge, vertex[i] -> num_nextedge,
							bal_node, bal_node -> num_lastedge, bal_node -> num_nextedge);
						exit(0);
					}
				} else if(bal_node)	{
					edge1 -> bal_edge = edge1;
				}
			} else	{
				i ++;
			}
		}
	} while(num_vertex_old > num_vertex);
	return(num_vertex);
}

EDGE *merge_vertex(NODES *vertex)
{
	int	i, j, k, l, m, n, k1, k2, n1, n2, c, q, p, num;
	NODES	*begin, *end, *b0, *e0;
	int	*del_path;
	EDGE	*lastedge, *nextedge, *newedge, *new1edge;
	double	v, r;

	begin = vertex -> lastedge[0] -> begin;
	end = vertex -> nextedge[0] -> end;
	lastedge = vertex -> lastedge[0];
	nextedge = vertex -> nextedge[0];
	if(lastedge == nextedge)	{
		return(lastedge);
	}
	newedge = new_edge(vertex, begin, end, lastedge, nextedge);
	erasedge(lastedge);
	erasedge(nextedge);
	return(newedge);
}

EDGE *new_edge(NODES *vertex, NODES *begin, NODES *end, EDGE *edge1, EDGE *edge2)
{
	int	i, j, k, l, m;
	EDGE	**tmpedge, *edge;

	edge = (EDGE *) ckalloc(1 * sizeof(EDGE));
	tmpedge = (EDGE **) ckalloc(MAX_BRA * sizeof(EDGE *));

	for(i = 0; i < begin -> num_nextedge; i ++)	{
		tmpedge[i] = begin -> nextedge[i];
	}
	free((void **) begin -> nextedge);
	begin -> nextedge = (EDGE **) ckalloc((begin -> num_nextedge + MAX_BRA) * sizeof(EDGE *));
	for(i = 0; i < begin -> num_nextedge; i ++)	{
		begin -> nextedge[i] = tmpedge[i];
	}
	begin -> nextedge[i] = edge;
	begin -> num_nextedge ++;

	for(i = 0; i < end -> num_lastedge; i ++)	{
		tmpedge[i] = end -> lastedge[i];
	}
	free((void **) end -> lastedge);
	end -> lastedge = (EDGE **) ckalloc((end -> num_lastedge + MAX_BRA) * sizeof(EDGE *));
	for(i = 0; i < end -> num_lastedge; i ++)	{
		end -> lastedge[i] = tmpedge[i];
	}
	end -> lastedge[i] = edge;
	end -> num_lastedge ++;

	edge -> begin = begin;
	edge -> end = end;
	edge -> length = edge1 -> length + edge2 -> length - 1;
	if(edge1 -> seq && edge2 -> seq)	{
		edge -> seq = (char *) ckalloc(edge -> length * sizeof(char));
		for(i = 0; i < edge1 -> length; i ++)	{
			edge -> seq[i] = edge1 -> seq[i];
		}
		for(i = 1; i < edge2 -> length; i ++)	{
			edge -> seq[edge1 -> length - 1 + i] = edge2 -> seq[i];
		}
	}
	combine_readinterval(edge1, edge2, edge);

	free((void **) tmpedge);
	return(edge);
}

void combine_readinterval(EDGE *edge1, EDGE *edge2, EDGE *newedge)
{
	int	i, j, k, l, m, n, c, num1;
	int	reads;
	NODES	*vertex0;
	EDGE	*edge01, *edge02;

	l = edge1 -> length - 1;
	k = edge1 -> multip + edge2 -> multip;
	if(!(edge1 -> readinterval) || !(edge2 -> readinterval))	{
		newedge -> multip = min(edge1 -> multip, edge2 -> multip);
		return;
	}
	newedge -> readinterval = (READINTERVAL *) ckalloc(k * sizeof(READINTERVAL));
	n = 0;
	for(i = 0; i < edge1 -> multip; i ++)	{
		newedge -> readinterval[n ++] = edge1 -> readinterval[i];
	}
	num1 = n;
	for(i = 0; i < edge2 -> multip; i ++)	{
		c = findposition(newedge -> readinterval, num1, edge2 -> readinterval[i].eq_read,
				 edge2 -> readinterval[i].begin);
		if(c >= 0)	{
			newedge -> readinterval[c].length += edge2 -> readinterval[i].length - 1;
		} else	{
			newedge -> readinterval[n] = edge2 -> readinterval[i];
			newedge -> readinterval[n].offset += l;
			n ++;
		}
	}
	newedge -> multip = n;
}

int findposition(READINTERVAL *readinterval, int num, int readindex, int position)
{
	int	i, j, k, l;

	for(i = 0; i < num; i ++)	{
		if(readinterval[i].eq_read == readindex &&
		   readinterval[i].begin + readinterval[i].length - 1 == position)	{
			return(i);
		}
	}
	return(-1);
}
