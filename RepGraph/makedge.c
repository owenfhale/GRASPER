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

int makedge(NODES *node, EDGE **edge, int num_edge, LIST **list);
EDGE *newedge(EDGE **midedge, int num_midedge, LIST **list);
int copyreadinterval(READINTERVAL *readinterval, READINTERVAL *readcov);
void sortreadinterval(READINTERVAL *readinterval, int multip);
int readintervalcompar(READINTERVAL *a, READINTERVAL *b);
void sortreadinterval_index(READINTERVAL *readinterval, int multip);
int readintervalcompar_index(READINTERVAL *a, READINTERVAL *b);
void insert_readcov(READINTERVAL **readcov, POSITION *position, int offset);
void sort_nodepos(NODES *node);
int poscompar(POSPOINT *a, POSPOINT *b);

int makedge(NODES *node, EDGE **edge, int num_edge, LIST **list)
{
	int	i, j, k, l;
	int	num_midedge;
	EDGE	**midedge, *nedge;
	NODES	*node1;

	if(node -> num_lastedge == 1 && node -> num_nextedge == 1)	{
		return(num_edge);
	}
	node -> visit = 1;
	for(j = 0; j < node -> num_nextedge; j ++)	{
		num_midedge = 0;
		midedge = (EDGE **) ckalloc(MAX_TMP_LEG * sizeof(EDGE *));
		midedge[num_midedge ++] = node -> nextedge[j];
		node1 = midedge[num_midedge - 1] -> end;
		while(node1 -> num_lastedge == 1 && node1 -> num_nextedge == 1)	{
			midedge[num_midedge ++] = node1 -> nextedge[0];
			node1 = node1 -> nextedge[0] -> end;
		}
		if(num_midedge > 1)	{
			nedge = newedge(midedge, num_midedge, list);
			if(nedge)	{
				nedge -> visit = 1;
				edge[num_edge ++] = nedge;
			} else	{
				printf("No returned edge\n");
				exit(0);
			}
		} else if(midedge[0] -> visit == 0)	{
			nedge = midedge[0];
			nedge -> visit = 1;
			edge[num_edge ++] = nedge;
		}
		free((void **) midedge);
		if(node1 -> visit == 0)	{
			num_edge = makedge(node1, edge, num_edge, list);
		}
	}
	return(num_edge);
}

EDGE *newedge(EDGE **midedge, int num_midedge, LIST **list)
{
	int	i, j, k, l, n;
	READINTERVAL	*readcov, readinterval;
	NODES	*node, *node1;
	EDGE	*tmpedge;
	char	*label;
	POSITION *position, *position1;

/*	Setup the beginning and ending vertics of new edge	*/

	tmpedge = midedge[0];
	if(num_midedge > 1)	{
		tmpedge -> end = midedge[num_midedge - 1] -> end;
		add_lastedge(tmpedge -> end, tmpedge);
	}

/*	Set up the length of the edge	*/

	label = (char *) ckalloc((tmpedge -> multip + 10) * sizeof(char));
	for(i = 1; i < num_midedge; i ++)	{
		for(j = 0; j < midedge[i] -> multip; j ++)	label[j] = 0;
		for(k = 0; k < tmpedge -> multip; k ++)	{
			for(j = 0; j < midedge[i] -> multip; j ++)	{
				if(midedge[i] -> readinterval[j].eq_read == tmpedge -> readinterval[k].eq_read &&
				   midedge[i] -> readinterval[j].begin == tmpedge -> readinterval[k].begin +
						 tmpedge -> length - tmpedge -> readinterval[k].offset - 1)	{
					tmpedge -> readinterval[k].length = tmpedge -> length -
						 tmpedge -> readinterval[k].offset;
					label[j] = 1;
					break;
				}
			}
		}
		tmpedge -> length += midedge[i] -> length - 1;
		for(j = 0; j < midedge[i] -> multip; j ++)	{
			if(label[j] == 0)	{
				readinterval = midedge[i] -> readinterval[j];
				readinterval.offset = i;
				insert_readinterval(tmpedge, readinterval);
			}
		}
	}

/*	Define read coverage of the edge	*/

	for(i = 0; i < tmpedge -> multip; i ++)	{
		tmpedge -> readinterval[i].length = tmpedge -> length;
	}

	if(tmpedge -> multip == 0)	{
		printf("No read cover: tmpedge %d length %d.\n", tmpedge, num_midedge);
		exit(-1);
	}
	sortreadinterval(tmpedge -> readinterval, tmpedge -> multip);

/*	Delete nodes on the path	*/

	for(i = 1; i < num_midedge; i ++)	{
		node = midedge[i] -> begin;
		position = node -> position;
		while(position)	{
			k = position -> readindex;
			l = position -> position;
			if(l >= 0)	{
				list[k][l].node = (NODES *) NULL;
			}
			position = position -> next;
		}
	}
	if(i > 1)	{
		k = searcherase(midedge[i - 1] -> end -> lastedge, midedge[i - 1],
				 midedge[i - 1] -> end -> num_lastedge);
		eraselast(midedge[i - 1] -> end, k);
	}
	for(i = 1; i < num_midedge; i ++)	{
		node = midedge[i] -> begin;
		if(node)	 {
			free_nodes(node);
		}
		free((void *) midedge[i]);
	}
	free((void *) label);
	return(tmpedge);
}

void insert_readcov(READINTERVAL **readcov, POSITION *position, int offset)
{
	int	i, j, k, l;
	READINTERVAL	*readcov0, *readcov_last, *readcov_new;
	int	readindex, pos;

	if(position -> position < 0)	{
		return;
	}

	readcov0 = readcov[0];
	readcov_last = (READINTERVAL *) NULL;
	while(readcov0)	{
		if(readcov0 -> eq_read == position -> readindex)	{
			if(offset != readcov0 -> offset)	{
				readcov0 -> cov = 1;
			}
			if(position -> position < readcov0 -> begin)	{
				readcov0 -> length += readcov0 -> begin - position -> position;
				readcov0 -> offset = offset;
				readcov0 -> begin = position -> position;
			} else if(position -> position > readcov0 -> begin + readcov0 -> length - 1)	{
				readcov0 -> length = position -> position - readcov0 -> begin + 1;
			}
			readcov[0] = readcov0;
			return;
		} else if(readcov0 -> eq_read > position -> readindex)	{
			break;
		}
		readcov_last = readcov0;
		readcov0 = readcov0 -> next;
	}
	readcov_new = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
	readcov_new -> eq_read = position -> readindex;
	readcov_new -> begin = position -> position;
	readcov_new -> length = 1;
	readcov_new -> offset = offset;
	if(readcov_last) 	{
		readcov_last -> next = readcov_new;
	} else	{
		readcov[1] = readcov_new;
	}
	readcov_new -> next = readcov0;
	readcov[0] = readcov_new;
}

int copyreadinterval(READINTERVAL *readinterval, READINTERVAL *readcov)
{
	int	i, n;

	i = 0;
	while(readcov)	{
		if(readcov -> cov == 1)	{
			readinterval[i] = *readcov;
			readinterval[i ++].next = (READINTERVAL *) NULL;
		}
		readcov = readcov -> next;
	}
	return(i);
}

void sortreadinterval(READINTERVAL *readinterval, int multip)
{
	int	n;

	qsort((void *) readinterval, multip, sizeof(READINTERVAL), (void *) readintervalcompar);
}

int readintervalcompar(READINTERVAL *a, READINTERVAL *b)
{
	if(a -> offset > b -> offset)	return(1);
	if(a -> offset < b -> offset)	return(-1);
	return(0);
}

void sortreadinterval_index(READINTERVAL *readinterval, int multip)
{
	int	n;

	qsort((void *) readinterval, multip, sizeof(READINTERVAL), (void *) readintervalcompar_index);
}

int readintervalcompar_index(READINTERVAL *a, READINTERVAL *b)
{
	if(a -> eq_read > b -> eq_read)	return(1);
	if(a -> eq_read < b -> eq_read)	return(-1);
	return(0);
}

void sort_nodepos(NODES *node)
{
	int	i, j, k;
	POSITION *position;
	POSPOINT *tmppos;

	tmppos = (POSPOINT *) ckalloc(node -> npos * sizeof(POSPOINT));
	i = 0;
	position = node -> position;
	while(position)	{
		tmppos[i ++].position = position;
		position = position -> next;
	}
	if(i != node -> npos)	{
		printf("npos not equal %d %d\n", node -> npos, i);
		exit(0);
	}
	qsort((void *) tmppos, node -> npos, sizeof(POSPOINT), (void *) poscompar);
	node -> position = tmppos[0].position;
	for(i = 0; i < node -> npos - 1; i ++)	{
		tmppos[i].position -> next = tmppos[i + 1].position;
	}
	tmppos[i].position -> next = (POSITION *) NULL;
	free((void *) tmppos);
}

int poscompar(POSPOINT *a, POSPOINT *b)
{
	if(a -> position -> readindex > b -> position -> readindex)	 return(1);
	if(a -> position -> readindex < b -> position -> readindex)	 return(-1);
	return(0);
}
