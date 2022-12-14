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

EDGE *insert_edge(NODES *node, NODES *node_next, int read, int pos1, int pos2);
void  insert_readinterval(EDGE *edge, READINTERVAL readinterval_new);
void add_nextedge(NODES *node, EDGE *edge);
void add_lastedge(NODES *node, EDGE *edge);

EDGE *insert_edge(NODES *node, NODES *node_next, int read, int pos1, int pos2)
{
        int     i, j;
        EDGE *edge;
	READINTERVAL	readinterval, *readinterval0;

	readinterval.eq_read = read;
	readinterval.offset = 0;
	readinterval.begin = pos1;
	readinterval.length = pos2 - pos1 + 1;
	readinterval.next = (READINTERVAL *) NULL;

        for(i = 0; i < node -> num_nextedge; i ++)      {
                if(node -> nextedge[i] -> end == node_next && node -> nextedge[i] -> length ==
			readinterval.length)     break;
        }
        if(i < node -> num_nextedge)    {
		insert_readinterval(node -> nextedge[i], readinterval);
		if(readinterval.length > node -> nextedge[i] -> length)	{
			node -> nextedge[i] -> length = readinterval.length;
		}
		return((EDGE *) NULL);
        }

        edge = (EDGE *) ckalloc(1 * sizeof(EDGE));
        edge -> begin = node;
        edge -> end = node_next;
        edge -> multip = 1;
	edge -> readinterval = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
	edge -> readinterval[0] = readinterval;
	edge -> length = pos2 - pos1 + 1;
        add_nextedge(node, edge);
        add_lastedge(node_next, edge);
        return(edge);
}

void  insert_readinterval(EDGE *edge, READINTERVAL readinterval_new)
{
	int	i, j, k, l;
	READINTERVAL	*readinterval0;

	readinterval0 = (READINTERVAL *) ckalloc((edge -> multip + 1) * sizeof(READINTERVAL));
	for(j = 0; j < edge -> multip; j ++)	{
		readinterval0[j] = edge -> readinterval[j];
	}
	free((void *) edge -> readinterval);
	edge -> multip ++;
	edge -> readinterval = (READINTERVAL *) ckalloc(edge -> multip * sizeof(READINTERVAL));
	for(j = 0; j < edge -> multip - 1; j ++)	{
		edge -> readinterval[j] = readinterval0[j];
	}
	edge -> readinterval[j] = readinterval_new;
	free((void *) readinterval0);
}

void add_nextedge(NODES *node, EDGE *edge)
{
        int     i, j, k, l;
        EDGE    **tmpedge;

        tmpedge = (EDGE **) ckalloc((node -> num_nextedge + 1) * sizeof(EDGE *));
        for(i = 0; i < node -> num_nextedge; i ++)      {
                tmpedge[i] = node -> nextedge[i];
                if(edge == tmpedge[i])  {
                        break;
                }
        }
        if(i == node -> num_nextedge)   {
                if(node -> num_nextedge > 0)    free((void **) node -> nextedge);
                node -> num_nextedge ++;
                node -> nextedge = (EDGE **) ckalloc(node -> num_nextedge * sizeof(EDGE *));
                for(i = 0; i < node -> num_nextedge - 1; i ++)  {
                        node -> nextedge[i] = tmpedge[i];
                }
                node -> nextedge[i] = edge;
        }
        free((void **) tmpedge);
}

void add_lastedge(NODES *node, EDGE *edge)
{
        int     i, j, k, l;
        EDGE    **tmpedge;

        tmpedge = (EDGE **) ckalloc((node -> num_lastedge + 1) * sizeof(EDGE *));
        for(i = 0; i < node -> num_lastedge; i ++)      {
                tmpedge[i] = node -> lastedge[i];
                if(edge == tmpedge[i])  {
                        break;
                }
        }
        if(i == node -> num_lastedge)   {
                if(node -> num_lastedge > 0)    free((void **) node -> lastedge);
                node -> num_lastedge ++;
                node -> lastedge = (EDGE **) ckalloc(node -> num_lastedge * sizeof(EDGE *));
                for(i = 0; i < node -> num_lastedge - 1; i ++)  {
                        node -> lastedge[i] = tmpedge[i];
                }
                node -> lastedge[i] = edge;
        }
        free((void **) tmpedge);
}
