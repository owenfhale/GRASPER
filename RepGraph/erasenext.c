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

void erasenext(NODES *vertex, int n);
void eraselast(NODES *vertex, int n);
int searcherase(EDGE **edge, EDGE *e, int num);

void erasenext(NODES *vertex, int n)
{
	int	i;

	for(i = n; i < vertex -> num_nextedge - 1; i ++)	{
		vertex -> nextedge[i] = vertex -> nextedge[i + 1];
	}
	vertex -> num_nextedge --;
}

void eraselast(NODES *vertex, int n)
{
	int	i;

	for(i = n; i < vertex -> num_lastedge - 1; i ++)	{
		vertex -> lastedge[i] = vertex -> lastedge[i + 1];
	}
	vertex -> num_lastedge --;
}

int searcherase(EDGE **edge, EDGE *e, int num)
{
	int	i, n;

	n = 0;
	for(i = 0; i < num; i ++)	{
		if(edge[i] == e)	{
			n = i;
			break;
		}
	}

	if(i == num)	{
		printf("e %d %d %d %d %d\n", e, e -> begin, e -> end, e -> length, e -> bal_edge);
		for(i = 0; i < num; i ++)	{
			printf("%d %d %d %d %d\n", i, edge[i], edge[i] -> length, edge[i] -> begin, edge[i] -> end);
		}
		printf("Not found\n");
		getchar();
		n = -1;
	}


	return(n);
}
