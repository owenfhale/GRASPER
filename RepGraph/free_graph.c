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

void free_graph(NODES **vertex, int num_vertex);

void free_graph(NODES **vertex, int num_vertex)
{
	int	i, j;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			if(!(vertex[i] -> nextedge[j] -> seq))	{
				free((void *) vertex[i] -> nextedge[j] -> seq);
			}
			free((void *) vertex[i] -> nextedge[j] -> readinterval);
			free((void *) vertex[i] -> nextedge[j]);
		}
		free((void **) vertex[i] -> lastedge);
		free((void **) vertex[i] -> nextedge);
		free((void *) vertex[i] -> position);
		free((void *) vertex[i]);
	}
}
