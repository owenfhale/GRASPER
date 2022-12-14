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

void initialize(LIST **list, int *len_seq, int num_seq);

void initialize(LIST **list, int *len_seq, int num_seq)
{
	int	i, j, k, l;
	NODES	*vertex1, *vertex2;

	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			vertex1 = (NODES *) ckalloc(1 * sizeof(NODES));
			vertex1 -> npos = 1;
			vertex1 -> position = (POSITION *) ckalloc(1 * sizeof(POSITION));
			vertex1 -> position[0].readindex = i;
			vertex1 -> position[0].position = j;
			list[i][j].node = vertex1;
			vertex2 = (NODES *) ckalloc(1 * sizeof(NODES));
			vertex2 -> npos = 1;
			vertex2 -> position = (POSITION *) ckalloc(1 * sizeof(POSITION));
			vertex2 -> position[0].readindex = i + num_seq;
			vertex2 -> position[0].position = len_seq[i] - j - 1;
			list[num_seq + i][len_seq[i] - j - 1].node = vertex2;
			vertex1 -> bal_node = vertex2;
			vertex2 -> bal_node = vertex1;
		}
	}
}
