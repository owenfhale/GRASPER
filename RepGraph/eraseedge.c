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

void erasedge(EDGE *edge);

void erasedge(EDGE *edge)
{
	int	i, j, k, l, n, m;
	NODES	*begin, *end;

	begin = edge -> begin;
	end = edge -> end;
	k = searcherase(begin -> nextedge, edge, begin -> num_nextedge);
	erasenext(begin, k);
	k = searcherase(end -> lastedge, edge, end -> num_lastedge);
	eraselast(end, k);

	free((void *) edge -> readinterval);
	free((void *) edge);
}
