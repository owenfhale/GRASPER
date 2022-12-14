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

ALIGN *free_align(ALIGN *align);
int size_align(ALIGN *align);

ALIGN *free_align(ALIGN *align)
{
	int	i;
	ALIGN	*cl;

	cl = align -> next;
	for(i = 0; i < 2; i ++)
		free((void *) align -> pos[i]);
	if(align -> last)	{
		align -> last -> next = align -> next;
	}
	if(align -> next)	{
		align -> next -> last = align -> last;
	}
	free((void *) align);
	return(cl);
}

int size_align(ALIGN *align)
{
	int	n;
	ALIGN	*cl;

	n = 0;
	cl = align;
	while(cl)	{
		n ++;
		cl = cl -> next;
	}
	return(n);
}
