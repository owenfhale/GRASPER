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

char char2int(char c);
char char2intgen(char c);

char char2int(char c)
{
	int	i, k;

	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = ran_number(4, &idum);
	}

	if(k > 3)	{
		k = ran_number(4, &idum);
	}

	return(k);
}

char char2intgen(char c)
{
	int	i, k;

	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = ran_number(4, &idum);
	}

	return(k);
}
