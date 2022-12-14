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

#include <stdio.h>

int parsestr(char *str, char **items, char *delimiter);

int parsestr(char *str, char **items, char *delimiter)
{
	int	i, j, k, len;
	int	num_item, l_del;

	num_item = 0;
	len = strlen(str);
	l_del = strlen(delimiter);
	i = k = 0;
	while(i < len)	{
		if(strncmp(&str[i], delimiter, l_del))	{ 
			items[num_item][k ++] = str[i];
			i ++;
		} else	{
			for(i ++; i < len; i ++)	{
				if(strncmp(&str[i], delimiter, l_del))	break;
			}
			items[num_item][k] = '\0';
			num_item ++;
			k = 0;
		}
	}
	return(num_item + 1);
}
