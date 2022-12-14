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
#include <extfunc.h>

int readins(int *insertpos, int *insertlen, FILE *fp);

int readins(int *insertpos, int *insertlen, FILE *fp)
{
	int	i, j, k, n;
	char	str[300];

	n = 0;
	while(fgets(str, 290, fp))	{
		if(str[0] == '#')	continue;
		sscanf(str, "%d%d%*s", &insertpos[n], &insertlen[n]);
		n ++;
	}
	return(n);
}
