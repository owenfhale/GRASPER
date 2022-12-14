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
#include <perdef.h>
#include <extfunc.h>

int input_reg(FILE *fp, SEGMENT *segment);

int input_reg(FILE *fp, SEGMENT *segment)
{
	int	i, j, k, l, n;
	char	str[400];

	n = 0;
	while(fgets(str, 390, fp))	{
		if(str[0] != '#')	{
			sscanf(str, "%*s%d%*s%d%d%d%d\n", &k, &(segment[n].pos[0]), &(segment[n].pos[1]),
				&(segment[n].eq_pos[0]), &(segment[n].eq_pos[1]));
			segment[n].chro = k - 1;
			n ++;
		}
	}
	return(n);
}
