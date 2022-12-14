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

int input_segment(SEGMENT *segment, char **chrname, int num_chro, FILE *fp);

int input_segment(SEGMENT *segment, char **chrname, int num_chro, FILE *fp)
{
	int	k;
	int	num_segment;
	char	str[500], name[100], temp[100];

	num_segment = 0;
	while(fgets(str, 400, fp))	{
		if(!strncmp(str, "chromosome", 10))	continue;
//		sscanf(str, "%s%d%d%*d%*s%s%d%d%d", name, &(segment[num_segment].pos[0]),
		sscanf(str, "%s%d%d%*s%s%d%d%d", name, &(segment[num_segment].pos[0]),
			&(segment[num_segment].pos[1]), temp,
			&(segment[num_segment].src_pos[0]), &(segment[num_segment].src_pos[1]), 
			&(segment[num_segment].length));
		segment[num_segment].pos[0] --;
		segment[num_segment].pos[1] --;
		if(segment[num_segment].src_pos[0] > segment[num_segment].src_pos[1])	{
			k = segment[num_segment].src_pos[1];
			segment[num_segment].src_pos[1] = segment[num_segment].src_pos[0];
			segment[num_segment].src_pos[0] = k;
			segment[num_segment].eq_pos[1] = 1;
		}
		segment[num_segment].chro = findgenname(name, chrname, num_chro);
		segment[num_segment].eq_pos[0] = atoi(&temp[1]);
		num_segment ++;
	}
	return(num_segment);
}
