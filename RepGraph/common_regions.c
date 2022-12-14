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

#define MIN_SEG_LENGTH 1000

char comp_segments(SEGMENT segment1, SEGMENT segment2);
char common_regions(SEGMENT *segment1, int len1, SEGMENT *segment2, int len2, int *coord1, int *coord2);

char common_regions(SEGMENT *segment1, int len1, SEGMENT *segment2, int len2, int *coord1, int *coord2)
{
	int	i, j, k, l, end, min_l, m;
	char	c;

	min_l = 0;
	for(i = 0; i < len1; i ++)	{
		for(j = 0; j < len2; j ++)	{
			c = comp_segments(segment1[i], segment2[j]);
			if(c == 1)	{
				for(k = j + 1, m = i + 1; k < len2 && m < len1; k ++, m ++)	{
					c = comp_segments(segment1[m], segment2[k]);
					if(c != 1)	{
						break;
					}
				}
				l = min(segment1[m - 1].pos[1] - segment1[i].pos[0] + 1,
					segment2[k - 1].pos[1] - segment2[j].pos[0] + 1);
				if(l > min_l)	{
					coord1[0] = segment1[i].chro;
					coord1[1] = segment1[i].pos[0] + 1;
					coord1[2] = segment1[m - 1].pos[1] + 1;
					coord2[0] = segment2[j].chro;
					coord2[1] = segment2[j].pos[0] + 1;
					coord2[2] = segment2[k - 1].pos[1] + 1;
					min_l = l;
				}
			}
		}
		for(j = 0; j < len2; j ++)	{
			c = comp_segments(segment1[i], segment2[j]);
			if(c == -1)	{
				for(k = j - 1, m = i + 1; k >= 0 && m < len2; k --, m ++)	{
					c = comp_segments(segment1[m], segment2[k]);
					if(c != -1)	{
						break;
					}
				}
				l = min(segment1[m - 1].pos[1] - segment1[i].pos[0] + 1,
					segment2[j].pos[1] - segment2[k + 1].pos[0] + 1);
				if(l > min_l)	{
					coord1[0] = segment1[i].chro;
					coord1[1] = segment1[i].pos[0] + 1;
					coord1[2] = segment1[m - 1].pos[1] + 1;
					coord2[0] = segment2[j].chro;
					coord2[1] = segment2[j].pos[1] + 1;
					coord2[2] = segment2[k + 1].pos[0] + 1;
					min_l = l;
				}
			}
		}
	}
	if(min_l > MIN_SEG_LENGTH)	{
		return(1);
	} else	{
		return(0);
	}
}

char comp_segments(SEGMENT segment1, SEGMENT segment2)
{
	int	i;

	if(segment1.eq_pos[1] == segment2.eq_pos[1] &&
	   segment1.eq_pos[0] == segment2.eq_pos[0])	{
		return(1);
	} else if(segment1.eq_pos[1] != segment2.eq_pos[1] &&
	   segment1.eq_pos[0] == segment2.eq_pos[0])	{
		return(-1);
	} else	{
		return(0);
	}
}
