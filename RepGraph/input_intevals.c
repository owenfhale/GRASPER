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

extern int min_length;

int input_intevals(FILE *fp, SEGMENT *segment, int *len, int num_seq);
int readchrolist(int *len_seq, char **chrname, int *reallen, FILE *fp);
int read_chro_centro(char **chrname, int *reallen, int *centro, FILE *fp);

int input_intevals(FILE *fp, SEGMENT *segment, int *len, int num_seq)
{
	int	i, j, k, l, dir, n, m, multip, length;
	char	str[400];

	k = n = 0;
	while(fgets(str, 390, fp))	{
		if(!strncmp(str, "EDGE", 4))	{
			sscanf(str, "%*s%d%*s%d%*s%d", &k, &length, &multip);
		} else if(!strncmp(str, "INTV", 4))	{
			sscanf(str, "%*s%d%d%d%d", &dir, &i, &j, &m);
			segment[n].length = length;
			if(dir < num_seq)	{
				segment[n].chro = dir;
				segment[n].src_pos[0] = m;
				segment[n].src_pos[1] = min(m + j - 1, segment[n].length - 1);
				segment[n].pos[0] = i;
				segment[n].pos[1] = i + j - 1;
				segment[n].eq_pos[0] = k;
				segment[n].eq_pos[1] = 0;
			} else	{
				segment[n].chro = dir - num_seq;
				segment[n].src_pos[0] = m;
				segment[n].src_pos[1] = min(m + j - 1, segment[n].length - 1);
				segment[n].pos[0] = len[dir - num_seq] - (i + j);
				segment[n].pos[1] = segment[n].pos[0] + j - 1;
				segment[n].eq_pos[0] = k;
				segment[n].eq_pos[1] = 1;
			}
/*	ignore single multiplicity edges	*/
			if(multip > 1 && length > min_length)	{
				n ++;
			}
		}
	}
	return(n);
}

int readchrolist(int *len_seq, char **chrname, int *reallen, FILE *fp)
{
	int	i;
	char	str[400];

	i = 0;
	while(fgets(str, 395, fp))	{
		sscanf(str, "%*d%d%s%d", &len_seq[i], chrname[i], &reallen[i]);
		i ++;
	}
	return(i);
}

int read_chro_centro(char **chrname, int *reallen, int *centro, FILE *fp)
{
	int	i;
	char	str[400];

	i = 0;
	while(fgets(str, 395, fp))	{
		sscanf(str, "%*d%d%s%d%*d%d", &reallen[i], chrname[i], &centro[2 * i], &centro[2 * i + 1]);
		i ++;
	}
	return(i);
}
