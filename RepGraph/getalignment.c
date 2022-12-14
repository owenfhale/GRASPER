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

int getalignment(int *coord1, int *coord2, int num_seq, int *len_seq, char **seq, int *pos, char **alnseq);
int getal(char *seq1, char *seq2, int stpos, char *alnseq1, char *alnseq2, int lt);
int getpos(char *seq1, char *seq2, int p1, int p2, int length, int *pos, int len_seq);
char chk_cont(int *coord, int *ov_coord);

int getalignment(int *coord1, int *coord2, int num_seq, int *len_seq, char **seq, int *pos, char **alnseq)
{
	int	i, j, k, l, lt;
	char	c1, c2;
	int	length;

	for(i = 0; i < num_seq; i ++)	{
		c1 = chk_cont(coord1, &pos[8 * i + 1]);
		c2 = chk_cont(coord2, &pos[8 * i + 5]);
		if(c1 && c2)	{
			l = coord1[2] - coord1[1] + 1;
			lt = getpos(seq[2 * i], seq[2 * i + 1], coord1[1] - pos[8 * i + 2], coord2[2] - pos[8 * i + 6], l, pos, len_seq[i]);
printf("1 lt %d\n", lt);
			length = getal(seq[2 * i], seq[2 * i + 1], pos[0], alnseq[0], alnseq[1], lt);
			return(length);
		} else	{
			c1 = chk_cont(coord1, &pos[8 * i + 5]);
			c2 = chk_cont(coord2, &pos[8 * i + 1]);
			if(c1 && c2)	{
				l = coord1[2] - coord1[1] + 1;
				lt = getpos(seq[2 * i], seq[2 * i + 1], coord2[1] - pos[8 * i + 2],
					coord1[2] - pos[8 * i + 6], l, pos, len_seq[i]);
printf("2 lt %d\n", lt);
				length = getal(seq[2 * i], seq[2 * i + 1], pos[0], alnseq[1], alnseq[0], lt);
				return(length);
			}
		}
	}
	printf("Not found.\n");
	exit(0);
}

int getal(char *seq1, char *seq2, int stpos, char *alnseq1, char *alnseq2, int lt)
{
	int	i, j, k;
	int	length;

	length = 0;
	for(i = stpos; i < stpos + lt; i ++)	{
		if(seq1[i] != 9 && seq2[i] != 9)	{
			alnseq1[length] = seq1[i];
			alnseq2[length ++] = seq2[i];
		}
	}
	return(length);
}

int getpos(char *seq1, char *seq2, int p1, int p2, int length, int *pos, int len_seq)
{
	int	i, j, k, l, lt, k1, k2;

	k1 = k2 = j = 0;
	for(i = 0; i < len_seq; i ++)	{
		if(seq1[i] != 9)	{
			k1 ++;
			j ++;
		}
		if(seq2[i] != 9)	{
			k2 ++;
		}
		if(k1 == p1)	{
			if(k2 == p2)	{
				pos[0] = i;
				pos[1] = i;
			} else	{
				printf("Alignment error\n");
				exit(0);
			}
			j = 0;
		}
		if(j == l)	{
			if(seq1[i] != 9 && seq2[i] != 9)	{
				lt = i - pos[0];
			} 
		}
	}
}

char chk_cont(int *coord, int *ov_coord)
{
	int	i;

	if(coord[0] == ov_coord[0] - 1 && coord[1] >= ov_coord[1] &&
	   coord[2] <= ov_coord[2])	{
		return(1);
	} else	{
		return(0);
	}
}
