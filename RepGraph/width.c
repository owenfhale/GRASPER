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

int countwidth(NODES *node);
int countthickness(NODES *node);
void countallnode(LIST **list, int *len_seq, int num_seq);
char check_width(NODES *node1, NODES *node2);

int countwidth(NODES *node)
{
	int	i, j, k, l, n, width;
	int	read, readold;
	int	num_read, *readall, *widthall;
	POSITION	*position;

	sort_nodepos(node);
	readall = (int *) ckalloc(node -> npos * sizeof(int));
	widthall = (int *) ckalloc(node -> npos * sizeof(int));
	position = node -> position;
	num_read = 0;
	k = 1;
	readold = -1;
	while(position)	{
		read = position -> readindex;
		if(read != readold)	{
			readall[num_read] = read;
			widthall[num_read] = k;
			num_read ++;
			if(position -> position >= 0)	{
				k = 1;
			} else	{
				k = 0;
			}
			readold = read;
		} else if(position -> position >= 0)	{
			k ++;
		}
		position = position -> next;
	}
	width = 0;
	for(i = 0; i < num_read; i ++)	{
		if(widthall[i] > width)	{
			width = widthall[i];
		}
	}
	free((void *) readall);
	free((void *) widthall);
	return(width);
}

int countthickness(NODES *node)
{
	int	i, j, k, l, n, thickness;
	int	read, readold;
	int	mini, maxi;
	POSITION	*position, *pos1;

	sort_nodepos(node);
	position = node -> position;
	mini = MAX_TMP_LEG;
	while(position)	{
		read = position -> readindex;
		pos1 = position -> next;
		while(pos1 && pos1 -> readindex == read)	{
			if(mini > abs(position -> position - pos1 -> position) + 1)	{
				mini = abs(position -> position - pos1 -> position) + 1;
			}
			pos1 = pos1 -> next;
		}
		position = position -> next;
	}
	if(mini > MIN_WHIRL_SIZE)	{
		thickness = 1;
	} else	{
		thickness = mini;
	}
	return(thickness);
}

char check_width(NODES *node1, NODES *node2)
{
	int	i, j, k, l, n, thickness;
	int	read1, read2;
	int	mini, maxi;
	POSITION	*position, *pos1, *pos2;

	pos1 = node1 -> position;
	while(pos1)	{
		read1 = pos1 -> readindex;
		pos2 = node2 -> position;
		while(pos2)	{
			read2 = pos2 -> readindex;
			if(read1 == read2 && abs(pos1 -> position - pos2 -> position) < MIN_WHIRL_SIZE)	{
				return(1);
			}
			pos2 = pos2 -> next;
		}
		pos1 = pos1 -> next;
	}
	return(0);
}

void countallnode(LIST **list, int *len_seq, int num_seq)
{
	int	i, j, k, l, n;

	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			if(list[i][j].node && list[i][j].node -> visit == 0)	{
				list[i][j].node -> visit = 1;
				n = countwidth(list[i][j].node);
				list[i][j].node -> num_path = n;
			}
		}
	}
	cleannode(list, len_seq, num_seq);
}
