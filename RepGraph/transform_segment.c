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

void transform_segment(SEGMENT *segment, int num_segment, SEGMENT *insertreg, int num_insertreg);
void new_segment(int chro, long *pos, SEGMENT *insertreg, int num);

void transform_segment(SEGMENT *segment, int num_segment, SEGMENT *insertreg, int num_insertreg)
{
	int	i, j, k, l;

	for(i = 0; i < num_segment; i ++)	{
		new_segment(segment[i].chro, segment[i].pos, insertreg, num_insertreg);
	}
}

void new_segment(int chro, long *pos, SEGMENT *insertreg, int num)
{
	int	i, j, k, l, n, n1, n2;

/*
	for(k = 0; k < num; k ++)	{
		if(insertreg[k].chro != chro)	continue;
		printf("%d %d %d %d\n", chro, insertreg[k].chro, insertreg[k].pos[0], insertreg[k].pos[1]);
		getchar();
	}
*/

	n1 = n2 = 0;
	for(k = 0; k < num; k ++)	{
		if(insertreg[k].chro != chro)	continue;
		if(k == 0 || insertreg[k - 1].chro != chro)	{
			if(pos[0] < insertreg[k].pos[0])	{
				n1 = 1;
				n2 = 1;
			}
			if(n1 == 1 && n2 == 1)	{
				return;
			}
		}
		if(k < num - 1 && pos[0] > insertreg[k].pos[1] && pos[0] < insertreg[k + 1].pos[0])	{
			pos[0] = pos[0] - insertreg[k].pos[1] + insertreg[k].eq_pos[1];
			pos[1] = pos[1] - insertreg[k].pos[1] + insertreg[k].eq_pos[1];
			n2 = 1;
			n1 = 1;
		} else if(k > 0 && pos[0] > insertreg[k - 1].pos[1] && pos[0] < insertreg[k].pos[0])	{
			pos[0] = pos[0] - insertreg[k - 1].pos[1] + insertreg[k - 1].eq_pos[1];
			pos[1] = pos[1] - insertreg[k - 1].pos[1] + insertreg[k - 1].eq_pos[1];
			n2 = 1;
			n1 = 1;
		} else if(pos[0] >= insertreg[k].pos[0] && pos[0] <= insertreg[k].pos[1])	{
			pos[0] = insertreg[k].eq_pos[1];
			pos[1] = pos[1] - insertreg[k].pos[1] + insertreg[k].eq_pos[1];
			n2 = 1;
			n1 = 1;
		}
/*
		if(k < num - 1 && pos[1] > insertreg[k].pos[1] && pos[1] <= insertreg[k + 1].pos[0])	{
			pos[1] = pos[1] - insertreg[k].pos[1] + insertreg[k].eq_pos[1];
			n2 = 1;
		} else if(k > 0 && pos[1] > insertreg[k - 1].pos[1] && pos[1] <= insertreg[k].pos[0])	{
			pos[1] = pos[1] - insertreg[k - 1].pos[1] + insertreg[k - 1].eq_pos[1];
			n2 = 1;
		} else if(pos[1] >= insertreg[k].pos[0] && pos[1] <= insertreg[k].pos[1])	{
			pos[1] = insertreg[k].eq_pos[0];
			n2 = 1;
		}
*/
		if(n1 == 1 && n2 == 1)	{
			return;
		}
		if(k == num - 1 || insertreg[k + 1].chro != chro)	{
			if(pos[0] >= insertreg[k].pos[1])	{
				pos[0] = pos[0] - insertreg[k].pos[1] + insertreg[k].eq_pos[1];
				n1 = 1;
			}
			if(pos[1] >= insertreg[k].pos[1])	{
				pos[1] = pos[1] - insertreg[k].pos[1] + insertreg[k].eq_pos[1];
				n2 = 1;
			}
			if(n1 == 1 && n2 == 1 && pos[0] <= pos[1])	{
				return;
			}
		}
	}
	printf("num %d\n", num);
	printf("k %d chro %d    %d %d\n %d %d\n %d %d\n", k, insertreg[k - 1].chro,insertreg[k].pos[0], insertreg[k].pos[1],
		insertreg[k - 1].pos[0], insertreg[k - 1].pos[1],
		insertreg[k + 1].pos[0], insertreg[k + 1].pos[1]);
	printf("pos chro %d %d %d n1 %d n2 %d\n", chro, pos[0], pos[1], n1, n2);
	getchar();
}
