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
#include <color.h>
#include <perdef.h>
#include <extfunc.h>

#define MAX_LEG 2000
#define MIN_LEG 200

int length_intv(int length);
int decompcolor(int colorindex, double *label);
int output_region(FILE *fp, int *base, char *name, char *chrname, int *range, int scale,
	 int *breaks, int length, int *colorindex, int *dir);

int output_region(FILE *fp, int *base, char *name, char *chrname, int *range, int scale,
	 int *breaks, int length, int *colorindex, int *dir)
{
	int	i, j, k, l, intv;
	int	coord[2], point[4][2];
	LINE	line[2];
	char	temp[100];
	int	mark;
	double	label[5];

/*	print name	*/
	prt_Text(fp, name, base, 6, 0);

/*	print range	*/
	coord[0] = base[0] + 15;
	coord[1] = base[1] + 12;
	prt_Text(fp, chrname, coord, 6, 0);
	coord[0] += 40; 
	sprintf(temp, "%d", range[0]);
	prt_Text(fp, temp, coord, 6, 0);
	sprintf(temp, "%d", range[1] - range[0] + 1);
	coord[0] += 35;
	prt_Text(fp, temp, coord, 6, 0);

/*	print boxes	*/
	point[0][0] = base[0] + 15;
	point[0][1] = base[1];
	point[3][0] = base[0] + 15;
	point[3][1] = base[1] + 5;
	for(i = 1; i < length; i ++)	{
		intv = length_intv(breaks[i] - breaks[i - 1] + 1);
/*
printf("breaks %d %d intv %d\n", breaks[i], breaks[i - 1], intv);
printf("range %d %d\n", range[0], range[1]);
getchar();
*/
		point[1][0] = point[0][0] + intv;
		point[1][1] = point[0][1];
		point[2][0] = point[3][0] + intv;
		point[2][1] = point[3][1];
		mark = decompcolor(colorindex[i - 1], label);
		prt_box(fp, point, label, (colorindex[i - 1] - 1) % total_color, mark, colorindex[i - 1]);
		if(dir[i - 1] == 0)	{
			line[0].bgnpoint[0] = point[1][0];
			line[0].bgnpoint[1] = point[1][1];
			line[0].endpoint[0] = point[1][0] - 3;
			line[0].endpoint[1] = point[1][1] - 3;
			line[1].bgnpoint[0] = point[2][0];
			line[1].bgnpoint[1] = point[2][1];
			line[1].endpoint[0] = point[2][0] - 3;
			line[1].endpoint[1] = point[2][1] + 3;
		} else	{
			line[0].bgnpoint[0] = point[0][0];
			line[0].bgnpoint[1] = point[0][1];
			line[0].endpoint[0] = point[0][0] + 3;
			line[0].endpoint[1] = point[0][1] - 3;
			line[1].bgnpoint[0] = point[3][0];
			line[1].bgnpoint[1] = point[3][1];
			line[1].endpoint[0] = point[3][0] + 3;
			line[1].endpoint[1] = point[3][1] + 3;
		}
		prt_line(fp, line, 2, 2, label);
		if(point[1][0] < 450)	{
			point[0][0] = point[1][0];
			point[0][1] = point[1][1];
			point[3][0] = point[2][0];
			point[3][1] = point[2][1];
		} else	{
			base[1] -= 25;
			point[0][0] = base[0] + 15;
			point[0][1] = base[1];
			point[3][0] = base[0] + 15;
			point[3][1] = base[1] + 5;
		}
	}

	return(point[1][0] + 10);
}

int decompcolor(int colorindex, double *label)
{
	int	i, j, k;
	int	mark;

	if(colorindex == 0)	{
		for(j = 0; j < 3; j ++)
			label[j] = 0;
		return(0);
	}
	i = (colorindex - 1) % total_color;
	mark = (colorindex - 1) / total_color + 1;
	for(j = 0; j < 3; j ++)
		label[j] = color_box[i][j];
	return(mark);
}

int length_intv(int length)
{
	if(length <= 200)	{
		return(10);
	} else if(length <= 500)	{
		return(20);
	} else if(length <= 1000)	{
		return(30);
	} else if(length <= 5000)	{
		return(40);
	} else {
		return(50);
	} 
}
