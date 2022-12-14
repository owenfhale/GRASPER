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
#include "extfunc.h"

void prt_line(FILE *fp, LINE *line, int num_line, int linewidth, double *label);
void prt_box(FILE *fp, int point[][2], double *label, int color, char mark, int index);

void prt_line(FILE *fp, LINE *line, int num_line, int linewidth, double *label)
{
	int	i, j;

	fprintf(fp, "gsave\n");
	fprintf(fp, "newpath\n");
	for(i = 0; i < num_line / 4; i ++)	{
		for(j = 0; j < 4; j ++)	{
			fprintf(fp, "%4d%4d m", line[4 * i + j].bgnpoint[0],
				line[4 * i + j].bgnpoint[1]);
			fprintf(fp, "%4d%4d l", line[4 * i + j].endpoint[0],
				line[4 * i + j].endpoint[1]);
		}
		fprintf(fp, "\n");
		if(i % 50 == 0)	{
			fprintf(fp, "\n");
			fprintf(fp, "%d setlinewidth\n", linewidth);
			if(label[0] >= 0.0)	{
				fprintf(fp, "%.2f %.2f %.2f setrgbcolor\n", label[0], label[1], label[2]);
			}
			fprintf(fp, "stroke\n\n");
		}
	}
	for(j = 0; j < num_line % 4; j ++)	{
		fprintf(fp, "%4d%4d m", line[4 * i + j].bgnpoint[0],
			line[4 * i + j].bgnpoint[1]);
		fprintf(fp, "%4d%4d l", line[4 * i + j].endpoint[0],
			line[4 * i + j].endpoint[1]);
	}
	fprintf(fp, "\n");
	if(label[0] >= 0)	{
		fprintf(fp, "%f %f %f setrgbcolor\n", label[0], label[1], label[2]);
	}
	fprintf(fp, "%d setlinewidth\n", linewidth);
	fprintf(fp, "stroke\n");
	fprintf(fp, "grestore\n\n");
}

void prt_box(FILE *fp, int point[][2], double *label, int color, char mark, int index)
{
	int	i, j;
	int	coord[2];
	char	temp[10];
	
	fprintf(fp, "gsave\n");
	if(label[0] < 0.0)	{
		fprintf(fp, "1 setlinewidth\n");
		fprintf(fp, "newpath\n");
		fprintf(fp, "%4d%4d m\n", point[0][0], point[0][1]);
		fprintf(fp, "%4d%4d l\n", point[1][0], point[1][1]);
		fprintf(fp, "%4d%4d m\n", point[2][0], point[2][1]);
		fprintf(fp, "%4d%4d l\n", point[3][0], point[3][1]);
		fprintf(fp, "stroke\n");
	} else	{
		if(label[0] >= 0.0)	{
			fprintf(fp, "%.2f %.2f %.2f setrgbcolor\n", label[0], label[1], label[2]);
/*
			fprintf(fp, "%.2f setgray\n", label[0]);
*/
		}
		fprintf(fp, "1 setlinewidth\n");
		fprintf(fp, "newpath\n");
		fprintf(fp, "%4d%4d m\n", point[0][0], point[0][1]);
		fprintf(fp, "%4d%4d l\n", point[1][0], point[1][1]);
		fprintf(fp, "%4d%4d l\n", point[2][0], point[2][1]);
		fprintf(fp, "%4d%4d l\n", point[3][0], point[3][1]);
		fprintf(fp, "closepath\n");
		if(label[0] >= 0.0)	{
			fprintf(fp, "fill\n");
		} else	{
			fprintf(fp, "stroke\n");
		}
		if(mark > 0)	{
			coord[0] = (point[0][0] + point[1][0]) / 2 - 2;
			coord[1] = point[3][1] + 1;
			sprintf(temp, "s%d", index);
			prt_Text(fp, temp, coord, 6, 0);
		}
	}
	fprintf(fp, "grestore\n\n");
}
