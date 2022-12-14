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


void prt_Text(FILE *fp, char *text, int position[2], int size, int angle);

void prt_Text(FILE *fp, char *text, int position[2], int size, int angle)
{
	fprintf(fp, "gsave\n");
	fprintf(fp, "0 0 0 setcolor\n");
	fprintf(fp, "/Times-Roman findfont %d scalefont setfont\n", size);
	fprintf(fp, "%4d%4d m\n", position[0], position[1]);
	if(angle != 0)	{
/*
		fprintf(fp, "%4d%4d translate\n", position[0], position[1]);
*/
		fprintf(fp, "%d rotate\n", angle);
	}
	fprintf(fp, "(%s) show\n", text);
	fprintf(fp, "grestore\n");
}
