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


void prt_Header(FILE *fp, char *title, char *creator, char *creationdate, int pages);

void prt_Header(FILE *fp, char *title, char *creator, char *creationdate, int pages)
{
	fprintf(fp, "%%!PSADOBE-1.0\n");
	fprintf(fp, "%%%%Creator: %s\n", creator);
	fprintf(fp, "%%%%Title: %s\n", title);
	fprintf(fp, "%%%%Creation date: %s\n", creationdate);
	if(pages != 0)
		fprintf(fp, "%%%%Pages: %d\n", pages);
	else
		fprintf(fp, "%%%%Pages: (atend)\n");
	fprintf(fp, "%%%%DocumentFonts: Times-Roman Times-Italic Times-Bold\n");
	fprintf(fp, "%%%%BoundingBox: 0 0 900 900\n");
	fprintf(fp, "%%%%EndComments\n\n");
}
