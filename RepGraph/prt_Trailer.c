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


void restore_scale(FILE *fp);
void prt_Trailer(FILE *fp, int pages);
void showpage(FILE *fp);
void prt_Page(FILE *fp, int page);
void change_scale(FILE *fp, double *scale);
void page_trailer(FILE *fp);

void page_trailer(FILE *fp)
{
	fprintf(fp, "%%%%PageTrailer\n");
}


void prt_Trailer(FILE *fp, int pages)
{
	fprintf(fp, "%%%%Trailer\n");
}



void showpage(FILE *fp)
{
	fprintf(fp, "\nshowpage\n");
}



/*
void prt_Page(FILE *fp, int page, char *pagefont)
*/
void prt_Page(FILE *fp, int page)
{
/*
	fprintf(fp, "%%%%Pagefont: %s\n", pagefont);
*/
	fprintf(fp, "%%%%Page %d %d\n", page, page);
	fprintf(fp, "%%%%PageBoundingBox %d %d %d %d\n", 0, 0, 900, 900);
	fprintf(fp, "%%%%EndPageComments\n\n");
	
}


void change_scale(FILE *fp, double *scale)
{
	fprintf(fp, "gsave\n");
	fprintf(fp, "%2.1f %2.1f scale\n", scale[0], scale[1]);
}

void restore_scale(FILE *fp)
{
	fprintf(fp, "grestore\n");
}
