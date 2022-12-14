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

void outputmatrix(FILE *fp, double *EvoDist, int num_seq, char **repnames);

void outputmatrix(FILE *fp, double *EvoDist, int num_seq, char **repnames)
{
	int	i, j, k, l;

	for(i = 0; i < num_seq; i ++)	{
		fprintf(fp, "Segment %d %s\n", i + 1, repnames[i]);
	}
	fprintf(fp, "--------------------------------------------------\n");
	fprintf(fp, "       ");
	for(i = 0; i < num_seq; i ++)	{
		fprintf(fp, "%-6d ", i + 1);
	}
	fprintf(fp, "\n");
	for(i = 0; i < num_seq; i ++)	{
		fprintf(fp, "%6d ", i + 1);
		for(j = 0; j < num_seq; j ++)	{
			k = numc(i, j);
			if(i == j)	{
				fprintf(fp, "--    ", EvoDist[k]);
			} else	{
				fprintf(fp, "%-6.3f ", EvoDist[k]);
			}
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "--------------------------------------------------\n");
}
