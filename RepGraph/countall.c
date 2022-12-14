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
#include <math.h>
#include <perdef.h>
#include <param.h>
#include <extvab.h>
#include <extfunc.h>

char resfile[100];
char chriname[100];

int min_length, nchro;

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	char	name[1000], temp[1000], c, str[100];
	int	nv, ne, nst, nt, nseg, nsub;
	FILE	*fp;

	sprintf(resfile, "split/watch.%s.simp", argv[1]);
/*	Input chromsomal information	*/
	fp = ckopen(resfile, "r");
	while(fgets(str, 300, fp))	{
		if(!strncmp(str, "Vertices", 8))	{
			fgets(str, 300, fp);
			fgets(str, 300, fp);
			sscanf(str, "%d%d%*d%*d%d%d", &nv, &ne, &nt, &nst);
			break;
		}
	}
	fclose(fp);

	sprintf(name, "count_rep -l chromosome.length <split/%s.txt >tmpfile", argv[1]);
	system(name);

	fp = ckopen("tmpfile", "r");
	fgets(str, 300, fp);
	sscanf(str, "%d%*s%*s%*s%d", &nseg, &nsub);
	fclose(fp);

	printf("%s	%d	%d	%d	%d	%d	%d\n",
		argv[1], nv, ne, nt, nst, nsub, nseg);
}
