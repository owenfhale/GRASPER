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
#include <extfunc.h>

int getchseq(char *seq, int *coord, int dir, char **chrname);

int getchseq(char *seq, int *coord, int dir, char **chrname)
{
	int	i, j, k, l, m, n;
	char	**src_seq, **src_name;
	int	num_seq, len_seq[10];
	char	temp[100];
	FILE	*fp;

	src_seq = (char **) ckalloc(2 * sizeof(char *));
	src_name = (char **) ckalloc(2 * sizeof(char *));
	src_name[0] = (char *) ckalloc(100 * sizeof(char));
	src_name[1] = (char *) ckalloc(100 * sizeof(char));
	sprintf(temp, "../hgs/%s.fa", chrname[coord[0]]);
	fp = ckopen(temp, "r");
	num_seq = readseq1by1(src_seq, src_name, len_seq, fp);
	if(num_seq != 1)	{
		printf("Sequence file error: %s\n", temp);
		exit(0);
	}
	fclose(fp);
	if(coord[1] < coord[2])	{
		for(i = coord[1] - 1, m = 0; i < coord[2]; i ++, m ++)	{
			seq[m] = src_seq[0][i];
		}
	} else	{
		for(i = coord[1] - 1, m = 0; i > coord[2]; i --, m ++)	{
			seq[m] = rev(src_seq[0][i]);
		}
	}
	free((void *) src_seq[0]);
	free((void **) src_seq);
	free((void *) src_name[0]);
	free((void *) src_name[1]);
	free((void **) src_name);
	return(m);
}
