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
#include <param.h>
#include <extfunc.h>

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	*len_seq, num_seq;
	int	len1, len2, pos1, pos2;
	int	**num_pa;
	char	**src_seq, **src_name;
	char	str[300], dir[3];
	int	num_ins, *insertpos, *insertlen;
	FILE	*fp, *fp1;

	if(argc < 7)	{
		printf("Usage: insencode seq_file reput_file insert_reg_file out_seq_file out_reput_file out_reg_file\n");
		exit(-1);
	}

	len_seq = (int *) ckalloc(2 * sizeof(int));
	src_seq = (char **) ckalloc(2 * sizeof(char *));
	src_name = (char **) ckalloc(1 * sizeof(char *));
	src_name[0] = (char *) ckalloc(100 * sizeof(char));
	fp = ckopen(argv[1], "r");
	num_seq = readseq1by1(src_seq, src_name, len_seq, fp);
	fclose(fp);

	insertpos = (int *) ckalloc(10000 * sizeof(int));
	insertlen = (int *) ckalloc(10000 * sizeof(int));
	fp = ckopen(argv[3], "r");
	num_ins = readins(insertpos, insertlen, fp);
	fclose(fp);

	fp = ckopen(argv[2], "r");
	fp1 = ckopen(argv[5], "w");
	while(fgets(str, 290, fp))	{
		if(str[0] != '#')	{
			sscanf(str, "%d%d%s%d%d", &len1, &pos1, dir, &len2, &pos2);
			k = len1 + pos1 - 1;
			pos1 = reculate_pos(pos1, insertpos, insertlen);
			k = reculate_pos(k, insertpos, insertlen, num_ins);
			len1 = k - pos1 + 1;
			k = len2 + pos2 - 1;
			pos2 = reculate_pos(pos2, insertpos, insertlen);
			k = reculate_pos(k, insertpos, insertlen, num_ins);
			len2 = k - pos2 + 1;
			fprintf(fp1, "%d %d %s %d %d\n", len1, pos1, len2, pos2, dir);
		}
	}
	fclose(fp);
	fclose(fp1);

	fp = ckopen(argv[6], "w");
	l = 0;
	for(i = 0; i < num_ins; i ++)	{
		fprintf(fp, "%d %d\n", insertpos[i] - l, insertlen[i]);
		l += insertlen[i];
	}
	fclose(fp);

	fp = ckopen(argv[4], "w");
	fprintf(fp, ">seq_no_common_repeat\n");
	k = n = 0;
	for(i = 0; i < len_seq[0]; i ++)	{
		while(n < num_ins && i == insertpos[n] + 1)	{
			i += insertlen[n];
			n ++;
		}
		fprintf(fp, "%c", na_name[src_seq[0][i]]);
		if(k % 50 == 49)	{
			fprintf(fp, "\n");
		}
		k ++;
	}
	if(k % 50 != 0)	{
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("Genome length after removal: %d\n", k);
	free((void *) insertpos);
	free((void *) insertlen);
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) src_seq[i]);
	}
	for(i = 0; i < num_seq; i ++)	{
		free((void *) src_name[i]);
	}
	free((void **) src_seq);
	free((void **) src_name);
	free((void *) len_seq);
}

int reculate_pos(int pos, int *insertpos, int *insertlen, int num_ins)
{
	int	i, j, k, l;

	k = l = 0;
	for(i = 0; i < num_ins; i ++)	{
		if(pos < insertpos[i])	{
			return(pos - l);
		} else if(pos < insertpos[i] + insertlen[i])	{
			l += (pos - insertpos[i] + 1);
		} else	{
			l += insertlen[i];
		}
	}
}
