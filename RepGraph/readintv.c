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

int readintv(char **seq, FILE *fp, int *len_seq, char **seqname);
int readtargintv(FILE *fp, int **intvlist, char **regname, int *len_chro, char **chrname, int num_chr);

int readintv(char **seq, FILE *fp, int *len_seq, char **seqname)
{
	int	i, j, k, l, m, n;
	int	len1, len2;
	char	*seq1, *seq2;
	int	num_seq;
	char	str[500];

	seq1 = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));
	seq2 = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));
	num_seq = k = 0;
	len1 = len2 = 0;
	while(fgets(str, 400, fp))	{
		if(str[0] == '>')	{
			if(k == 0 && len1 > 0 && len2 > 0)	{
				if(len1 != len2)	{
					printf("Length no equal: %d %d\n", len1, len2);
					printf("str %s\n", str);
					exit(0);
				}
				seq[2 * num_seq] = (char *) ckalloc(len1 * sizeof(char));
				seq[2 * num_seq + 1] = (char *) ckalloc(len1 * sizeof(char));
				for(i = 0; i < len1; i ++)	{
					seq[2 * num_seq][i] = seq1[i];
					seq[2 * num_seq + 1][i] = seq2[i];
				}
				len_seq[num_seq ++] = len1;
				len1 = len2 = 0;
			}
			l = strlen(str);;
			for(i = 1; i < l - 1; i ++)	{
				seqname[2 * num_seq + k][i - 1] = str[i];
			}
			seqname[2 * num_seq + k][l - 1] = '\0';
			if(k == 1)	k = 0;
			else		k = 1;
		} else	{
			if(k == 1)	{
				len1 = getseq(str, seq1, len1);
			} else	{
				len2 = getseq(str, seq2, len2);
			}
		}
	}

	if(len1 != len2)	{
		printf("Length no equal: %d %d\n", len1, len2);
		exit(0);
	}
	seq[2 * num_seq] = (char *) ckalloc(len1 * sizeof(char));
	seq[2 * num_seq + 1] = (char *) ckalloc(len1 * sizeof(char));
	for(i = 0; i < len1; i ++)	{
		seq[2 * num_seq][i] = seq1[i];
		seq[2 * num_seq + 1][i] = seq2[i];
	}
	len_seq[num_seq ++] = len1;

	free((void *) seq1);
	free((void *) seq2);
	return(num_seq);
}

int readtargintv(FILE *fp, int **intvlist, char **regname, int *len_chro, char **chrname, int num_chro)
{
	int	i, j, k, n;
	char	str[100], namestr[100];

	i = 0;
	while(fgets(str, 90, fp))	{
		sscanf(str, "%s%s%d%d", regname[i], namestr, &intvlist[i][1], &intvlist[i][2]);
		intvlist[i][0] = findgenname(namestr, chrname, num_chro) + 1;
		i ++;
	}
	return(i);
}
