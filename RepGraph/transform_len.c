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

void outputseq(FILE *fp, char *seq, int len_seq, char *seqname);
int transform_len(char *seqname, char **chrname, int num_chro, int *len_chro, int *posindex);
int getstr_new(char *seqname, int *range, int *len_chro, char **chrname, int num_chro);

void outputseq(FILE *fp, char *seq, int len_seq, char *seqname)
{
	int	i, j, k, l;

	fprintf(fp, ">%s\n", seqname);
	for(i = 0; i < len_seq; i ++)	{
		if(seq[i] < 4)	{
			fprintf(fp, "%c", na_name[seq[i]]);
		} else	{
			fprintf(fp, "-");
		}
		if(i % 60 == 59)	{
			fprintf(fp, "\n");
		}
	}
	if(i % 60 != 0)	{
		fprintf(fp, "\n");
	}
}

int transform_len(char *seqname, char **chrname, int num_chro, int *len_chro, int *posindex)
{
	posindex[0] = getstr_new(seqname, &posindex[1], len_chro, chrname, num_chro);
	if(posindex[0] < 0)	{
		posindex[0] = -posindex[0];
		return(-1);
	} else	{
		return(1);
	}
}


int getstr_new(char *str, int *range, int *len_chro, char **chrname, int num_chro)
{
	int	i, j, k, l;
	char	name[200];

	l = strlen(str);
	for(i = 0; i < l; i ++)	{
		if(str[i] == ':' || str[i] == '-')	{
			str[i] = ' ';
		}
	}
	sscanf(str, "%s %d%d", name, &range[0], &range[1]);
	k = findgenname(name, chrname, num_chro);
	k ++;
	if(range[0] > range[1])	{
		i = range[1];
		range[1] = range[0];
		range[0] = i;
		k = -k;
	}

	return(k);
}
