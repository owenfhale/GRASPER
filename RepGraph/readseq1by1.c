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

#define MAXNUM 250000000

int readseq1by1(char **src_seq, char **src_name, int *len_seq, FILE *fp);
int readseq1by1gen(char **src_seq, char **src_name, int *len_seq, FILE *fp);

int readseq1by1(char **src_seq, char **src_name, int *len_seq, FILE *fp)
{
	int	i, j, k, n;
	char	*seq, c;
	char	str[500];

	seq = (char *) ckalloc(MAXNUM * sizeof(char));

	n = 0;
	k = -1;
	while(fgets(str, 450, fp))	{
		if(str[0] == '#')	continue;
		if(str[0] == '>')	{
			if(k >= 0)	{
				len_seq[k] = n;
				src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
				for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
			}
			n = 0;
			k ++;
			sscanf(&str[1], "%s", src_name[k]);
		} else {
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					c = char2int(str[i]);
					seq[n ++] = c;
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					c = char2int(str[i] - 'A' + 'a');
					seq[n ++] = c;
				}
			}
		}
	}

	if(k >= 0)	{
		len_seq[k] = n;
		src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
		for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
	}
	k ++;

	free((void *) seq);
	return(k);
}

int readseq1by1gen(char **src_seq, char **src_name, int *len_seq, FILE *fp)
{
	int	i, j, k, n;
	char	*seq, c;
	char	str[500];

	seq = (char *) ckalloc(MAXNUM * sizeof(char));

	n = 0;
	k = -1;
	while(fgets(str, 450, fp))	{
		if(str[0] == '#')	continue;
		if(str[0] == '>')	{
			if(k >= 0)	{
				len_seq[k] = n;
				src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
				for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
			}
			n = 0;
			k ++;
			sscanf(&str[1], "%s", src_name[k]);
		} else {
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					c = char2intgen(str[i]);
					seq[n ++] = c;
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					c = char2intgen(str[i] - 'A' + 'a');
					seq[n ++] = c;
				}
			}
		}
	}

	if(k >= 0)	{
		len_seq[k] = n;
		src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
		for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
	}
	k ++;

	free((void *) seq);
	return(k);
}
