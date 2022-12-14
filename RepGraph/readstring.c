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

int getstr(char *str, int *range, int *len_seq, char **chrname, int num_chro);
int getseq(char *str, char *seq, int len);
double seq2sapp(char *seq1, char *seq2, int len1, int len2, int *sapp, int *len);
int readstring(ALIGN **align, FILE *fp, int MIN_LEG, double MIN_ID, int *len_seq, char **chrname, int num_chro);
int readlen(FILE *fp, int *len_seq, char **chrname);
int readlenhetero(FILE *fp, int *len_seq, char **chrname, int *hpos);
char **alloc_name(int n1, int n2);
char **free_name(char **chrname, int n1);
int findgenname(char *name, char **chrname, int num_chro);

char **alloc_name(int n1, int n2)
{
	int	i;
	char	**chrname;

	chrname = (char **) ckalloc(n1 * sizeof(char *));
	for(i = 0; i < n1; i ++)	{
		chrname[i] = (char *) ckalloc(n2 * sizeof(char));
	}
	return(chrname);
}

char **free_name(char **chrname, int n1)
{
	int	i;

	for(i = 0; i < n1; i ++)	{
		free((void *) chrname[i]);
	}
	free((void **) chrname);
	return((char **) NULL);
}

int readlen(FILE *fp, int *len_seq, char **chrname)
{
	int	i, j, k, l, n;
	char	str[500];

	n = 0;
	while(fgets(str, 400, fp))	{
		sscanf(str, "%d%d%s", &j, &k, chrname[n]);
		len_seq[j - 1] = k;
		n ++;
	}
	return(n);
}

int readlenhetero(FILE *fp, int *len_seq, char **chrname, int *hpos)
{
	int	i, j, k, l, n;
	char	str[500];

	n = 0;
	while(fgets(str, 400, fp))	{
		sscanf(str, "%d%d%s%d%*d%d", &j, &k, chrname[n], &hpos[2 * n], &hpos[2 * n + 1]);
		len_seq[j - 1] = k;
		n ++;
	}
	return(n);
}

int readstring(ALIGN **align, FILE *fp, int MIN_LEG, double MIN_ID, int *len_seq, char **chrname, int num_chro)
{
	int	i, j, k, l, m, n, s1, s2, pos, l1, l2;
	int	nchro, range[2], q1, q2;
	int	len1, len2, pos1, pos2, len[2];
	char	*seq1, *seq2;
	int	nall;
	int	*sapp;
	double	id;
	char	str[500];

	seq1 = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));
	seq2 = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));
	sapp = (int *) ckalloc(2 * MAX_TMP_LEG * sizeof(int));
	nall = 0;
	k = 0;
	len1 = len2 = 0;
	while(fgets(str, 400, fp))	{
		pos = 0;
		if(str[0] == '>')	{
			if(k == 0 && len1 > 0 && len2 > 0)	{
				if(len1 != len2)	{
					printf("Length no equal: %d %d\n", len1, len2);
					printf("str %s\n", str);
					printf("pos1 %d pos2 %d\n", pos1, pos2);
					exit(0);
				}
				id = seq2sapp(seq1, seq2, len1, len2, sapp, len);
				if(q1 * q2 > 0)	{
					s1 = abs(q1) - 1;
					s2 = abs(q2) - 1;
					if(len1 >= MIN_LEG && id >= MIN_ID)	{
						align[s1] = new_align(sapp, align[s1], s1, s2,
							    pos1, pos2, len[0], len[1]);
						nall ++;
					}
				} else if(q1 * q2 < 0)	{
					if(len1 >= MIN_LEG && id >= MIN_ID)	{
						if(q2 > 0)	{
							s1 = -q1 - 1 + num_chro;
							s2 = q2 - 1;
							align[s1] = new_align(sapp, align[s1], s1, s2,
								    pos1, pos2, len[0], len[1]);
						} else	{
							s1 = q1 - 1;
							s2 = -q2 - 1 + num_chro;
							align[s1] = new_align(sapp, align[s1], s1, s2,
								    pos1, pos2, len[0], len[1]);
						}
						nall ++;
					}
				}
				len1 = len2 = 0;
				n ++;
			}
			nchro = getstr(&str[1], range, len_seq, chrname, num_chro);
			if(k == 1)	{
				q2 = nchro;
				pos2 = range[0];
				l2 = range[1];
				k = 0;
			} else	{
				q1 = nchro;
				pos1 = range[0];
				l1 = range[1];
				k = 1;
			}
		} else	{
			if(k == 1)	{
				len1 = getseq(str, seq1, len1);
			} else	{
				len2 = getseq(str, seq2, len2);
			}
		}
	}

	pos2 = range[0];
	if(len1 != len2)	{
		printf("Length no equal: %d %d\n", len1, len2);
		exit(0);
	}
	id = seq2sapp(seq1, seq2, len1, len2, sapp, len);
	if(q1 * q2 > 0)	{
		s1 = abs(q1) - 1;
		s2 = abs(q2) - 1;
		if(len1 >= MIN_LEG && id >= MIN_ID)	{
			align[s1] = new_align(sapp, align[s1], s1, s2,
				    pos1, pos2, len[0], len[1]);
			nall ++;
		}
	} else if(q1 * q2 < 0)	{
		if(len1 >= MIN_LEG && id >= MIN_ID)	{
			if(q2 > 0)	{
				s1 = -q1 - 1 + num_chro;
				s2 = q2 - 1;
				align[s1] = new_align(sapp, align[s1], s1, s2,
					    pos1, pos2, len[0], len[1]);
			} else	{
				s1 = q1 - 1;
				s2 = -q2 - 1 + num_chro;
				align[s1] = new_align(sapp, align[s1], s1, s2,
					    pos1, pos2, len[0], len[1]);
			}
			nall ++;
		}
	}

	free((void *) sapp);
	free((void *) seq1);
	free((void *) seq2);
	return(nall);
}

int getstr(char *str, int *range, int *len_seq, char **chrname, int num_chro)
{
	int	i, j, k, l;
	char	name[100];

	l = strlen(str);
	for(i = 0; i < l; i ++)	{
		if(str[i] == ':' || str[i] == '-')	{
			str[i] = ' ';
		}
	}
	sscanf(str, "%s %d%d", name, &range[0], &range[1]);
	k = findgenname(name, chrname, num_chro);
	k ++;

	range[0] --;
	range[1] --;
	if(range[0] < range[1])	{
		return(k);
	} else {
		range[0] = len_seq[k - 1] - 1 - range[0];
		range[1] = len_seq[k - 1] - 1 - range[1];
/*
		j = range[0];
		range[0] = range[1];
		range[1] = j;
*/
		return(-k);
	}
}

double seq2sapp(char *seq1, char *seq2, int len1, int len2, int *sapp, int *len)
{
	int	i, j, k, l, n;

	n = 0;
	k = j = 0;
	len[0] = len[1] = 0;
	for(i = 0; i < len1; i ++)	{
		if(seq1[i] == seq2[i])	{
			if(j != 0)	{
				sapp[k ++] = j;
			}
			n ++;
			sapp[k ++] = 0;
			len[0] ++;
			len[1] ++;
			j = 0;
		} else if(seq1[i] == 9)	{
			len[1] ++;
			j ++;
		} else if(seq2[i] == 9)	{
			len[0] ++;
			j --;
		} else	{
			if(j != 0)	{
				sapp[k ++] = j;
			}
			sapp[k ++] = 0;
			len[0] ++;
			len[1] ++;
			j = 0;
		}
	}
	if(j != 0)	{
		sapp[k ++] = j;
	}
	return(((double) n) / len1);
}

int getseq(char *str, char *seq, int len)
{
	int	i, j, k, l;

	l = strlen(str);
	for(i = 0; i < l; i ++)	{
		if(str[i] >= 'A' && str[i] <= 'Z')	{
			seq[len ++] = char2int(str[i] - 'A' + 'a');
		} else if(str[i] >= 'a' && str[i] <= 'z')	{
			seq[len ++] = char2int(str[i]);
		} else if(str[i] == '-')	{
			seq[len ++] = 9;
		}
	}
	return(len);
}

int findgenname(char *name, char **chrname, int num_chro)
{
	int	i, j;

	for(i = 0; i < num_chro; i ++)	{
		if(!strcmp(name, chrname[i]))	{
			return(i);
		}
	}
	printf("%s not found.\n", name);
	exit(0);
}
