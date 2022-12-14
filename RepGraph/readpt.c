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

int readblast(ALIGN **align, FILE *fp, int MIN_LEG, double MIN_ID);

int readblast(ALIGN **align, FILE *fp, int MIN_LEG, double MIN_ID)
{
	int	i, j, k, l, m, n, s1, s2;
	int	length, *pos[2], len;
	int	nall;
	int	ma, tot_len;
	int	beg1, beg2, mis;
	double	id;
	ALIGN	*align0;
	char	str[500], temp[100], temp1[100], temp2[100],
		seq1[1000], seq2[1000];

	pos[0] = (int *) ckalloc(10000 * sizeof(int));
	pos[1] = (int *) ckalloc(10000 * sizeof(int));

	nall = 0;
	while(fgets(str, 400, fp))	{
		sscanf(str, "%s", temp);
		if(temp[0] == '(')	{
			len = atoi(&temp[1]);
			break;
		}
	}
	while(fgets(str, 400, fp))	{
		if(!strncmp(&str[1], "Identities", 10))	{
			sscanf(str, "%*s%*s%s", temp);
			for(i = 0; i < strlen(temp); i ++)	{
				if(temp[i] == '/')	temp[i] = ' ';
			}
			sscanf(temp, "%d%d", &ma, &tot_len);
			id = ((double) ma) / tot_len;
			if(tot_len >= MIN_LEG && id >= MIN_ID)	{
				fgets(str, 400, fp);
				sscanf(str, "%*s%*s%s%*s%s", temp1, temp2);
				if(!strncmp(temp1, "Plus", 4))	{
					s1 = 0;
				} else	{
					s1 = 1;
				}
				if(!strncmp(temp2, "Plus", 4))	{
					s2 = 0;
				} else	{
					s2 = 1;
				}

/*	Don't consider alignments both from '-' strands	*/
				if(s1 == 1 && s2 == 1)	continue;

				while(fgets(str, 400, fp))	{
					if(!strncmp(str, "Query", 5))	{
						sscanf(str, "%*s%d%s", &beg1, seq1);
					}
					if(!strncmp(str, "Sbjct", 5))	{
						sscanf(str, "%*s%d%s", &beg2, seq2);
						break;
					}
				}
/*	If s1 == s2 && beg1 == beg2, alignments start from the same position;
	if s1 != s2, choose only the half of the all alignments	*/
				if(s1 == s2 && beg1 >= beg2 || s1 == 0 && s2 == 1 && beg1 <= len / 2 ||
				   s1 == 1 && s2 == 0 && beg2 <= len / 2)	continue;
				mis = 0;
				align0 = (ALIGN *) ckalloc(1 * sizeof(ALIGN));
				align0 -> reads[0] = s1;
				align0 -> reads[1] = s2;
				pos[0][0] = beg1;
				pos[1][0] = beg2;
				length = 1;
				do	{
					l = strlen(seq1);
					for(i = 0; i < l; i ++)	{
						if(seq1[i] == '-')	{
							beg2 ++;
							if(seq1[i + 1] != '-')	{
								pos[0][length] = beg1;
								pos[1][length] = beg2;
								length ++;
							}
						} else if(seq2[i] == '-')	{
							beg1 ++;
							if(seq2[i + 1] != '-')	{
								pos[0][length] = beg1;
								pos[1][length] = beg2;
								length ++;
							}
						} else	{
							if(seq1[i] != seq2[i])	mis ++;
							beg1 ++;
							beg2 ++;
						}
					}
					while(fgets(str, 400, fp))	{
						if(!strncmp(str, " Score", 6))	break;
						if(!strncmp(str, "Query", 5))	{
							sscanf(str, "%*s%*d%s", seq1);
						}
						if(!strncmp(str, "Sbjct", 5))	{
							sscanf(str, "%*s%*d%s", seq2);
							break;
						}
					}
					if(!strncmp(str, " Score", 6))	break;
				} while(l > 0);
				pos[0][length] = beg1 + 1;
				pos[1][length] = beg2 + 1;
				length ++;
				align0 -> pos[0] = (int *) ckalloc(length * sizeof(int));
				align0 -> pos[1] = (int *) ckalloc(length * sizeof(int));
				align0 -> length = length;
				for(j = 0; j < length; j ++)	{
					align0 -> pos[0][j] = pos[0][j];
					align0 -> pos[1][j] = pos[1][j];
				}
				align0 -> mis_match = mis;
				nall ++;
				if(s1 == 0)	{
					align0 -> next = align[0];
					align[0] = align0;
				} else	{
					align0 -> next = align[1];
					align[1] = align0;
				}
			}
		}
	}
	free((void *) pos[0]);
	free((void *) pos[1]);
	return(nall);
}
