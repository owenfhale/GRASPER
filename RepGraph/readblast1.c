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

#define MAX_LEN 1000000000
#define MAX_F 2000

void *ckalloc(int amount);
FILE *ckopen(char *name, char *mode);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, s1, s2, nf;
	int	length, targ_pos[2 * MAX_F], src_pos[2 * MAX_F], alnlen[MAX_F], dir[MAX_F];
	int	ma, tot_len, last;
	double	id[MAX_F];
	FILE	*fp;
	char	str[500], temp[100], temp1[100], temp2[100],
		seq1[1000], seq2[1000], name[1000];

	if(argc < 2)	{
		printf("No input file.\n");
		exit(0);
	}
	fp = ckopen(argv[1], "r");
	n = 0;
	while(fgets(str, 400, fp))	{
		if(str[0] == '>')	{
			if(n > 0)	{
				for(i = 0; i <= nf; i ++)	{
					if(id[i] < 0.85 || alnlen[i] < 1000)	continue;
					printf("%d %d %s(%d) %d %d %d %d/%d %f\n", targ_pos[2 * i],
						targ_pos[2 * i + 1], name, i + 1, src_pos[2 * i], src_pos[2 * i + 1],
						dir[i], alnlen[i], length, id[i]);
				}
			}
			sscanf(&str[1], "%s", name);
			fgets(str, 400, fp);
			sscanf(str, "%*s%*s%d", &length);
			for(i = 0; i < MAX_F; i ++)	{
				src_pos[2 * i] = MAX_LEN;
				src_pos[2 * i + 1] = 0;
				targ_pos[2 * i] = MAX_LEN;
				targ_pos[2 * i + 1] = 0;
			}
			nf = -1;
			n ++;
		} else if(!strncmp(&str[1], "Strand", 6))	{
			sscanf(str, "%*s%*s%s%*s%s", temp1, temp2);
			if(!strcmp(temp2,"Minus"))	{
				dir[nf] = 1;
			} else	{
				dir[nf] = 0;
			}
		} else if(!strncmp(&str[1], "Identities", 10))	{
			sscanf(str, "%*s%*s%s", temp);
			for(i = 0; i < strlen(temp); i ++)	{
				if(temp[i] == '/')	temp[i] = ' ';
			}
			sscanf(temp, "%d%d", &ma, &tot_len);
			id[nf ++] = ((double) ma) / tot_len;
			alnlen[nf] = tot_len;
		} else if(!strncmp(str, "Query:", 6))	{
			sscanf(str, "%*s%d", &k);
			if(k < targ_pos[2 * nf])	{
				targ_pos[2 * nf] = k;
			}
			if(k > targ_pos[2 * nf + 1])	{
				targ_pos[2 * nf + 1] = k;
			}
		} else if(!strncmp(str, "Sbjct:", 6))	{
			sscanf(str, "%*s%d", &k);
			if(k < src_pos[2 * nf])	{
				src_pos[2 * nf] = k;
			}
			if(k > src_pos[2 * nf + 1])	{
				src_pos[2 * nf + 1] = k;
			}
		}
	}
	fclose(fp);
}

FILE *ckopen(char *name, char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)	{
		printf("Cannot open %s.\n", name);
		exit(-1);
	}
	return(fp);
}


/* ckalloc - allocate space; check for success */

void *ckalloc(int amount)
{
	void *p;

	if(amount == 0)	{
		amount = (unsigned) 100;
	}
	if ((p = (void *) calloc( (unsigned) amount, 1)) == NULL)	{
		printf("Ran out of memory.\n");
                printf("There may be errors as follows:\n");
                printf("1) Not enough memory.\n");
                printf("2) The ARRAY may be overrode.\n");
                printf("3) The wild pointers.\n");
                exit(-1);
	}
	return(p);
}
