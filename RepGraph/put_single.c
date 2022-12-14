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
#include <unistd.h>
#include <math.h>
#include <perdef.h>
#include <param.h>
#include <extvab.h>
#include <extfunc.h>

#define MAX_NUM_ALN 40000

#define INTV -1000

int min_length;
int tel_dist, centro_dist;
int MIN_SEG;

char inpfile[100], lenfile[100], alnfile[100], outfile[100], targfile[100];

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);
double cal_id(char **alnseq, int length);
int chk_ovp_n(int *posindex1, int **posindex2, int num_intv);
char chk_ovp(int *posindex1, int *posindex2);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, n1, centro[200], len_chroseq[100];
	int	coord1[3], coord2[3], posindex[3];
	char	**alnseq;
	char	**chrname;
	int	num_chro;
	int	p1, p2;
	int	length;
	SEGMENT	*segment;
	int	**intvlist, num_intv;
	char	**regname;
	int	num_segment;
	int	*repeats, num_repeats;
	char	name[1000], temp[1000], c;
	double	id;
	FILE	*fp;

	readpar();
	initenv(argc, argv);

/*	Input chromsomal information	*/
	chrname = alloc_name(100, 100);
	fp = ckopen(lenfile, "r");
	num_chro = read_chro_centro(chrname, len_chroseq, centro, fp);
	fclose(fp);

/*	input segments	*/

	segment = (SEGMENT *) ckalloc(20000 * sizeof(SEGMENT));
	fp = ckopen(inpfile, "r");
	num_segment = input_segment(segment, chrname, num_chro, fp);
	fclose(fp);

/*	Input target intevals	*/

	fp = ckopen(targfile, "r");
	intvlist = (int **) ckalloc(100 * sizeof(int *));
	regname = (char **) ckalloc(100 * sizeof(char *));
	for(i = 0; i < 100; i ++)	{
		regname[i] = (char *) ckalloc(20 * sizeof(char));
		intvlist[i] = (int *) ckalloc(3 * sizeof(int));
	}
	num_intv = readtargintv(fp, intvlist, regname, len_chroseq, chrname, num_chro);
	fclose(fp);
	printf("number intevals: %d\n", num_intv);
	for(i = 0; i < num_intv; i ++)	{
		intvlist[i][0] --;
		printf("%d %d %d\n", intvlist[i][0], intvlist[i][1], intvlist[i][2]);
	}

/*	Define repeats from sub-repeats	*/

	repeats = (int *) ckalloc(num_segment * sizeof(int));
	num_repeats = 0;
	for(i = 0; i < num_segment; i ++)	{
		if(i == num_segment - 1 || segment[i + 1].pos[0] > segment[i].pos[1] + min_length ||
		   segment[i + 1].chro != segment[i].chro)	{
			repeats[num_repeats ++] = i;
		}
	}

/*	Output the list of non-centromere/telemere repeats	*/

/*
	k = 0;
	for(i = 0; i < num_repeats; i ++)	{
		p1 = segment[k].pos[0];
		p2 = segment[repeats[i]].pos[1];
		n = repeats[i] - k + 1;
		posindex[0] = segment[k].chro;
		posindex[1] = p1;
		posindex[2] = p2;
		n1 = chk_ovp_n(posindex, intvlist, num_intv);
		if(n1 >= 0)	{
			intvlist[n1][1] = p1;
			intvlist[n1][2] = p2;
		}
		k = repeats[i] + 1;
	}
	for(i = 0; i < num_intv; i ++)	{
		printf("%d %d %d\n", intvlist[i][0], intvlist[i][1], intvlist[i][2]);
	}
*/

	k = 0;
	for(i = 0; i < num_repeats; i ++)	{
		p1 = segment[k].pos[0];
		p2 = segment[repeats[i]].pos[1];
		n = repeats[i] - k + 1;
		if((p1 > centro[2 * segment[k].chro + 1] + centro_dist && p2 < len_chroseq[segment[k].chro] - tel_dist ||
		   p2 < centro[2 * segment[k].chro] - centro_dist && p1 > tel_dist) && n <= MIN_SEG)	{
			for(j = k; j <= repeats[i]; j ++)	{
				for(m = 0; m < num_segment; m ++)	{
					posindex[0] = segment[m].chro;
					posindex[1] = segment[m].pos[0];
					posindex[2] = segment[m].pos[1];
					n1 = chk_ovp_n(posindex, intvlist, num_intv);
					if(segment[m].eq_pos[0] == segment[j].eq_pos[0] && n1 >= 0 &&
					   posindex[2] - posindex[1] > 500)	{
					  if(segment[m].eq_pos[1] == 0)	{
						printf("%s	%d	%d	%s-s%d	:	",
						  chrname[segment[m].chro], segment[m].pos[0] + 1, segment[m].pos[1] + 1,
						  inpfile, segment[m].eq_pos[0]);
				  	  } else	{
						printf("%s	%d	%d	%s-s%d	:	",
						  chrname[segment[m].chro], segment[m].pos[0] + 1, segment[m].pos[1] + 1,
						  inpfile, segment[m].eq_pos[0]);
					  }
					  if(segment[j].eq_pos[1] == 0)	{
						printf("%s	%d	%d\n",
						  chrname[segment[j].chro], segment[j].pos[0] + 1, segment[j].pos[1] + 1);
					  } else	{
						printf("%s	%d	%d\n",
						  chrname[segment[j].chro], segment[j].pos[0] + 1, segment[j].pos[1] + 1);
					  }
					}
				}
/*
				for(m = 0; m < num_segment; m ++)	{
					if(segment[m].chro == 21 && segment[m].eq_pos[0] == segment[j].eq_pos[0])	{
					  if(segment[m].eq_pos[1] == 0)	{
						printf("%s	%d	%d	%d	%s-s%d	%d	%d	%d ID	:	",
						  chrname[segment[m].chro], segment[m].pos[0] + 1, segment[m].pos[1] + 1,
						  len_chroseq[segment[m].chro], inpfile, segment[m].eq_pos[0], segment[m].src_pos[0],
						  segment[m].src_pos[1], segment[m].length);
				  	  } else	{
						printf("%s	%d	%d	%d	%s-s%d	%d	%d	%d ID	:	",
						  chrname[segment[m].chro], segment[m].pos[0] + 1, segment[m].pos[1] + 1,
						  len_chroseq[segment[m].chro], inpfile, segment[m].eq_pos[0], segment[m].src_pos[1],
						  segment[m].src_pos[0], segment[m].length);
					  }
					}
				}
				if(segment[j].eq_pos[1] == 0)	{
					printf("%s	%d	%d	%d	%s-s%d	%d	%d	%d\n",
					  chrname[segment[j].chro], segment[j].pos[0] + 1, segment[j].pos[1] + 1,
					  len_chroseq[segment[j].chro], inpfile, segment[j].eq_pos[0], segment[j].src_pos[0],
					  segment[j].src_pos[1], segment[j].length);
				} else	{
					printf("%s	%d	%d	%d	%s-s%d	%d	%d	%d\n",
					  chrname[segment[j].chro], segment[j].pos[0] + 1, segment[j].pos[1] + 1,
					  len_chroseq[segment[j].chro], inpfile, segment[j].eq_pos[0], segment[j].src_pos[1],
					  segment[j].src_pos[0], segment[j].length);
				}
*/
			}
		}
		k = repeats[i] + 1;
	}
/*
*/

	for(i = 0; i < 100; i ++)	{
		free((void *) intvlist[i]);
		free((void *) regname[i]);
	}
	free((void **) intvlist);
	free((void **) regname);
	free((void **) alnseq);
	chrname = free_name(chrname, 100);
	free((void *) repeats);
	free((void *) segment);
}

int chk_ovp_n(int *posindex1, int **posindex2, int num_intv)
{
	int	i;
	char	c;

	for(i = 0; i < num_intv; i ++)	{
		c = chk_ovp(posindex1, posindex2[i]);
		if(c)	{
			return(i);
		}
	}
	return(-1);
}

char chk_ovp(int *posindex1, int *posindex2)
{
	int	i, j, k, l;
	int	range1[2], range2[2], nchro1, nchro2;

	nchro1 = posindex1[0];
	nchro2 = posindex2[0];
	if(nchro1 != nchro2)	{
		return(0);
	}
	range1[0] = posindex1[1];
	range1[1] = posindex1[2];
	range2[0] = posindex2[1];
	range2[1] = posindex2[2];
	if(range2[0] >= range1[0])	{
		if(range2[0] <= range1[1] - INTV)	{
			return(1);
		} else	{
			return(0);
		}
	} else 	{
		if(range1[0] <= range2[1] - INTV)	{
			return(1);
		} else	{
			return(0);
		}
	}
}


int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> chro > b -> chro)	return(1);
	else if(a -> chro == b -> chro && a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	MIN_SEG = 1;
	min_length = 5000;
	centro_dist = 5000000;
	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:a:L:l:t:c:T:s:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'L':
			  sscanf(optarg,"%d", &min_length);
			  continue;
			case 'T':
			  sscanf(optarg,"%s", targfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 't':
			  sscanf(optarg,"%d", &tel_dist);
			  continue;
			case 's':
			  sscanf(optarg,"%d", &MIN_SEG);
			  continue;
			case 'c':
			  sscanf(optarg,"%d", &centro_dist);
			  continue;
			default:
			  printf("put_single -i InpFile -l LenFile -T TargFile [-L min_leg -t Tel_dist -c Centro_dist]\n");
			  printf("-i InpFile: The input file name of segments\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  printf("-T TargFile: inpput file for ranges.\n");
			  printf("-L min_length: minimum interval length between repeats\n");
			  printf("-t Tel_dist: distance from the telemere\n");
			  printf("-c Centro_dist: distance from the centromere\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0)	{
		printf("put_single -i InpFile -l LenFile -T TargFile [-L min_leg -t Tel_dist -c Centro_dist]\n");
		printf("-i InpFile: The input file name of segments\n");
		printf("-l lenfile: inpput file for chromosome length.\n");
		printf("-T TargFile: inpput file for ranges.\n");
		printf("-L min_length: minimum interval length between repeats\n");
		printf("-t Tel_dist: distance from the telemere\n");
		printf("-c Centro_dist: distance from the centromere\n");
		exit(-1);
	}
}
