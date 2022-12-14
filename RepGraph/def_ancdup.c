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

#define par 30
#define MAX_NUM_ALN 40000
#define MIN_SEG_LEG 600
#define MAX_SUBREP 10000
#define MAX_COPY 700
#define MIN_SEG 2
#define MIN_LEG 100000

int min_length;
int seglength;
int simtimes;
double probthresh;

char inpfile[100], lenfile[100], outfile[100];

void initenv(int argc, char **argv);
int segcompar(SEGMENT *a, SEGMENT *b);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, centro[200], len_chroseq[100], k1, k2;
	int	min_l, min_n;
	int	num_seq;
	int	p1, p2, l1, l2, n1;
	char	**chrname, **repnames;
	int	num_chro;
	char	*visit;
	int	**subrep, **rep, *index_subrep, *repindex, num_subrep;
	double	*distchr;
	int	*num_chrseg, num_seg;
	int	*len_subrep;
	int	*ncomp;
	SEGMENT	*segment;
	int	num_segment;
	int	*repeats, num_repeats;
	char	name[1000], temp[1000], c;
	double	*prob, *probi, *probj, probi0, probj0, nprob;
	FILE	*fp;

	readpar();
	initenv(argc, argv);
	random1(&idum);

/*	Input chromsomal information	*/
	chrname = alloc_name(100, 100);
	fp = ckopen(lenfile, "r");
	num_chro = read_chro_centro(chrname, len_chroseq, centro, fp);
	fclose(fp);

/*	input segments	*/

	segment = (SEGMENT *) ckalloc(70000 * sizeof(SEGMENT));
	fp = ckopen(inpfile, "r");
	num_segment = input_segment(segment, chrname, num_chro, fp);
	fclose(fp);
	i = 0;
	while(i < num_segment)	{
		if(segment[i].pos[1] - segment[i].pos[0] < MIN_SEG_LEG)	{
			segment[i] = segment[num_segment - 1];
			num_segment --;
		} else	{
			i ++;
		}
	}
	printf("num_segment %d\n", num_segment);

/*	sort the segments	*/
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);
	printf("sorted\n");

/*	Define repeats from sub-repeats	*/

	repeats = (int *) ckalloc(num_segment * sizeof(int));
	num_repeats = 0;
	for(i = 0; i < num_segment; i ++)	{
		if(i == num_segment - 1 || segment[i + 1].pos[0] > segment[i].pos[1] + min_length ||
		   segment[i + 1].chro != segment[i].chro)	{
			repeats[num_repeats ++] = i;
		}
	}
	printf("%d repeats found.\n", num_repeats);

/*	get the sequences of each repeat	*/

	subrep = (int **) ckalloc(MAX_SUBREP * sizeof(int *));
	rep = (int **) ckalloc(MAX_SUBREP * sizeof(int *));
	len_subrep = (int *) ckalloc(num_segment * sizeof(int));
	for(i = 0; i < MAX_SUBREP; i ++)	{
		subrep[i] = (int *) ckalloc(MAX_COPY * sizeof(int));
		rep[i] = (int *) ckalloc(MAX_COPY * sizeof(int));
	}
	k = -1;
	index_subrep = (int *) ckalloc(MAX_SUBREP * sizeof(int));
	repindex = (int *) ckalloc(MAX_SUBREP * sizeof(int));
	num_subrep = 0;
	for(i = 0; i < num_repeats; i ++)	{
		m = k + 1;
		k = repeats[i];
		for(j = m; j <= k; j ++)	{
			n = findsegment(segment[j].eq_pos[0], repindex, num_subrep);
			if(n < 0)	{
				repindex[num_subrep] = segment[j].eq_pos[0];
				subrep[num_subrep][index_subrep[num_subrep]] = j;
				rep[num_subrep][index_subrep[num_subrep]] = i;
				index_subrep[num_subrep] ++;
				if(segment[j].pos[1] - segment[j].pos[0] > len_subrep[num_subrep] - 1)	{
					len_subrep[num_subrep] = segment[j].pos[1] - segment[j].pos[0] + 1;
				}
				num_subrep ++;
			} else	{
				subrep[num_subrep][index_subrep[num_subrep]] = j;
				rep[num_subrep][index_subrep[num_subrep]] = i;
				index_subrep[num_subrep] ++;
			}
		}
	}
	k = 0;
	for(i = 0; i < num_subrep; i ++)	{
		if(index_subrep[i] > k)	{
			k = index_subrep[i];
		}
	}
	printf("# subrepeats: %d, maximal copy: %d\n", num_subrep, k);

/*	Statistics of the frequecy distribution over the chromosomes	*/

	num_chrseg = (int *) ckalloc(num_chro * sizeof(int));
	for(i = 0; i < num_chro; i ++)	{
		num_chrseg[i] = len_chroseq[i] / seglength + 1;
	}
	for(i = 1; i < num_chro; i ++)	{
		num_chrseg[i] += num_chrseg[i - 1];
	}
	distchr = (double *) ckalloc(num_chrseg[num_chro - 1] * sizeof(double));
	for(i = 0; i < num_segment; i ++)	{
		n = segment[i].pos[0] / seglength;
		if(segment[i].chro == 0)	{
			distchr[n] += 1;
		} else	{
			distchr[num_chrseg[segment[i].chro - 1] + n] += 1;
		}
	}
	for(i = 0; i < num_chrseg[num_chro - 1]; i ++)	{
		distchr[i] = distchr[i] / num_segment;
	}
	printf("distribution done %d \n", num_chrseg[num_chro - 1]);

/*	Simulation of P-value of a pair of subrepeats to both be ancestral duplication units */

	len_subrep = (int *) ckalloc(num_subrep * sizeof(int));
	prob = (double *) ckalloc(pow1(num_subrep) * sizeof(double));
	for(i = 0; i < num_subrep; i ++)	{
		n = 0;
		for(j = 0; j < index_subrep[i]; j ++)	{
			n += segment[rep[i][j]].pos[1] - segment[rep[i][j]].pos[0] + 1;
		}
		len_subrep[i] = n / index_subrep[i];
	}

	p1 = 0;
	for(i = 0; i < num_subrep; i ++)	{
		for(j = i + 1; j < num_subrep; j ++)	{
			n = numc(i, j);
			min_l = seglength;
			nprob = 0;
			probi = (double *) ckalloc(index_subrep[i] * sizeof(double));
			probj = (double *) ckalloc(index_subrep[j] * sizeof(double));
			for(k = 0; k < index_subrep[i]; k ++)	{
				probi[k] = 0;
			}
			for(l = 0; l < index_subrep[j]; l ++)	{
				probj[l] = 0;
			}
			probi0 = probj0 = 1;
			for(k = 0; k < index_subrep[i]; k ++)	{
				for(l = 0; l < index_subrep[j]; l ++)	{
					if(segment[rep[i][k]].chro == segment[rep[j][l]].chro &&
					   abs(segment[rep[i][k]].pos[0] - segment[rep[j][l]].pos[0]) < 2 * seglength) {
/*
						nprob = simulate_pair_all(segment[rep[i][k]].pos[1] - segment[rep[i][k]].pos[0] + 1,
						    segment[rep[j][l]].pos[1] - segment[rep[j][l]].pos[0] + 1,
							index_subrep[i], index_subrep[j],
							abs(segment[rep[i][k]].pos[0] - segment[rep[j][l]].pos[0]),
							 distchr, num_chrseg, num_chro, num_chrseg[num_chro - 1], seglength, simtimes);
*/
						nprob = 0;
						if(segment[rep[i][k]].pos[0] < segment[rep[j][l]].pos[0])	{
							if(segment[rep[j][l]].pos[0] - segment[rep[i][k]].pos[1] + 1 < MIN_LEG)	{
								nprob = 1;
							}
						} else {
							if(segment[rep[i][k]].pos[0] - segment[rep[j][l]].pos[1] + 1 < MIN_LEG)	{
								nprob = 1;
							}
						}
						if(probi[k] < nprob)	{
							probi[k] = nprob;
						}
						if(probj[l] < nprob)	{
							probj[l] = nprob;
						}
					}
				}
			}
			for(l = 0; l < index_subrep[j]; l ++)	{
				if(probi0 > probj[l])	{
					probi0 = probj[l];
				}
			}
			for(k = 0; k < index_subrep[i]; k ++)	{
				if(probj0 > probi[k])	{
					probj0 = probi[k];
				}
			}
			free((void *) probi);
			free((void *) probj);
			prob[n] = max(probi0, probj0);
			if(prob[n] > probthresh)	{
				p1 ++;
			}
/*
printf("%d %d %f\n", i, j, prob[n]);
*/
		}
	}
	free((void *) distchr);
	free((void *) num_chrseg);
	printf("prob edges computed: %d\n", p1);

/*	Compute the minimal number of independent duplication units -- the largest vertex cover */

	min_n = 0;
	visit = (char *) ckalloc(num_subrep * (num_subrep - 1) / 2 * sizeof(char));
	ncomp = (int *) ckalloc(num_subrep * sizeof(int));
	for(k = 0; k < num_subrep; k ++)	{
		for(j = k + 1; j < num_subrep; j ++)	{
			n = numc(k, j);
			if(prob[n] > 0)	{
				ncomp[k] ++;
				ncomp[j] ++;
			}
		}
	}
	do	{
		n1 = 0;
		for(k = 0; k < num_subrep; k ++)	{
			if(ncomp[k] > n1)	{
				n1 = ncomp[k];
				l1 = k;
			}
		}
		if(n1 > 0)	{
			ncomp[l1] = 0;
			for(l = 0; l < num_subrep; l ++)	{
				if(l == l1)	continue;
				n = numc(l1, l);
				if(prob[n] > 0)	{
					ncomp[l] --;
				}
			}
			min_n += 1;
		}
	} while(n1 > 0);

	free((void *) visit);
	free((void *) ncomp);

	printf("The minimal vertex cover has %d duplication units.\n", min_n);

	free((void *) prob);
	free((void *) len_subrep);
	for(i = 0; i < MAX_SUBREP; i ++)	{
		free((void *) subrep[i]);
		free((void *) rep[i]);
	}
	free((void **) subrep);
	free((void **) rep);
	free((void *) index_subrep);
	free((void *) repindex);
	free((void *) repeats);
	free((void *) segment);
}

int findsegment(int eq_pos, int *repindex, int num_subrep)
{
	int	i, j;

	for(i = 0; i < num_subrep; i ++)	{
		if(repindex[i] == eq_pos)	{
			return(i);
		}
	}
	return(-1);
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

	min_length = 5000;
	simtimes = 1000;
	inpseq = outseq = 0;
	probthresh = 0.001;
	seglength = 10000000;

	while ((copt=getopt(argc,argv,"i:o:l:L:s:p:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 'L':
			  sscanf(optarg,"%d", &min_length);
			  continue;
			case 's':
			  sscanf(optarg,"%d", &simtimes);
			  continue;
			case 'p':
			  sscanf(optarg,"%lf", &probthresh);
			  continue;
			default:
			  printf("def_ancdup -i InpFile -l LenFile [-L min_length]\n");
			  printf("-i InpFile: The input file name of segments\n");
			  printf("-l lenfile: inpput file for chromosome length.\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0)	{
		printf("def_ancdup -i InpFile -l LenFile [-L min_length]\n");
		printf("-i InpFile: The input file name of segments\n");
		printf("-l lenfile: inpput file for chromosome length.\n");
		exit(-1);
	}
}
