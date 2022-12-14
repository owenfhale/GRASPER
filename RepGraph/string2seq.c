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
#include <param.h>
#include <extfunc.h>

char	inpfile[100], outfile[100], lenfile[100], segfile[100], newlenfile[100];
int	qualinp;
int	min_leg;
double	min_id;

int segcompar(SEGMENT *a, SEGMENT *b);

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	last, front;
	int	*num_segment, num_chro, num_class;
	int	len_seq[100];
	int	read1, read2;
	char	temp[100];
	char	**chrname;
	int	**range, **newrange, length[100], num[100];
	ALIGN	**eq_class, *align, *aln, *aln0;
	SEGMENT **segment;
	FILE	*fp;

	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	input the length of the genomes	*/
	chrname = alloc_name(100, 100);
	fp = ckopen(lenfile, "r");
	num_chro = readlen(fp, len_seq, chrname);
	fclose(fp);

/*      read in pairwise alignments 	*/
	fp = ckopen(inpfile, "r");
	eq_class = (ALIGN **) ckalloc(2 * num_chro * sizeof(ALIGN *));
	num_class = readclass(eq_class, num_chro, fp);
	fclose(fp);
	printf("# equivalent readintervales input: %d\n", num_class);

	segment = (SEGMENT **) ckalloc(num_chro * sizeof(SEGMENT *));
	num_segment = (int *) ckalloc(num_chro * sizeof(int));
	for(i = 0; i < num_chro; i ++)	{
		segment[i] = (SEGMENT *) ckalloc(2 * num_class * sizeof(SEGMENT));
	}
	k = 0;
	for(i = 0; i < 2 * num_chro; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			if(read1 < num_chro)	{
				segment[read1][num_segment[read1]].pos[0] = align -> pos[0][0];
				segment[read1][num_segment[read1]].pos[1] = align -> pos[0][align -> length - 1] - 1;
				num_segment[read1] ++;
			} else	{
				segment[read1 - num_chro][num_segment[read1 - num_chro]].pos[0] =
					len_seq[read1 - num_chro] - align -> pos[0][align -> length - 1];
				segment[read1 - num_chro][num_segment[read1 - num_chro]].pos[1] =
					len_seq[read1 - num_chro] - 1 - align -> pos[0][0];
				num_segment[read1 - num_chro] ++;
			}
			if(read2 < num_chro)	{
				segment[read2][num_segment[read2]].pos[0] = align -> pos[1][0];
				segment[read2][num_segment[read2]].pos[1] = align -> pos[1][align -> length - 1] - 1;
				num_segment[read2] ++;
			} else	{
				segment[read2 - num_chro][num_segment[read2 - num_chro]].pos[0] =
					len_seq[read2 - num_chro] - align -> pos[1][align -> length - 1];
				segment[read2 - num_chro][num_segment[read2 - num_chro]].pos[1] =
					len_seq[read2 - num_chro] - 1 - align -> pos[1][0];
				num_segment[read2 - num_chro] ++;
			}
			k ++;
			align = align -> next;
		}
	}
	printf("Processed # alignments: %d\n", k);

	range = (int **) ckalloc(num_chro * sizeof(int *));
	newrange = (int **) ckalloc(num_chro * sizeof(int *));
	for(i = 0; i < num_chro; i ++)	{
		range[i] = (int *) ckalloc(2 * (num_segment[i] + 1) * sizeof(int));
		newrange[i] = (int *) ckalloc(2 * (num_segment[i] + 1) * sizeof(int));
	}
	fp = ckopen(segfile, "w");
	for(i = 0; i < num_chro; i ++)	{
		qsort((void *) segment[i], num_segment[i], sizeof(SEGMENT), (void *) segcompar);

		n = 0;
		last = 0;
		l = 0;
		if(num_segment[i] > 0 && segment[i][0].pos[0] > 0)	{
			range[i][0] = 0;
			newrange[i][0] = 0;
			last = range[i][1] = segment[i][0].pos[0] - 1;
			l = SEGLEN;
			newrange[i][1] = l - 1;
			n ++;
		}
		front = 0;
		for(j = 1; j < num_segment[i]; j ++)	{
			front = max(front, segment[i][j - 1].pos[1] + 100);
			if(segment[i][j].pos[0] > front)	{
				range[i][2 * n] = front + 1;
				l += range[i][2 * n] - last;
				newrange[i][2 * n] = l - 1;
				last = range[i][2 * n + 1] = segment[i][j].pos[0] - 1;
				l += SEGLEN;
				newrange[i][2 * n + 1] = l - 1;
				n ++;
			}
		}
		if(num_segment[i] > 0)	{
			front = max(front, segment[i][j - 1].pos[1]);
			range[i][2 * n] = front + 1;
			l += range[i][2 * n] - last;
			newrange[i][2 * n] = l - 1;
			if(front < len_seq[i] - 1)	{
				range[i][2 * n + 1] = len_seq[i] - 1;
				l += SEGLEN;
				newrange[i][2 * n + 1] = l - 1;
				n ++;
			}
		}
		num[i] = n;
		for(j = 0; j < num[i]; j ++)	{
			fprintf(fp, "chromosome %d Range %d %d %d %d\n", i + 1, newrange[i][2 * j], newrange[i][2 * j + 1],
				range[i][2 * j], range[i][2 * j + 1]);
		}
		length[i] = l;
	}
	fclose(fp);
	printf("code done\n");

/*	Write alignments	*/

	fp = ckopen(outfile, "w");
	for(m = 0; m < 2 * num_chro; m ++)	{
		n = size_align(eq_class[m]);
		fwrite(&n, sizeof(int), 1, fp);
		aln = eq_class[m];
		while(aln)	{
			for(i = 0; i < aln -> length; i ++)	{
				if(aln -> reads[0] < num_chro)	{
					aln -> pos[0][i] = trans_loc(aln -> pos[0][i], range[aln -> reads[0]],
						 num[aln -> reads[0]]);
				} else	{
					j = trans_loc(len_seq[aln -> reads[0] - num_chro] - 1 - aln -> pos[0][i],
						 range[aln -> reads[0] - num_chro], num[aln -> reads[0] - num_chro]);
					aln -> pos[0][i] = length[aln -> reads[0] - num_chro] - 1 - j;
				}
				if(aln -> reads[1] < num_chro)	{
					aln -> pos[1][i] = trans_loc(aln -> pos[1][i], range[aln -> reads[1]],
						 num[aln -> reads[1]]);
				} else	{
					j = trans_loc(len_seq[aln -> reads[1] - num_chro] - 1 - aln -> pos[1][i],
						 range[aln -> reads[1] - num_chro], num[aln -> reads[1] - num_chro]);
					aln -> pos[1][i] = length[aln -> reads[1] - num_chro] - 1 - j;
				}
				if(aln -> pos[0][i] < 0 || aln -> pos[1][i] < 0)	{
					printf("reads %d %d pos %d %d\n", aln -> reads[0], aln -> reads[1],
						aln -> pos[0][i], aln -> pos[1][i]);
				}
			}
			fwrite(&(aln -> reads[1]), sizeof(int), 1, fp);
			fwrite(&(aln -> mis_match), sizeof(int), 1, fp);
			fwrite(&(aln -> length), sizeof(int), 1, fp);
			fwrite(aln -> pos[0], sizeof(int), aln -> length, fp);
			fwrite(aln -> pos[1], sizeof(int), aln -> length, fp);
			aln0 = aln -> next;
			free((void *) aln -> pos[0]);
			free((void *) aln -> pos[1]);
			free((void *) aln);
			aln = aln0;
		}
	}
	fclose(fp);
	printf("Done...\n");

	for(i = 0; i < num_chro; i ++)	{
		free((void *) range[i]);
		free((void *) newrange[i]);
	}
	free((void **) range);
	free((void **) newrange);

/*	Output encoded genome length	*/
	fp = ckopen(newlenfile, "w");
	for(i = 0; i < num_chro; i ++)	{
		fprintf(fp, "%d %d %s %d\n", i + 1, length[i], chrname[i], len_seq[i]);
	}
	fclose(fp);

	chrname = free_name(chrname, 100);
	free((void *) num_segment);
	for(i = 0; i < num_chro; i ++)	{
		free((void *) segment[i]);
	}
	free((void **) segment);
	free((void **) eq_class);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	min_leg = 500;
	min_id = 0.99;
	inpseq = qualinp = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:L:c:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 'L':
			  sscanf(optarg,"%s", newlenfile);
			  continue;
			case 'c':
			  sscanf(optarg,"%s", segfile);
			  continue;
			default:
			  printf("string2seq -i InpFile -o outfile -l lenfile -L newlenfile -c segfile\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-o OutFile: output alignment file\n");
			  printf("-l lenfile: the input file of the length of the genome.\n");
			  printf("-L lenfile: the output file of the length of the encoded genome.\n");
			  printf("-c segfile: the output file of encoded segments.\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("string2seq -i InpFile -o outfile -l lenfile -L newlenfile -c segfile\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-o OutFile: output alignment file\n");
		printf("-l lenfile: the input file of the length of the genome.\n");
		printf("-L lenfile: the output file of the length of the encoded genome.\n");
		printf("-c segfile: the output file of encoded segments.\n");
		exit(-1);
	}
}

int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}
