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

char	inpfile[100], outfile[100], seqfile[100];

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	dist[20];
	int	reads;
	int	num_vertex, num_class, num_edge;
	int	*len_seq, num_seq, num_remain;
	int	**num_pa;
	char	**src_seq, **src_name;
	char	temp[100];
	int 	*range, *newrange;
	int	read1, read2;
	ALIGN	**eq_class, *align, *aln, *aln0;
	char	*mark;
	PATH	*path;
	int	num_path;
	int	num, len;
	READINTERVAL	*readinterval;
	FILE	*fp, *fp1;


	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	Input the length of the reads (required) */

	len_seq = (int *) ckalloc(2 * sizeof(int));
	src_seq = (char **) ckalloc(2 * sizeof(char *));
	src_name = (char **) ckalloc(1 * sizeof(char *));
	src_name[0] = (char *) ckalloc(100 * sizeof(char));

/*	for each component, build-repeat graph	*/

	fp = ckopen(seqfile, "r");
	num_seq = readseq1by1(src_seq, src_name, len_seq, fp);
	fclose(fp);
	printf("Genome length: %d\n", len_seq[0]);

/*	Make reverse complements of input sequences rev(i) --> i + num_seq	*/

	len_seq[1] = len_seq[0];
	src_seq[1] = (char *) ckalloc(len_seq[0] * sizeof(char));
	for(j = 0; j < len_seq[0]; j ++)	{
		src_seq[1][j] = rev(src_seq[0][len_seq[0] - j - 1]);
	}

/*	Input equivalent readintervales between reads --
	see the format of the equivalent readinterval files	*/

	printf("Read equivalent readintervales...\n");
	eq_class = (ALIGN **) ckalloc(2 * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	num_class = readclass(eq_class, num_seq, fp);
	fclose(fp);
	printf("# equivalent readintervales input: %d\n", num_class);

	k = 0;
	mark = (char *) ckalloc(len_seq[0] * sizeof(char));
	for(i = 0; i < 2 * num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			if(read1 == 0)	{
				for(m = align -> pos[0][0]; m <= align -> pos[0][align -> length - 1]; m ++)	{
					mark[m] = 1;
				}
			} else	{
				for(m = align -> pos[0][0]; m <= align -> pos[0][align -> length - 1]; m ++)	{
					mark[len_seq[0] - m - 1] = 1;
				}
			}
			if(read2 == 0)	{
				for(m = align -> pos[1][0]; m <= align -> pos[1][align -> length - 1]; m ++)	{
					mark[m] = 1;
				}
			} else	{
				for(m = align -> pos[1][0]; m <= align -> pos[1][align -> length - 1]; m ++)	{
					mark[len_seq[0] - m - 1] = 1;
				}
			}
			k ++;
			align = align -> next;
		}
	}
	printf("# alignments: %d\n", k);
	range = (int *) ckalloc(8 * k * sizeof(int));
	n = 0;
	for(i = 0; i < len_seq[0]; i ++)	{
		if(mark[i] == 0)	{
			if(i == 0 || mark[i - 1] != 0)
				range[2 * n] = i;
		} else if(i > 0 && mark[i - 1] == 0)	{
			range[2 * n + 1] = i - 1;
			n ++;
		}
	}
	if(mark[i - 1] == 0)	range[2 * n + 1] = i - 1;
	n ++;
	printf("# ranges: %d\n", n);
	num = n;

/*	Output coded sequence	*/

/*	ranges after encoding	*/
	newrange = (int *) ckalloc(2 * num * sizeof(int));
	fp = ckopen(outfile, "w");
	k = 0;
	fprintf(fp, ">codeseq\n");
	for(i = 0; i < range[0]; i ++)	{
		fprintf(fp, "%c", na_name[src_seq[0][i]]);
		if(k % 50 == 49)	{
			fprintf(fp, "\n");
		}
		k ++;
	}
	newrange[0] = k;
	for(i = 0; i < n - 1; i ++)	{
		for(j = 0; j < SEGLEN; j ++)	{
			fprintf(fp, "%c", na_name[0]);
			if(k % 50 == 49)	{
				fprintf(fp, "\n");
			}
			k ++;
		}
		newrange[2 * i + 1] = k - 1; 
		for(j = range[2 * i + 1] + 1; j < range[2 * i + 2]; j ++)	{
			fprintf(fp, "%c", na_name[src_seq[0][j]]);
			if(k % 50 == 49)	{
				fprintf(fp, "\n");
			}
			k ++;
		}
		newrange[2 * i + 2] = k;
	}
	for(j = 0; j < SEGLEN; j ++)	{
		fprintf(fp, "%c", na_name[0]);
		if(k % 50 == 49)	{
			fprintf(fp, "\n");
		}
		k ++;
	}
	newrange[2 * i + 1] = k - 1;
	for(i = range[2 * n - 1] + 1; i < len_seq[0]; i ++)	{
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
	printf("Coded sequence length: %d.\n", k);
	len = k;

/*	Write alignments	*/

	sprintf(temp, "%s.aln", outfile);
	fp = ckopen(temp, "w");
	for(m = 0; m < 2 * num_seq; m ++)	{
		n = size_align(eq_class[m]);
		fwrite(&n, sizeof(int), 1, fp);
		aln = eq_class[m];
		while(aln)	{
			for(i = 0; i < aln -> length; i ++)	{
				if(aln -> reads[0] == 0)	{
					aln -> pos[0][i] = trans_loc(aln -> pos[0][i], range, num);
				} else	{
					j = trans_loc(len_seq[0] - 1 - aln -> pos[0][i], range, num);
					aln -> pos[0][i] = len - 1 - j;
				}
				if(aln -> reads[1] == 0)	{
					aln -> pos[1][i] = trans_loc(aln -> pos[1][i], range, num);
				} else	{
					j = trans_loc(len_seq[0] - 1 - aln -> pos[1][i], range, num);
					aln -> pos[1][i] = len - 1 - j;
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

	for(i = 0; i < num; i ++)	{
		printf("Range %d %d %d %d\n", newrange[2 * i], newrange[2 * i + 1], range[2 * i], range[2 * i + 1]);
	}

	free((void **) eq_class);

	free((void *) range);
	free((void *) newrange);
	free((void *) mark);
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

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:s:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", seqfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			default:
			  printf("longencode -i InpFile -s SeqFile -o outfile\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-s SeqFile: The input file name of reads\n");
			  printf("-o Outfile: The output file name of contigs\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("longencode -i InpFile -s SeqFile -o outfile\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-s SeqFile: The input file name of reads\n");
		printf("-o Outfile: The output file name of contigs\n");
		exit(-1);
	}
}
