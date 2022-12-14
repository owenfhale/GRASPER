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

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100], seqfile[100], lenfile[100], segfile[100];
int	qualinp;
int	min_leg;
double	min_id;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	char	temp[100];
	int	len_seq[100], num_chro;
	char	**chrname;
	ALIGN	**align, *aln, *aln0;
	FILE	*fp;

	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	read in chromosome length	*/
	chrname = alloc_name(100, 100);
	fp = ckopen(lenfile, "r");
	num_chro = readlen(fp, len_seq, chrname);
	fclose(fp);
	printf("# chromosome %d\n", num_chro);

/*      read in pairwise alignments 	*/

	align = (ALIGN **) ckalloc(2 * num_chro * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	n = readstring(align, fp, min_leg, min_id, len_seq, chrname, num_chro);
	fclose(fp);
	printf("# alignments input: %d.\n", n);

/*	Write alignments	*/

	fp = ckopen(outfile, "w");
	for(m = 0; m < 2 * num_chro; m ++)	{
		n = size_align(align[m]);
		fwrite(&n, sizeof(int), 1, fp);
		aln = align[m];
		while(aln)	{
			fwrite(&(aln -> reads[1]), sizeof(int), 1, fp);
			fwrite(&(aln -> mis_match), sizeof(int), 1, fp);
			fwrite(&(aln -> length), sizeof(int), 1, fp);
			fwrite(aln -> pos[0], sizeof(int), aln -> length, fp);
			fwrite(aln -> pos[1], sizeof(int), aln -> length, fp);
/*
printf("aln %d %d length %d %d %d %d %d\n", aln -> reads[0], aln -> reads[1], aln -> length,
 aln -> pos[0][aln -> length - 1], aln -> pos[0][0],
 aln -> pos[1][aln -> length - 1], aln -> pos[1][0]);
getchar();
*/
			aln0 = aln -> next;
			free((void *) aln -> pos[0]);
			free((void *) aln -> pos[1]);
			free((void *) aln);
			aln = aln0;
		}
	}
	fclose(fp);
	printf("Done...\n");

	chrname = free_name(chrname, 100);
	free((void **) align);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	min_leg = 500;
	min_id = 0.99;
	inpseq = qualinp = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:d:c:")) != EOF)	{
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
			  sscanf(optarg,"%d", &min_leg);
			  continue;
			case 'c':
			  sscanf(optarg,"%s", lenfile);
			  continue;
			case 'd':
			  sscanf(optarg,"%lf", &min_id);
			  continue;
			default:
			  printf("string2aln -i InpFile -c LenFile -o outfile [-l min_leg -d min_id]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-c lenfile: inpput file for chromosome length.\n");
			  printf("-o OutFile: output alignment file\n");
			  printf("-l min_leg: minimum length of repeats.\n");
			  printf("-d min_id: minimum identity of repeats.\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("string2aln -i InpFile -c LenFile -o outfile [-l min_leg -d min_id]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-c lenfile: inpput file for chromosome length.\n");
		printf("-o OutFile: output alignment file\n");
		printf("-l min_leg: minimum length of repeats.\n");
		printf("-d min_id: minimum identity of repeats.\n");
		exit(-1);
	}
}
