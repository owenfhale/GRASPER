
#include <stdinc.h>
#include <param.h>
#include <extfunc.h>

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100];
int	qualinp;
int	min_leg;
double	min_id;

void initenv(int argc, char **argv);

int main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	char	**src_seq, **src_name;
	char	temp[100];
	ALIGN	**align, *aln, *aln0;
	FILE	*fp;

	readpar();
	random1(&idum);
	initenv(argc, argv);

/*      read in pairwise alignments by BLAST	*/

	align = (ALIGN **) ckalloc(2 * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	n = readblast(align, fp, min_leg, min_id);
	fclose(fp);
	printf("# alignments input: %d.\n", n);

/*	Write alignments	*/

	fp = ckopen(outfile, "w");
	for(m = 0; m < 2; m ++)	{
		n = size_align(align[m]);
		fwrite(&n, sizeof(int), 1, fp);
		aln = align[m];
		while(aln)	{
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

	free((void **) align);
	
	return 0;
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	min_leg = 500;
	min_id = 0.99;
	inpseq = qualinp = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:l:d:")) != EOF)	{
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
			case 'd':
			  sscanf(optarg,"%lf", &min_id);
			  continue;
			default:
			  printf("bl2aln -i InpFile -o outfile [-l min_leg -d min_id]\n");
			  printf("-i InpFile: The input file name of reads\n");
			  printf("-o OutFile: output alignment file\n");
			  printf("-l min_leg: minimum length of repeats.\n");
			  printf("-d min_id: minimum identity of repeats.\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("bl2aln -i InpFile -o outfile [-l min_leg -d min_id]\n");
		printf("-i InpFile: The input file name of reads\n");
		printf("-o OutFile: output alignment file\n");
		printf("-l min_leg: minimum length of repeats.\n");
		printf("-d min_id: minimum identity of repeats.\n");
		exit(-1);
	}
}
