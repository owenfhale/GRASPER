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
#include <color.h>
#include <param.h>
#include <extfunc.h>

char	inpfile[100], graphfile[100], edgefile[100], outfile[100], seqfile[100];

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	dist[20];
	int	reads;
	int	num_vertex, num_edge;
	int	**num_pa;
	char	temp[100];
	EDGE	**edge, *edge1, *edge2, *bal_edge1, *bal_edge2;
	NODES	**vertex, *begin, *node, *node_next, *start_node;
	READINTERVAL	*readinterval;
	POSITION	*position;
	FILE	*fp, *fp1;

	readpar();
	random1(&idum);
	initenv(argc, argv);

/*	input graph and edge files	*/
	read_graph_file(edgefile, graphfile,
		  &num_vertex, &vertex,
		  &num_edge, &edge);

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}
	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("%d vertices %d edges (%d source %d sinks) in the beginning\n", num_vertex, num_edge,
		num_pa[0][1], num_pa[1][0]);

/*	Remove cycles shorter than some threshold in the graph	*/
	printf("Shaving graph...\n");
	num_vertex = rem_cycle(vertex, num_vertex);
	printf("done.\n%d vertices remained.\n", num_vertex);
	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("%d vertices %d edges (%d source %d sinks) remained.\n", num_vertex, num_edge,
		num_pa[0][1], num_pa[1][0]);

/*	Remove short edges	*/
/*
	printf("Remove short edges...\n");
	num_vertex = rem_short_edge(vertex, num_vertex);
	printf("done.\n%d vertices remained.\n", num_vertex);
*/

/*	Remove short end edges	*/
	printf("Remove shorte end dges...\n");
	num_vertex = rem_short_end(vertex, num_vertex);
	printf("done.\n%d vertices remained.\n", num_vertex);

	num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	printf("%d vertices %d edges (%d source %d sinks) remained.\n", num_vertex, num_edge,
		num_pa[0][1], num_pa[1][0]);

/*	Output graph & contigs	*/
	sprintf(temp, "%s.new", edgefile);
	fp = ckopen(temp, "w");
	sprintf(temp, "%s.new", graphfile);
	fp1 = ckopen(temp, "w");
	write_graph(vertex, num_vertex, fp, fp1);
	fclose(fp);
	fclose(fp1);

/*	Output graphviz format graph	*/

	sprintf(temp, "%s", outfile);
	fp = ckopen(temp, "w");
	output_graph(vertex, num_vertex, fp);
	fclose(fp);

	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
	free_graph(vertex, num_vertex);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	inpseq = outseq = 0;

	while ((copt=getopt(argc,argv,"i:o:e:g:l:c:d:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  sprintf(edgefile, "%s.edge", inpfile);
			  sprintf(graphfile, "%s.graph", inpfile);
			  continue;
			case 'e':
			  inpseq = 1;
			  sscanf(optarg,"%s", edgefile);
			  continue;
			case 'g':
			  sscanf(optarg,"%s", graphfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'l':
			  sscanf(optarg, "%d", &SHORTLEG);
			  continue;
			case 'c':
			  sscanf(optarg, "%d", &SHORTCYC);
			  continue;
			case 'd':
			  sscanf(optarg, "%d", &SHORTEND);
			  continue;
			default:
			  printf("simplifygraph -i InpFile -o outfile [-e edgefile -g graphfile -l shorleg -c min_cycle_length -d minimal_end_length]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-o Outfile: The output file name of contigs\n");
			  printf("-e edgefile (optional): default using inpfile.edge\n");
			  printf("-g graphfile (optional): default using inpfile.graph\n");
			  printf("-l shortleg (optional): default 50\n");
			  printf("-c min_cycle_length (optional): default 100\n");
			  printf("-d min_end_length (optional): default 50\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("simplifygraph -i InpFile -o outfile [-e edgefile -g graphfile -l shorleg -c min_cycle_length -d minimal_end_length]\n");
		printf("-i InpFile: The input file name of alignments\n");
		printf("-o Outfile: The output file name of contigs\n");
		printf("-e edgefile (optional): default using inpfile.edge\n");
		printf("-g graphfile (optional): default using inpfile.graph\n");
		printf("-l shortleg (optional): default 50\n");
		printf("-c min_cycle_length (optional): default 100\n");
		printf("-d min_end_length (optional): default 50\n");
		exit(-1);
	}
}
