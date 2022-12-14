#include <stdio.h>
#include <math.h>
#include <perdef.h>
#include <param.h>
#include <extvab.h>
#include <extfunc.h>

#define MIN_NUM 6
#define BIGGAP 5000

int min_length;

int npages = 1;

FILE	*fp1;
double	scale[2];

void outputline(FILE *fp, char **str, int num, int *base, int len, int size);
void output_one_length(FILE *fp, int base0, int length, char *name);
void output_repeat_legend(FILE *fp, int *startpos, int **breaks, int *numsegment, int **colorindex, int num);
void output_length_legend(FILE *fp);
int segcompar(SEGMENT *a, SEGMENT *b);
void output_legend(FILE *fp, int *repeat_length, int *multip, int num_repeat, int **repeats);
char ckname(char *name, char *str);
char concastr(char *str, char *newstr);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, len_seq[100], reallen[100];
	char	**chrname;
	int	*dir, **breaks, **colorindex, base[2];
	int	*multip, *multip2, *repeat_length, num_repeat, *startpos, *numsegment;
	int	num_ins, *insertpos, *insertlen;
	int	global_scale;
	SEGMENT	*segment, *insertreg, segment0;
	long int START_POS;
	int	num_segment, num_insertreg;
	int	**repeats;
	int	range[2];
	int	num_chro;
	char	name[1000], temp[1000], c, str[500], newstr[500];
	char	*label, pt;
	FILE	*fp;

	if(argc < 5)	{
		printf("Usage: makefig subrepeatfile graphvizgraph psblockshowfile min_length\n");
		exit(-1);
	}

	min_length = atoi(argv[4]);

/*	Input subrepeats	*/
	fp = ckopen(argv[2], "r");
	num_region = input_segment(segment, num_segment, chroname, fp);
	fclose(fp);

/*	Input graph	*/
	fp = ckopen(argv[1], "r");
	num_edges = input_graph(edges, sizenote, fp);
	fclose(fp);

	fp = ckopen(argv[3], "r");
	fprintf(fp, "<html>\n<head>\n<title>Subrepeats</title>\n</head>\n<body>\n");
	fprintf(fp, "<table border="1" cellspacing="0">\n");
	fprintf(fp, "<tr>\n");
	fprintf(fp, "<td width=\"80\" valign=\"top\" >\n");
	fprintf(fp, "<p> Region </p>\n");
	fprintf(fp, "</td>\n");
	fprintf(fp, "<td width=\"80\" valign=\"top\" >\n");
	fprintf(fp, "<p> Chromosome </p>\n");
	fprintf(fp, "</td>\n");
	fprintf(fp, "<td width=\"160\" valign=\"top\" >\n");
	fprintf(fp, "<p> Range(from-to)<sup>1</sup> </p>\n");
	fprintf(fp, "</td>\n");
	fprintf(fp, "</tr>\n");
	for(i = 0; i < num_region; i ++)	{
		fprintf(fp, "<tr>\n");
		fprintf(fp, "<td width=\"80\" valign=\"top\" >\n");
		fprintf(fp, "<p>%d</p>", i + 1);
		fprintf(fp, "</td>\n");
		fprintf(fp, "<td width=\"80\" valign=\"top\" >\n");
		fprintf(fp, "<p>%s</p>\n", chroname[i]);
		fprintf(fp, "</td>\n");
		fprintf(fp, "<td width=\"160\" valign=\"top\" >\n");
		fprintf(fp, "<p><a href=\"region-%d.html\">%d-%d</a></p>\n", i + 1,
			 segment[i][0].pos[0], segment[i][num_segment[i] - 1].pos[1]);
		fprintf(fp, "</td>\n");
		fprintf(fp, "</tr>\n");
		sprintf(temp, "tmp/region-%d.html", i + 1);
		fp1 = ckopen(temp, "w");
		sprintf(temp, "tmp/region-%d.gvz", i + 1);
		fp2 = ckopen(temp, "w");
		for(j = 0; j < num_segment[i]; j ++)	{
		}
		fclose(fp1);
		fclose(fp2);
	}
	fprintf(fp, "<sup>1</sup> Click to links see the sub-repeats composition of each region\n");
	fprintf(fp, "</table>\n");
	fprintf(fp, "</body>\n</html>\n");
	fclose(fp);
}

char concastr(char *str, char *newstr)
{
	int	i, l;

	l = strlen(str);
	for(i = l - 1; i >= 0; i --)	{
		if(str[i] == '\"')	{
			break;
		}
	}
	if(i >= 0)	{
		strncpy(newstr, str, i + 1);
		newstr[i + 1] = '\0';
		strcat(newstr, ", color=red, style=bold");
		strcat(newstr, &str[i + 1]);
	} else	{
		strcpy(newstr, str);
	}
}

char ckname(char *name, char *str)
{
	int	i, j, k, l, n;

	l = strlen(str);
	n = strlen(name);
	for(i = 0; i < l; i ++)	{
		if(!strncmp(&str[i], "[label = ", 9))	{
			i += 10;
			break;
		}
	}
	if(strncmp(name, &str[i], n))	{
		return(0);
	} else	{
		if(str[i + n] == ' ')	{
			return(1);
		} else	{
			return(0);
		}
	}
}

void output_repeat_legend(FILE *fp, int *startpos, int **breaks, int *numsegment, int **colorindex, int num)
{
	int	i, j, k, l, m, n;
	int	base[50], point[4][2];
	char	c, name[100], temp[100];
	char	**str;
	double	label[3];

	str = (char **) ckalloc(10 * sizeof(char *));
	for(i = 0; i < 10; i ++)	{
		str[i] = (char *) ckalloc(100 * sizeof(char));
	}
	sprintf(str[0], "Repeat ID");
	sprintf(str[1], "Length");
	sprintf(str[2], "From-to");
	sprintf(str[3], "Subrepeats");
	sprintf(str[4], "From-to");
	base[0] = 50;
	base[1] = 800;
	fprintf(fp, "Repeat ID	Length	From-to	Subrepeats	From-to\n");
	outputline(fp, str, 5, base, 90, 10);
	base[1] -= 20;
	for(i = 0; i < num; i ++)	{
		for(j = 0; j < numsegment[i]; j ++)	{
			if(j == 0)	{
				sprintf(str[0], "R%d", i + 1);
				sprintf(str[1], "%d", breaks[i][numsegment[i]]);
				sprintf(str[2], "%d-%d", startpos[i], startpos[i] + breaks[i][numsegment[i]] - 1);
				fprintf(fp1, "R%d	%d	%d-%d	", i + 1, breaks[i][numsegment[i]],
					startpos[i], startpos[i] + breaks[i][numsegment[i]] - 1);
			} else	{
				str[0][0] = str[1][0] = str[2][0] = '\0';
				fprintf(fp1, "			");
			}
			sprintf(str[3], "%d", colorindex[i][j]);
			m = (colorindex[i][j] - 1) % total_color;
			c = m + 'a';
			k = (colorindex[i][j] - 1) / total_color + 1;
			sprintf(str[3], "%c%d", c, k);
			sprintf(str[4], "%d-%d", startpos[i] + breaks[i][j], startpos[i] + breaks[i][j + 1] - 1);
			fprintf(fp1, "%c%d	%d-%d\n", c, k, startpos[i] + breaks[i][j], startpos[i] + breaks[i][j + 1] - 1);
			base[0] = 50;
			base[1] -= 10;
			outputline(fp, str, 5, base, 90, 8);
			if(base[1] <= 50)	{
				base[0] = 50;
				base[1] = 800;
				showpage(fp);
				fprintf(fp, "%%%%Page: P%d\n", npages);
				npages ++;
				change_scale(fp, scale);
			}
		}
	}
	for(i = 0; i < 10; i ++)	{
		free((void *) str[i]);
	}
	free((void *) str);
}

void output_legend(FILE *fp, int *repeat_length, int *multip, int num_repeat, int **repeats)
{
	int	i, j, k, l, m, n;
	int	base[50], point[4][2];
	char	c, name[100], temp[100];
	char	**str;
	double	label[3];

	str = (char **) ckalloc(10 * sizeof(char *));
	for(i = 0; i < 10; i ++)	{
		str[i] = (char *) ckalloc(100 * sizeof(char));
	}
	sprintf(str[0], "Subrepeat ID");
	sprintf(str[1], "Length");
	sprintf(str[2], "Multiplicity");
	sprintf(str[3], "Repeats");
	sprintf(str[4], "Color coding");
	base[0] = 50;
	base[1] = 680;
	outputline(fp, str, 5, base, 60, 10);
	for(i = 0; i < num_repeat; i ++)	{
		m = i % total_color;
		c = m + 'a';
		k = i / total_color + 1;
		sprintf(name, "%c%d", c, k);
		strcpy(str[0], name);
		sprintf(str[1], "%d", repeat_length[i]);
		sprintf(str[2], "%d", multip[i]);
		str[3][0] = '\0';
		for(j = 0; j < min(multip[i], MIN_NUM); j ++)	{
			sprintf(temp, "R%d ", repeats[i][j] + 1);
			strcat(str[3], temp);
		}
		base[0] = 50;
		base[1] -= 10;
		outputline(fp, str, 4, base, 60, 8);
		point[0][0] = 370;
		point[0][1] = base[1];
		point[1][0] = 370 + length_intv(repeat_length[i]);
		point[1][1] = base[1];
		point[2][0] = 370 + length_intv(repeat_length[i]);
		point[2][1] = base[1] + 5;
		point[3][0] = 370;
		point[3][1] = base[1] + 5;
		for(j = 0; j < 3; j ++)	{
			label[j] = color_box[m][j];
		}
		prt_box(fp, point, label, 0);
		base[0] = 370 + length_intv(repeat_length[i]) / 2;
		base[1] += 2;
		prt_Text(fp, name, base, 8, 0);
		if(base[1] <= 50 && i < num_repeat - 1)	{
			base[0] = 50;
			base[1] = 800;
			showpage(fp);
			fprintf(fp, "%%%%Page: P%d\n", npages);
			npages ++;
			change_scale(fp, scale);
		}
	}

	for(i = 0; i < 10; i ++)	{
		free((void *) str[i]);
	}
	free((void *) str);
}

void outputline(FILE *fp, char **str, int num, int *base, int len, int size)
{
	int	i, j, k, l;

	for(i = 0; i < num; i ++)	{
		prt_Text(fp, str[i], base, size, 0);
		base[0] += len;
	}
}

void output_length_legend(FILE *fp)
{
	int	i, j, k, l;
	int	base[2];

	base[0] = 50;
	base[1] = 835;
	prt_Text(fp, "Scale of repeat length", base, 10, 0);
	output_one_length(fp, 800, 10, "50-200");
	output_one_length(fp, 785, 20, "201-500");
	output_one_length(fp, 770, 30, "501-1000");
	output_one_length(fp, 755, 40, "1001-5000");
	output_one_length(fp, 740, 50, ">5000");
}

void output_one_length(FILE *fp, int base0, int length, char *name)
{
	int	point[4][2], base[2];
	double	label[3];

	base[0] = 50;
	base[1] = base0;
	prt_Text(fp, name, base, 8, 0);
	point[0][0] = 95;
	point[0][1] = base0;
	point[1][0] = 95 + length;
	point[1][1] = base0;
	point[2][0] = 95 + length;
	point[2][1] = base0 + 5;
	point[3][0] = 95;
	point[3][1] = base0 + 5;
	label[0] = label[1] = label[2] = 0;
	prt_box(fp, point, label, 0);
}

int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}
