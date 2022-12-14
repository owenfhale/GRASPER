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
#include <math.h>
#include <perdef.h>
#include <param.h>
#include <extvab.h>
#include <extfunc.h>

#define MIN_NUM 6
#define BIGGAP 5000

int min_length = 10000;

int npages = 1;

double	scale[2];

void outputline(FILE *fp, char **str, int num, int *base, int len, int size);
void output_one_length(FILE *fp, int base0, int length, char *name);
void output_repeat_legend(FILE *fp, int *startpos, int **breaks, int *numsegment, int **colorindex, char **chrname, int num);
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
	char	**brchrname;
	int	*multip, *multip2, *repeat_length, num_repeat, *startpos, *numsegment;
	int	num_ins, *insertpos, *insertlen;
	int	global_scale;
	SEGMENT	*segment, segment0;
	long int START_POS;
	int	num_segment;
	int	**repeats;
	int	range[2];
	int	num_chro;
	char	name[1000], temp[1000], c;
	char	pt;
	FILE	*fp, *fp1;

	if(argc < 4)	{
		printf("Usage: makefig intv_file output_file(ps format) chrolistfile\n");
		exit(-1);
	}

	chrname = alloc_name(100, 100);
	fp = ckopen(argv[3], "r");
	num_chro = readchrolist(len_seq, chrname, reallen, fp);
	fclose(fp);
	segment = (SEGMENT *) ckalloc(20000 * sizeof(SEGMENT));
	fp = ckopen(argv[1], "r");
	num_segment = input_segment(segment, chrname, num_chro, fp);
	fclose(fp);

/*	Sort the segments	*/
	qsort((void *) segment, num_segment, sizeof(SEGMENT), (void *) segcompar);
/*
for(i = 0; i < num_segment; i ++)	{
	printf("aftersort segment %d %d %d\n", segment[i].pos[0], segment[i].pos[1], segment[i].eq_pos[1]);
}
getchar();
*/

/*	Compute multiplicity of each sub-repeat	*/
	repeat_length = (int *) ckalloc(20000 * sizeof(int));
	multip = (int *) ckalloc(20000 * sizeof(int));
	num_repeat = 0;
	k = 0;
	for(i = 0; i < num_segment; i ++)	{
		repeat_length[segment[i].eq_pos[0] - 1] = segment[i].length;
		multip[segment[i].eq_pos[0] - 1] ++;
		if(segment[i].eq_pos[0] > num_repeat)	{
			num_repeat = segment[i].eq_pos[0];
		}
	}
	printf("num_repeat %d\n", num_repeat);

	global_scale = 25;
	multip2 = (int *) ckalloc(num_repeat * sizeof(int));
	repeats = (int **) ckalloc(num_repeat * sizeof(int *));
	for(i = 0; i < num_repeat; i ++)	{
		repeats[i] = (int *) ckalloc((1 + multip[i]) * sizeof(int));
	}
	startpos = (int *) ckalloc(num_segment * sizeof(int));
	numsegment = (int *) ckalloc(num_segment * sizeof(int));
	breaks = (int **) ckalloc(num_segment * sizeof(int *));
	brchrname = (char **) ckalloc(num_segment * sizeof(char *));
	colorindex = (int **) ckalloc(num_segment * sizeof(int *));
	for(i = 0; i < num_segment; i ++)	{
		breaks[i] = (int *) ckalloc(num_segment * sizeof(int));
		brchrname[i] = (char *) ckalloc(100 * sizeof(char));
		colorindex[i] = (int *) ckalloc(num_segment * sizeof(int));
	}
	dir = (int *) ckalloc(num_segment * sizeof(int));
	k = m = 0;
	fp = ckopen(argv[2], "w");
	prt_Header(fp, "Tang Hai-xu", "OUTPUT", "2-21-2004", 0);
	prt_Macro(fp);
	scale[0] = 0.7;
	scale[1] = 0.8;
	change_scale(fp, scale);
	base[0] = 50;
	base[1] = 850;
	for(i = 0; i < num_segment; i ++)	{
		if(i == num_segment - 1 || segment[i + 1].pos[0] > segment[i].pos[1] + min_length ||
		   segment[i + 1].chro != segment[i].chro)	{
			for(j = k; j <= i; j ++)	{
				breaks[m][j - k + 1] += segment[j].pos[1] - segment[k].pos[0] + 1;
				colorindex[m][j - k] = segment[j].eq_pos[0];
				dir[j - k] = segment[j].eq_pos[1];
				repeats[segment[j].eq_pos[0] - 1][multip2[segment[j].eq_pos[0] - 1] ++] = m;
			}
			startpos[m] = range[0] = segment[k].pos[0];
			numsegment[m] = i - k + 1;
			strcpy(brchrname[m], chrname[segment[i].chro]);
			range[1] = segment[i].pos[1];
			sprintf(name, "R%d", m + 1);
			l = output_region(fp, base, name, chrname[segment[k].chro], range,
				 global_scale, breaks[m], i - k + 2, colorindex[m], dir);
			m ++;
			base[0] = 50;
			base[1] -= 25;
			if(base[1] < 150)	{
				showpage(fp);
				fprintf(fp, "%%%%Page: P%d\n", npages);
				npages ++;
				change_scale(fp, scale);
				base[0] = 50;
				base[1] = 850;
			}
			k = i + 1;
		}
	}
	showpage(fp);
	fprintf(fp, "%%%%Page: P%d\n", npages);
	npages ++;
	change_scale(fp, scale);
	fflush(fp);
	for(i = 0; i < num_repeat; i ++)	{
		if(multip[i] != multip2[i])	{
			printf("i %d multip %d %d\n", i, multip[i], multip2[i]);
			exit(0);
		}
	}
	output_length_legend(fp);
	output_legend(fp, repeat_length, multip, num_repeat, repeats);
	showpage(fp);
	fprintf(fp, "%%%%Page: P%d\n", npages);
	npages ++;
	change_scale(fp, scale);
	output_repeat_legend(fp, startpos, breaks, numsegment, colorindex, brchrname, m);
	showpage(fp);
	page_trailer(fp);
	restore_scale(fp);
	prt_Trailer(fp, k + 1);
	fclose(fp);
	printf("%d regions output.\n", m);

	chrname = free_name(chrname, 100);
	for(i = 0; i < num_repeat; i ++)	{
		free((void *) repeats[i]);
	}
	free((void **) repeats);
	free((void *) repeat_length);
	for(i = 0; i < num_segment; i ++)	{
		free((void *) brchrname[i]);
		free((void *) breaks[i]);
		free((void *) colorindex[i]);
	}
	free((void *) multip);
	free((void *) multip2);
	free((void **) brchrname);
	free((void **) breaks);
	free((void **) colorindex);
	free((void *) startpos);
	free((void *) numsegment);
	free((void *) dir);
	free((void *) segment);
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

void output_repeat_legend(FILE *fp, int *startpos, int **breaks, int *numsegment, int **colorindex, char **chrname, int num)
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
	sprintf(str[2], "Chromosome");
	sprintf(str[3], "From-to");
	sprintf(str[4], "Subrepeats");
	sprintf(str[5], "From-to");
	base[0] = 50;
	base[1] = 800;
	outputline(fp, str, 6, base, 90, 10);
	base[1] -= 20;
	for(i = 0; i < num; i ++)	{
		for(j = 0; j < numsegment[i]; j ++)	{
			if(j == 0)	{
				sprintf(str[0], "R%d", i + 1);
				sprintf(str[1], "%d", breaks[i][numsegment[i]]);
				sprintf(str[2], "%s", chrname[i]);
				sprintf(str[3], "%d-%d", startpos[i], startpos[i] + breaks[i][numsegment[i]] - 1);
			} else	{
				str[0][0] = str[1][0] = str[2][0] = str[3][0] = '\0';
			}
			m = (colorindex[i][j] - 1) % total_color;
			sprintf(str[4], "s%d", colorindex[i][j]);
			sprintf(str[5], "%d-%d", startpos[i] + breaks[i][j], startpos[i] + breaks[i][j + 1] - 1);
			base[0] = 50;
			base[1] -= 10;
			outputline(fp, str, 6, base, 90, 8);
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
		if(multip[i] == 0)	continue;
		sprintf(name, "s%d", i + 1);
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
		m = i % total_color;
		for(j = 0; j < 3; j ++)	{
			label[j] = color_box[m][j];
		}
		prt_box(fp, point, label, 0, 0, 0);
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
	output_one_length(fp, 800, 10, "100-200");
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
	prt_box(fp, point, label, 0, 0, 0);
}

int segcompar(SEGMENT *a, SEGMENT *b)
{
	if(a -> chro > b -> chro)	return(1);
	else if(a -> chro == b -> chro && a -> pos[0] > b -> pos[0])	return(1);
	else				return(-1);
}
