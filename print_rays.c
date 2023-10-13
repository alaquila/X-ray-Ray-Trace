/*
* print_rays
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include "ray.h"


static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Converts binary Ray format to a human readable CSV string format.\n"
"It also converst text files to Ray format using the --text flag.\n"
"Note: All values are metric.\n"
"\n"
"  -h, --help               Display this help message.\n"
"  -i, --input=<file>       Input filename. Default: stdin.\n"
"  -o, --output=<file>      Output filename. Default: stdout.\n"
"  -n, --precision=<int>    Specifies the number of digits to be printed after\n"
"                             the decimal point (defalut = 15).\n"
"  -t, --text               Flag to convert a CSV text input file to a Ray\n"
"                             output file.\n"
"\n");
}

int main(int argc, char *argv[])
{
	/* initial variable */
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;
	char line[RAY_STRING_LENGTH];	
	Ray myRay;
	int c;
	int precision = 15;
	int textIN = 0;
	
	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"precision",         1, NULL,               'n'},
		{"text",              0, NULL,               't'},
		{0, 0, NULL, 0}
	};


	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:n:t",
                longopts, NULL)) != -1)
        {

                switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			inFileName = strdup(optarg);
			break;

			case 'o' :
			outFileName = strdup(optarg);
			break;

			case 'n' :
			precision = atoi(optarg);
			break;

			case 't' :
			textIN = 1;
			break;

			default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;
		}
	}
	/* open files if necessary */
	if ( inFileName != NULL ) {
		if (textIN == 0) in = fopen(inFileName, "rb");
		else             in = fopen(inFileName, "r" );
		if ( in == NULL ) {
			fprintf(stderr,"Failed to open input file: %s!\n", inFileName);
			return 1;
		}
		free(inFileName);
	}
	

	if ( outFileName != NULL ) {
		if (textIN == 0) out = fopen(outFileName, "w" );
		else             out = fopen(outFileName, "wb");	
		if ( out == NULL ) {
			fprintf(stderr,"Failed to open output file: %s!\n", outFileName);
			return 1;
		}
		free(outFileName);
	}

	if (textIN == 0) {
		write_header(line);
		fprintf(out, "%s", line);
		/* loop over vectors */
		if (feof(in)) return 0;
		do {
			c = read_ray(in, &myRay);
			if (c != 0) {
				if (feof(in)) break;
				fprintf(stderr,"Unexpected error reading input file!\n");
				return 1;
			}
			write_ray_string(line, myRay, precision);
			fprintf(out, "%s", line);
		} while (!feof(in));
	} else {
		fgets(line, RAY_STRING_LENGTH - 1, in); //read header line
		while (fgets(line, RAY_STRING_LENGTH - 1, in) != NULL) {	
			myRay = read_ray_string(line);
			c = write_ray(out, myRay);
						if (c != 0) {
				fprintf(stderr,"Unexpected error writing output file!\n");
				return 1;
			}
		}
	}
	fclose(in);
	fclose(out);
	return 0;
}
		
