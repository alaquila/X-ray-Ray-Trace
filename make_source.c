/*
* make_source
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include "ray.h"
#include "rng-double.h"

enum {
	GAUSSIAN,
	UNIFORM,
	NONE
};

enum {
	R_SEED,
	R_DOUBLE,
	R_FREE
};

double* rand_double(int c, int n)
{
/* 
   An encapsulation object that returns a pointer to an array of double precision
   random numbers.  The idea is to decouple the random number generator from the
   rest of the code, making a change of random number generators easier.
   if c = R_SEED, then n = seed value
   if c = R_DOUBLE, then n = # of returned doubles
   if c = R_FREE, then n = void
*/	
	static int length = 1024; //default length
	static int index = -1;
	static double *rand_array = NULL;
	
	switch (c) {

		case R_SEED:
		ranf_start(n);
		index = 0;
		return NULL;
		break;

		case R_DOUBLE:
		if (index < 0) //initialize 
			rand_double(R_SEED, time(NULL));
		if (rand_array == NULL) 
			rand_array = malloc(length * sizeof(double));
			// NOTE: should really check return value from malloc
		if (n >= length) {
			free(rand_array);
			length = n;
			rand_array = malloc(length * sizeof(double));
			index = 0;
		}			
		if ((index == 0) || (n + index >= length)) {
			ranf_array(rand_array,length);
			index = 0;
		}
		index = index + n;		
		return rand_array + index - n; 
		//ok I could add a variable to make math easier
		break;

		case R_FREE:
		free(rand_array);
		return NULL;
		
		default:
		return NULL;
	}
	return NULL;	
}


static double gauss_rand()
{
	//creats 0 mean & unit sigma gaussian random numbers	
	//using Box-Muller transfomre (there are better methods)
	double dist, mag;
	double *var;
	do {
		var = rand_double(R_DOUBLE,2);		
		var[0] = 2.0*var[0]-1.0;
		var[1] = 2.0*var[1]-1.0;
		dist = (var[0]*var[0]) + (var[1]*var[1]);
	} while (dist >= 1.0 || dist == 0.0);
	mag = sqrt(-2.0*log(dist)/dist);
	dist = mag*var[0];
	return dist;
}

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Generates random rays for a source point.\n"
"Note: all values are metric.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"  -p, --position='[x,y,z]'  Center location for rays (default = [0,0,0]).\n"
"  -v, --direction='[x,y,z]' Center direction of travel (default = [0,0,1]).\n"
"  -d, --delta='[x,y,z]'     Width of position distribution vector.\n"
"  -s, --sigma='[x,y,z]'     Width of angular distribution vector.\n"
"  -w, --wavelength=<num>  Center wavelength of radiation (defalut is 1.0E-9 m).\n"
"  -b, --bandwidth=<num>   Bandwidth of wavelength.\n"
"  -n, --Nrays=n           Number of rays to generate (default is 1000).\n"
"  -r, --random=<type>     Use 'type' to select type of random number generator.\n"
"                          Choose from:\n"
"                            gaussian    : gaussian distributed values 1 sigam(default)\n"
"                            uniform     : unifrom random number generator\n"
"                            none        : uniformally spaced between position\n"
"                                          and position+delta (& direction+sigma)\n"
"  -i, --intensity=<num>   Intensity of the produced rays. (Constant value)\n"
"  --polarization=<num>    Value for polarization S = 1, P = -1. (Constant value)\n"
"  --seed=<num>            Seed for random number generator (default is current time)\n"
"  --start=<num>           Starting index number for ray\n"
"\n");
}

int main(int argc, char *argv[])
{

	FILE *fh = NULL;
	char *filename = NULL;
	Ray base, delta, current;
	int N = 1000;
	int i, c;
	int seed = time(NULL);
	int random = GAUSSIAN;

	/* Default settings */
	base.pos = make_vector(0.0, 0.0, 0.0);
	base.v  = make_vector(0.0, 0.0, 1.0);
	base.w     = 1.0E-9;
	base.l     = 0.0;
	base.i     = 1.0;
	base.p     = 0.0;
	base.ray   = 0;

	delta.pos = make_vector(0.0, 0.0, 0.0);
	delta.v = make_vector(0.0, 0.0, 0.0);
	delta.w     = 0.0;
	delta.l     = 0.0;
	delta.i     = 0.0;
	delta.p     = 0.0;
	delta.ray   = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"position",           1, NULL,               'p'},
		{"direction",          1, NULL,               'v'},
		{"delta",              1, NULL,               'd'},
		{"sigma",              1, NULL,               's'},
		{"wavelength",         1, NULL,               'w'},
		{"bandwidth",          1, NULL,               'b'},
		{"Nrays",              1, NULL,               'n'},
		{"random",             1, NULL,               'r'},
		{"intensity",          1, NULL,               'i'},
		{"seed",               1, NULL,               'a'},
		{"start",              1, NULL,               'y'},
		{"polarization",       1, NULL,               'c'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:p:v:d:s:w:b:n:r:i:a:b:y:",
	        longopts, NULL)) != -1)
	{

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'o' :
			filename = strdup(optarg);
			break;

			case 'p' :
			sscanf(optarg,"[%le,%le,%le]", &base.pos.x, &base.pos.y,
			       &base.pos.z);
			break;

			case 'v' :
			sscanf(optarg,"[%le,%le,%le]", &base.v.x, &base.v.y,
			       &base.v.z);
			break;

			case 'd' :
			sscanf(optarg,"[%le,%le,%le]", &delta.pos.x, &delta.pos.y,
			       &delta.pos.z);
			break;

			case 's' :
			sscanf(optarg,"[%le,%le,%le]", &delta.v.x, &delta.v.y,
			       &delta.v.z);
			break;

			case 'w' :
			base.w = atof(optarg);
			break;

			case 'b' :
			delta.w = atof(optarg);
			break;

			case 'n' :
			N = atoi(optarg);
			break;

			case 'y' :
			base.ray = atoi(optarg);
			break;

			case 'r' :
			if (strcmp(optarg, "gaussian") == 0 ) {
				random = GAUSSIAN;
			} else if (strcmp(optarg, "uniform") == 0 ) {
				random = UNIFORM;
			} else if (strcmp(optarg, "none") == 0) {
				random = NONE;
			} else {
				fprintf(stderr,"Invalid random type: '%s'\n", optarg);
				return 1;
			}
			break;

			case 'i' :
			base.i = atof(optarg);
			break;

			case 'a' :
			seed = atoi(optarg);
			break;

			case 'c' :
			base.p = atof(optarg);
			break;

			case 0 :
			break;

			default :
			fprintf(stderr,"Unexpected arguement \n");
			show_help(argv[0]);
			return 1;
		}
	}

	if ( delta.w >= base.w) {
		fprintf(stderr,"Bandwidth must be less then the central wavelength\n");
		return 1;
	}

	/* set random seed */
	rand_double(R_SEED, seed);

	/* open file or stdout */
	if ( filename != NULL ) {
		fh = fopen(filename, "w");
		if ( fh == NULL ) {
			fprintf(stderr,"Failed to open output file\n");
			return 1;
		}
		free(filename);
	} else {
		fh = stdout;
	}

	/* Ray 0 is the optical axis from starting point */
	write_ray(fh, base);

	current.i = base.i;
	current.p = base.p;
	current.l = 0.0;
	/* loop over all other rays */
	// increment 2
	for(i=2; i<=2*N; i+=2) {
		int n_rand = 7;
		double *d_rand;
		current.ray = i + base.ray;
		switch (random) {
			case GAUSSIAN:
			current.pos.x = base.pos.x + (delta.pos.x*gauss_rand());
			current.pos.y = base.pos.y + (delta.pos.y*gauss_rand());
			current.pos.z = base.pos.z + (delta.pos.z*gauss_rand());
			current.v.x = base.v.x + (delta.v.x*gauss_rand());
			current.v.y = base.v.y + (delta.v.y*gauss_rand());
			current.v.z = base.v.z + (delta.v.z*gauss_rand());
			current.w = base.w + (delta.w*gauss_rand());
			break;

			case UNIFORM:

			d_rand = rand_double(R_DOUBLE,n_rand);	
			current.pos.x = base.pos.x + (delta.pos.x*(2.0*d_rand[0]-1.0));
			current.pos.y = base.pos.y + (delta.pos.y*(2.0*d_rand[1]-1.0));
			current.pos.z = base.pos.z + (delta.pos.z*(2.0*d_rand[2]-1.0));
			current.v.x = base.v.x + (delta.v.x*(2.0*d_rand[3]-1.0));
			current.v.y = base.v.y + (delta.v.y*(2.0*d_rand[4]-1.0));
			current.v.z = base.v.z + (delta.v.z*(2.0*d_rand[5]-1.0));
			current.w = base.w + (delta.w*(2.0*d_rand[6]-1.0));
			break;

			case NONE:
			current.pos.x = base.pos.x + (delta.pos.x*i/N);
			current.pos.y = base.pos.y + (delta.pos.y*i/N);
			current.pos.z = base.pos.z + (delta.pos.z*i/N);
			current.v.x = base.v.x + (delta.v.x*i/N);
			current.v.y = base.v.y + (delta.v.y*i/N);
			current.v.z = base.v.z + (delta.v.z*i/N);
			current.w = base.w + (delta.w*i/N);
			break;

			default :
			fprintf(stderr,"Unexpected random number type.\n");
			return 1;
		}
		
		if (current.w < 0.0) current.w = current.w *  -1.0;
		current.v = unit(current.v);
		write_ray(fh, current);	
	}

	fclose(fh);
	rand_double(R_FREE,0);
	return 0;	
}
