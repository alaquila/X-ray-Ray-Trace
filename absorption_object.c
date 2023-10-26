/*
* object
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include "ray.h"
#include "rng-double.h"

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
        }
        return NULL;
}

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Creates an object and propagates all rays to the aperture.\n"
"Sets all rays outside object to 0 intensity. object defined by the \n"
"center position C, the radius  R (on the object), and the normal\n"
"(perpendicular) vector to the object N.\n"
"Note: all values are metric.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --C='[x,y,z]'       Center of object (default = [0,0,0]).\n"
"      --R=r               Radius of object (default = 1.0).\n"
"      --N='[x,y,z]'       Vector defining normal (default = [0,1,0]).\n"
"      --seed=<num>        Seed for random number generator (default is current time)\n"
"      --invert            Inverts the aperture making a beam stop.\n"
"  -p,    --absorb            absorption factor (default = 0.1).\n"
"  -z, --zero              Enables propagation of 0 intensity rays\n"
"                            By default not propaged.\n"
"\n");
}

int main(int argc, char *argv[])
{

        /* initial variable */
        //file stuff
        FILE *in = stdin;
        FILE *out = stdout;
        char *inFileName = NULL;
        char *outFileName = NULL;
        // input options stuff
        int c;
	int seed = time(NULL);
        int invert = 0;
 	double absorb_factor = 0.1;
	int useZero = 0;
        Vector Center = make_vector(0.0, 0.0, 0.0);
        Vector Normal = make_vector(0.0, 1.0, 0.0);
        double Radius = 1.0;

        /* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
                {"output",            1, NULL,               'o'},
                {"N",                 1, NULL,               'a'},
                {"R",                 1, NULL,               'b'},
                {"C",                 1, NULL,               'c'},
                {"invert",            0, NULL,               'f'},
		{"absorb",            1, NULL,               'p'},
		{"seed",              1, NULL,               's'},
		{"zero",              0, NULL,               'z'},
                {0, 0, NULL, 0}
        };

         /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:c:f:p:s:z",
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

                        case 'a' :
                        sscanf(optarg,"[%le,%le,%le]", &Normal.x, &Normal.y, &Normal.z);
                        break;

                        case 'b' :
                        Radius = atof(optarg);
                        break;

                        case 'c' :
                        sscanf(optarg,"[%le,%le,%le]", &Center.x, &Center.y, &Center.z);
                        break;

                        case 'f' :
                        invert = 1;
                        break;
			
			case 'p' :
                        absorb_factor = atof(optarg);
                        break;
			
			case 's' :
                        seed = atoi(optarg);
                        break;

                        case 'z' :
                        useZero = 1;
                        break;

                        default :
                        fprintf(stderr,"Unexpected arguement \n");
                        show_help(argv[0]);
                        return 1;

                }
        }


         /* set random seed */
        rand_double(R_SEED, seed);	


	/* open files if necessary */
	if ( inFileName != NULL ) {
                in = fopen(inFileName, "rb");
                if ( in == NULL ) {
                        fprintf(stderr,"Failed to open input file: %s!\n", inFileName);
                        return 1;
                }
                free(inFileName);
        }

        if ( outFileName != NULL ) {
                out = fopen(outFileName, "wb");
                if ( out == NULL ) {
                        fprintf(stderr,"Failed to open output file: %s!\n", outFileName);
                        return 1;
                }
                free(outFileName);
        }

        Normal = unit(Normal); // make positive N is a unit normal vector

        /* loop over vectors */
        while (!feof(in)) {
                double distance, *d_rand;
                Ray myRay;
                c = read_ray(in, &myRay);
                if (c != 0) {
                        if (feof(in)) break;
                        fprintf(stderr,"Unexpected error reading input file!\n");
                        return 1;
                }
                myRay.v = unit(myRay.v); //check that direction is a unit vector

                /* check if using 0 or skip */
                if (useZero == 0 && myRay.i == 0.0) {
                        write_ray(out, myRay);
                        continue;
                }

                /* check ray is not parallel to aperture */
                if (dot(myRay.v,Normal) == 0.0) {
                        fprintf(stderr,"Warning: Ray %lu propagates parallel to"
                         " object!\n Skipping ray, but copying it to output.\n", myRay.ray);
                        write_ray(out, myRay);
                        continue;
                }

                //caluclate distance from the location point in myRay to 
                //aperture plane along propagation direction
                distance = dot(sub(Center,myRay.pos),Normal)/dot(myRay.v,Normal);
                if (distance < 0.0) {
                        fprintf(stderr,"Warning Ray %lu is furthern then the object"
                         " back propagating to object\n", myRay.ray);
                }
                myRay.pos = add(myRay.pos,mult(myRay.v,distance));
                myRay.l = myRay.l + distance;

	        //calculate where myRay.pos is with respect to aperture center and if inside radius
                if (Radius >= magnitude(sub(myRay.pos,Center))) {
			//in aperture
			if (invert == 0) { //is a beam block & inside beam block
		
		                d_rand =  rand_double(R_DOUBLE,1);
                                if (d_rand[0] < absorb_factor)
                                       myRay.i = 0.0; // Absorb the ray 
			}
	
		} else {
		
			if (invert == 1) { // is an aperture & outsid aperture
		
				d_rand =  rand_double(R_DOUBLE,1);
                                if (d_rand[0] < absorb_factor)
                                       myRay.i = 0.0; // Absorb the ray 
			}
		}

		//print ray
		write_ray(out, myRay);	}

	fclose(in);
	fclose(out);
	return 0;
}
