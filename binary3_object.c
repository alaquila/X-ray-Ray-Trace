/*
* Aperture
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
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
"Creates a parallelogram aperture and propagates all rays to the aperture.\n"
"Sets all rays outside aperture to 0 intensity. Aperture defined by the \n"
"parallelogram with vertices p0, p1, p2, and (p0 + p1 + p2).\n"
"Note: all values are metric.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --P0='[x,y,z]'      First vertex (default = [0,0,0].\n"
"      --P1='[x,y,z]'      Second vertex (default = [1,0,0].\n"
"      --p2='[x,y,z]'      Third vertex (default = [0,1,0].\n"
"      --seed=<num>        Seed for random number generator (default is current time)\n"
"      --invert            Inverts the object making a beam stop.\n"
"  -p, --absorb            absorption factor (default = 0.1).\n"
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
        double absorb_factor = 0.1;
	int invert = 0;
	int useZero = 0;
	Vector P0 = make_vector(0.0, 0.0, 0.0);
	Vector P1 = make_vector(1.0, 0.0, 0.0);
	Vector P2 = make_vector(0.0, 1.0, 0.0);
	// math stuff
	Vector B1, B2, B1unitV, B2unitV, N; //basis vectors & normal
	double B1mag, B2mag;

	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"P0",                1, NULL,               'a'},
		{"P1",                1, NULL,               'b'},
		{"P2",                1, NULL,               'c'},
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
			sscanf(optarg,"[%le,%le,%le]", &P0.x, &P0.y, &P0.z);
			break;

			case 'b' :
			sscanf(optarg,"[%le,%le,%le]", &P1.x, &P1.y, &P1.z);
			break;

			case 'c' :
			sscanf(optarg,"[%le,%le,%le]", &P2.x, &P2.y, &P2.z);
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
	
	//calculate cross product sanity check
	N = cross(sub(P1,P0),sub(P2,P0));
	if (magnitude(N) == 0.0) {  //see if points are colinear
		fprintf(stderr,"Error points are colinear!\n");
		return 1;
	}
	N=unit(N);

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

	/* calculate unit basis vectors */
	B1 = sub(P1,P0);
	B1unitV = unit(B1);
	B1mag = magnitude(B1);
	B2 = sub(P2,P0);
	B2unitV = unit(B2);
	B2mag = magnitude(B2);

//	int image[256][256];
//	int i;
//	int j;
//	for (i=0 ; i<256; i++){
//		for(j=0; j<256; j++){
//			image[i][j] = 1;
//		}	
//	}

	const char *filename = "image";
	FILE *infile = fopen(filename, "rb");
	if (infile == NULL) {
		fprintf(stderr, "Failed to open input file: %s!\n", filename);
		return 1;
	}
	unsigned char object[65536];
	int bytes_read = fread(object, 1, sizeof(object), infile);
	
	if (bytes_read == 0) {
		fprintf(stderr, "Failed to read from the file.\n");
		fclose(infile);
		return 1;
	}

	fclose(infile);

//	int j;
//	for (j = 0; j< 65536; j++){     
//		fprintf(stderr, "object, bytes_read:  %x\n", object[j]);
//	}


	/* loop over vectors */
	while (!feof(in)) {				
		double distance, m, n, point_B1, point_B2, B1_B2, *d_rand;
		int pixel1, pixel2, index;
		Ray myRay;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		Vector point;
		myRay.v = unit(myRay.v); //check that direction is a unit vector
		
		/* check if using 0 or skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		/* check ray is not parallel to aperture */
		if (dot(myRay.v,N) == 0.0) {
			fprintf(stderr,"Warning: Ray %lu propagates parallel to"
			 " aperture!\n Skipping ray, but copying it to output.\n", myRay.ray);
			write_ray(out, myRay);
			continue;
		} 

		//caluclate distance from the location point in myRay to 
		//aperture plane along propagation direction
		distance = dot(sub(P0,myRay.pos),N)/dot(myRay.v,N);
		if (distance < 0.0) {
			fprintf(stderr,"Warning Ray %lu is furthern then the aperture"
			 " back propagating to aperture\n", myRay.ray);
			continue; // remove the rays are further than the object and back propagating to the object
		}
		myRay.pos = add(myRay.pos,mult(myRay.v,distance));
		myRay.l = myRay.l + distance;

		//calculate where myRay.pos is with respect to P0 & basis vectors
		//B1=(P1-P0), B2=(P2-P0).
		point = sub(myRay.pos, P0);
		point_B1 = dot(point,B1unitV);
		point_B2 = dot(point,B2unitV);
		B1_B2    = dot(B1unitV,B2unitV);		
		m = (point_B1-(point_B2*B1_B2))/(1-(B1_B2*B1_B2)); //Basis 1 length
		n = point_B2-(m*B1_B2); //Basis 2 length
		m = m/B1mag; // was using unit lenght basis vectors to make math eaier
		n = n/B2mag; // convert back to regular size vector
	
                if ( m >= 0.0 && m <= 1.0 && n >= 0.0 && n <= 1.0 ) {
                       //in object

                        pixel1 = floor(m*256);
                        pixel2 = floor(n*256);

                        index = (256 * pixel2) + pixel1;

                        //fprintf(stderr,"point_B1, point_B2, pixel1, pixel2, index: %f %f %d %d %d\n", point_B1, point_B2, pixel1, pixel2, index); // TEST

                        if (invert == 0) {//is a beam block & inside beam bloc

                                if (object[index] == 1) {
                                        d_rand =  rand_double(R_DOUBLE,1);
                                        if (d_rand[0] < absorb_factor)
                                                myRay.i = 0.0; // Absorb the ray 
                                }
                	}
		} else {

                        pixel1 = floor(m*256);
                        pixel2 = floor(n*256);

                        index = (256 * pixel2) + pixel1;

                        if (invert == 1) { //is a beam block & inside beam block

                                if (object[index] == 1) {
                                        d_rand =  rand_double(R_DOUBLE,1);
                                        if (d_rand[0] < absorb_factor)
                                                myRay.i = 0.0; // Absorb the ray 
                                }
                        }
                }
                //print ray
                write_ray(out, myRay);
        }

        fclose(in);
        fclose(out);
        return 0;
}
	
