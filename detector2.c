/*
* detector
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <complex.h>
#include "ray.h"

#define PI 3.14159265358979323846264338327

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Propagates input rays to a parallelogram detector.  The detector corrners are\n"
"Loacated at P0, P1, P2, and ( P0 + P1 + P2) The detector consists of n pixels\n"
"along the (P1-P0) direction and m pixels along the (P2-P0) direction.  The\n"
"input rays can be either added coherently or incoherently.  All rays outside\n"
"the detectors boundry are ignored in calculation.\n"
"The output is an n by m array of intensity values.\n"
"Notes: All values are metric.\n"
"       Skips rays of 0 intensity.\n"
"       Skips rays propagating parallel to detector.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"  -e, --E-range='[min,max]' Minimum and maximum Energy.\n"                    
"      --P0='[x,y,z]'      First vertex (default = [0,0,0].\n"
"      --P1='[x,y,z]'      Second vertex (default = [1,0,0].\n"
"      --p2='[x,y,z]'      Third vertex (default = [0,1,0].\n"
"  -n,                     Number of pixels between P0 and P1\n"
"                           (default 512)\n"
"  -m,                     Number of pixels between P0 and P2\n"
"                           (default 512)\n"
"      --coherent          Enables coherent addition of vectors\n"
"                            (default is to add the intensity of\n"
"                             the rays, ignoring their phase)\n"
"\n");
}

int main(int argc, char *argv[])
{
	/* initial variable */
	FILE *in = stdin;
	FILE *out = stdout;
	char *inFileName = NULL;
	char *outFileName = NULL;
	double Energy_min = 0.0;  //  TEST
	double Energy_max = 999999.0; //  TEST
	int c, n, m, k, j;
	int coherent = 0;
	Vector P0 = make_vector(0.0, 0.0, 0.0);
	Vector P1 = make_vector(1.0, 0.0, 0.0);
	Vector P2 = make_vector(0.0, 1.0, 0.0);
	Vector B1, B2, B1unitV, B2unitV, Normal;
	double B1mag, B2mag;
	double complex *detector_c; //for coherent case
	double *detector; //for incoherent case
	n = 512;
	m = 512;
	
	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"E-range",           1, NULL,               'e'},
		{"P0",                1, NULL,               'a'},
		{"P1",                1, NULL,               'b'},
		{"P2",                1, NULL,               'c'},
		{"coherent",          0, NULL,               'z'},
		{0, 0, NULL, 0}
	};


	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:e:a:b:c:n:m:z",
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

			case 'f' :
			outFileName = strdup(optarg);
			break;

			case 'e':
		        sscanf(optarg, "[%le,%le]", &Energy_min, &Energy_max);
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

			case 'n' :
			n = atoi(optarg);
			break;

			case 'm' :
			m = atoi(optarg);
			break;

			case 'z' :
			coherent = 1;
			break;

                        default :
                        fprintf(stderr,"Unexpected arguement \n");
                        return 1;
		}
	}

	//calculate cross product sanity check
	Normal = cross(sub(P1,P0),sub(P2,P0));
	if (magnitude(Normal) == 0.0) {  //see if points are colinear
		fprintf(stderr,"Error points are colinear!\n");
		return 1;
	}
	Normal = unit(Normal);

	if (dot(sub(P1,P0),sub(P2,P0)) != 0.0) {
		fprintf(stderr,"Warning: detector axes are not perpendicular!\n"
		   "Calculation will, however, continue.\n");
	}

	/* open file if necessary */
	if ( inFileName != NULL ) {
		in = fopen(inFileName, "rb");
		if ( in == NULL ) {
			fprintf(stderr,"Failed to open input file: %s!\n", inFileName);
			return 1;
		}
		free(inFileName);
	}

	if ( outFileName != NULL ) {
		out = fopen(outFileName, "w");
		if ( out == NULL ) {
			fprintf(stderr,"Failed to open output file: %s!\n", outFileName);
			return 1;
		}
		free(outFileName);
	}

	/* allocate detector array */
	if (coherent == 1) {
		detector_c = (double complex *) calloc( n * m , sizeof(double complex));
	} else {
		detector = (double *) calloc( n * m, sizeof(double));
	}
	
	/* calculate unit basis vectors */
	B1 = sub(P1,P0);
	B1unitV = unit(B1);
	B1mag = magnitude(B1);
	B2 = sub(P2,P0);
	B2unitV = unit(B2);
	B2mag = magnitude(B2);

	/* loop over vectors */
	while (!feof(in)) {				
		Ray myRay;		
		double distance, a, b, point_B1, point_B2, B1_B2;
		c = read_ray(in, &myRay);		
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}

		// Calculate energy
                double E = (12.4*1e-10)/myRay.w;
                
		// Check if energy is within the desired interval
                //if (E >= 6 && E <= 9) {
		if (E >= Energy_min && E <= Energy_max) {

			//fprintf(stderr,"E = %g\n",E); // TEST
                        Vector point;
                        myRay.v = unit(myRay.v); //check that direction is a unit vector

                        /* check if intensity is 0 if so skip */
                        if (myRay.i == 0.0) continue;

                        /* check ray is not parallel to plane */
                        if (dot(myRay.v,Normal) == 0.0) continue;

			//caluclate distance from the location point in myRay to 
                        //detector plane along propagation direction
                        distance = dot(sub(P0,myRay.pos),Normal)/dot(myRay.v,Normal);
			if (distance < 0.0) continue;
			if (dot(myRay.v,Normal) > 0.0) continue;
			myRay.pos = add(myRay.pos,mult(myRay.v,distance));
                  	myRay.l = myRay.l + distance;

                        //calculate where myRay.pos is with respect to P0 & basis vectors
                        //B1=(P1-P0), B2=(P2-P0).
                        point = sub(myRay.pos, P0);
                        point_B1 = dot(point,B1unitV);
                        point_B2 = dot(point,B2unitV);
                        B1_B2    = dot(B1unitV,B2unitV);
                        a = (point_B1-(point_B2*B1_B2))/(1-(B1_B2*B1_B2)); //Basis 1 length
                        b = point_B2-(a*B1_B2); //Basis 2 length
                        a = a/B1mag; // was using unit lenght basis vectors to make math eaier
                        b = b/B2mag; // convert back to regular size vector

                        //check if outside of detector
                       	if ( a <= 0.0 || a >= 1.0 || b <= 0.0 || b >= 1.0 ) continue;
			
	                //get detector pixel index
	                k = floor( a * m );
        	        j = floor( b * n );
                	if (coherent == 1) {
                      		double complex value;
                        	double amplitude, phase;
                                //unsigned long int n_wavelengths;
                                amplitude = sqrt(myRay.i);
                                //n_wavelengths = floor(myRay.l/myRay.w);
                                phase = 2 * PI * myRay.l / (myRay.w);
                                //phase = 2 * PI * myRay.l / (n_wavelengths * myRay.w);
                                value = (amplitude * cos(phase)) + (amplitude * sin(phase))*I;
                               	detector_c[(k*n)+j] += value;
                       	} else {
                    	       detector[(k*n)+j] += myRay.i;
                       	}
		}

        }


	/* print results */
	for( k = 0; k < m; k++ ) {
		for ( j = 0; j < n; j++ ) {
			double value;
			if (coherent == 1) {
				value = cabs(detector_c[(k*n)+j] * conj(detector_c[(k*n)+j]));
			} else {
				value = detector[(k*n)+j];
			}
			if (j != n-1) {
				fprintf(out,"%f, ", value);
			} else {
				fprintf(out,"%f\n", value);
			}
		}
	}
	
	fclose(in);
	fclose(out);
	if (coherent == 1) {
		free(detector_c);
	} else {
		free(detector);
	}
	
	return 0;
}
