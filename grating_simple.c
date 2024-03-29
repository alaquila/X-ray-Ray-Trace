/*
* grating_simple
NOTE: Check phase change upon reflection!!!!
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "ray.h"
#define PI 3.14159265358979323846264338327

static void show_help(const char *s)
{
printf("Syntax: %s [options]\n\n", s);
printf(
"Creates a rectangular refletive grating with a constant d spacing.\n"
"Diffracts all rays in grating aperture with order m, and it sets all\n"
"rays outside aperture to 0 intensity. The grating is defined by a\n"
"rectangle with with vertices P0, P1, P2, and (P0 + P1 + P2). With\n"
"P1-P0 the grating direction (perindicular to grating grooves).\n"
"Note: all values are metric.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --P0='[x,y,z]'      First vertex (default = [0,0,0]).\n"
"      --P1='[x,y,z]'      Second vertex (default = [1,0,0]).\n"
"      --p2='[x,y,z]'      Third vertex (default = [0,1,0]).\n"
"  -d, --dspacing=<num>    D-spacing [m] along the P1-P0\n"
"                            direction (default = 1000.0).\n"
"  -m, --order=<int>       Grating order (default = 1)\n"
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
	int useZero = 0;
	int order = 1;
	double dSpacing = 1000.0;
	Vector P0 = make_vector(0.0, 0.0, 0.0);
	Vector P1 = make_vector(1.0, 0.0, 0.0);
	Vector P2 = make_vector(0.0, 1.0, 0.0);

	// math stuff
	Vector B1, B2, B1unitV, B2unitV, Normal; //basis vectors & normal
	double B1mag, B2mag;

	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"P0",                1, NULL,               'a'},
		{"P1",                1, NULL,               'b'},
		{"P2",                1, NULL,               'c'},
		{"dspacing",          1, NULL,               'd'},
		{"order",             1, NULL,               'm'},
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:c:d:m:z",
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

			case 'd' :
			dSpacing = atof(optarg);
			break;

			case 'm' :
			order = atoi(optarg);
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
	
	//calculate cross product sanity check
	Normal = cross(sub(P1,P0),sub(P2,P0));

	if (magnitude(Normal) == 0.0) {  //see if points are colinear
		fprintf(stderr,"Error points are colinear!\n");
		return 1;
	}
	Normal = unit(Normal);

	if (dot(sub(P1,P0),sub(P2,P0)) != 0.0) { //points are not perpendicular
		fprintf(stderr,"Error points are not perpendicular!\n");
		// NOTE: suggest direction to set P2?
		return 1;
	}

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
	
	/* loop over vectors */
	while (!feof(in)) {				
		double distance, m, n, point_B1, point_B2, B1_B2;
		Ray myRay;
		Vector point;
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		
		/* check if using 0 or skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		/* check ray is not parallel to mirror */
		if (dot(myRay.v,Normal) == 0.0) {
			fprintf(stderr,"Warning: Ray %lu propagates parallel to"
			 " grating surface!\n Skipping ray, but copying it to output.\n", myRay.ray);
			write_ray(out, myRay);
			continue;
		} 

		myRay.v = unit(myRay.v); //make the direction a unit vector
		//caluclate distance from the location point in myRay to 
		//aperture plane along propagation direction
		distance = dot(sub(P0,myRay.pos),Normal)/dot(myRay.v,Normal);
		if (distance < 0.0) {
			fprintf(stderr,"Warning Ray %lu is furthern then the grating's surface"
			 " back propagating to grating.\n", myRay.ray);
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
			//in grating
			double i,j,k, theta, alpha, beta;
			Vector Out;
		
			myRay.v = mult(myRay.v, -1.0);  //Note math easier if I flip my incident ray
			if(acos(dot(Normal,myRay.v)) > PI/2.0) 
				Normal = mult(Normal, -1.0); //Normal facing wrong way

			//Calculate angles			
			theta = atan(dot(B2unitV, myRay.v) / dot(Normal, myRay.v)); //incidence angle perp to grating
			alpha = atan(dot(B1unitV, myRay.v) / dot(Normal, myRay.v)); //incidence angle parallel to grating
			beta = asin((order*myRay.w/dSpacing) - sin(alpha)); //diffraction angle

			i = cos((PI/2.0)-beta);
			j = cos((PI/2.0)-theta);
			k = sqrt(1.0 - (i*i) - (j*j)); 
			//Note there are two solutions (-k works too)! 
			//I am taking the positive "reflective" one as Normal is defined outside surface is positive

			// i*B1unitV + j*B2unitV + k*Normal = Out_vector
			Out = add(mult(B1unitV,i),mult(B2unitV,j));
			Out = add(Out,mult(Normal, k));
			Out = unit(Out); // just checking, but needless if math above is correct
			myRay.v = Out;
//fprintf(stderr,"alpha = %f, beta = %f, theta = %f\n", alpha*180/PI, beta*180/PI, theta*180/PI);
//fprintf(stderr,"i = %f, j = %f, k = %f\n",i,j,k);
//fprintf(stderr,"Out = (%f, %f, %f)\n",Out.x, Out.y, Out.z);
		} else {
			myRay.i = 0.0; // ray is outside mirror aperture
		}

		//print ray
		write_ray(out, myRay);
	}

	fclose(in);
	fclose(out);
	return 0;
}
