/*
* flat_mirror
NOTE: NO phase change upon reflection included yet!!!!
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
"Creates a parallelogram flat mirror and propagates all rays to the mirror.\n"
"All rays in mirror aperture are reflected, and it sets all rays outside\n"
"mirror aperture to 0 intensity. The mirror is defined by the parallelogram\n"
"with vertices p0, p1, p2, and (p0 + p1 + p2).\n"
"Note: all values are metric.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<file>      Input filename. Default: stdin.\n"
"  -o, --output=<file>     Output filename. Default: stdout.\n"
"      --P0='[x,y,z]'      First vertex (default = [0,0,0].\n"
"      --P1='[x,y,z]'      Second vertex (default = [1,0,0].\n"
"      --p2='[x,y,z]'      Third vertex (default = [0,1,0].\n"
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
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:c:z",
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
		c = read_ray(in, &myRay);
		if (c != 0) {
			if (feof(in)) break;
			fprintf(stderr,"Unexpected error reading input file!\n");
			return 1;
		}
		Vector point;
		myRay.v = unit(myRay.v); //make the direction a unit vector
		
		/* check if using 0 or skip */
		if (useZero == 0 && myRay.i == 0.0) {
			write_ray(out, myRay);
			continue;
		}

		/* check ray is not parallel to mirror */
		if (dot(myRay.v,Normal) == 0.0) {
			fprintf(stderr,"Warning: Ray %lu propagates parallel to"
			 " mirror surface!\n Skipping ray, but copying it to output.\n", myRay.ray);
			write_ray(out, myRay);
			continue;
		} 

		//caluclate distance from the location point in myRay to 
		//aperture plane along propagation direction
		distance = dot(sub(P0,myRay.pos),Normal)/dot(myRay.v,Normal);
		if (distance < 0.0) {
			fprintf(stderr,"Warning Ray %lu is furthern then the mirror's surface"
			 " back propagating to mirror.\n", myRay.ray);
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
			//in mirror
			double k, j, theta, VdotN;
			Vector Out;
			myRay.v = mult(myRay.v, -1.0);  //Note math easier if I flip my incident ray
			if(acos(dot(Normal,myRay.v)) > PI/2.0) 
				Normal = mult(Normal, -1.0); //Normal facing wrong way

			theta = acos(dot(Normal, myRay.v));

			VdotN = dot(myRay.v, Normal); //note both unit vectors
			if (VdotN == 1.0) {
				//vectors are colinear! Normal reflection!
				k = 1;
				j = 0;
			} else {
				k = ((cos(2.0*theta)*VdotN)-cos(theta))/(pow(VdotN,2.0) - 1.0);
				j = cos(2.0*theta) - (k*VdotN);

			} 
			Out = add(mult(myRay.v,j),mult(Normal,k));
			Out = unit(Out); // just checking, but needless otherwise math is wrong
			myRay.v = Out;		
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
