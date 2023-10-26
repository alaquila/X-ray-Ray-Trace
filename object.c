/*
* Aperture
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include "ray.h"

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
"      --N='[x,y,z]'       Vector defining normal (default = [0,1,0]).\n"
"      --invert            Inverts the aperture making a beam stop.\n"
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
	int invert = 0;
	int useZero = 0;
	int letterB[5][5]={
		{1,1,1,1,0},
		{1,0,0,0,1},
		{1,1,1,1,0},
                {1,0,0,0,1},
		{1,1,1,1,0}
	};
	Vector Pr = make_vector(letterB[1][1], letterB[1][2], letterB[1][3]);
	Vector Pc = make_vector(letterB[1][1], letterB[2][1], letterB[3][1]);
	Vector Normal = make_vector(0.0, 0.1, 0.0);
	// math stuff
	Vector B1, B2, B1unitV, B2unitV, N; //basis vectors & normal
	double B1mag, B2mag;

	/* Long options */
        const struct option longopts[] = {
                {"help",              0, NULL,               'h'},
                {"input",             1, NULL,               'i'},
		{"output",            1, NULL,               'o'},
		{"N",                 1, NULL,               'a'},
		{"invert",            0, NULL,               'f'},
		{"zero",              0, NULL,               'z'},
		{0, 0, NULL, 0}
	};

	 /* Short options */
        while ((c = getopt_long(argc, argv, "hi:o:a:b:c:fz",
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

			case 'f' :
			invert = 1;
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
	N = cross(Pr, Pc);
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
	Bx = letterB[x][y]-letterB[x-1][y];
	BxunitV = unit(Bx);
	Bxmag = magnitude(Bx);
	By = letterB[x][y]-letterB[x][y-1];
	ByunitV = unit(By);
	Bymag = magnitude(By);

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
			//in aperture
			if (invert == 1) myRay.i = 0.0; // is a beam block & inside beam block
		} else {
			if (invert == 0) myRay.i = 0.0; // is an aperture & outsid aperture
		}

		//print ray
		write_ray(out, myRay);
	}

	fclose(in);
	fclose(out);
	return 0;
}
