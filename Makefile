all: make_source propagate_by propagate_to optic_baffle aperture detector reflective_optic transform print_rays flat_mirror grating_simple grating_flat circ_aperture grep_rays set_rays

print_rays : vector.o ray.o
	gcc -lm -Wall print_rays.c vector.o ray.o -o print_rays

set_rays : vector.o ray.o
	gcc -lm -Wall set_rays.c vector.o ray.o -o set_rays

grep_rays : vector.o ray.o
	gcc -lm -Wall grep_rays.c vector.o ray.o -o grep_rays

reflective_shape : vector.o ray.o
	gcc -lm -Wall reflective_shape.c vector.o ray.o -o reflective_shape

reflective_optic : vector.o ray.o
	gcc -lm -Wall reflective_optic.c vector.o ray.o -o reflective_optic

flat_mirror : vector.o ray.o
	gcc -lm -Wall flat_mirror.c vector.o ray.o -o flat_mirror

grating_flat : vector.o ray.o
	gcc -lm -Wall grating_flat.c vector.o ray.o -o grating_flat

grating_simple : vector.o ray.o
	gcc -lm -Wall grating_simple.c vector.o ray.o -o grating_simple

detector : vector.o ray.o
	gcc -lm -Wall detector.c vector.o ray.o -o detector

transform : vector.o ray.o
	gcc -lm -Wall transform.c vector.o ray.o -o transform

aperture : vector.o ray.o
	gcc -lm -Wall aperture.c vector.o ray.o -o aperture

circ_aperture : vector.o ray.o
	gcc -lm -Wall circ_aperture.c vector.o ray.o -o circ_aperture

optic_baffle: vector.o ray.o
	gcc -lm -Wall optic_baffle.c vector.o ray.o -o optic_baffle

propagate_to: vector.o ray.o
	gcc -lm -Wall propagate_to.c vector.o ray.o -o propagate_to

propagate_by: vector.o ray.o
	gcc -lm -Wall propagate_by.c vector.o ray.o -o propagate_by

make_source: vector.o ray.o rng-double.o
	gcc -lm -Wall make_source.c rng-double.o vector.o ray.o -o make_source

test: vector.o
	gcc -lm -Wall testing.c vector.o -o test

random:
	gcc -lm -Wall -c rng-double.c rng-double.o

ray: vector.o
	gcc -lm -Wall -c ray.c vector.o

vector: 
	gcc -lm -Wall -c vector.c
clean:
	rm -rf *.o 
	rm -rf test make_source propagate_by propagte_to optic_baffle aperture
	rm -rf detector reflective_optic print_rays transform propagate_to grep_rays
	rm -rf flat_mirror grating_simple grating_flat reflective_shape circ_aperture
	rm -rf set_rays	
