all: strassen matrix_generator

strassen: strassen.c
	gcc -std=c11 strassen.c -o strassen

matrix_generator: matrix_generator.c
	gcc -std=c11 matrix_generator.c -o matrix_generator