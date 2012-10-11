CC=gcc
CFLAGS=-Wall -Wextra -g -lgmp
CKATTIS = -O2 -std=gnu99 -static 

all: build

build: main.c
	$(CC) $(CFLAGS) main.c -o factoring

test: matrix_generator.c
	$(CC) $(CFLAGS) matrix_generator.c -o matrix_generator

generate_matrix_test:
	./matrix_generator

run_test:
	time --output time_${SIZE}.time ./factoring < ../data/data_in_${SIZE}.dat | sort > answer.tmp 
	sort ../data/data_out_${SIZE}.dat > reference.tmp
	diff answer.tmp reference.tmp

submit:
	./submit.py -f --problem=oldkattis:factoring main.c primes.h

clean:
	rm -rf factoring generator matrix_generator
