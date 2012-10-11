# Copyright 2005 Paul Zimmermann, INRIA Lorraine/LORIA.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

GMP=/usr/local
CC=gcc
CFLAGS=-g -O2 -Wall -W -Wmissing-prototypes -ansi -pedantic
N=
VERSION=0.1

##############################################################################
# this part concerns the compilation of the package

all: sieve combine buildmatrix gauss lanczos

sieve: sieve.c aux.o multiplier.o
	$(CC) $(CFLAGS) -I$(GMP)/include $< multiplier.o aux.o $(GMP)/lib/libgmp.a -lm -o $@

combine: combine.c aux.o
	$(CC) $(CFLAGS) -I$(GMP)/include $< aux.o $(GMP)/lib/libgmp.a -lm -o $@

aux.o: aux.c
	$(CC) $(CFLAGS) -I$(GMP)/include $< -c

multiplier.o: multiplier.c
	$(CC) $(CFLAGS) -I$(GMP)/include $< -c

buildmatrix: buildmatrix.c aux.o
	$(CC) $(CFLAGS) -I$(GMP)/include $< aux.o $(GMP)/lib/libgmp.a -lm -o $@

gauss: gauss.c aux.o
	$(CC) $(CFLAGS) -I$(GMP)/include $< aux.o $(GMP)/lib/libgmp.a -lm -o $@

lanczos: lanczos.c aux.o
	$(CC) $(CFLAGS) -I$(GMP)/include $< aux.o $(GMP)/lib/libgmp.a -lm -o $@

distclean:
	rm sieve combine buildmatrix gauss lanczos *.o

dist:
	mkdir mpqs-$(VERSION)
	cp COPYING Makefile README README.orig buildmatrix.c combine.c sieve.c aux.c mpqs.h gauss.c lanczos.c multiplier.c c[6-9]0 c[6-9]0.magma TIMINGS mpqs-$(VERSION)
	tar cf mpqs-$(VERSION).tar mpqs-$(VERSION)
	gzip --best mpqs-$(VERSION).tar
	rm -fr mpqs-$(VERSION)

###############################################################################
# this part concerns the use of the package to factor different numbers

factor: $(C).matrix
	./gauss -fb $(C).fbase -p $(C) -m $(C).matrix

$(C).matrix: $(C).cycles
	./buildmatrix -fb $(C).fbase -f $(C).fulls -c $(C).cycles > $(C).matrix

$(C).cycles: $(C).partials
	./combine -in $(C).partials -out $(C).cycles

$(C).partials: $(C)
	./sieve -p $(C) -n $(N)

clean:
	rm -f $(C).a_val $(C).fbase $(C).pps $(C).status $(C).fulls $(C).partials $(C).adata $(C).cycles $(C).matrix



