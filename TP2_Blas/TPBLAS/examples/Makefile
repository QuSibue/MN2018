all: mntest mnswap mnasum mnnrm2 mndot mnaxpy

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas -lm
LIBDYN=-lmnblasdyn -lm

OPTIONS_COMPIL  =-O2 -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  = -L$(LIB_DIR) $(LIBST) -lgsl -lgslcblas -lm
OPTIONS_LINK_DYN  = -L$(LIB_DIR) $(LIBDYN) -lgsl -lgslcblas -lm

mntest: mntest.o $(LIB_DIR)/libmnblas.a fonctions_test.o
	gcc -o mntest mntest.o fonctions_test.o $(OPTIONS_LINK_STATIC)

mnswap: mnswap.o $(LIB_DIR)/libmnblas.a
	gcc -o mnswap mnswap.o $(OPTIONS_LINK_STATIC)

mnasum: mnasum.o $(LIB_DIR)/libmnblas.a fonctions_test.o
	gcc -o mnasum mnasum.o fonctions_test.o $(OPTIONS_LINK_STATIC)

mnaxpy: mnaxpy.o $(LIB_DIR)/libmnblas.a fonctions_test.o
	gcc -o mnaxpy mnaxpy.o fonctions_test.o $(OPTIONS_LINK_STATIC)

mnnrm2: mnnrm2.o $(LIB_DIR)/libmnblas.a fonctions_test.o
	gcc -o mnnrm2 mnnrm2.o fonctions_test.o $(OPTIONS_LINK_STATIC)

mnswapdyn: mnswap.o $(LIB_DIR)/libmnblasdyn.so
	gcc -o mnswapdyn mnswap.o $(OPTIONS_LINK_DYN)

mndot: mndot.o $(LIB_DIR)/libmnblas.a fonctions_test.o
	gcc -o mndot mndot.o fonctions_test.o $(OPTIONS_LINK_STATIC)

mntestdyn: mntest.o $(LIB_DIR)/libmnblasdyn.so fonctions_test.o
	gcc -o mntestdyn mntest.o fonctions_test.o $(OPTIONS_LINK_DYN)

mnnrm2dyn: mnnrm2.o $(LIB_DIR)/libmnblasdyn.so fonctions_test.o
	gcc -o mnnrm2dyn mnnrm2.o fonctions_test.o$(OPTIONS_LINK_DYN)

mntest.o: mntest.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h $(INC_DIR)/fonctions_test.h
	gcc $(OPTIONS_COMPIL) -c mntest.c

mnasum.o: mnasum.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h $(INC_DIR)/fonctions_test.h
	gcc $(OPTIONS_COMPIL) -c mnasum.c

mnaxpy.o: mnaxpy.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h $(INC_DIR)/fonctions_test.h
	gcc $(OPTIONS_COMPIL) -c mnaxpy.c

mnnrm2.o: mnnrm2.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h $(INC_DIR)/fonctions_test.h
	gcc $(OPTIONS_COMPIL) -c mnnrm2.c

mnswap.o: mnswap.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h $(INC_DIR)/fonctions_test.h
	gcc $(OPTIONS_COMPIL) -c mnswap.c

mndot.o: mndot.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h $(INC_DIR)/fonctions_test.h
	gcc $(OPTIONS_COMPIL) -c mndot.c

fonctions_test.o: fonctions_test.c $(INC_DIR)/complex.h $(INC_DIR)/fonctions_test.h
	gcc $(OPTIONS_COMPIL) -c fonctions_test.c

clean:
	rm -f *.o mntest mntestdyn *~
