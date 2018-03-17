all: mntest mntestdyn

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-O2 -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  =   -L$(LIB_DIR) $(LIBST)  -lgsl -lgslcblas
OPTIONS_LINK_DYN  = -L$(LIB_DIR) $(LIBDYN)  -lgsl -lgslcblas

mntest: mntest.o $(LIB_DIR)/libmnblas.a
	gcc -o mntest mntest.o $(OPTIONS_LINK_STATIC)

mntestdyn: mntest.o $(LIB_DIR)/libmnblasdyn.so
	gcc -o mntestdyn mntest.o $(OPTIONS_LINK_DYN)

mntest.o: mntest.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c mntest.c

clean:
	rm -f *.o mntest mntestdyn *~
