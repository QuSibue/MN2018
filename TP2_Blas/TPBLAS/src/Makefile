# une bibliotheque statique (libmnblas.a)
# une bibliotheque dynamique (libmnblasdyn.so)

all: libmnblas.a libmnblasdyn.so

LIB_DIR=../lib
INC_DIR=../include

OPTIM=-O2
OPTIONS=-fPIC $(OPTIM) -I$(INC_DIR)

# Bibliothèque statique

libmnblas.a: copy.o swap.o dot.o nrm2.o asum.o gemv.o axpy.o gemm.o complex.o
	rm -f libmnblas.a $(LIB_DIR)/libmnblas.a
	ar -r libmnblas.a copy.o swap.o dot.o nrm2.o asum.o axpy.o gemm.o complex.o gemv.o
	cp libmnblas.a $(LIB_DIR)

# Bibliothèque dynamique

libmnblasdyn.so:copy.o swap.o dot.o nrm2.o asum.o gemv.o axpy.o gemm.o complex.o
	rm -f libmnblasdyn.so $(LIB_DIR)/libmnblasdyn.so
	gcc -shared -o libmnblasdyn.so copy.o swap.o dot.o nrm2.o asum.o gemv.o axpy.o gemm.o complex.o
	cp libmnblasdyn.so $(LIB_DIR)

# BLAS 1

copy.o : copy.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc $(OPTIONS) -c copy.c

swap.o : swap.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c swap.c

dot.o : dot.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c dot.c

nrm2.o : nrm2.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c nrm2.c -lm

asum.o : asum.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c asum.c

axpy.o : axpy.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c axpy.c -lm

complex.o : complex.c $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c complex.c -lm

# BLAS 2
gemv.o : gemv.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c gemv.c

# BLAS 3
gemm.o : gemm.c $(INC_DIR)/mnblas.h $(INC_DIR)/complex.h
	gcc  $(OPTIONS) -c gemm.c


clean:
	rm -f *.o *.a *.so *~  $(INC_DIR)/*~ $(LIB_DIR)/libmnblasdyn.so  $(LIB_DIR)/libmnblas.a
