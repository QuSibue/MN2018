#ifndef COMPLEX_H
#define COMPLEX_H

struct complex_simple{
  float real;
  float imaginary;
};

struct complex_double{
  double real;
  double imaginary;
};


/*
 * ===========================================================================
 * Opération Sur Des Complex Simple
 * ===========================================================================
 */

struct complex_simple addition_cs(struct complex_simple c1, struct complex_simple c2);

struct complex_simple soustraction_cs(struct complex_simple c1,struct complex_simple c2);

struct complex_simple conjugue_cs(struct complex_simple c1);

struct complex_simple multiplication_cs(struct complex_simple c1, struct complex_simple c2);

/*
 * ===========================================================================
 * Opération Sur Des complex Double
 * ===========================================================================
 */

 struct complex_double addition_cd(struct complex_double c1, struct complex_double c2);

 struct complex_double soustraction_cd(struct complex_double c1,struct complex_double c2);

 struct complex_double conjugue_cd(struct complex_double c1);

 struct complex_double multiplication_cd(struct complex_double c1, struct complex_double c2);

#endif
