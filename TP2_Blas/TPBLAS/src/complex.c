#include "complex.h"


/*
 * ===========================================================================
 * Opération complex Simple
 * ===========================================================================
 */

struct complex_simple addition_cs(struct complex_simple c1,
                                   struct complex_simple c2) {
  struct complex_simple c;
  c.real = c1.real + c2.real;
  c.imaginary = c1.imaginary + c2.imaginary;
  return c;
}

struct complex_simple soustraction_cs(struct complex_simple c1,
                                       struct complex_simple c2) {

  struct complex_simple c;
  c.real = c1.real - c2.real;
  c.imaginary = c1.imaginary - c2.imaginary;
  return c;
}

struct complex_simple conjugue_cs(struct complex_simple c1) {

  struct complex_simple c;
  c.real = c1.real;
  c.imaginary = -c1.imaginary;
  return c;
}

struct complex_simple multiplication_cs(struct complex_simple c1,
                                         struct complex_simple c2) {

  struct complex_simple c;
  // partie real = x.real*y.rel - x.img*y.img

  c.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
  // partie img = x.real*y.img + x.img*y.real

  c.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;
  return c;
}

/*
 * ===========================================================================
 * Opération complex Double
 * ===========================================================================
 */

struct complex_double addition_cd(struct complex_double c1,
                                   struct complex_double c2) {
  struct complex_double c;
  c.real = c1.real + c2.real;
  c.imaginary = c1.imaginary + c2.imaginary;
  return c;
}

struct complex_double soustraction_cd(struct complex_double c1,
                                       struct complex_double c2) {

  struct complex_double c;
  c.real = c1.real - c2.real;
  c.imaginary = c1.imaginary - c2.imaginary;
  return c;
}

struct complex_double conjugue_cd(struct complex_double c1) {

  struct complex_double c;
  c.real = c1.real;
  c.imaginary = -c1.imaginary;
  return c;
}

struct complex_double multiplication_cd(struct complex_double c1,
                                         struct complex_double c2) {

  struct complex_double c;
  // partie real = x.real*y.rel - x.img*y.img

  c.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
  // partie img = x.real*y.img + x.img*y.real

  c.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;
  return c;
}
