#ifndef RAN_H
#define RAN_H 

double std_normal(int* idum);

// Linear congruential generator
// Returns a uniform random deviate
// between 0.0 and 1.0. Set idum to
// any integer value to initialize
// the sequence, but do not alter its
// value between successive calls.
double ran(int* idum);

// Inverse Transform Method
// for standard normal variable
double normal(int* idum, int* count);

// Acceptance-Rejection Method
// for standard normal variable
double arnormal(int* idum, int* count);

// Box-Muller Method
// for standard normal variable
double bmnormal(int* idum, int* count);

#endif /* RAN_H */
