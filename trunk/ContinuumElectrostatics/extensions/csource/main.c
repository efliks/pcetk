# include "StateVector.h"

int main (void ) {
  StateVector *v;

  v = StateVector_Allocate (100);

  StateVector_Reset (v);

  StateVector_Deallocate (v);

  return 1;
}
