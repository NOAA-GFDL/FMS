#include <stdio.h>

#define GETPEAKRSS_FC FC_FUNC (getpeakrss, GETPEAKRSS)
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif

size_t GETPEAKRSS_FC();

int main()
{
  size_t maxrss = GETPEAKRSS_FC();
  unsigned short return_val = 0;

  /* This test is mostly dubious, as maxrss can never be less than zero */
  if (maxrss <= 0) {
    return_val = 1;
  }

  return return_val;
}
