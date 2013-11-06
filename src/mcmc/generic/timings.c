#include <time.h>
#include <sys/timeb.h>
#include "timings.h"

double toMilliseconds(struct timeb t_start, struct timeb t_end)
{
  int t_diff = (int) (1000.0 * (t_end.time - t_start.time)  + (t_end.millitm - t_start.millitm));
  return((double) t_diff);
}
