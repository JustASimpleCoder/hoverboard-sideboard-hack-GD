#ifndef MADGWICK_FILTER_H
#define MADGWICK_FILTER_H

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#include "defines.h"



void madgwick_init(Quaternion* q);

void madgwick_update(Quaternion* q, double ax, double ay, double az, 
                     double gx, double gy, double gz, double dt);

#endif