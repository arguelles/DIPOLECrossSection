#ifndef __F2_H
#define __F2_H

#include "maddip.h"
#include "gbw.h"
#include "cgc.h"

#include "wavefunction.h"
#include "vegas.h"

#ifdef __cplusplus
extern "C" {
#endif

double CalculateF2(double x, double q2, double dd);

#ifdef __cplusplus
}
#endif

#endif
