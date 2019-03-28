#ifndef __DIPOLE_H
#define __DIPOLE_H

class Dipole{
    public:
      Dipole(){};
      virtual double SigmaD(double, double) { return 0.0;};
};

#endif  // __DIPOLE_H
