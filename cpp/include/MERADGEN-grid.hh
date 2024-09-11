#ifndef __MERADGEN_GRID_HH
#define __MERADGEN_GRID_HH

#include <vector>

//TODO: WHAT THE HECK ARE THESE???
const int nv = /* define nv here */;
const int nt1 = /* define nt1 here */;
const int nz = /* define nz here */;

struct MERADGENgr {
    std::vector<double> grv(nv);
    std::vector<double> grt1(nt1);
    std::vector<double> grz(nz);
};

#endif
