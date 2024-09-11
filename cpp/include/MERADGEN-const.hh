#ifndef __MERADGEN_HH
#define __MERADGEN_HH

/************************************************************
 *  TODO: 
 * M: MASS OF ELECTRON
 * M2: MASS OF ELECTRON SQUARED
 * PI: NOT SURE ABOUT THIS ONE
 * ALFA: FINE STRUCTURE CONSTANT
 * BARN: BARN CONVERSION
 * 
 * EN: ?
 * S: HOLDER FOR CENTER OF MASS ENERGY???
 * COEB: ?
 * COER: ? 
 * EGMIN: ?
 * 
 * SHOULD I SET VALUES HERE -- CONSTANTS SURE... WHAT ABOUT OTHERS???
 ***********************************************************/

struct MERADGENgr {
	int iy;

    double m;
    double m2;
    double als;
    double pi;
    double alfa;
    double barn;
    double coeb;
    double coer;
    double egmin;
    double en;
    double s;
};

#endif
