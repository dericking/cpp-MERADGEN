#include <cmath>
#include <cstdio>
#include <iostream>

#include "meradgen_globals.hpp"
#include "meradgen_constants.hpp"
#include "meradgen_api.hpp"
#include "meradgen_parity_trace.hpp"

namespace meradgen {

double fsir(double t_in, double t1, double v, double z, double pl_in, int& nn, int ikey) {
    // Map Fortran argument names to more direct variables to follow original formulas.
    const double t  = t_in;        // original Mandelstam t
    const double pl = pl_in;

    // Local kinematic variables (from Fortran fsir.f)
    double u, uu, vv, v1, vv1, z1, z2, zz, zz1, zz2;
    double dz, dz1, vt, pls, tt, tt1, dt, sd, ds, fir;
    double sr1, sr2, sr3, sr4, sr5, sr6, sr7, sr8, sr9, sr10;
    double aj1, aj2, aj3, aj4, aj5, aj6, aj7, aj8, aj9, aj10;
    double aj11, aj12, aj13, aj14, aj15, aj16, aj17, aj18, aj19, aj20;
    double aj21, aj22, aj23, aj24, aj25, aj26, aj27, aj28, aj29, aj30;
    double aj31, aj32, aj33, aj34, aj35, aj36, aj37, aj38, aj39, aj40;
    double aj41, aj42, aj43, aj44, aj45, aj46, aj47, aj48, aj49, aj50;
    double aj51;

    u   = v - s - t + 4.0 * m2;
    uu  = 1.0 / u;
    v1  = s + u + t1 - 4.0 * m2;
    z1  = z - t1 + t;
    z2  = z - s - t1 + 4.0 * m2;
    zz  = 1.0 / z;
    zz1 = 1.0 / z1;
    zz2 = 1.0 / z2;
    vv1 = 1.0 / v1;
    vv  = 1.0 / v;
    pls = -pl / als;
    tt  = 1.0 / t;
    tt1 = 1.0 / t1;
    dt  = 1.0 / (t - t1);
    dz  = 1.0 / (s + t1 - 4.0 * m2);
    dz1 = 1.0 / (s + t - 4.0 * m2);
    ds  = 1.0 / (s - 4.0 * m2);
    vt  = 1.0 / (v - t);

    if (ikey == 0) {
        // Explicit (z*z): same association as FORTRAN z**2. Left-associative
        // ((-az)*z)*z is what meradgen-cpp/ uses; under ~5e8:1 cancellation it
        // was the 30363 hex residual. This is the conventional grouping, not a
        // FORTRAN-matching rewrite of the formula.
        sd = 1.0 / std::sqrt(-az * (z * z) - 2.0 * bz * z - cz);

        aj1  = sd;
        aj2  = z * sd;
        aj3  = z * z * sd;
        aj4  = zz * sd;
        aj5  = m2 * zz * zz * sd;
        aj6  = zz1 * sd;
        aj7  = m2 * zz1 * zz1 * sd;
        aj8  = zz2 * sd;
        aj9  = m2 * zz2 * zz2 * sd;
        aj10 = t1 * aj1;
        aj11 = tt1 * aj1;
        aj12 = m2 * tt1 * tt1 * aj1;
        aj13 = vv1 * aj1;
        aj14 = m2 * vv1 * vv1 * aj1;
        aj15 = tt1 * aj2;
        aj16 = vv1 * aj2;
        aj17 = m2 * vv1 * vv1 * aj2;
        aj18 = vv1 * aj3;
        aj19 = m2 * vv1 * vv1 * aj3;
        aj20 = t1 * aj4;
        aj21 = tt1 * aj4;
        aj22 = m2 * tt1 * tt1 * aj4;
        aj23 = dt * (aj4 - aj6);
        aj24 = vv1 * aj4;
        aj25 = dz * aj4;
        aj26 = m2 * dz * dz * aj4;
        aj27 = m2 * dz * dz * dz * aj4;
        aj28 = tt1 * aj5;
        aj29 = tt1 * tt1 * aj5;
        aj30 = dz * aj5;
        aj31 = dz * dz * aj5;
        aj32 = t1 * aj6;
        aj33 = t1 * t1 * aj6;
        aj34 = tt1 * aj6;
        aj35 = m2 * tt1 * tt1 * aj6;
        aj36 = vv1 * aj6;
        aj37 = t1 * aj7;
        aj38 = t1 * t1 * aj7;
        aj39 = tt1 * aj7;
        aj40 = tt1 * tt1 * aj7;
        aj41 = t1 * aj8;
        aj42 = tt1 * aj8;
        aj43 = vv1 * aj8;
        aj44 = m2 * vv1 * vv1 * aj8;
        aj45 = dz * aj8;
        aj46 = m2 * dz * dz * aj8;
        aj47 = m2 * dz * dz * dz * aj8;
        aj48 = vv1 * aj9;
        aj49 = vv1 * vv1 * aj9;
        aj50 = dz * aj9;
        aj51 = dz * dz * aj9;

    } else if (ikey == 1) {
        aj1  = pi / std::sqrt(az);
        aj2  = -pi * bz / std::pow(az, 1.5);
        aj3  = pi * (3.0 * bz * bz - az * cz) / 2.0 / std::pow(az, 2.5);
        aj4  = pi / std::sqrt(cz);
        aj5  = -m2 * pi * bz / std::pow(cz, 1.5);
        aj6  = pi / std::sqrt(cz1);
        aj7  = -m2 * pi * bz1 / std::pow(cz1, 1.5);
        aj8  = -pi / std::sqrt(cz2);
        aj9  = m2 * pi * bz2 / std::pow(cz2, 1.5);
        aj10 = t1 * aj1;
        aj11 = tt1 * aj1;
        aj12 = m2 * tt1 * tt1 * aj1;
        aj13 = vv1 * aj1;
        aj14 = m2 * vv1 * vv1 * aj1;
        aj15 = tt1 * aj2;
        aj16 = vv1 * aj2;
        aj17 = m2 * vv1 * vv1 * aj2;
        aj18 = vv1 * aj3;
        aj19 = m2 * vv1 * vv1 * aj3;
        aj20 = t1 * aj4;
        aj21 = tt1 * aj4;
        aj22 = m2 * tt1 * tt1 * aj4;
        aj23 = dt * (aj4 - aj6);
        aj24 = vv1 * aj4;
        aj25 = dz * aj4;
        aj26 = m2 * dz * dz * aj4;
        aj27 = m2 * dz * dz * dz * aj4;
        aj28 = tt1 * aj5;
        aj29 = tt1 * tt1 * aj5;
        aj30 = dz * aj5;
        aj31 = dz * dz * aj5;
        aj32 = t1 * aj6;
        aj33 = t1 * t1 * aj6;
        aj34 = tt1 * aj6;
        aj35 = m2 * tt1 * tt1 * aj6;
        aj36 = vv1 * aj6;
        aj37 = t1 * aj7;
        aj38 = t1 * t1 * aj7;
        aj39 = tt1 * aj7;
        aj40 = tt1 * tt1 * aj7;
        aj41 = t1 * aj8;
        aj42 = tt1 * aj8;
        aj43 = vv1 * aj8;
        aj44 = m2 * vv1 * vv1 * aj8;
        aj45 = dz * aj8;
        aj46 = m2 * dz * dz * aj8;
        aj47 = m2 * dz * dz * dz * aj8;
        aj48 = vv1 * aj9;
        aj49 = vv1 * vv1 * aj9;
        aj50 = dz * aj9;
        aj51 = dz * dz * aj9;

    } else if (ikey == 2 || ikey == -1) {
        aj1  = pi*v/(m2+v);
        aj2  = -pi*std::pow(v,2)*(v-s+2.0*m2)/2.0/std::pow((m2+v),2);
        aj3  = 0.0;
        aj4  = pi/std::sqrt(std::pow((s-v),2)-4.0*m2*s)*std::log(std::pow((s-v-2.0*m2+std::sqrt(std::pow((s-v),2)-4.0*m2*s)),2)/4.0/m2/(m2+v));
        aj5  = pi/v;
        aj6  = pi/std::sqrt(std::pow((v-u),2)-4.0*m2*u)*std::log(std::pow((v-u+2.0*m2+std::sqrt(std::pow((v-u),2)-4.0*m2*u)),2)/4.0/m2/(m2+v));
        aj7  = pi/v;
        aj8  = pi*std::log(4.0*m2*std::pow(u,2)*(m2+v)/std::pow((v*(v-u+std::sqrt(std::pow((v-u),2)-4.0*m2*u))-2.0*m2*u),2))/std::sqrt(std::pow((v-u),2)-4.0*m2*u);
        aj9  = pi*v/std::pow(u,2);
        aj10 = pi*v*(2.0*m2*t+(t-v)*v)/2.0/std::pow((m2+v),2);
        aj11 = pi*std::log(4.0*m2*std::pow(t,2)*(m2+v)/std::pow((v*(v-t+std::sqrt(std::pow((v-t),2)-4.0*m2*t))-2.0*m2*t),2))/std::sqrt(std::pow((v-t),2)-4.0*m2*t);
        aj12 = pi*v/std::pow(t,2);
        aj13 = pi*std::log(std::pow((v-t+4.0*m2+std::sqrt(std::pow((v-t),2)-4.0*m2*t)),2)/4.0/m2/(m2+v))/std::sqrt(std::pow((v-t),2)-4.0*m2*t);
        aj14 = pi/v;
        aj15 = pi*v*(v*(2.0*m2-t+v)-s*(t+v))/(std::pow((v-t),2)-4.0*m2*t)/(m2+v)+pi*t*(s*(t-v)+2.0*m2*v)/std::pow((std::pow((v-t),2)-4.0*m2*t),(1.5))*std::log(4.0*m2*std::pow(t,2)*(m2+v)/std::pow((v*(v-t+std::sqrt(std::pow((v-t),2)-4.0*m2*t))-2.0*m2*t),2));
        aj16 = pi*v*(v*(2.0*m2-t+v)-s*(t+v))/(std::pow((v-t),2)-4.0*m2*t)/(m2+v)-pi*v*(u*(t-v)+2.0*m2*v)/std::pow((std::pow((v-t),2)-4.0*m2*t),(1.5))*std::log(4.0*m2*(m2+v)/std::pow((v-t+std::sqrt(std::pow((v-t),2)-4.0*m2*t)-2.0*m2),2));
        aj17 = -pi*(1.0+(s*(t-v)+2.0*m2*v)/(std::pow((v-t),2)-4.0*m2*t))+m2*pi*(s*(t+v)+v*(t-v-2.0*m2))/std::pow((std::pow((v-t),2)-4.0*m2*t),(1.5))*std::log(4.0*m2*(m2+v)/std::pow((v-t+std::sqrt(std::pow((v-t),2)-4.0*m2*t)-2.0*m2),2));
        aj18 = pi*((std::pow(v,2)*(4.0*ipow3(m2)*(4.0*s*t+(8.0*t-3.0*v)*v) -2.0*std::pow(m2,2)*(2.0*s*t*(s+7.0*t)-2.0*(s - 6.0*t)*t*v -(4.0*s+31.0*t)*std::pow(v,2)+13*ipow3(v))+(t-v)*(3.0*std::pow((t-v),2)*std::pow(v,2)+2.0*s*(t-v)*v*(2.0*t+3.0*v)+std::pow(s,2)*(-std::pow(t,2)+4.0*t*v+3.0*std::pow(v,2)))+2.0*m2*(std::pow(s,2)*(4.0*std::pow(t,2)-t*v-std::pow(v,2))+s*(t+v)*(3.0*std::pow(t,2)-15.0*t*v+8.0*std::pow(v,2))+(t-v)*v*(2.0*std::pow(t,2)-13.0*t*v+8.0*std::pow(v,2)))))/(2.0*std::pow((-4.0*m2*t+std::pow((t-v),2)),2)*std::pow((m2+v),2))-(std::pow(v,2)*(std::pow(u,2)*std::pow((t-v),2)+6.0*std::pow(m2,2)*std::pow(v,2)-2.0*m2*u*(s*t+2.0*v*(v-t)))*std::log((4*m2*(m2+v))/std::pow((2.0*m2-t+std::sqrt(-4.0*m2*t+std::pow((t-v),2))+v),2)))/std::pow((-4.0*m2*t+std::pow((t-v),2)),2.5));
        aj19 = pi*((2.0*v*(std::pow((u-4.0*m2),2)*std::pow((t-v),2)*v-4.0*ipow3(m2)*(4.0*(s-t)*t+4.0*t*v-3.0*std::pow(v,2))+4.0*std::pow(m2,2)*((s-2.0*t)*t*(s+t)-3.0*(s-3.0*t)*t*v-(2.0*s+9.0*t)*std::pow(v,2)+4.0*ipow3(v))+m2*(2.0*std::pow(s,2)*std::pow((t+v),2)+2.0*s*(t-v)*(std::pow(t,2)-3.0*t*v+4.0*std::pow(v,2))+std::pow((t-v),2)*(std::pow(t,2)-10*t*v+6*std::pow(v,2)))))/(std::pow((std::pow((t-v),2)-4.0*m2*t),2)*(m2+v))+m2*(4.0*v*(-2.0*std::pow(m2,2)*(4.0*s*t+(4.0*t-3.0*v)*v)+m2*(2.0*s*t*(s+t)+2.0*t*(s+t)*v-3.0*t*std::pow(v,2)+ipow3(v))+u*(t-v)*((t-v)*v+s*(2.0*t+v)))*std::log((4.0*m2*(m2+v))/std::pow((2.0*m2-t+std::sqrt(-4*m2*t+std::pow((t-v),2))+v),2)))/std::pow((std::pow((t-v),2)-4.0*m2*t),2.5))/2.0;
        aj20 = -((s*(t+v)-v*(v-t+2.0*m2))*aj1+(s*t*(v-s)+m2*(4.0*s*t-2.0*std::pow(v,2)))*aj4)/(std::pow((s-v),2)-4.0*m2*s);
        aj21 = pi/t*std::log(((s-2.0*m2)*(s-2.0*m2+std::sqrt(als)))/2.0/std::pow(m2,2)-1.0)/std::sqrt(als);
        aj22 = -pi/ipow3(t)/als*((s*(v-t)-2.0*m2*v)*v+std::log(((s-2.0*m2)*((s-2.0*m2)+std::sqrt(als)))/2.0/std::pow(m2,2)-1.0)*(t*(s*v-als)-2.0*m2*std::pow(v,2))*m2/std::sqrt(als));
        aj23 = 2.0*pi*std::log((std::sqrt(std::pow(t,2)-4.0*m2*t)-t+2.0*m2)/2.0/m2)/v/std::sqrt(std::pow(t,2)-4.0*m2*t);
        aj24 = 2.0*pi*std::log((std::sqrt(std::pow(u,2)-4.0*m2*u)-u+2.0*m2)/2.0/m2)/v/std::sqrt(std::pow(u,2)-4.0*m2*u);
        aj25 = -pi*std::log((std::pow(m2,2)*(std::pow(v,2)-2.0*std::pow(u,2))+std::pow((m2*v-(s-2.0*m2)*u),2)+std::sqrt(als)*u*(s*u-2.0*m2*(u+v)))/2.0/m2/(m2*std::pow((u+v),2)-s*u*v))/u/std::sqrt(als);
        aj26 = (-pi/als/std::pow(u,2)*(v*(s*u*(s+u-2.0*v)+2.0*m2*(v*(u+v)-2.0*s*u))/(s*u*v-m2*std::pow((u+v),2)))+(s*u*(v-s)+m2*(4.0*s*u-2.0*std::pow(v,2)))*aj25/std::pow(u,2)/als)*m2;
        // ikey=2 aj27: gfortran multiply tree, not std::pow. Event 138994.
        // Original (pow) kept for examination:
        // aj27 = (-pi*v*(ipow3(s)*ipow3(u)*v*(std::pow(u,2)-3.0*std::pow(s,2)-2.0*s*u+2.0*(5.0*s+u)*v-6.0*std::pow(v,2))+2.0*m2*std::pow(s,2)*std::pow(u,2)*(2.0*s*std::pow(u,2)*(s+u)-u*(s*u-15.0*std::pow(s,2)+4.0*std::pow(u,2))*v+(std::pow(s,2)-34.0*s*u+2.0*std::pow(u,2))*std::pow(v,2)-(11.0*s - 5.0*u)*ipow3(v) + 12.0*ipow4(v)) + 2.0*std::pow(m2,2)*s*u*(2.0*std::pow(u,2)*(13.0*s*u-24.0*std::pow(s,2)+std::pow(u,2))*v-8.0*s*ipow3(u)*(2.0*s + u)-2.0*s*(4.0*s-37.0*u)*u*std::pow(v,2)+(62.0*s-11.0*u)*u*ipow3(v) +6.0*(s-4.0*u)*ipow4(v)-15.0*ipow5(v))+4.0*ipow3(m2)*(u+v)*(3.0*ipow3(v)*std::pow((u+v),2)-12.0*s*u*v*std::pow((u+v),2)+8.0*std::pow(s,2)*std::pow(u,2)*(2.0*u+v))))*m2/(2.0*ipow4(u)*std::pow((s*u*v-m2*std::pow((u+v),2)),2)*std::pow(als,2)) +(aj25*(2.0*m2*s*u*(3.0*s+u-3.0*v)*std::pow(v,2)+6.0*std::pow(m2,2)*std::pow(v,2)*(std::pow(v,2)-4.0*s*u)+std::pow(u,2)*std::pow((als-s*v),2)))*m2/(ipow4(u)*std::pow(als,2));
        const double s2 = s * s;
        const double s3 = s2 * s;
        const double u2 = u * u;
        const double u3 = u2 * u;
        const double u4 = u2 * u2;
        const double v2 = v * v;
        const double v3 = v2 * v;
        const double v4 = v2 * v2;
        const double v5 = v4 * v;
        const double m2_2 = m2 * m2;
        const double m2_3 = m2_2 * m2;
        const double upv = u + v;
        const double upv2 = upv * upv;
        const double als2 = als * als;
        const double suv_minus = s * u * v - m2 * upv2;
        const double suv_minus2 = suv_minus * suv_minus;
        const double als_minus_sv = als - s * v;
        const double als_minus_sv2 = als_minus_sv * als_minus_sv;

        const double aj27_poly =
            s3 * u3 * v * (u2 - 3.0 * s2 - 2.0 * s * u + 2.0 * (5.0 * s + u) * v - 6.0 * v2)
            + 2.0 * m2 * s2 * u2 * (2.0 * s * u2 * (s + u) - u * (s * u - 15.0 * s2 + 4.0 * u2) * v + (s2 - 34.0 * s * u + 2.0 * u2) * v2 - (11.0 * s - 5.0 * u) * v3 + 12.0 * v4)
            + 2.0 * m2_2 * s * u * (2.0 * u2 * (13.0 * s * u - 24.0 * s2 + u2) * v - 8.0 * s * u3 * (2.0 * s + u) - 2.0 * s * (4.0 * s - 37.0 * u) * u * v2 + (62.0 * s - 11.0 * u) * u * v3 + 6.0 * (s - 4.0 * u) * v4 - 15.0 * v5)
            + 4.0 * m2_3 * upv * (3.0 * v3 * upv2 - 12.0 * s * u * v * upv2 + 8.0 * s2 * u2 * (2.0 * u + v));
        const double aj27_den = 2.0 * u4 * suv_minus2 * als2;
        const double aj27_term1 = (-pi * v * aj27_poly) * m2 / aj27_den;
        const double aj27_aj25_poly =
            2.0 * m2 * s * u * (3.0 * s + u - 3.0 * v) * v2
            + 6.0 * m2_2 * v2 * (v2 - 4.0 * s * u)
            + u2 * als_minus_sv2;
        const double aj27_term2 = (aj25 * aj27_aj25_poly) * m2 / (u4 * als2);
        aj27 = aj27_term1 + aj27_term2;
        // FORTRAN fsir.f:231,238 have /2/m2 on this log; same log in aj21/aj22 uses /2/m2**2. See README_MERADGEN_TYPO.md.
        aj28 = m2*pi/std::pow(t,2)/als*((s*(v-t)-2.0*m2*v)*std::log((s-2.0*m2)*(s-2.0*m2+std::sqrt(als))/2.0/std::pow(m2,2)-1.0)/std::sqrt(als)-(t*(s*v-als)-2.0*m2*std::pow(v,2))/m2/v);
        aj29 = -pi/ipow4(t)/std::pow(als,(2.5))*(std::sqrt(als)*(als*std::pow(t,2)*(2.0*s*v-als)+std::pow(v,2)*(2.0*s*t*v*(s+2.0*m2)-std::pow(v,2)*(std::pow((s-2.0*m2),2)+8.0*m2)-2.0*s*t*(s*t+2.0*(2.0*s+t)*m2-16.0*std::pow(m2,2))))/v+2.0*m2*(2.0*std::pow(m2,2)*v*(8.0*s*t-3.0*std::pow(v,2))+s*v*m2*(3.0*std::pow(v,2)+2.0*t*(v-2.0*s-t))-s*t*(2.0*s*u*v+t*als))*std::log((s-2.0*m2)*(s-2.0*m2+std::sqrt(als))/2.0/std::pow(m2,2)-1.0));
        aj30 = -(pi*(als*u-s*u*v+2.0*m2*std::pow(v,2))/v/als-m2*(m2*(4.0*s*u-2.0*v*(u+v))-s*u*(s+u-2.0*v))*aj25)/std::pow(u,2);
        aj31 = pi*(u*(std::pow(u,2)/v+(std::pow(u,2)*(v-4.0*m2-2.0*s))/als+12.0*v*(-s*t*u+m2*std::pow(v,2))*m2/std::pow(als,2)+(std::pow(u,2)*(std::pow((s-u),2)*v+4.0*m2*(std::pow(u,2)+(u-s)*v)))*m2/(als*(-s*u*v+m2*std::pow((u+v),2))))+(2.0*(3.0*s*u*(s*u-8.0*std::pow(m2,2)+2.0*m2*s)*std::pow(v,2)-9.0*m2*s*u*ipow3(v)+s*als*std::pow(u,2)*(v-t)+6.0*std::pow(m2,2)*ipow3(v)*(u+v) -2.0*s*(s-m2)*std::pow(u,2)*v*(2.0*(v-t)-u))*std::log(std::pow((std::sqrt(als)*u+s*u-2.0*m2*(u+v)),2)/(4.0*m2*(-s*u*v+m2*std::pow((u+v),2)))))*m2/std::pow(als,2.5))/ipow5(u);
        aj32 = ((u*(t+v)-v*(v-t+2.0*m2))*aj1-(u*t*(v-u)+m2*(4.0*u*t-2.0*std::pow(v,2)))*aj6)/(std::pow((u-v),2)-4.0*m2*u);
        aj33 = (-(pi*v*(4.0*m2*(4.0*m2-u)*std::pow(u,2)*std::pow((t-v),2)+(t-v)*u*(32.0*ipow3(m2)-4.0*std::pow(m2,2)*(5.0*s-24.0*u)+3.0*std::pow(u,2)*(s+u)-2.0*m2*u*(6.0*s+19.0*u))*v+u*(128.0*ipow3(m2)-4.0*std::pow(m2,2)*(33.0*s-40.0*u)-(s-9.0*u)*u*(s+u)+2.0*m2*(11.0*std::pow(s,2)-6.0*s*u-42.0*std::pow(u,2)))*std::pow(v,2)-(20.0*ipow3(m2)-2.0*std::pow(m2,2)*(8.0*s+5.0*u)+2.0*m2*(std::pow(s,2)+19.0*s*u-21.0*std::pow(u,2))+u*(9.0*std::pow(u,2)-5.0*std::pow(s,2)-2.0*s*u))*ipow3(v)-(14.0*std::pow(m2,2)+std::pow(s,2)+4.0*s*u-3.0*std::pow(u,2)+6.0*m2*(u-2.0*s))*ipow4(v)+2.0*m2*ipow5(v)))/(2.0*std::pow((m2+v),2))+aj6*(std::pow(v,2)*(6.0*std::pow(m2,2)*std::pow(v,2)-(s*t*u*(2.0*m2+u)))+t*u*(u-4.0*m2)*(ipow3(v)-2.0*t*u*v+(u-4.0*m2)*(t*u - std::pow(v,2)))))/std::pow((std::pow((u - v),2)-4.0*m2*u),2);
        aj34 = -pi*2.0*std::log(2.0*m2/(2.0*m2-u+std::sqrt(u*(u-4.0*m2))))/t/std::sqrt(u*(u-4.0*m2));
        aj35 = -pi*(v*(u*(t-v)+2.0*m2*v)+2.0*m2*(t*u*(s+t)-2.0*m2*std::pow(v,2))*std::log((std::sqrt(u*(u-4.0*m2))-u)/2.0/m2-1.0)/std::sqrt(u*(u-4.0*m2)))/ipow3(t)/u/(u-4.0*m2);
        aj36 = pi*2.0*std::log((s-2.0*m2+std::sqrt(als))/2.0/m2)/v/std::sqrt(als);
        aj37 = pi*((t/v+(s*(u-v)+2.0*m2*v)/(4.0*m2*u-std::pow((u-v),2))-1.0)+m2*(t*u-s*v+2.0*m2*v)/std::pow((std::pow((u-v),2)-4.0*m2*u),(1.5))*std::log(std::pow((v-u+2.0*m2+std::sqrt(std::pow((v-u),2)-4.0*m2*u)),2)/4.0/m2/(m2+v)));
        aj38 = pi*((2.0*(4.0*m2*u-std::pow((u-v),2))*(m2*std::pow(t,2)*std::pow(u,2)*std::pow((u-4.0*m2),2)+std::pow(t,2)*std::pow(u,2)*(24.0*std::pow(m2,2)-10.0*m2*u+std::pow(u,2))*v+2*t*u*(m2*(22.0*m2-5.0*s)*u-2.0*std::pow(m2,2)*(4.0*m2+s)+(s-9.0*m2)*std::pow(u,2)+ipow3(u))*std::pow(v,2)-t*u*(16.0*std::pow(m2,2)+6.0*m2*(s-3.0*u)+u*(s+3.0*u))*ipow3(v)+(12.0*ipow3(m2)-std::pow(u,2)*(s+u)-4.0*std::pow(m2,2)*(s+4.0*u)+m2*(std::pow(s,2)+4.0*s*u+8*std::pow(u,2)))*ipow4(v)+(8.0*std::pow(m2,2)-4.0*m2*u+std::pow(u,2))*ipow5(v)))/(v*(m2+v))+4.0*m2*std::sqrt(std::pow((u-v),2)-4.0*m2*u)*(std::pow(t,2)*(4.0*m2-u)*std::pow(u,2)+t*u*(8.0*std::pow(m2,2)+(s-u)*u+m2*(2.0*u-6.0*s))*v+t*u*(u-2.0*s)*std::pow(v,2)-m2*(2.0*(m2-s)+u)*ipow3(v)+m2*ipow4(v))*std::log(std::pow((2.0*m2-u+std::sqrt(std::pow((u-v),2)-4.0*m2*u)+v),2)/(4.0*m2*(m2+v))))/(2.0*ipow3((4.0*m2*u-std::pow((u-v),2))));
        aj39 = pi/std::pow(t,2)/(u*(u-4.0*m2))*((u*(v-t)-2.0*m2*v)*2.0*std::log(2.0*m2/(2.0*m2-u+std::sqrt(u*(u-4.0*m2))))*m2/std::sqrt(u*(u-4.0*m2))-(u*t*(v-u)+m2*(4.0*u*t-2.0*std::pow(v,2)))/v);
        aj40 = pi*((u*(std::pow(t,2)*std::pow(u,2)*(u-4.0*m2)*(u-4.0*m2-2.0*v)-2.0*t*u*(8.0*std::pow(m2,2)+2.0*m2*(s-3.0*u)+u*(s+u))*std::pow(v,2)+(12.0*std::pow(m2,2)-4.0*m2*u+std::pow(u,2))*ipow4(v)))/(std::pow((u-4.0*m2),2)*v)+m2*(2.0*std::sqrt(-u)*(-(std::pow(t,2)*(4.0*m2-u)*std::pow(u,2))+2.0*t*u*(-m2*(4.0*m2+s)+(m2+s)*u)*v+3.0*m2*(2.0*m2-u)*ipow3(v))*std::log(((2.0*m2-u)*(2.0*m2+std::sqrt(4.0*m2-u)*std::sqrt(-u) -u))/(2.0*std::pow(m2,2))-1.0))/std::pow((4.0*m2-u),2.5))/(ipow4(t)*ipow3(u));
        aj41 = -pi*(v*(s*v-t*u-2.0*m2*v)/(m2+v)-((std::pow(v,2)-4.0*u*m2)*(s-4.0*m2)-u*(s*v-2.0*m2*(3.0*v-2.0*u)))*std::log(std::pow((v*(v-u+std::sqrt(std::pow((v-u),2)-4.0*m2*u))-2.0*m2*u),2)/4.0/m2/std::pow(u,2)/(m2+v))/std::sqrt(std::pow((v-u),2)-4.0*m2*u))/(std::pow((v-u),2)-4.0*m2*u);
        aj42 = 2.0*pi*std::log(std::pow((v*std::sqrt(s-4.0*m2)+std::sqrt(4.0*u*t*m2+std::pow(v,2)*(s-4.0*m2))),2)/(4.0*m2*t*u))/std::sqrt((s-4.0*m2)*(4.0*u*t*m2+std::pow(v,2)*(s-4.0*m2)));
        aj43 = 2.0*pi*std::log((std::sqrt(t*(t-4.0*m2))+2.0*m2-t)/(2.0*m2))/u/std::sqrt(t*(t-4.0*m2));
        aj44 = pi*(std::sqrt(t*(t-4.0*m2))*(2.0*m2*std::pow(v,2)-t*u*(s+u))/v+2.0*m2*(t*(v-u)-2.0*m2*v)*std::log((std::sqrt(t*(t-4.0*m2))+2.0*m2-t)/(2.0*m2)))/std::pow((t*(t-4.0*m2)),1.5)/std::pow(u,2);
        aj45 = pi*std::log(std::pow((v*(std::sqrt(als)+s)-2.0*m2*(u+v)),2)/4.0/m2/(m2*std::pow((u+v),2)-s*u*v))/std::sqrt(als)/u;
        aj46 = m2*(-pi*v*(s*u*(t+v)-2.0*m2*v*(u+v))/u/(m2*std::pow((u+v),2)-s*u*v)+(s*(u-v)+2.0*m2*v)*aj45)/u/als;
        // FORTRAN fsir.f:310 has v*(v-v) (identically zero). Intended v*(v-u); see README_MERADGEN_TYPO.md.
        aj47 = m2*(2.0*aj45*(std::pow(s,2)*std::pow((u-v),2)+6.0*std::pow(m2,2)*std::pow(v,2)-2.0*m2*s*(t*u+2.0*v*(v-u)))-pi*(v*(-(ipow3(s)*std::pow(u,2)*v*(std::pow(s,2)-3.0*std::pow(u,2)+10.0*u*v-6.0*std::pow(v,2)+2.0*s*(v-u)))+4.0*ipow3(m2)*(u+v)*(8.0*std::pow(s,2)*std::pow(u,2) +std::pow((u+v),2)*(4.0*s*u-3.0*std::pow(v,2)))+2.0*m2*std::pow(s,2)*u*(2.0*ipow3(u)*v-ipow4(u)+13.0*std::pow(u,2)*std::pow(v,2)-7.0*u*ipow3(v)-6.0*ipow4(v)+std::pow(s,2)*u*(u+5.0*v) +s*v*(8.0*u*v-5.0*std::pow(u,2)+std::pow(v,2)))-2.0*std::pow(m2,2)*s*(8.0*s*std::pow(u,2)*v*(s+2.0*v)+(u+v)*(8.0*std::pow(s,2)*std::pow(u,2)+2.0*ipow3(u)*(s+u)+6.0*ipow3(u)*v + 6.0*s*u*std::pow(v,2) +2.0*std::pow(u,2)*std::pow(v,2)-15.0*u*ipow3(v)-3.0*ipow4(v)))))/(u*std::pow((s*u*v-m2*std::pow((u+v),2)),2)))/(2.0*std::pow(als,2)*std::pow(u,2));
        aj48 = pi*(t*v*(4.0*m2-t)*(t*(u-v)+2.0*m2*v)-2.0*m2*std::sqrt(t*(t-4.0*m2))*(t*u*(s+u)-2.0*m2*std::pow(v,2))*std::log((2.0*m2-t+std::sqrt(t*(t-4.0*m2)))/2.0/m2))/std::pow((t*(t-4.0*m2)),2)/ipow3(u);
        aj49 = pi*((t*(t-4.0*m2)*(-2.0*m2*(t-6.0*m2)*ipow4(v)+t*((-64.0*ipow3(m2)+t*(std::pow((s+t),2)-(s+3.0*t)*u)+2.0*m2*(2.0*std::pow(s,2)+5.0*t*(2.0*u-v)))*std::pow(v,2)+(t*u-8.0*std::pow(m2,2)-2.0*m2*(3.0*s-2.0*t+u))*ipow3(v))-std::pow(t,2)*(t-4.0*m2)*(2.0*(2.0*m2-u)*std::pow(v,2)+std::pow(u,2)*(s+u+v))))/v-4.0*m2*std::sqrt(t*(t-4.0*m2))*(-3.0*m2*(t-2.0*m2)*ipow3(v)+t*u*(2.0*s*(t-m2)*v+(t-4.0*m2)*(t*u+2.0*m2*v)))*std::log((2.0*m2-t+std::sqrt(t*(t-4.0*m2)))/(2.0*m2)))/(ipow3(t)*ipow3((t-4.0*m2))*ipow4(u));
        aj50 = -pi/ipow3(u)/als*(v*(s*(v-u)-2.0*m2*v)+m2*(s*u*(t+v)-2.0*m2*v*(v+u))/std::sqrt(als)*std::log(std::pow((v*std::sqrt(als)+(s*v-2.0*m2*(u+v))),2)/(4.0*m2*(m2*std::pow((u+v),2)-s*u*v))));
        aj51 = pi*((v*(-(ipow3(s)*u*std::pow((u-v),2)*v)-4.0*std::pow(m2,2)*s*((2.0*s-u)*std::pow(u,2)*(s+u)-std::pow(u,2)*(11.0*s+3.0*u)*v+u*(u-2.0*s)*std::pow(v,2)+6.0*u*ipow3(v)+ipow4(v))+m2*std::pow(s,2)*(std::pow(u,2)*(std::pow(u,2)+std::pow((s+u),2)-8.0*(s+u)*v)+std::pow(v,2)*(2.0*std::pow(u,2)+4.0*u*v+std::pow(v,2)))+4.0*ipow3(m2)*(std::pow((u+v),2)*(3.0*std::pow(v,2)-4.0*s*u)+4.0*s*u*(s*u-v*(u+v)))))/(std::pow(als,2)*(m2*std::pow((u+v),2)-s*u*v))-(2.0*(-3.0*(std::pow(m2,2)+std::pow((s-m2),2))*u*std::pow(v,2)+3.0*m2*(s-2.0*m2)*ipow3(v)-s*(2.0*m2+s)*std::pow(u,2)*(v-t)+2.0*s*u*v*(s*u+(s-2.0*m2)*(v-t)))*std::log(std::pow((std::sqrt(als)*v+s*v-2.0*m2*(u+v)),2)/(4.0*m2*(m2*std::pow((u+v),2)-s*u*v)))*m2)/std::pow(als,2.5))/ipow4(u);
    // The following sr1..sr10 expressions are direct translations of the Fortran code,
    sr1 = 4.0*(2.0*(2.0*(2.0*t-3.0*v+5.0*s-2.0*m2)*m2*pls-((3.0*(t-v)+2.0*s)*pls*s+1.0))*aj7-(6.0*(2.0*m2-s)*m2*pls+pls*std::pow(s,2)+1.0)*aj6+2.0*(pls*std::pow(s,2)-1.0)*aj5+4.0*((4.0*(t-v+s)+pls*s*std::pow(t,2))*m2-(2.0*(pls*std::pow(t,2)+2.0)*std::pow(m2,2)+std::pow((s+t-v),2)))*aj40-(pls*std::pow(s,2)-1.0)*aj4-2.0*(2.0*(((t-v)*(t-2.0*v)+6.0*std::pow(s,2)+2.0*(3.0*t-2.0*v)*s)*pls+4.0-4.0*(2.0*s+t)*m2*pls)*m2-((2.0*(std::pow((t-v),2)+std::pow(s,2))+(3.0*t-2.0*v)*s)*pls*s+2.0*(t-v+s)))*aj39-2.0*(2.0*m2-s)*aj37*pls-2.0*(2.0*(4.0*(tt*v-1.0-2.0*s*tt)-(2.0*t-v)*pls*t+8.0*m2*tt)*m2+4.0*std::pow(s,2)*tt-t-4.0*(tt*v-1.0)*s+pls*s*std::pow(t,2))*aj35-(2.0*(((2.0*tt*v-9.0-4.0*s*tt)*s+(2.0*tt*v-3.0)*v)*pls*s-(4.0*(tt*v-3.0-s*tt)*s*tt+6.0*tt*v-5.0)-16.0*(2.0*pls*s-tt)*std::pow(m2,2)*tt-2.0*(((2.0*tt*v-5.0)*v+t+4.0*(tt*v-2.0-3.0*s*tt)*s)*pls-4.0*(tt*v-1.0-2.0*s*tt)*tt)*m2)*m2+2.0*(2.0*tt*v-3.0-2.0*s*tt)*s-(2.0*(tt*v-2.0)*v+3.0*t)+(2.0*s+t)*pls*std::pow(s,2))*aj34-4.0*std::pow((2.0*m2-s),2)*aj29+4.0*(2.0*(4.0*m2-3.0*s)*m2*pls+pls*std::pow(s,2)-1.0)*aj28*s+2.0*((4.0*(((6.0*std::pow(s,2)-std::pow(v,2))*tt-2.0*(tt*v-2.0)*s)*pls+2.0*(tt*v-1.0-2.0*s*tt)*tt-4.0*(2.0*pls*s-tt)*m2*tt)*m2-(2.0*(2.0*(tt*v-3.0-s*tt)*s*tt+3.0*tt*v-1.0)-(2.0*(tt*v-7.0-2.0*s*tt)*s+(2.0*tt*v+3.0)*v)*pls*s))*m2-((tt*v-1.0)*v+t-2.0*(tt*v-1.0-s*tt)*s-(t-v+2.0*s)*pls*std::pow(s,2)))*aj23+2.0*(8.0*(tt*v-1.0-2.0*s*tt+2.0*m2*tt)*m2+4.0*std::pow(s,2)*tt-t-4.0*(tt*v-1.0)*s)*aj22+(2.0*(4.0*(((6.0*std::pow(s,2)-std::pow(v,2))*tt-2.0*(tt*v-2.0)*s)*pls+2.0*(tt*v-1.0-2.0*s*tt)*tt-4.0*(2.0*pls*s-tt)*m2*tt)*m2+((2.0*tt*v-9.0-4.0*s*tt)*s+(2.0*tt*v+1.0)*v)*pls*s-(4.0*(tt*v-3.0-s*tt)*s*tt+6.0*tt*v-1.0))*m2+2.0*(2.0*tt*v-1.0-2.0*s*tt)*s-(2.0*(tt*v-1.0)*v+t)+(t-2.0*v+2.0*s)*pls*std::pow(s,2))*aj21+4.0*(pls*s*t-1.0-2.0*m2*pls*t)*aj12-2.0*(4.0*std::pow(m2,2)*pls-2.0*m2*pls*s+1.0)*aj11);
    sr2 = -4.0*((2.0*(tt*v-1.0+t*vt-2.0*(tt-vt-6.0*s*tt*vt)*s-2.0*(t*vt+tt*v-4.0*(tt-vt)*s)*pls*std::pow(s,2)-4.0*(5.0*(tt-vt)*pls*s+6.0*tt*vt)*m2*s)*m2-(tt*v-1.0+t*vt-2.0*(tt-vt-2.0*s*tt*vt)*s-(tt*v-1.0+t*vt-2.0*(tt-vt)*s)*pls*std::pow(s,2))*s+16.0*((tt*v+1.0+t*vt+2.0*(tt-vt)*s)*pls+2.0*tt*vt)*ipow3(m2))*aj36+(2.0*(((2.0*tt-11.0*vv-8.0*s*tt*vv)*s-(3.0*t*vv-4.0*tt*v-2.0))*pls*s-(tt+6.0*vv+6.0*s*tt*vv))*m2-(((tt-4.0*vv-2.0*s*tt*vv)*s-2.0*(t*vv-1.0))*(pls*std::pow(s,2)+1.0)+32.0*(s*tt+1.0)*ipow3(m2)*pls*vv)-8.0*((2.0*tt*v+1.0-2.0*t*vv-(5.0*s*tt+4.0)*s*vv)*pls-2.0*tt*vv)*std::pow(m2,2))*aj6+(8.0*(2.0*(tt+vv)-5.0*s*tt*vv+4.0*m2*tt*vv)*std::pow(m2,2)*pls*s-((2.0*s*vv-3.0)*s*tt+tt*v-1.0)*(pls*std::pow(s,2)-1.0)+2.0*((8.0*std::pow(s,2)*vv-7.0*s+2.0*v)*pls*s-(2.0*s*vv+1.0))*m2*tt)*aj4+(2.0*(4.0*(((2.0*(t*vv+1.0)-5.0*(vt-vv)*s)*s-(std::pow(t,2)*vv-2.0*t+v))*pls-2.0*(tt+2.0*vv-3.0*(vt-vv)*s*tt)-2.0*(((vt-vv)*t+1.0-2.0*(vt-vv)*s)*pls+2.0*(vt-vv)*tt)*m2)*m2-(t*vt-7.0*t*vv+4.0*tt*v+2.0+12.0*(vt-vv)*std::pow(s,2)*tt+2.0*(vt-9.0*vv-2.0*tt)*s+(3.0*t-2.0*v-8.0*(vt-vv)*std::pow(s,2)-((2.0*vt-7.0*vv)*t-2.0)*s)*pls*s))*m2+2.0*(2.0*(vt-vv)*s*tt+vt-5.0*vv)*std::pow(s,2)-(3.0*std::pow(t,2)*vv-4.0*t+2.0*v)+((vt-9.0*vv)*t+6.0)*s-(s*vt-s*vv-t*vv)*(2.0*s+t)*pls*std::pow(s,2))*aj34-(((pls*std::pow(s,2)+1.0)*(s*tt+1.0)-16.0*ipow3(m2)*pls*tt)*vv+8.0*(tt+vv+2.0*s*tt*vv)*std::pow(m2,2)*pls-2.0*((2.0*tt+3.0*vv+4.0*s*tt*vv)*pls*s+tt*vv)*m2)*aj32-(2.0*(4.0*(2.0*((tt*v+1.0+t*vt+2.0*(tt-vt)*s)*pls+2.0*tt*vt)*m2+(2.0*(tt*v-3.0-t*vt)*s-(5.0*(tt-vt)*std::pow(s,2)-tt*std::pow(v,2)))*pls+2.0*(2.0*tt+vt-3.0*s*tt*vt))*m2+7.0*tt*v-1.0+t*vt-6.0*(3.0*tt+vt-2.0*s*tt*vt)*s-((7.0*tt*v-18.0-4.0*t*vt)*s-(8.0*(tt-vt)*std::pow(s,2)-3.0*v))*pls*s)*m2+(3.0*tt*v-4.0)*v+2.0*t+2.0*(5.0*tt+vt-2.0*s*tt*vt)*std::pow(s,2)-(9.0*tt*v-5.0+t*vt)*s-((tt*v-4.0)*v+2.0*t+2.0*(tt-vt)*std::pow(s,2)-(3.0*tt*v-7.0-t*vt)*s)*pls*std::pow(s,2))*aj24+((2.0*(vt+vv-2.0*(vt-vv)*s*tt)*s-((vt-vv)*t+2.0)-((vt+vv)*t-2.0*(vt-vv)*s)*pls*std::pow(s,2))*s+16.0*((t*vt+1.0-2.0*(vt-vv)*s)*pls+2.0*(vt-vv)*tt)*ipow3(m2)-8.0*(((2.0*vt+vv)*t+2.0-5.0*(vt-vv)*s)*pls*s-2.0*(tt+vt-3.0*(vt-vv)*s*tt))*std::pow(m2,2)+2.0*((vt-vv)*t+2.0+12.0*(vt-vv)*std::pow(s,2)*tt-2.0*(3.0*vt+vv+2.0*tt)*s+((4.0*vt+vv)*t+2.0-8.0*(vt-vv)*s)*pls*std::pow(s,2))*m2)*aj21-((s*vv-1.0-2.0*m2*vv)*(pls*std::pow(s,2)-1.0)+8.0*std::pow(m2,2)*pls*s*vv)*aj20*tt-2.0*(4.0*(tt-2.0*vt)*std::pow(m2,2)*pls+tt-2.0*((tt-2.0*vt)*pls*s+2.0*tt*vt)*m2)*aj16-2.0*(4.0*(2.0*vt+vv)*std::pow(m2,2)*pls+vv-2.0*((2.0*vt+vv)*pls*s-2.0*tt*vt)*m2)*aj15-2.0*(2.0*tt*v-1.0-(3.0*tt+vt)*s+pls*std::pow(s,2)+8.0*(tt+vt)*ipow3(m2)*pls+4.0*((2.0*tt*v-1.0-t*vt-2.0*s*tt)*pls-2.0*tt*vt)*std::pow(m2,2)-((3.0*tt*v-1.0-2.0*t*vt-2.0*s*tt)*pls*s-2.0*(2.0*s*vt+1.0)*tt)*m2)*aj13-2.0*(2.0*(2.0*(((vt-vv)*t+3.0)*pls+2.0*(vt-vv)*tt-2.0*(vt-vv)*m2*pls)*m2-(2.0*((vt-vv)*s*tt-2.0*vv)+(t*vt+2.0)*pls*s))*m2-(t*vv-1.0-(vt-3.0*vv)*s-pls*std::pow(s,2)))*aj11-2.0*(4.0*(tt+vv+2.0*s*tt*vv-2.0*m2*tt*vv)*std::pow(m2,2)*pls+(s*vv-1.0)*tt+(tt+vv)*pls*std::pow(s,2)-((tt+3.0*vv+2.0*s*tt*vv)*pls*s+2.0*tt*vv)*m2)*aj1);
    sr3 = -4.0*(2.0*((((2.0*m2*tt-s*vv)*m2*pls+(pls*std::pow(s,2)-1.0)*vv)*aj2-2.0*(pls*s-tt-2.0*m2*pls)*aj19+(4.0*std::pow(m2,2)*pls-2.0*m2*pls*s+1.0)*aj18*vv)*tt-((4.0*tt*v-1.0-2.0*s*tt)*pls*s-2.0*(2.0*tt*v-1.0-2.0*s*tt)*tt+2.0*((4.0*s-3.0*v)*pls-4.0*tt)*m2*tt)*aj17+(2.0*tt-vv-2.0*s*tt*vv+pls*std::pow(s,2)*vv-(2.0*(2.0*tt-3.0*vv-2.0*s*tt*vv)*tt-(tt-3.0*vv-2.0*s*tt*vv)*pls*s)*m2-2.0*((tt*v-3.0-4.0*s*vv)*pls+4.0*tt*vv)*std::pow(m2,2)*tt)*aj16+(2.0*((6.0*std::pow(s,2)+std::pow(v,2)-(tt*v+4.0)*s*v)*pls-4.0*(s-v)*tt-2.0*((4.0*s-tt*std::pow(v,2))*pls-2.0*tt)*m2)*m2*tt-(((2.0*tt*v-1.0)*v+2.0*std::pow(s,2)*tt-(3.0*tt*v-1.0)*s)*pls*s+2.0*(2.0*tt*v-1.0-s*tt)*s*tt-(2.0*std::pow(tt,2)*std::pow(v,2)-2.0*tt*v+1.0)))*aj14)-(2.0*((2.0*(6.0*tt-7.0*vv-2.0*s*tt*vv)*s-(tt*v+3.0)*(tt*v-1.0))*pls*s-(4.0*(2.0*tt-3.0*vv-s*tt*vv)*s*tt-(3.0*std::pow(tt,2)*v-6.0*tt+2.0*vv))-8.0*((4.0*s*vv-tt*v)*pls-2.0*tt*vv)*std::pow(m2,2)*tt+2.0*((((2.0*tt*v-3.0)*v+12.0*std::pow(s,2)*vv)*tt-2.0*(std::pow(tt,2)*v+6.0*tt-4.0*vv)*s)*pls+4.0*(2.0*tt-vv-2.0*s*tt*vv)*tt)*m2)*m2-(3.0*tt*v-4.0+2.0*t*vv-2.0*(3.0*tt-2.0*vv-2.0*s*tt*vv)*s-(tt*v-4.0+2.0*t*vv-2.0*(tt-2.0*vv)*s)*pls*std::pow(s,2)))*aj13-((pls*std::pow(s,2)-1.0)*vv+2.0*(pls*s-vv)*m2*tt)*aj10*tt+(4.0*((2.0*(3.0*std::pow(s,2)*std::pow(vv,2)+tt*v)-(tt+8.0*vv)*s)*pls-4.0*(s*vv-1.0)*tt*vv)*std::pow(m2,2)*tt-(8.0*((4.0*s*std::pow(vv,2)-tt)*pls-2.0*tt*std::pow(vv,2))*ipow3(m2)*tt-(tt-vv-2.0*s*tt*vv)*(pls*std::pow(s,2)-1.0))+2.0*(((9.0*tt-vv-2.0*s*tt*vv)*s*vv-(tt*v+3.0)*tt)*pls*s-(2.0*(2.0*tt-vv-s*tt*vv)*s*tt*vv-(3.0*std::pow(tt,2)-tt*vv+std::pow(vv,2))))*m2)*aj1);
    sr4 = -4.0*(((2.0*(2.0*t*vv-3.0+s*vv)*s+2.0*std::pow(t,2)*vv-5.0*t+4.0*v)*(pls*std::pow(s,2)+1.0)-16.0*(2.0*t*vv-3.0+6.0*s*vv)*ipow3(m2)*pls-2.0*(((11.0*s*vv+18.0*t*vv-21.0)*s+9.0*std::pow(t,2)*vv-20.0*t+12.0*v)*pls*s+9.0*t*vv-10.0+7.0*s*vv)*m2+4.0*(((22.0*t*vv-25.0+20.0*s*vv)*s+4.0*std::pow(t,2)*vv-18.0*t+13.0*v)*pls+8.0*vv)*std::pow(m2,2))*aj6-2.0*(((11.0*t-8.0*v+5.0*s)*s+2.0*(4.0*std::pow(t,2)-6.0*t*v+3.0*std::pow(v,2)))*pls*s+4.0*(t-v+s)+16.0*(2.0*t-v+4.0*s)*std::pow(m2,2)*pls-2.0*(((26.0*t-17.0*v+18.0*s)*s+2.0*(4.0*std::pow(t,2)-6.0*t*v+3.0*std::pow(v,2)))*pls+8.0)*m2)*aj7-2.0*((2.0*(t-v)-3.0*std::pow(s,2)*vv+2.0*(t*vv+1.0)*s)*pls*s-(2.0*(t*vv-1.0)-s*vv)+4.0*(5.0*(s*vv-1.0)*pls*s+vv-8.0*m2*pls*s*vv)*m2)*aj4*m2-2.0*(2.0*((2.0*(7.0*t-6.0*v+4.0*s)*std::pow(s,2)+(3.0*std::pow(t,2)-3.0*t*v+2.0*std::pow(v,2))*(t-v)+(13.0*std::pow(t,2)-17.0*t*v+8.0*std::pow(v,2))*s)*pls+8.0*(t-v+s))*m2-(8.0*(((4.0*(t-v)+5.0*s)*s+2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2))*pls+3.0)*std::pow(m2,2)+((2.0*std::pow(s,2)+3.0*s*t-2.0*s*v+4.0*std::pow(t,2)-4.0*t*v+2.0*std::pow(v,2))*pls*s+2.0*(s+t-v))*(s+t-v)-32.0*ipow3(m2)*pls*s))*aj39+4.0*(2.0*m2-s)*aj38*pls-4.0*((7.0*t-6.0*v+9.0*s-8.0*m2)*m2*pls-((3.0*(t-v)+2.0*s)*pls*s+1.0))*aj37-(((2.0*(2.0*(t*vv-1.0)+s*vv)*s+3.0*(std::pow(t,2)*vv-2.0*t+v))*s+ipow3(t)*vv-3.0*std::pow(t,2)+3.0*t*v-std::pow(v,2))*(pls*std::pow(s,2)+1.0)+64.0*(tt+vv)*ipow4(m2)*pls*s+16.0*(((2.0*(tt*v+1.0-2.0*t*vv)-(5.0*tt+7.0*vv)*s)*s-((t*vv-3.0)*t+(tt*v+1.0)*v))*pls-3.0*(tt+vv))*ipow3(m2)-4.0*((6.0*(tt*v+2.0-3.0*t*vv)*std::pow(s,2)-(8.0*ipow3(s)*tt+18.0*ipow3(s)*vv+2.0*ipow3(t)*vv-10.0*std::pow(t,2)+13.0*t*v-5.0*std::pow(v,2))-(3.0*(2.0*t*vv-7.0)*t+2.0*(2.0*tt*v+3.0)*v)*s)*pls+2.0*(tt*v+3.0-7.0*t*vv-(4.0*tt+7.0*vv)*s))*std::pow(m2,2)-2.0*((2.0*ipow3(t)*vv-10.0*std::pow(t,2)+11.0*t*v-3.0*std::pow(v,2)+2.0*(tt+5.0*vv)*ipow3(s)-(2.0*tt*v+13.0-15.0*t*vv)*std::pow(s,2)+((7.0*t*vv-18.0)*t+2.0*(tt*v+3.0)*v)*s)*pls*s+7.0*std::pow(t,2)*vv-10.0*t+3.0*v+2.0*(tt+5.0*vv)*std::pow(s,2)-(2.0*(tt*v+5.0)-15.0*t*vv)*s)*m2)*aj34-(8.0*(t*vv-4.0+2.0*s*vv-2.0*m2*vv)*std::pow(m2,2)*pls+(t*vv-2.0+s*vv)*(pls*std::pow(s,2)+1.0)-2.0*(3.0*(t*vv-2.0+s*vv)*pls*s+2.0*vv)*m2)*aj32+(8.0*((4.0*ipow3(s)*tt+std::pow(v,2)-3.0*(tt*v-1.0)*std::pow(s,2)+((2.0*tt*v-1.0)*v-6.0*t)*s)*pls-(tt*v-3.0-4.0*s*tt))*std::pow(m2,2)-(16.0*(((5.0*std::pow(s,2)+std::pow(v,2))*tt-2.0*(tt*v+1.0)*s)*pls+3.0*tt)*ipow3(m2)-((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)-2.0*s*v+std::pow(v,2))+64.0*ipow4(m2)*pls*s*tt))+2.0*((2.0*(tt*v-4.0-s*tt)*std::pow(s,2)+4.0*std::pow(t,2)-3.0*t*v-3.0*std::pow(v,2)-2.0*((tt*v-3.0)*v-2.0*t)*s)*pls*s-(2.0*std::pow(s,2)*tt-5.0*v-2.0*(tt*v-4.0)*s))*m2)*aj23-(16.0*((((5.0*tt+7.0*vv)*s-2.0*tt*v)*s+(tt*v+1.0)*v)*pls+3.0*(tt+vv))*ipow3(m2)-((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)*vv-2.0*s+v)+64.0*(tt+vv)*ipow4(m2)*pls)*s-2.0*((3.0*std::pow(t,2)*vv-5.0*t-2.0*tt*std::pow(v,2)-2.0*(tt+5.0*vv)*std::pow(s,2)+(2.0*tt*v+5.0+t*vv)*s)*pls*std::pow(s,2)-(std::pow(t,2)*vv+v+2.0*(tt+5.0*vv)*std::pow(s,2)-(2.0*(tt*v+1.0)+t*vv)*s))*m2+8.0*(((3.0*tt*v+1.0+t*vv-(4.0*tt+9.0*vv)*s)*s+(t*vv-1.0)*t-(2.0*tt*v+1.0)*v)*pls*s-(t*vv-tt*v+(4.0*tt+7.0*vv)*s))*std::pow(m2,2))*aj21+2.0*((3.0*s*vv-2.0)*pls*s-vv-4.0*m2*pls*s*vv)*aj20*m2-(4.0*(2.0*m2-s)*(t*vv-2.0)*m2*pls+(pls*std::pow(s,2)+1.0)*(t*vv-1.0))*aj15+(4.0*((2.0*(4.0*(t*vv-1.0)+3.0*s*vv)*s+2.0*std::pow(t,2)*vv-7.0*t+7.0*v)*pls+6.0*vv)*std::pow(m2,2)+(std::pow(t,2)*vv-3.0*t+2.0*v+(3.0*t*vv-2.0)*s)*(pls*std::pow(s,2)+1.0)-16.0*(2.0*s+t)*ipow3(m2)*pls*vv-2.0*(((9.0*t*vv-8.0+2.0*s*vv)*s+4.0*std::pow(t,2)*vv-11.0*t+9.0*v)*pls*s+7.0*t*vv-8.0+2.0*s*vv)*m2)*aj11-(8.0*(t*vv-3.0+3.0*s*vv-2.0*m2*vv)*std::pow(m2,2)*pls+(t*vv-2.0+s*vv)*(pls*std::pow(s,2)+1.0)-2.0*((5.0*t*vv-7.0+4.0*s*vv)*pls*s+4.0*vv)*m2)*aj1)*uu;
    sr5 = -4.0*(((16.0*(((2.0*(t*uu+3.0)-7.0*s*uu)*s-2.0*(std::pow(t,2)*uu+t+v))*pls-3.0*uu)*ipow3(m2)-((2.0*(t*uu+1.0-s*uu)*s-(std::pow(t,2)*uu+t+v))*(pls*std::pow(s,2)+1.0)-64.0*ipow4(m2)*pls*uu)*s-8.0*(((5.0*t*uu+11.0-9.0*s*uu)*s-(3.0*std::pow(t,2)*uu+5.0*t+4.0*v))*pls*s-7.0*(s*uu-1.0))*std::pow(m2,2)+2.0*((5.0*t*uu+12.0-10.0*s*uu)*s-(2.0*std::pow(t,2)*uu+2.0*t+v)-(2.0*(5.0*s*uu-4.0*t*uu-6.0)*s+4.0*std::pow(t,2)*uu+6.0*t+5.0*v)*pls*std::pow(s,2))*m2-4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s+v)*dz1*m2)*ds+2.0*(2.0*((t*uu+3.0)*t+3.0*std::pow(s,2)*uu-(3.0*t*uu+2.0)*s)*pls*s-3.0*(t*uu+1.0-2.0*s*uu))*m2+(2.0*(t*uu+1.0-s*uu)*s-(std::pow(t,2)*uu+t+v))*(pls*std::pow(s,2)+1.0)-8.0*(2.0*(s-t)*s*uu+std::pow(t,2)*uu+t+v)*std::pow(m2,2)*pls+((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)-2.0*s*v+std::pow(v,2))-32.0*ipow3(m2)*pls*s-2.0*(2.0*(4.0*std::pow(s,2)-3.0*s*v+std::pow(v,2))*pls*s+8.0*s-3.0*v)*m2+8.0*((5.0*std::pow(s,2)-2.0*s*v+std::pow(v,2))*pls+3.0)*std::pow(m2,2))*dz1+4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s)*std::pow(ds,2)*m2)*aj45+(((pls*std::pow(s,2)+1.0)*(s+t)-48.0*ipow3(m2)*pls+4.0*(7.0*t-4.0*v+21.0*s)*std::pow(m2,2)*pls-2.0*((11.0*t-4.0*v+11.0*s)*pls*s+4.0)*m2)*dz1-2.0*(2.0*(8.0*m2-5.0*s)*m2*pls+pls*std::pow(s,2)+1.0))*aj8-((pls*std::pow(s,2)+1.0)*(s+t)-48.0*ipow3(m2)*pls+4.0*(7.0*t-4.0*v+13.0*s)*std::pow(m2,2)*pls-2.0*((7.0*t-4.0*v+7.0*s)*pls*s+4.0)*m2)*aj6*dz1-16.0*aj5*m2*pls*s+4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s)*aj46*ds-(2.0*(4.0*((3.0*(3.0*s*vt+2.0)*s+std::pow(t,2)*vt+3.0*t+v)*pls*s+t*vt+1.0+7.0*s*vt)*std::pow(m2,2)+((pls*std::pow(s,2)+1.0)*std::pow(s,2)+32.0*ipow4(m2)*pls)*s*vt-8.0*(((7.0*s*vt+4.0)*s+std::pow(t,2)*vt+t+v)*pls+3.0*vt)*ipow3(m2)+((t*vt+1.0-10.0*s*vt)*s-(std::pow(t,2)*vt+t+v)-2.0*(5.0*std::pow(s,2)*vt+2.0*s+2.0*t)*pls*std::pow(s,2))*m2)*dz1+16.0*(6.0*t*vt-1.0+6.0*s*vt-8.0*m2*vt)*ipow3(m2)*pls-((2.0*std::pow(s,2)+std::pow(t,2))*uu-(2.0*uu-vt)*s*t)*(pls*std::pow(s,2)+1.0)+4.0*((t-2.0*v-2.0*(uu+vt)*std::pow(t,2)+((4.0*uu-15.0*vt)*t-4.0*(uu+2.0*vt)*s)*s)*pls-2.0*vt)*std::pow(m2,2)-2.0*((2.0*(3.0*(uu-vt)*t-1.0-(3.0*uu+vt)*s)*s+2.0*t-v-2.0*(uu+vt)*std::pow(t,2))*pls*s+(3.0*uu+vt)*t-1.0-6.0*s*uu)*m2+(((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)-2.0*s*t+std::pow(t,2))+64.0*ipow4(m2)*pls)*s-16.0*((7.0*std::pow(s,2)-2.0*s*t+2.0*std::pow(t,2))*pls+3.0)*ipow3(m2)+8.0*((9.0*std::pow(s,2)-5.0*s*t+3.0*std::pow(t,2))*pls+7.0)*std::pow(m2,2)*s-2.0*(10.0*std::pow(s,2)-5.0*s*t+2.0*std::pow(t,2)+2.0*(5.0*std::pow(s,2)-4.0*s*t+2.0*std::pow(t,2))*pls*std::pow(s,2))*m2)*(uu-vt)*ds)*aj43-((((t-v+3.0*s)*(t-v)-2.0*(s*vt-2.0)*std::pow(s,2))*(pls*std::pow(s,2)+1.0)-64.0*ipow4(m2)*pls*s*vt+16.0*(((7.0*s*vt-8.0)*s+std::pow(t,2)*vt-2.0*t+2.0*v)*pls+3.0*vt)*ipow3(m2)+2.0*(((10.0*(s*vt-2.0)*s-(14.0*t-13.0*v))*s-(6.0*std::pow(t,2)-9.0*t*v+4.0*std::pow(v,2)))*pls*s-((t*vt+20.0-10.0*s*vt)*s-(std::pow(t,2)*vt-8.0*t+9.0*v)))*m2-4.0*((2.0*(9.0*s*vt-16.0)*std::pow(s,2)-(5.0*std::pow(t,2)-7.0*t*v+4.0*std::pow(v,2))+(2.0*std::pow(t,2)*vt-17.0*t+16.0*v)*s)*pls+2.0*(t*vt-10.0+7.0*s*vt))*std::pow(m2,2))*dz1-(2.0*(((6.0*t*vt-1.0+2.0*s*vt)*s+2.0*std::pow(t,2)*vt-6.0*t+9.0*v)*pls*s-t*vt)*m2-(((t*vt-1.0)*s-(t-2.0*v))*(pls*std::pow(s,2)+1.0)+128.0*ipow4(m2)*pls*vt-16.0*(6.0*t*vt+5.0+6.0*s*vt)*ipow3(m2)*pls)-4.0*(((15.0*t*vt+7.0+8.0*s*vt)*s+2.0*std::pow(t,2)*vt-3.0*t+7.0*v)*pls+2.0*vt)*std::pow(m2,2))+4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s)*std::pow(ds,2)*m2+(16.0*(((2.0*(t*vt+3.0)-7.0*s*vt)*s-2.0*(std::pow(t,2)*vt+t+v))*pls-3.0*vt)*ipow3(m2)-((2.0*(t*vt+1.0-s*vt)*s-(std::pow(t,2)*vt+t+v))*(pls*std::pow(s,2)+1.0)-64.0*ipow4(m2)*pls*vt)*s-8.0*(((5.0*t*vt+11.0-9.0*s*vt)*s-(3.0*std::pow(t,2)*vt+5.0*t+4.0*v))*pls*s-7.0*(s*vt-1.0))*std::pow(m2,2)+2.0*((5.0*t*vt+12.0-10.0*s*vt)*s-(2.0*std::pow(t,2)*vt+2.0*t+v)-(2.0*(5.0*s*vt-4.0*t*vt-6.0)*s+4.0*std::pow(t,2)*vt+6.0*t+5.0*v)*pls*std::pow(s,2))*m2-4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s+v)*dz1*m2)*ds)*aj42-4.0*(2.0*m2-s)*aj41*dz1*m2*pls+8.0*(s+t-4.0*m2)*aj4*dz1*m2*pls*s+2.0*(4.0*((3.0*(3.0*s*vt+2.0)*s+std::pow(t,2)*vt+3.0*t+v)*pls*s+t*vt+1.0+7.0*s*vt)*std::pow(m2,2)+((pls*std::pow(s,2)+1.0)*std::pow(s,2)+32.0*ipow4(m2)*pls)*s*vt-8.0*(((7.0*s*vt+4.0)*s+std::pow(t,2)*vt+t+v)*pls+3.0*vt)*ipow3(m2)+((t*vt+1.0-10.0*s*vt)*s-(std::pow(t,2)*vt+t+v)-2.0*(5.0*std::pow(s,2)*vt+2.0*s+2.0*t)*pls*std::pow(s,2))*m2)*aj36*dz1+(((t-v+3.0*s)*(t-v)-2.0*(s*vt-2.0)*std::pow(s,2))*(pls*std::pow(s,2)+1.0)+64.0*(tt-vt)*ipow4(m2)*pls*s+16.0*(((2.0*(tt*v-4.0)-(5.0*tt-7.0*vt)*s)*s+std::pow(t,2)*vt-2.0*t+2.0*v)*pls-3.0*(tt-vt))*ipow3(m2)+2.0*(((2.0*(tt*v-10.0-(tt-5.0*vt)*s)*s-(14.0*t-13.0*v))*s-(6.0*std::pow(t,2)-9.0*t*v+4.0*std::pow(v,2)))*pls*s+std::pow(t,2)*vt-8.0*t+9.0*v-2.0*(tt-5.0*vt)*std::pow(s,2)+(2.0*(tt*v-10.0)-t*vt)*s)*m2-4.0*((2.0*(3.0*tt*v-16.0-(4.0*tt-9.0*vt)*s)*std::pow(s,2)-(5.0*std::pow(t,2)-7.0*t*v+4.0*std::pow(v,2))+(2.0*std::pow(t,2)*vt-17.0*t+16.0*v)*s)*pls+2.0*(3.0*tt*v-10.0+t*vt-(4.0*tt-7.0*vt)*s))*std::pow(m2,2))*aj34*dz1+4.0*(2.0*m2-s)*aj32*dz1*m2*pls-4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s)*aj30*ds+4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s)*aj28*ds-4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s)*aj26*ds-((16.0*(((2.0*(t*uu+3.0)-7.0*s*uu)*s-2.0*(std::pow(t,2)*uu+t+v))*pls-3.0*uu)*ipow3(m2)-((2.0*(t*uu+1.0-s*uu)*s-(std::pow(t,2)*uu+t+v))*(pls*std::pow(s,2)+1.0)-64.0*ipow4(m2)*pls*uu)*s-8.0*(((5.0*t*uu+11.0-9.0*s*uu)*s-(3.0*std::pow(t,2)*uu+5.0*t+4.0*v))*pls*s-7.0*(s*uu-1.0))*std::pow(m2,2)+2.0*((5.0*t*uu+12.0-10.0*s*uu)*s-(2.0*std::pow(t,2)*uu+2.0*t+v)-(2.0*(5.0*s*uu-4.0*t*uu-6.0)*s+4.0*std::pow(t,2)*uu+6.0*t+5.0*v)*pls*std::pow(s,2))*m2-4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s+v)*dz1*m2)*ds+2.0*(2.0*((t*uu+3.0)*t+3.0*std::pow(s,2)*uu-(3.0*t*uu+2.0)*s)*pls*s-3.0*(t*uu+1.0-2.0*s*uu))*m2+(2.0*(t*uu+1.0-s*uu)*s-(std::pow(t,2)*uu+t+v))*(pls*std::pow(s,2)+1.0)-8.0*(2.0*(s-t)*s*uu+std::pow(t,2)*uu+t+v)*std::pow(m2,2)*pls+((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)-2.0*s*v+std::pow(v,2))-32.0*ipow3(m2)*pls*s-2.0*(2.0*(4.0*std::pow(s,2)-3.0*s*v+std::pow(v,2))*pls*s+8.0*s-3.0*v)*m2+8.0*((5.0*std::pow(s,2)-2.0*s*v+std::pow(v,2))*pls+3.0)*std::pow(m2,2))*dz1+4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s)*std::pow(ds,2)*m2)*aj25+((((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)-2.0*s*t+std::pow(t,2))+64.0*ipow4(m2)*pls)*s-16.0*((7.0*std::pow(s,2)-2.0*s*t+2.0*std::pow(t,2))*pls+3.0)*ipow3(m2)+8.0*((9.0*std::pow(s,2)-5.0*s*t+3.0*std::pow(t,2))*pls+7.0)*std::pow(m2,2)*s-2.0*(10.0*std::pow(s,2)-5.0*s*t+2.0*std::pow(t,2)+2.0*(5.0*std::pow(s,2)-4.0*s*t+2.0*std::pow(t,2))*pls*std::pow(s,2))*m2)*(uu-vt)*ds-(2.0*((2.0*(3.0*t*uu+1.0-3.0*s*uu)*s-(2.0*std::pow(t,2)*uu+6.0*t-v))*pls*s+3.0*(t*uu+1.0-2.0*s*uu))*m2-((2.0*(t*uu+1.0-s*uu)*s-(std::pow(t,2)*uu+t+v))*(pls*std::pow(s,2)+1.0)-8.0*(2.0*(s-t)*s*uu+std::pow(t,2)*uu+t+v)*std::pow(m2,2)*pls)))*aj24-(8.0*((4.0*ipow3(s)*tt+std::pow(v,2)-3.0*(tt*v-1.0)*std::pow(s,2)-(6.0*t+v)*s)*pls-(3.0*(tt*v-1.0)-4.0*s*tt))*std::pow(m2,2)+16.0*((2.0*(tt*v+1.0)-5.0*s*tt)*pls*s-3.0*tt)*ipow3(m2)+(pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)-2.0*s*v+std::pow(v,2))+64.0*ipow4(m2)*pls*s*tt+2.0*(((2.0*(tt*v-4.0-s*tt)*s+4.0*t+5.0*v)*s+4.0*std::pow(t,2)-t*v-2.0*std::pow(v,2))*pls*s-(2.0*std::pow(s,2)*tt-3.0*v-2.0*(tt*v-4.0)*s))*m2)*aj23*dz1+((16.0*(((2.0*(t*vt+3.0)-7.0*s*vt)*s-2.0*(std::pow(t,2)*vt+t+v))*pls-3.0*vt)*ipow3(m2)-((2.0*(t*vt+1.0-s*vt)*s-(std::pow(t,2)*vt+t+v))*(pls*std::pow(s,2)+1.0)-64.0*ipow4(m2)*pls*vt)*s-8.0*(((5.0*t*vt+11.0-9.0*s*vt)*s-(3.0*std::pow(t,2)*vt+5.0*t+4.0*v))*pls*s-7.0*(s*vt-1.0))*std::pow(m2,2)+2.0*((5.0*t*vt+12.0-10.0*s*vt)*s-(2.0*std::pow(t,2)*vt+2.0*t+v)-(2.0*(5.0*s*vt-4.0*t*vt-6.0)*s+4.0*std::pow(t,2)*vt+6.0*t+5.0*v)*pls*std::pow(s,2))*m2-4.0*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*(2.0*m2-s+v)*dz1*m2)*ds-4.0*((2.0*m2-s+v)*dz1*tt-(2.0*m2-s)*std::pow(ds,2))*(8.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)-6.0*m2+pls*ipow3(s)+s)*m2)*aj21+4.0*(2.0*m2-s)*aj16*m2*pls*vt-4.0*(2.0*m2-s)*aj15*m2*pls*vt-(32.0*ipow3(m2)*pls*vt-8.0*std::pow(m2,2)*pls*s*vt-16.0*std::pow(m2,2)*pls*t*vt+4.0*std::pow(m2,2)*pls+8.0*m2*pls*s*t*vt-2.0*m2*pls*s-4.0*m2*vt+pls*std::pow(s,2)+1.0)*aj13+(32.0*ipow3(m2)*pls*vt-8.0*std::pow(m2,2)*pls*s*vt-16.0*std::pow(m2,2)*pls*t*vt-20.0*std::pow(m2,2)*pls+8.0*m2*pls*s*t*vt+6.0*m2*pls*s-4.0*m2*vt-pls*std::pow(s,2)-1.0)*aj11);
    sr6 = 4.0*((((2.0*(s*tt+2.0)*s*vv+3.0*(t*vv-1.0))*s+std::pow(t,2)*vv-2.0*t+v)*(pls*std::pow(s,2)+1.0)+64.0*ipow4(m2)*pls*s*tt*vv-16.0*((tt*v-1.0+2.0*t*vv+(7.0*s*vv+2.0)*s*tt)*pls+3.0*tt*vv)*ipow3(m2)+4.0*(((2.0*(3.0*tt+4.0*vv+9.0*s*tt*vv)*s-(2.0*tt*v-1.0))*s+(2.0*t*vv-1.0)*t+(2.0*tt*v-1.0)*v)*pls+2.0*(2.0*tt+3.0*vv+7.0*s*tt*vv))*std::pow(m2,2)+2.0*(((2.0*tt*v+3.0-3.0*t*vv-2.0*(tt+6.0*vv+5.0*s*tt*vv)*s)*s+(t*vv-1.0)*t-(2.0*tt*v-1.0)*v)*pls*s+tt*v+3.0-5.0*t*vv-(3.0*tt+13.0*vv+10.0*s*tt*vv)*s)*m2)*aj6-2.0*(((pls*std::pow(s,2)+1.0)*std::pow(s,2)+32.0*ipow4(m2)*pls)*s*tt-8.0*((7.0*std::pow(s,2)*tt-4.0*s+tt*std::pow(v,2))*pls+3.0*tt)*ipow3(m2)-(2.0*((5.0*s*tt-2.0)*s-2.0*(t-v))*pls*std::pow(s,2)+(10.0*std::pow(s,2)+s*v+std::pow(v,2))*tt)*m2+4.0*(((tt*v+2.0)*v-2.0*t+3.0*(3.0*s*tt-2.0)*s)*pls*s+(7.0*s-v)*tt)*std::pow(m2,2))*aj36+(((pls*std::pow(s,2)+1.0)*(s*tt+1.0)-32.0*ipow3(m2)*pls*tt)*vv+8.0*(tt+vv+3.0*s*tt*vv)*std::pow(m2,2)*pls-4.0*((tt+vv+2.0*s*tt*vv)*pls*s+tt*vv)*m2)*aj33+(((tt-4.0*vv-2.0*s*tt*vv)*s-2.0*(t*vv-1.0))*(pls*std::pow(s,2)+1.0)-16.0*(tt-4.0*vv-2.0*s*tt*vv)*ipow3(m2)*pls+4.0*(4.0*tt*v-1.0-4.0*t*vv+(tt-18.0*vv-10.0*s*tt*vv)*s)*std::pow(m2,2)*pls-2.0*(((2.0*tt-15.0*vv-8.0*s*tt*vv)*s-(5.0*t*vv-4.0*tt*v-2.0))*pls*s-(2.0*tt+7.0*vv+5.0*s*tt*vv))*m2)*aj32+(8.0*(4.0*s*tt+1.0-4.0*m2*tt)*std::pow(m2,2)*pls+(pls*std::pow(s,2)+1.0)*(s*tt+1.0)-2.0*(5.0*pls*std::pow(s,2)+1.0)*m2*tt)*aj2*vv+(4.0*(2.0*m2-s)*(2.0*tt-vv)*m2*pls+(pls*std::pow(s,2)+1.0)*(tt-vv))*aj18-(4.0*((2.0*(8.0*tt-vv-3.0*s*tt*vv)*s-(7.0*tt*v-3.0))*pls-6.0*tt*vv)*std::pow(m2,2)-((2.0*(tt*v-1.0)-(3.0*tt-2.0*vv)*s)*(pls*std::pow(s,2)+1.0)-32.0*(s*vv-1.0)*ipow3(m2)*pls*tt)+2.0*((9.0*tt*v-4.0-t*vv-(13.0*tt-5.0*vv-2.0*s*tt*vv)*s)*pls*s-(3.0*(3.0*tt-vv)-2.0*s*tt*vv))*m2)*aj16+(((tt*v-1.0)*v+2.0*(2.0*tt-vv)*std::pow(s,2)-(3.0*tt*v-2.0)*s)*(pls*std::pow(s,2)+1.0)-64.0*ipow4(m2)*pls*s*tt*vv+16.0*(((5.0*std::pow(s,2)*vv+2.0*v)*tt-2.0*(4.0*tt+vv)*s)*pls+3.0*tt*vv)*ipow3(m2)+4.0*((2.0*(16.0*tt-3.0*vv-4.0*s*tt*vv)*std::pow(s,2)+(5.0*tt*v-3.0)*v-(17.0*tt*v+6.0-12.0*t*vv)*s)*pls+2.0*(10.0*tt-3.0*vv-4.0*s*tt*vv))*std::pow(m2,2)+2.0*(((14.0*tt*v-3.0-4.0*t*vv+2.0*(s*tt*vv-10.0*tt+4.0*vv)*s)*s-((4.0*t*vv-7.0)*t+(6.0*tt*v-1.0)*v))*pls*s+2.0*(s*tt*vv-10.0*tt+4.0*vv)*s+8.0*tt*v-5.0)*m2)*aj13+(((pls*std::pow(s,2)+1.0)*(s*tt+1.0)-32.0*ipow3(m2)*pls*tt)*vv+8.0*(tt+vv+5.0*s*tt*vv)*std::pow(m2,2)*pls-4.0*(3.0*(s*tt+1.0)*pls*s+tt)*m2*vv)*aj10-(16.0*((2.0*(tt-vv)+5.0*std::pow(s,2)*tt*std::pow(vv,2)-2.0*(3.0*tt+2.0*vv)*s*vv)*pls+3.0*tt*std::pow(vv,2))*ipow3(m2)-(((tt-3.0*vv-2.0*s*tt*vv)*s-(t*vv-2.0))*(pls*std::pow(s,2)+1.0)+64.0*ipow4(m2)*pls*s*tt*std::pow(vv,2))+2.0*(((2.0*(s*vv-6.0)*s*tt*vv+8.0*tt-5.0*vv)*s+t*vv-3.0*tt*v)*pls*s+2.0*(s*vv-6.0)*s*tt*vv+6.0*tt-5.0*vv)*m2+4.0*((2.0*(11.0*tt+2.0*vv-4.0*s*tt*vv)*std::pow(s,2)*vv+2.0*t*vv-1.0-(13.0*tt+4.0*vv-4.0*t*std::pow(vv,2))*s)*pls-2.0*(4.0*s*vv-7.0)*tt*vv)*std::pow(m2,2))*aj1)*uu;
    sr7 = -4.0*((2.0*((2.0*uu-3.0*vv+7.0*tt-12.0*s*tt*vv)*s-(2.0*(tt*v-1.0)+t*vv)-((2.0*uu-5.0*vv-5.0*tt+2.0*(uu+4.0*vv)*s*tt)*s+4.0*tt*v+3.0+(2.0*uu-3.0*vv)*t)*pls*std::pow(s,2))*m2+((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)*vv-2.0*s+v)*s-256.0*ipow5(m2)*pls*uu)*tt+64.0*(tt+2.0*uu+(5.0*uu-vv)*s*tt)*ipow4(m2)*pls-16.0*(((4.0*(2.0*uu-vv)+5.0*tt+(10.0*uu-vv)*s*tt)*s+2.0*tt*v+1.0+t*uu)*pls+(uu+3.0*vv)*tt)*ipow3(m2)+8.0*(((5.0*uu-7.0*vv+tt+(5.0*uu+4.0*vv)*s*tt)*s+3.0*(tt*v+1.0)+(2.0*uu-vv)*t)*pls*s-(uu-vv+2.0*tt-(uu+9.0*vv)*s*tt))*std::pow(m2,2))*aj45-(16.0*(6.0*tt+5.0*vv+13.0*s*tt*vv-12.0*m2*tt*vv)*ipow3(m2)*pls+((s*vv-2.0)*s*tt-(t*vv+1.0))*(pls*std::pow(s,2)+1.0)-4.0*(2.0*tt*v+3.0+2.0*t*vv+(16.0*s*tt*vv+15.0*tt+18.0*vv)*s)*std::pow(m2,2)*pls+2.0*((2.0*tt*v+3.0+4.0*t*vv+(8.0*tt+7.0*vv+s*tt*vv)*s)*pls*s+tt-vv-3.0*s*tt*vv)*m2)*aj8-2.0*(2.0*(((2.0*tt*v-17.0-2.0*s*tt)*std::pow(s,2)-(3.0*std::pow(t,2)-3.0*t*v+2.0*std::pow(v,2))-((2.0*tt*v-13.0)*v+20.0*t)*s+16.0*(tt*v-5.0-3.0*s*tt+4.0*m2*tt)*std::pow(m2,2))*pls-4.0*(((3.0*tt*v-16.0-4.0*s*tt)*s-((tt*v-4.0)*v+7.0*t))*pls-tt)*m2)*m2+((5.0*t-3.0*v+3.0*s)*s+2.0*(2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2)))*pls*s-2.0*t)*aj44-(2.0*((2.0*std::pow(t,2)*vv+4.0*t-3.0*v-2.0*(uu-vv)*ipow3(s)*tt-(2.0*uu-7.0*vv+2.0*tt)*std::pow(s,2)-((2.0*uu-5.0*vv)*t-2.0)*s)*pls*s+(2.0*uu+vv+tt)*s-(t*vv-4.0))*m2-(16.0*((tt*v-6.0+(uu-5.0*vv)*t-(tt-8.0*uu+14.0*vv-10.0*(uu-vv)*s*tt)*s)*pls+(uu-vv)*tt)*ipow3(m2)-(64.0*(2.0*(uu-2.0*vv-tt)+5.0*(uu-vv)*s*tt)*ipow4(m2)*pls-((pls*std::pow(s,2)+1.0)*(std::pow(s,2)+std::pow(t,2))*vv+256.0*(uu-vv)*ipow5(m2)*pls*tt)))-4.0*((2.0*std::pow(t,2)*vv+4.0*t-v-10.0*(uu-vv)*ipow3(s)*tt-5.0*(2.0*(uu-2.0*vv)+tt)*std::pow(s,2)-(2.0*tt*v-11.0+4.0*(uu-3.0*vv)*t)*s)*pls+2.0*(uu-vv+tt-(uu-vv)*s*tt))*std::pow(m2,2))*aj43+(((pls*std::pow(s,2)+1.0)*(s*tt+1.0)-32.0*ipow3(m2)*pls*tt)*vv+8.0*(tt+vv+3.0*s*tt*vv)*std::pow(m2,2)*pls-4.0*((tt+vv+2.0*s*tt*vv)*pls*s+tt*vv)*m2)*aj41+4.0*((std::pow(s,2)*tt*vv+s*vv+1.0)*pls*s-(s*tt+1.0)*vv-2.0*(3.0*s*vv-1.0-4.0*m2*vv)*m2*pls*s*tt)*aj4*m2-(2.0*((2.0*uu-3.0*vv+7.0*tt-12.0*s*tt*vv)*s-(2.0*(tt*v-1.0)+t*vv)-((2.0*uu-5.0*vv-5.0*tt+2.0*(uu+4.0*vv)*s*tt)*s+4.0*tt*v+3.0+(2.0*uu-3.0*vv)*t)*pls*std::pow(s,2))*m2+((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)*vv-2.0*s+v)*s-256.0*ipow5(m2)*pls*uu)*tt+64.0*(tt+2.0*uu+(5.0*uu-vv)*s*tt)*ipow4(m2)*pls-16.0*(((4.0*(2.0*uu-vv)+5.0*tt+(10.0*uu-vv)*s*tt)*s+2.0*tt*v+1.0+t*uu)*pls+(uu+3.0*vv)*tt)*ipow3(m2)+8.0*(((5.0*uu-7.0*vv+tt+(5.0*uu+4.0*vv)*s*tt)*s+3.0*(tt*v+1.0)+(2.0*uu-vv)*t)*pls*s-(uu-vv+2.0*tt-(uu+9.0*vv)*s*tt))*std::pow(m2,2))*aj25+(((pls*std::pow(s,2)+1.0)*(2.0*std::pow(s,2)-2.0*s*v+std::pow(v,2))-256.0*ipow5(m2)*pls*uu)*tt+64.0*(tt+2.0*uu+5.0*s*tt*uu)*ipow4(m2)*pls+8.0*((5.0*ipow3(s)*tt*uu+2.0*s*t*uu+tt*std::pow(v,2)+(6.0*tt+5.0*uu)*std::pow(s,2))*pls+2.0*tt-uu+s*tt*uu)*std::pow(m2,2)-16.0*((tt*v+1.0+t*uu+(5.0*tt+8.0*uu+10.0*s*tt*uu)*s)*pls+tt*uu)*ipow3(m2)-2.0*(((2.0*tt*v+5.0)*v-4.0*t+2.0*(4.0*tt+uu+s*tt*uu)*std::pow(s,2)-(3.0*tt*v+2.0-2.0*t*uu)*s)*pls*s+2.0*(3.0*tt-uu)*s-5.0*tt*v)*m2)*aj24-2.0*((3.0*s*vv-2.0)*pls*s-vv-4.0*m2*pls*s*vv)*aj20*m2*tt+4.0*(2.0*m2-s)*aj19*pls*tt-4.0*((2.0*tt*v+1.0-s*tt)*pls*s-tt-(3.0*(tt*v+1.0)-2.0*s*tt-4.0*m2*tt)*m2*pls)*aj17+((pls*std::pow(s,2)+1.0)*(tt+vv)-16.0*ipow3(m2)*pls*tt*vv+2.0*(tt-3.0*vv-2.0*s*tt*vv)*m2*pls*s+8.0*(tt+vv+2.0*s*tt*vv)*std::pow(m2,2)*pls)*aj16+2.0*(((3.0*tt*v-1.0-2.0*s*tt)*s-2.0*((2.0*tt*v-1.0)*v+2.0*t))*pls*s+2.0*(tt*v+1.0-s*tt)-8.0*(3.0*s*tt+4.0-4.0*m2*tt)*std::pow(m2,2)*pls-2.0*(((4.0*tt*v-5.0-6.0*s*tt)*s-3.0*(t+tt*std::pow(v,2)))*pls-4.0*tt)*m2)*aj14-(2.0*(8.0*(3.0*(tt+vv+s*tt*vv)-4.0*m2*tt*vv)*std::pow(m2,2)*pls+tt-3.0*vv-2.0*s*tt*vv+(s+6.0*t)*pls*s*vv)*m2-(2.0*tt*v+1.0+t*vv-(2.0*tt+vv)*s)*(pls*std::pow(s,2)+1.0)-4.0*((2.0*tt*v+3.0+2.0*t*vv+(5.0*(tt+2.0*vv)+2.0*s*tt*vv)*s)*pls-2.0*tt*vv)*std::pow(m2,2))*aj13+(((pls*std::pow(s,2)+1.0)*(s*tt+1.0)-32.0*ipow3(m2)*pls*tt)*vv+8.0*(tt+vv+3.0*s*tt*vv)*std::pow(m2,2)*pls-2.0*((5.0*tt+4.0*vv+2.0*s*tt*vv)*pls*s+tt*vv)*m2)*aj1);
    sr8 = -4.0*((2.0*(2.0*(4.0*((2.0*(7.0*t*vv-4.0+10.0*s*vv)*s+2.0*std::pow(t,2)*vv-7.0*t+5.0*v)*pls+6.0*vv-16.0*m2*pls*s*vv)*m2-(((4.0*(9.0*s*vv+12.0*t*vv-8.0)*s+22.0*std::pow(t,2)*vv-37.0*t+16.0*v)*s+2.0*ipow3(t)*vv-9.0*std::pow(t,2)+12.0*t*v-4.0*std::pow(v,2))*pls+4.0*(5.0*t*vv-2.0+5.0*s*vv)))*m2+(((14.0*s*vv+25.0*t*vv-20.0)*s+16.0*std::pow(t,2)*vv-29.0*t+12.0*v)*s+5.0*ipow3(t)*vv-13.0*std::pow(t,2)+12.0*t*v-3.0*std::pow(v,2))*pls*s+2.0*(5.0*s*vv+9.0*t*vv-6.0)*s+12.0*std::pow(t,2)*vv-14.0*t+5.0*v)*m2-((2.0*(2.0*(t*vv-1.0)+s*vv)*s+4.0*std::pow(t,2)*vv-6.0*t+3.0*v)*s+2.0*ipow3(t)*vv-4.0*std::pow(t,2)+3.0*t*v-std::pow(v,2)-((2.0*t-v)*(t-v)-2.0*ipow3(s)*vv-4.0*(t*vv-1.0)*std::pow(s,2)-(2.0*std::pow(t,2)*vv-6.0*t+3.0*v)*s)*pls*std::pow(s,2)))*aj6-2.0*(2.0*(((2.0*t-v)*std::pow((t-v),2)+6.0*ipow3(s)+2.0*(3.0*t-4.0*v)*std::pow(s,2)+(6.0*std::pow(t,2)-11.0*t*v+6.0*std::pow(v,2))*s)*pls+4.0*(t-v+s))*m2-(4.0*((4.0*(2.0*(t-v)+3.0*s)*s+4.0*std::pow(t,2)-6.0*t*v+3.0*std::pow(v,2))*pls+6.0)*std::pow(m2,2)+((std::pow(s,2)-s*v+std::pow(t,2)-2.0*t*v+std::pow(v,2))*pls*s+s+t-v)*(s+t-v)-64.0*ipow3(m2)*pls*s))*aj7-4.0*(2.0*(t-v+3.0*s-4.0*m2)*m2*pls-((t-v+s)*pls*s+1.0))*aj38+2.0*(2.0*((4.0*t-3.0*v)*(t-v)+4.0*std::pow(s,2)+(12.0*t-11.0*v)*s)*m2*pls-(8.0*(4.0*t-3.0*v+2.0*s)*std::pow(m2,2)*pls+(pls*std::pow(s,2)+3.0*pls*s*t-3.0*pls*s*v+2.0)*(s+t-v)))*aj37+2.0*(4.0*std::pow(m2,2)*pls-2.0*m2*pls*s+1.0)*(4.0*m2*vv-s*vv-t*vv+1.0)*aj33+2.0*(((2.0*(t*vv-1.0)+s*vv)*s+std::pow(t,2)*vv-2.0*t+v)*(pls*std::pow(s,2)+1.0)-8.0*(4.0*t*vv-7.0+2.0*s*vv)*ipow3(m2)*pls+2.0*((20.0*t*vv-21.0+12.0*s*vv)*s+4.0*std::pow(t,2)*vv-11.0*t+6.0*v)*std::pow(m2,2)*pls-(((16.0*t*vv-15.0+9.0*s*vv)*s+7.0*std::pow(t,2)*vv-13.0*t+6.0*v)*pls*s+2.0*(3.0*(t*vv-1.0)+s*vv))*m2)*aj32-((t*vv-1.0+s*vv)*(pls*std::pow(s,2)+1.0)-32.0*ipow3(m2)*pls*vv+4.0*(2.0*t*vv-3.0+8.0*s*vv)*std::pow(m2,2)*pls-2.0*((3.0*t*vv-4.0+5.0*s*vv)*pls*s+3.0*vv)*m2)*aj2-((t*vv-1.0+s*vv)*(pls*std::pow(s,2)+1.0)-32.0*ipow3(m2)*pls*vv+8.0*(t*vv-1.0+3.0*s*vv)*std::pow(m2,2)*pls-2.0*((2.0*t*vv-1.0+4.0*s*vv)*pls*s+vv)*m2)*aj10+(((3.0*t*vv-4.0+2.0*s*vv)*s+std::pow(t,2)*vv-3.0*t+2.0*v)*(pls*std::pow(s,2)+1.0)-128.0*ipow4(m2)*pls*s*std::pow(vv,2)+8.0*((4.0*(t*vv-6.0+3.0*s*vv)*s*vv-(4.0*t*vv-11.0))*pls+6.0*std::pow(vv,2))*ipow3(m2)+2.0*(((s*vv-14.0)*std::pow(s,2)*vv-(5.0*std::pow(t,2)*vv-12.0*t+8.0*v)-(std::pow(t,2)*std::pow(vv,2)+15.0*t*vv-21.0)*s)*pls*s+(s*vv-10.0)*s*vv+std::pow(t,2)*std::pow(vv,2)-9.0*t*vv+11.0)*m2-4.0*((2.0*((t*vv-16.0+3.0*s*vv)*s*vv-(11.0*t*vv-16.0))*s-(2.0*std::pow(t,2)*vv-7.0*t+6.0*v))*pls+4.0*(t*vv-4.0+s*vv)*vv)*std::pow(m2,2))*aj1)*std::pow(uu,2);
    sr9 = 4.0*(2.0*(((s*vv-2.0)*s-(std::pow(t,2)*vv-v)-(t-2.0*v+3.0*s)*pls*std::pow(s,2)+8.0*(5.0*t*vv+4.0+11.0*s*vv-12.0*m2*vv)*ipow3(m2)*pls+((9.0*t-13.0*v+2.0*std::pow(s,2)*vv+(2.0*t*vv+21.0)*s)*pls*s+4.0*(t*vv+2.0))*m2-4.0*(((5.0*t*vv+12.0+6.0*s*vv)*s+std::pow(t,2)*vv+3.0*t-5.0*v)*pls+2.0*vv)*std::pow(m2,2)+(4.0*(((17.0*t-12.0*v+20.0*s)*s+3.0*std::pow(t,2)-5.0*t*v+std::pow(v,2))*pls+6.0)*std::pow(m2,2)+(t-v+s)*t+(2.0*s-v)*(s+t)*pls*std::pow(s,2)-8.0*(6.0*t-5.0*v+14.0*s)*ipow3(m2)*pls-(((20.0*t-13.0*v+21.0*s)*s+7.0*std::pow(t,2)-13.0*t*v+2.0*std::pow(v,2))*pls*s+4.0*(3.0*t-v+2.0*s))*m2)*dz1)*aj8-(std::pow(s,2)+s*t-t*v+(s+t-v)*(s+t)*pls*std::pow(s,2)-8.0*(6.0*t-5.0*v+6.0*s)*ipow3(m2)*pls+4.0*((13.0*t-9.0*v+10.0*s)*s+3.0*std::pow(t,2)-5.0*t*v+std::pow(v,2))*std::pow(m2,2)*pls-(((11.0*s+18.0*t-10.0*v)*s+7.0*std::pow(t,2)-10.0*t*v+2.0*std::pow(v,2))*pls*s+2.0*(t-2.0*v+3.0*s))*m2)*aj6*dz1)-(((3.0*s*vv-1.0)*s+std::pow(t,2)*vv-t+v+((s*vv-3.0)*s-(std::pow(t,2)*vv-t-v))*pls*std::pow(s,2))*s-16.0*((6.0*std::pow(s,2)*vv-2.0*s*t*vv+v)*pls+4.0*vv)*ipow3(m2)+8.0*(((t*vv-8.0+8.0*s*vv)*s-(std::pow(t,2)*vv-2.0*t-v))*pls*s+2.0*t*vv-1.0+8.0*s*vv)*std::pow(m2,2)-2.0*((4.0*t*vv-3.0+11.0*s*vv)*s+std::pow(t,2)*vv-t+v+((2.0*(t*vv-7.0)+7.0*s*vv)*s-(std::pow(t,2)*vv-5.0*v))*pls*std::pow(s,2))*m2+(8.0*((5.0*s-2.0*v)*pls*s*v-2.0*(3.0*s-v)-2.0*((2.0*s-v)*pls*v-2.0)*m2)*std::pow(m2,2)-(4.0*std::pow(s,2)-2.0*s*v+std::pow(v,2)-(2.0*s-v)*pls*std::pow(s,2)*v)*s+2.0*(12.0*std::pow(s,2)-6.0*s*v+std::pow(v,2)-4.0*(2.0*s-v)*pls*std::pow(s,2)*v)*m2)*dz1)*aj45+(16.0*(((11.0*s*vv+16.0*t*vv-8.0)*s+5.0*std::pow(t,2)*vv-2.0*t-2.0*v)*pls+4.0*vv-8.0*(3.0*s+2.0*t-2.0*m2)*m2*pls*vv)*ipow3(m2)-(((t*vv-2.0)*s+t)*s+3.0*ipow3(t)*vv-5.0*std::pow(t,2)+3.0*t*v-std::pow(v,2)+(((t*vv-2.0)*s+t)*s-(ipow3(t)*vv-3.0*std::pow(t,2)+3.0*t*v-std::pow(v,2)))*pls*std::pow(s,2))+2.0*(((3.0*t-v)*(t-v)+ipow3(s)*vv+(6.0*t*vv-11.0)*std::pow(s,2)+(std::pow(t,2)*vv+2.0*t+v)*s)*pls*s+(4.0*t*vv-7.0+s*vv)*s+11.0*std::pow(t,2)*vv-9.0*t+5.0*v)*m2-8.0*(((t*vv-1.0)*std::pow(t,2)+4.0*ipow3(s)*vv+(5.0*t*vv-2.0)*s*t+(10.0*t*vv-11.0)*std::pow(s,2))*pls+8.0*t*vv-5.0+2.0*s*vv)*std::pow(m2,2)-(8.0*(3.0*(2.0*t-v+6.0*s)*pls*std::pow(s,2)+2.0*(2.0*t-v+5.0*s)+32.0*std::pow(m2,2)*pls*s)*std::pow(m2,2)+(2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2)+2.0*std::pow(s,2)-(2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2)-2.0*std::pow(s,2))*pls*std::pow(s,2))*s-16.0*((20.0*std::pow(s,2)+std::pow(v,2)+2.0*(2.0*t-v)*s)*pls+6.0)*ipow3(m2)-2.0*(2.0*(2.0*(2.0*t-v)+5.0*s)*s+2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2)+2.0*((2.0*t-v+7.0*s)*s-(std::pow(t,2)-t*v+std::pow(v,2)))*pls*std::pow(s,2))*m2)*dz1)*aj43-2.0*(2.0*(2.0*(t*vv-2.0+3.0*s*vv-4.0*m2*vv)*m2*pls-((t*vv-2.0+s*vv)*pls*s+2.0*vv))*m2+t*vv-1.0+s*vv-pls*std::pow(s,2)+(((s+t)*std::pow(s,2)-48.0*ipow3(m2)+4.0*(3.0*t-2.0*v+9.0*s)*std::pow(m2,2))*pls-2.0*((3.0*t-2.0*v+5.0*s)*pls*s+2.0)*m2)*dz1)*aj41+((pls*std::pow(s,2)-1.0)*(std::pow(s,2)*vv+s-t)-64.0*ipow3(m2)*pls*s*vv+8.0*((5.0*s+2.0*t)*pls*s+1.0)*std::pow(m2,2)*vv-2.0*((5.0*std::pow(s,2)*vv+3.0*s-v)*pls*s+s*vv+1.0)*m2+2.0*(16.0*std::pow(m2,2)*pls*s-6.0*m2*pls*std::pow(s,2)+3.0*m2*pls*s*v-6.0*m2+pls*ipow3(s)-pls*std::pow(s,2)*t-s+t)*(4.0*m2-s-t)*dz1)*aj4+(8.0*(3.0*(2.0*t-v+6.0*s)*pls*std::pow(s,2)+2.0*(2.0*t-v+5.0*s)+32.0*std::pow(m2,2)*pls*s)*std::pow(m2,2)+(2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2)+2.0*std::pow(s,2)-(2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2)-2.0*std::pow(s,2))*pls*std::pow(s,2))*s-16.0*((20.0*std::pow(s,2)+std::pow(v,2)+2.0*(2.0*t-v)*s)*pls+6.0)*ipow3(m2)-2.0*(2.0*(2.0*(2.0*t-v)+5.0*s)*s+2.0*std::pow(t,2)-2.0*t*v+std::pow(v,2)+2.0*((2.0*t-v+7.0*s)*s-(std::pow(t,2)-t*v+std::pow(v,2)))*pls*std::pow(s,2))*m2)*aj36*dz1-2.0*(2.0*((3.0*t-2.0*v+3.0*s)*pls*s+4.0)*m2+48.0*ipow3(m2)*pls-36.0*std::pow(m2,2)*pls*s-12.0*std::pow(m2,2)*pls*t+8.0*std::pow(m2,2)*pls*v-s-t)*aj32*dz1+(((3.0*s*vv-1.0)*s+std::pow(t,2)*vv-t+v+((s*vv-3.0)*s-(std::pow(t,2)*vv-t-v))*pls*std::pow(s,2))*s-16.0*((6.0*std::pow(s,2)*vv-2.0*s*t*vv+v)*pls+4.0*vv)*ipow3(m2)+8.0*(((t*vv-8.0+8.0*s*vv)*s-(std::pow(t,2)*vv-2.0*t-v))*pls*s+2.0*t*vv-1.0+8.0*s*vv)*std::pow(m2,2)-2.0*((4.0*t*vv-3.0+11.0*s*vv)*s+std::pow(t,2)*vv-t+v+((2.0*(t*vv-7.0)+7.0*s*vv)*s-(std::pow(t,2)*vv-5.0*v))*pls*std::pow(s,2))*m2+(8.0*((5.0*s-2.0*v)*pls*s*v-2.0*(3.0*s-v)-2.0*((2.0*s-v)*pls*v-2.0)*m2)*std::pow(m2,2)-(4.0*std::pow(s,2)-2.0*s*v+std::pow(v,2)-(2.0*s-v)*pls*std::pow(s,2)*v)*s+2.0*(12.0*std::pow(s,2)-6.0*s*v+std::pow(v,2)-4.0*(2.0*s-v)*pls*std::pow(s,2)*v)*m2)*dz1)*aj25+((2.0*(s-v)*s+2.0*std::pow(t,2)+std::pow(v,2)+(2.0*(s-v)*s-(2.0*std::pow(t,2)-std::pow(v,2)))*pls*std::pow(s,2))*t+256.0*ipow4(m2)*pls*s-16.0*((12.0*std::pow(s,2)+std::pow(v,2)+4.0*(3.0*t-v)*s)*pls+6.0)*ipow3(m2)+8.0*((6.0*ipow3(s)+t*std::pow(v,2)+2.0*(2.0*t-v)*(t-v)*s+2.0*(7.0*t-2.0*v)*std::pow(s,2))*pls+2.0*(5.0*t-v+2.0*s))*std::pow(m2,2)+2.0*((3.0*(t-v)*t*v-2.0*ipow3(s)-2.0*(6.0*t-v)*std::pow(s,2)-(2.0*std::pow(t,2)-9.0*t*v+std::pow(v,2))*s)*pls*s-(2.0*(4.0*t-v+s)*s+10.0*std::pow(t,2)-4.0*t*v+std::pow(v,2)))*m2)*aj23*dz1-((s*vv+1.0-2.0*m2*vv)*(pls*std::pow(s,2)-1.0)+8.0*std::pow(m2,2)*pls*s*vv-2.0*(s+t-4.0*m2)*(pls*std::pow(s,2)-1.0)*dz1)*aj20-((pls*std::pow(s,2)+1.0)*(t*vv-1.0)-16.0*ipow3(m2)*pls*vv+8.0*(t*vv-2.0+2.0*s*vv)*std::pow(m2,2)*pls-2.0*((3.0*t*vv-5.0+s*vv)*pls*s+vv)*m2)*aj16+(16.0*(3.0*t*vv+1.0+s*vv-4.0*m2*vv)*ipow3(m2)*pls-(std::pow(t,2)*vv-v-(t*vv-2.0)*s)*(pls*std::pow(s,2)+1.0)+8.0*(((s*vv-6.0)*s-(std::pow(t,2)*vv-3.0*v))*pls+vv)*std::pow(m2,2)+2.0*((4.0*(t-2.0*v)-std::pow(s,2)*vv-3.0*(t*vv-3.0)*s)*pls*s-(t*vv-6.0+s*vv))*m2)*aj13-((t*vv-2.0+s*vv)*(pls*std::pow(s,2)+1.0)-32.0*ipow3(m2)*pls*vv+8.0*(t*vv-2.0+2.0*s*vv)*std::pow(m2,2)*pls-2.0*((2.0*t*vv-3.0+5.0*s*vv)*pls*s+2.0*vv)*m2)*aj1)*uu;
    sr10 = 4.0*(2.0*(2.0*(2.0*(t-v+3.0*s-4.0*m2)*m2*pls-((t-v+s)*pls*s+1.0))*aj9-(4.0*std::pow(m2,2)*pls-2.0*m2*pls*s+1.0)*aj8-2.0*std::pow((2.0*m2-s),2)*aj51+(4.0*(2.0*t*uu-1.0+2.0*s*uu-4.0*m2*uu)*m2-((4.0*t*uu-1.0)*s-(t-v)))*aj50+(pls*std::pow(s,2)-1.0)*aj5+2.0*(16.0*(t-v+2.0*s)*ipow3(m2)*pls-(32.0*ipow4(m2)*pls+std::pow(t,2))+(std::pow((s+t-v),2)*pls*s+4.0*t)*m2-2.0*((5.0*s+t-v)*(s+t-v)*pls+2.0)*std::pow(m2,2))*aj49-(2.0*(4.0*((4.0*t-3.0*v+6.0*s)*pls-2.0*uu-8.0*m2*pls)*m2+2.0*(2.0*t*uu-1.0+2.0*s*uu)-(6.0*s+2.0*t-v)*(s+t-v)*pls)*m2-((4.0*t*uu-1.0)*s-(t-v)-std::pow((s+t-v),2)*pls*s))*aj48+4.0*std::pow((2.0*m2-s),2)*aj47-(4.0*(2.0*t*uu-1.0+2.0*s*uu+3.0*pls*std::pow(s,2)-4.0*(pls*s+uu)*m2)*m2-((4.0*t*uu-3.0)*s-(t-v)+2.0*pls*ipow3(s)))*aj46)-(2.0*(2.0*(2.0*t*uu-1.0)*s*uu-(2.0*t*uu-3.0)-(6.0*s*t*uu-3.0*s+2.0*std::pow(t,2)*uu+2.0*t+3.0*v)*pls*s+64.0*ipow3(m2)*pls*uu-16.0*((2.0*t*uu+1.0+4.0*s*uu)*pls-std::pow(uu,2))*std::pow(m2,2))*m2+(2.0*s*uu-1.0)*s+2.0*std::pow(t,2)*uu+t+v+(t+v-s)*pls*std::pow(s,2)+8.0*(((8.0*t*uu+1.0+3.0*s*uu)*s+std::pow(t,2)*uu+t+v)*pls-2.0*(t*uu+1.0+s*uu)*uu)*std::pow(m2,2))*aj45+2.0*((std::pow(s,2)+2.0*std::pow(t,2)+(t+v)*s)*pls*s-2.0*t+2.0*(4.0*(t+2.0*v+s)*m2-(3.0*s+t)*(s+t+v))*m2*pls)*aj44+(2.0*(2.0*(2.0*t*uu-1.0)*s*uu-(2.0*t*uu-3.0)-(6.0*s*t*uu-3.0*s+2.0*std::pow(t,2)*uu+2.0*t-v)*pls*s+64.0*ipow3(m2)*pls*uu-8.0*((4.0*t*uu+1.0+8.0*s*uu)*pls-2.0*std::pow(uu,2))*std::pow(m2,2))*m2+(2.0*s*uu+1.0)*s+2.0*std::pow(t,2)*uu-t+v+(t-v-s)*pls*std::pow(s,2)+4.0*(((16.0*t*uu+1.0+6.0*s*uu)*s+2.0*std::pow(t,2)*uu+t-2.0*v)*pls-4.0*(t*uu+1.0+s*uu)*uu)*std::pow(m2,2))*aj43-(pls*std::pow(s,2)-1.0)*aj4-4.0*std::pow((2.0*m2-s),2)*aj31-4.0*(2.0*(4.0*m2-3.0*s)*m2*pls+pls*std::pow(s,2)-1.0)*aj30*s-8.0*std::pow((2.0*m2-s),2)*aj27+2.0*(4.0*(2.0*t*uu-1.0+2.0*s*uu+3.0*pls*std::pow(s,2)-4.0*(pls*s+uu)*m2)*m2-((4.0*t*uu-3.0)*s-(t-v)+2.0*pls*ipow3(s)))*aj26+(2.0*(2.0*(2.0*t*uu-1.0)*s*uu-(2.0*t*uu-3.0)-(6.0*s*t*uu-3.0*s+2.0*std::pow(t,2)*uu+2.0*t+3.0*v)*pls*s+64.0*ipow3(m2)*pls*uu-16.0*((2.0*t*uu+1.0+4.0*s*uu)*pls-std::pow(uu,2))*std::pow(m2,2))*m2+(2.0*s*uu-1.0)*s+2.0*std::pow(t,2)*uu+t+v+(t+v-s)*pls*std::pow(s,2)+8.0*(((8.0*t*uu+1.0+3.0*s*uu)*s+std::pow(t,2)*uu+t+v)*pls-2.0*(t*uu+1.0+s*uu)*uu)*std::pow(m2,2))*aj25+2.0*(16.0*((2.0*t*uu+1.0+4.0*s*uu)*pls-std::pow(uu,2)-4.0*m2*pls*uu)*ipow3(m2)-(std::pow(t,2)*uu+v+std::pow(s,2)*uu-(s-t)*pls*std::pow(s,2))+((2.0*std::pow(t,2)*uu+2.0*t+5.0*v+6.0*(t*uu-1.0)*s)*pls*s-2.0*((2.0*t*uu-1.0)*s*uu-(t*uu-2.0)))*m2-4.0*(((8.0*t*uu+1.0+3.0*s*uu)*s+std::pow(t,2)*uu+t+v)*pls-2.0*(t*uu+1.0+s*uu)*uu)*std::pow(m2,2))*aj24-2.0*(2.0*m2-s)*aj17*pls+2.0*((2.0*t+v)*pls*s-1.0-2.0*(t+2.0*v-2.0*m2)*m2*pls)*aj14-(4.0*(2.0*m2-s)*m2*pls+pls*std::pow(s,2)+1.0)*aj13);

    double result = coer * (sr1 + sr2 + sr3 + sr4 + sr5 + sr6 + sr7 + sr8 + sr9 + sr10);

    if (result < 0.0) {
        nn += 1;
    }

    if (ikey == -1) {
        fir = -4.0 *
              ((aj5 + aj7 + m2 * aj1 * vv * vv + aj14) +
               (t - 2.0 * m2) * (aj23 + aj13 * vv) +
               (s - 2.0 * m2) * (aj36 + vv * aj4) -
               (s + t - 2.0 * m2) * (vv * aj6 + aj24));

        result -= alfa / 4.0 / (pi * pi) * fir * sig(t, pl, 0);
    }

    // Grid calls only (ikey 0/1/2). Skip simpson ikey=-1 (thousands of samples).
    if (ikey >= 0) {
      ptrc_i("fsir.ikey", ikey);
      ptrc_d("fsir.t", t);
      ptrc_d("fsir.t1", t1);
      ptrc_d("fsir.v", v);
      ptrc_d("fsir.z", z);
      ptrc_d("fsir.az", az);
      ptrc_d("fsir.bz", bz);
      ptrc_d("fsir.cz", cz);
      const double ajs[51] = {
          aj1,  aj2,  aj3,  aj4,  aj5,  aj6,  aj7,  aj8,  aj9,  aj10,
          aj11, aj12, aj13, aj14, aj15, aj16, aj17, aj18, aj19, aj20,
          aj21, aj22, aj23, aj24, aj25, aj26, aj27, aj28, aj29, aj30,
          aj31, aj32, aj33, aj34, aj35, aj36, aj37, aj38, aj39, aj40,
          aj41, aj42, aj43, aj44, aj45, aj46, aj47, aj48, aj49, aj50,
          aj51};
      const double srs[10] = {sr1, sr2, sr3, sr4, sr5, sr6, sr7, sr8, sr9, sr10};
      char tag[8];
      for (int i = 0; i < 51; i++) {
        std::snprintf(tag, sizeof tag, "aj%02d", i + 1);
        ptrc_d(tag, ajs[i]);
      }
      for (int i = 0; i < 10; i++) {
        std::snprintf(tag, sizeof tag, "sr%d", i + 1);
        ptrc_d(tag, srs[i]);
      }
      ptrc_d("fsir.out", result);
    }

    return result;
}

}  // namespace meradgen

