#include <cmath>
#include <iostream>

double fsir(double t, double t1, double v, double z, double pl, int nn, int ikey) {
    // The cross section of real photon emission
    // t, t1, v, z are kinematic invariant
    // pl is degree of polarization
    // nn number of point where cross section is negative
    // ikey=0 radiative cross section ds/dt/dv/dt1/dz
    // ikey=1 radiative cross section integrated over z ds/dt/dv/dt1
    // ikey=2 radiative cross section integrated over z and t1 ds/dt/dv
    // ikey=-1 finite part of radiative cross section integrated over z and t1 ds/dt/dv

    const double m2 = 1.0; // Placeholder for m2, should be defined appropriately
    const double pi = 3.141592653589793;
    const double als = 1.0; // Placeholder for als, should be defined appropriately
    const double alfa = 1.0; // Placeholder for alfa, should be defined appropriately
    double u, uu, v1, z1, z2, zz, zz1, zz2, ds, dz, dz1, vt, pls, tt, tt1, dt, fir;
    double aj[51] = {0}; // Array to hold aj1 to aj51

    u = v - s - t + 4.0 * m2;
    uu = 1.0 / u;
    v1 = s + u + t1 - 4.0 * m2;
    z1 = z - t1 + t;
    z2 = z - s - t1 + 4.0 * m2;
    zz = 1.0 / z;
    zz1 = 1.0 / z1;
    zz2 = 1.0 / z2;
    vv1 = 1.0 / v1;
    vv = 1.0 / v;
    pls = -pl / als;
    tt = 1.0 / t;
    tt1 = 1.0 / t1;
    dt = 1.0 / (t - t1);
    dz = 1.0 / (s + t1 - 4 * m2);
    dz1 = 1.0 / (s + t - 4 * m2);
    ds = 1.0 / (s - 4 * m2);
    vt = 1.0 / (v - t);


    if (ikey == 0) {
        sd = 1.0 / sqrt(-az * z * z - 2.0 * bz * z - cz);
        aj1 = sd;
        aj2 = z * sd;
        aj3 = z * z * sd;
        aj4 = zz * sd;
        aj5 = m2 * zz * zz * sd;
        aj6 = zz1 * sd;
        aj7 = m2 * zz1 * zz1 * sd;
        aj8 = zz2 * sd;
        aj9 = m2 * zz2 * zz2 * sd;
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
        aj47 = m2 * dz * dz * aj8;
        aj48 = vv1 * aj9;
        aj49 = vv1 * vv1 * aj9;
        aj50 = dz * aj9;
        aj51 = dz * dz * aj9;
    } else if (ikey == 1) {
        aj1 = pi / sqrt(az);
        aj2 = -pi * bz / pow(az, 1.5);
        aj3 = pi * (3.0 * bz * bz - az * cz) / (2.0 * pow(az, 2.5));
        aj4 = pi / sqrt(cz);
        aj5 = -m2 * pi * bz / pow(cz, 1.5);
        aj6 = pi / sqrt(cz1);
        aj7 = -m2 * pi * bz1 / pow(cz1, 1.5);
        aj8 = -pi / sqrt(cz2);
        aj9 = m2 * pi * bz2 / pow(cz2, 1.5);
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
        aj27 = m2 * dz * dz * aj4;
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
        aj47 = m2 * dz * dz * aj8;
        aj48 = vv1 * aj9;
        aj49 = vv1 * vv1 * aj9;
        aj50 = dz * aj9;
        aj51 = dz * dz * aj9;
    } else if (ikey == 2 || ikey == -1) {
        aj1 = pi * v / (m2 + v);
        aj2 = -pi * v * v * (v - s + 2.0 * m2) / (2.0 * pow(m2 + v, 2));
        aj3 = 0.0; // Placeholder for actual calculation
        aj4 = pi / sqrt((s - v) * (s - v) - 4.0 * m2 * s) *
              log(pow((s - v - 2.0 * m2 + sqrt((s - v) * (s - v) - 4.0 * m2 * s)), 2) / (4.0 * m2 * (m2 + v)));
        aj5 = pi / v;
        aj6 = pi / sqrt((v - u) * (v - u) - 4.0 * m2 * u) *
              log(pow((v - u + 2.0 * m2 + sqrt((v - u) * (v - u) - 4.0 * m2 * u)), 2) / (4.0 * m2 * (m2 + v)));
        aj7 = pi / v;
        aj8 = pi * log(4.0 * m2 * u * u * (m2 + v) / (v * (v - u + sqrt((v - u) * (v - u) - 4.0 * m2 * u)) 
              - 2.0 * m2 * u) * (v - u) * (v - u)) / sqrt((v - u) * (v - u) - 4.0 * m2 * u);
        aj9 = pi * v / (u * u);
        aj10 = pi * v * (2.0 * m2 * t + (t - v) * v) / (2.0 * (m2 + v) * (m2 + v));
        aj11 = pi * log(4.0 * m2 * t * t * (m2 + v) / (v * (v - t + sqrt((v - t) * (v - t) - 4.0 * m2 * t)) 
             - 2.0 * m2 * t) * (v - t) * (v - t)) / sqrt((v - t) * (v - t) - 4.0 * m2 * t);
        aj12 = pi * v / (t * t);
        aj13 = pi * log(pow((v - t + 4.0 * m2 + sqrt((v - t) * (v - t) - 4.0 * m2 * t)), 2) 
             / (4.0 * m2 * (m2 + v))) / sqrt((v - t) * (v - t) - 4.0 * m2 * t);
        aj14 = pi / v;
        aj15 = pi * v * (v * (2.0 * m2 - t + v) - s * (t + v)) / ((v - t) * (v - t) - 4.0 * m2 * t) / (m2 + v) +
               pi * t * (s * (t - v) + 2.0 * m2 * v) / pow((v - t) * (v - t) - 4.0 * m2 * t, 1.5) *
               log(4.0 * m2 * t * t * (m2 + v) / (v * (v - t + sqrt((v - t) * (v - t) - 4.0 * m2 * t)) 
               - 2.0 * m2 * t) * (v - t) * (v - t));
        aj16 = pi * v * (v * (2.0 * m2 - t + v) - s * (t + v)) / ((v - t) * (v - t) - 4.0 * m2 * t) / (m2 + v) -
               pi * v * (u * (t - v) + 2.0 * m2 * v) / pow((v - t) * (v - t) - 4.0 * m2 * t, 1.5) *
               log(4.0 * m2 * (m2 + v) / (v - t + sqrt((v - t) * (v - t) - 4.0 * m2 * t) - 2.0 * m2) * (v - t) * (v - t));
        aj17 = -pi * (1.0 + (s * (t - v) + 2.0 * m2 * v) / ((v - t) * (v - t) - 4.0 * m2 * t)) +
               m2 * pi * (s * (t + v) + v * (t - v - 2.0 * m2)) / pow((v - t) * (v - t) - 4.0 * m2 * t, 1.5) *
               log(4.0 * m2 * (m2 + v) / (v - t + sqrt((v - t) * (v - t) - 4.0 * m2 * t) - 2.0 * m2) 
               * (v - t) * (v - t));
        aj18 = pi * ((v * v * (4.0 * m2 * m2 * (4.0 * s * t + (8.0 * t - 3.0 * v) * v) - 2.0 * m2 * m2 
               * (2.0 * s * t * (s + 7.0 * t) - 2.0 * (s - 6.0 * t) * t * v - (4.0 * s + 31.0 * t) * v * v 
               + 13.0 * v * v * v) + (t - v) * (3.0 * (t - v) * (t - v) * v * v - 2.0 * s * (t - v) * v * (2.0 * t 
               + 3.0 * v) + s * s * (-t * t + 4.0 * t * v + 3.0 * v * v)) - 2.0 * m2 * (s * s * (4.0 * t * t 
               - t * v - v * v) + s * (t + v) - (3.0 * t * t - 15.0 * t * v + 8.0 * v * v) + (t - v) * v * (2.0 
               * t * t - 13.0 * t * v + 8.0 * v * v)))) / (2.0 * (-4.0 * m2 * t + (t - v) * (t - v)) * (m2 + v) 
               * (m2 + v)) - (v * v * (u * u * (t - v) * (t - v) + 6.0 * m2 * m2 * v * v - 2.0 * m2 * u * (s * t 
               + 2.0 * v * (v - t))) * log((4.0 * m2 * (m2 + v)) / (2.0 * m2 - t + sqrt(-4.0 * m2 * t + (t - v) 
               * (t - v) + v)) * (2.0 * m2 - t + sqrt(-4.0 * m2 * t + (t - v) * (t - v) + v)) * (2.0 * m2 - t 
               + sqrt(-4.0 * m2 * t + (t - v) * (t - v) + v))) / (-4.0 * m2 * t + (t - v) * (t - v)) * (t - v) 
               * (t - v)) / (2.0 * (m2 + v) * (m2 + v)));
        aj19 = pi * ((2.0 * v * ((u - 4.0 * m2) * (u - 4.0 * m2) * (t - v) * (t - v) - 4.0 * m2 * m2 * (4.0 
               * (s - t) * t - 4.0 * t * v - 3.0 * v * v) + 4.0 * m2 * m2 * ((s - 2.0 * t) * t 
               * (s + t) - 3.0 * (s - 3.0 * t) * t * v - (2.0 * s + 9.0 * t) * v * v + 4.0 * v * v * v)) /(((t - v) 
               * (t - v) - 4.0 * m2 * t) * (m2 + v)) - m2 * (4.0 * v * (-2.0 * m2 * m2 * (4.0 * s * t 
               + (4.0 * t - 3.0 * v) * v) + m2 * (2.0 * s * t * (s + t) - 2.0 * t * (s + t) * v 
               - 3.0 * t * v * v + v * v) + u * (t - v) * ((t - v) * v - s * (2.0 * t + v))) 
               * log((4.0 * m2 * (m2 + v)) / (2.0 * m2 - t + sqrt(-4.0 * m2 * t + (t - v) * (t - v) + v)) 
               * (2.0 * m2 - t + sqrt(-4.0 * m2 * t + (t - v) * (t - v) + v))) / ((t - v) * (t - v) - 4.0 * m2 * t) 
               * (t - v) * (t - v)) / (2.0 * (m2 + v) * (m2 + v)));
        aj20 = -((s * (t + v) - v * (v - t + 2.0 * m2)) * aj1 + (s * t * (v - s) + m2 * (4.0 * s * t - 2.0 * v * v)) 
               * aj4) / ((s - v) * (s - v) - 4.0 * m2 * s);
        aj21 = pi / t * log(((s - 2.0 * m2) * (s - 2.0 * m2 + sqrt(als))) / (2.0 * m2 * m2) - 1.0) / sqrt(als);
        aj22 = -pi / (t * t * t) * ((s * (v - t) - 2.0 * m2 * v) * v + log(((s - 2.0 * m2) * ((s - 2.0 * m2) 
               + sqrt(als))) / (2.0 * m2 * m2) - 1.0) * (t * (s * v - als) - 2.0 * m2 * v * v) * m2 / sqrt(als));
        aj23 = 2.0 * pi * log((sqrt(t * t - 4.0 * m2 * t) - t + 2.0 * m2) / (2.0 * m2)) / v / sqrt(t * t - 4.0 * m2 * t);
        aj24 = 2.0 * pi * log((sqrt(u * u - 4.0 * m2 * u) - u + 2.0 * m2) / (2.0 * m2)) / v / sqrt(u * u - 4.0 * m2 * u);
        aj25 = -pi * log((m2 * m2 * (v * v - 2.0 * u * u) + (m2 * v - (s - 2.0 * m2) * u) * (s * u - 2.0 * m2 
               * (u + v)) / (2.0 * m2 * (m2 * (u + v) * (u + v) - s * u * v))) / u / sqrt(als);
        aj26 = (-pi / als / u / u * (v * (s * u * (s + u - 2.0 * v) + 2.0 * m2 * (v * (u + v) - 2.0 * s * u)) 
               / (s * u * v - m2 * (u + v) * (u + v)) + (s * u * (v - s) + m2 * (4.0 * s * u 
               - 2.0 * v * v)) * aj25 / (u * u) / als) * m2;
        aj27 = (-pi * v * (s * s * u * u * (u * u - 3.0 * s * s - 2.0 * s * u + 2.0 * (5.0 * s + u) * v 
               - 6.0 * v * v) + 2.0 * m2 * s * s * u * u * (2.0 * s * u * u * (s + u) - u * (s * u 
               - 15.0 * s * s - 4.0 * u * u) * v + (s * s - 34.0 * s * u + 2.0 * u * u) * v * v - (11.0 * s 
               - 5.0 * u) * v * v * v + 12.0 * v * v * v) + 4.0 * m2 * m2 * s * u * (2.0 * s * u * u * (s + u) 
               - 2.0 * s * (s - 2.0 * u) * u * v + (s - 2.0 * u) * (s - 2.0 * u) * v * v)) * m2 / (2.0 * u * u 
               * (s * u * v - m2 * (u + v) * (u + v)) * (s * u * v - m2 * (u + v) * (u + v)) * als * als) 
               + (-aj25 * (2.0 * m2 * s * u * (3.0 * s + u - 3.0 * v) * v * v + 6.0 * m2 * m2 * v * v) 
               - (v * v - 4.0 * s * u) + u * u * (als - s * v) * (s * u * v - m2 * (u + v) * (u + v))) 
               * m2 / (u * u * als * als);
        aj28 = m2 * pi / (t * t * als) * ((s * (v - t) - 2.0 * m2 * v) * log((s - 2.0 * m2) 
               * (s - 2.0 * m2 + sqrt(als)) / (2.0 * m2 * m2) - 1.0) / sqrt(als) 
               - (t * (s * v - als) - 2.0 * m2 * v * v) / m2 / v);
        aj29 = -pi / (t * t * t) / pow(als, 2.5) * (sqrt(als) * (als * t * t * (2.0 * s * v - als) + v * v 
               * (2.0 * s * t * v * (s + 2.0 * m2) - v * v * ((s - 2.0 * m2) * (s - 2.0 * m2) + 8.0 * m2) 
               - 2.0 * s * t * (s * t + 2.0 * (2.0 * s + t) * m2 - 16.0 * m2 * m2)) / v + 2.0 * m2 * (2.0 * m2 
               * m2 * v * (8.0 * s * t - 3.0 * v * v) + s * v * m2 * (3.0 * v * v + 2.0 * t * (v - 2.0 * s - t)) 
               - s * t * (2.0 * s * u * v + t * als)) * log((s - 2.0 * m2) * (s - 2.0 * m2 + sqrt(als)) 
               / (2.0 * m2 * m2) - 1.0));
        aj30 = -(pi * (als * u - s * u * v + 2.0 * m2 * v * v) / v / als - m2 * (m2 * (4.0 * s * u - 2.0 * v 
               * (u + v)) - s * u * (s + u - 2.0 * v)) * aj25) / (u * u);
        aj31 = pi * (u * (u * u / v + (u * u * (v - 4.0 * m2 - 2.0 * s)) / als + -12.0 * v * (-s * t * u + m2 
               * v * v) / (als * als) + (u * u * ((s - u) * (s - u) * v + -4.0 * m2 * (u * u + (u - s) * v))) 
               * m2 / (als * (-s * u * v + m2 * (u + v) * (u + v))) + -(2.0 * (3.0 * s * u * (s * u - 8.0 * m2 * m2 
               + 2.0 * m2 * s) * v * v - 9.0 * m2 * s * u * v * v + -s * als * u * u * (v - t) + 6.0 * m2 
               * m2 * v * v * (u + v) - 2.0 * s * (s - m2) * u * u * v * (2.0 * (v - t) - u)) * m2 
               / (als * als) / u / u);
        aj32 = ((u * (t + v) - v * (v - t + 2.0 * m2)) * aj1 - (u * t * (v - u) + m2 * (4.0 * u * t - 2.0 * v * v)) 
               * aj6) / ((u - v) * (u - v) - 4.0 * m2 * u);
        aj33 = (-(pi * v * (4.0 * m2 * (4.0 * m2 - u) * u * u * (t - v) * (t - v) + (t - v) * u 
               * -(32.0 * m2 * m2 - 4.0 * m2 * (5.0 * s - 24.0 * u) + 3.0 * u * u * (s + u) 
               - 2.0 * m2 * u * (6.0 * s + 19.0 * u)) * v + u * (128.0 * m2 * m2 - 4.0 * m2 * (33.0 * s 
               - 40.0 * u) - (s - 9.0 * u) * u * (s + u) + 2.0 * m2 * (11.0 * s * s - 6.0 * s * u - 42.0 * u * u)) 
               * v * v - (20.0 * m2 * m2 - 2.0 * m2 * (8.0 * s + 5.0 * u) - 2.0 * m2 * (s * s + 19.0 * s * u 
               - 21.0 * u * u) + u * (9.0 * u * u - 5.0 * s * s - 2.0 * s * u)) * v * v + 2.0 * m2 * v * v * v)) 
               / (2.0 * (m2 + v) * (m2 + v)) + -aj6 * (v * v * (6.0 * m2 * m2 * v * v - (s * t * u 
               * (2.0 * m2 + u))) - t * u * (u - 4.0 * m2) * (v * v - 2.0 * t * u * v + (u - 4.0 * m2) 
               * (t * u - v * v))) / ((u - v) * (u - v) - 4.0 * m2 * u) * (u - v) * (u - v));
        aj34 = -pi * 2.0 * log(2.0 * m2 / (2.0 * m2 - u + sqrt(u * (u - 4.0 * m2)))) / t / sqrt(u * (u - 4.0 * m2));
        aj35 = -pi * (v * (u * (t - v) + 2.0 * m2 * v) + 2.0 * m2 * (t * u * (s + u) - 2.0 * m2 * v * v) 
               * log((sqrt(u * (u - 4.0 * m2)) - u) / (2.0 * m2) - 1.0) / sqrt(u * (u - 4.0 * m2))) 
               / (t * t * u * (u - 4.0 * m2));
        aj36 = pi * 2.0 * log((s - 2.0 * m2 + sqrt(als)) / (2.0 * m2)) / v / sqrt(als);
        aj37 = pi * ((t / v + (s * (u - v) + 2.0 * m2 * v) / (4.0 * m2 * u - (u - v) * (u - v)) - 1.0) 
               + m2 * (t * u - s * v + 2.0 * m2 * v) / pow((u - v) * (u - v) - 4.0 * m2 * u, 1.5) 
               * log(pow((v - u + 2.0 * m2 + sqrt((v - u) * (v - u) - 4.0 * m2 * u)), 2) 
               / (4.0 * m2 * (m2 + v)));
        aj38 = pi * ((2.0 * (4.0 * m2 * u - (u - v) * (u - v)) * (m2 * t * t * u * u * (u - 4.0 * m2) 
               * (u - 4.0 * m2 - 2.0 * v) - t * t * u * u * (24.0 * m2 * m2 - 10.0 * m2 * u + u * u) 
               * v + 2.0 * t * u * -(m2 * (22.0 * m2 - 5.0 * s) * u - 2.0 * m2 * m2 * (4.0 * m2 + s) 
               + -(s - 9.0 * m2) * u * u + u * u * u) * v * v - t * u * (16.0 * m2 * m2 + 6.0 * m2 
               * (s - 3.0 * u) + -u * (s + 3.0 * u)) * v * v + (12.0 * m2 * m2 - u * u * (s + u) 
               - 4.0 * m2 * m2 * (s + 4.0 * u) + -m2 * (s * s + 4.0 * s * u + 8.0 * u * u)) * v * v 
               + (8.0 * m2 * m2 - 4.0 * m2 * u + u * u) * v * v)) / (v * (m2 + v)) + 4.0 * m2 * sqrt((u - v) 
               * (u - v) - 4.0 * m2 * u) * (t * t * (4.0 * m2 - u) * u * u - t * u * (8.0 * m2 * m2 + (s - u) 
               * u + m2 * (2.0 * u - 6.0 * s)) * v - t * u * (u - 2.0 * s) * v * v - m2 * (2.0 * (m2 - s) + u) 
               * v * v + m2 * v * v * v) * log((2.0 * m2 - u + sqrt((u - v) * (u - v) - 4.0 * m2 * u) + v) 
               / (2.0 * (4.0 * m2 * (m2 + v)))) / (2.0 * (4.0 * m2 * u - (u - v) * (u - v)) * (4.0 * m2 * u 
               - (u - v) * (u - v)) * (4.0 * m2 * u - (u - v) * (u - v)));
        aj39 = pi / (t * t) / (u * (u - 4.0 * m2)) * ((u * (v - t) - 2.0 * m2 * v) * 2.0 
               * log(2.0 * m2 / (2.0 * m2 - u + sqrt(u * (u - 4.0 * m2)))) * m2 / sqrt(u * (u - 4.0 * m2)) 
               - (u * t * (v - u) + m2 * (4.0 * u * t - 2.0 * v * v)) / v);
        aj40 = pi * ((u * (pow(t, 2) * pow(u, 2) * (u - 4.0 * m2) * (u - 4.0 * m2 - 2.0 * v) - 2.0 * t * u 
               * (8.0 * pow(m2, 2) + 2.0 * m2 * (s - 3.0 * u) + u * (s + u)) * pow(v, 2)
               - (12.0 * pow(m2, 2) - 4.0 * m2 * u + pow(u, 2)) * pow(v, 4)) / ((u - 4.0 * m2) * (u - 4.0 * m2) * v) 
               - m2 * (2.0 * sqrt(-u) * (-(pow(t, 2) * (4.0 * m2 - u) * pow(u, 2)) - 2.0 * t * u * (-m2 * 
               (4.0 * m2 + s) + (m2 + s) * u) * v + 3.0 * m2 * (2.0 * m2 - u) * pow(v, 3)) 
               * -log(((2.0 * m2 - u) * (2.0 * m2 + sqrt(4.0 * m2 - u) * sqrt(-u) - u) / (2.0 * pow(m2, 2)) - 1.0)) 
               / pow((4.0 * m2 - u), 2.5)) / (pow(t, 4) * pow(u, 3)));
        aj41 = -pi * (v * (s * v - t * u - 2.0 * m2 * v) / (m2 + v) -
               ((v * v - 4.0 * u * m2) * (s - 4.0 * m2) - u * (s * v - 2.0 * m2 * (3.0 * v - 2.0 * u))) *
               log((v * (v - u + sqrt((v - u) * (v - u) - 4.0 * m2 * u)) - 2.0 * m2 * u) * 
               (v * (v - u + sqrt((v - u) * (v - u) - 4.0 * m2 * u)) - 2.0 * m2 * u) / 
               (4.0 * m2 * u * u) / (m2 + v)) / sqrt((v - u) * (v - u) - 4.0 * m2 * u)) /
               ((v - u) * (v - u) - 4.0 * m2 * u));
        aj42 = 2.0 * pi * log((v * sqrt(s - 4.0 * m2) + sqrt(4.0 * u * t * m2 + v * v * (s - 4.0 * m2))) * 
               (v * sqrt(s - 4.0 * m2) + sqrt(4.0 * u * t * m2 + v * v * (s - 4.0 * m2))) / 
               (4.0 * m2 * t * u)) / 
               sqrt((s - 4.0 * m2) * (4.0 * u * t * m2 + v * v * (s - 4.0 * m2)));
        aj43 = 2.0 * pi * log((sqrt(t * (t - 4.0 * m2)) + 2.0 * m2 - t) / (2.0 * m2)) / (u * sqrt(t * (t - 4.0 * m2)));
        aj44 = pi * (sqrt(t * (t - 4.0 * m2)) * (2.0 * m2 * v * v - t * u * (s + u)) / v +
               2.0 * m2 * (t * (v - u) - 2.0 * m2 * v) * log((sqrt(t * (t - 4.0 * m2)) + 2.0 * m2 - t) / (2.0 * m2)) / 
               (t * (t - 4.0 * m2)) / (u * u));
        aj45 = pi * log((v * (sqrt(als) + s) - 2.0 * m2 * (u + v)) * (v * (sqrt(als) + s) - 2.0 * m2 * (u + v)) / 
               (4.0 * m2 * (m2 * (u + v) * (u + v) - s * u * v)) / sqrt(als) / u);
        aj45 = pi * log((v * (sqrt(als) + s) - 2.0 * m2 * (u + v)) * (v * (sqrt(als) + s) - 2.0 * m2 * (u + v)) 
               / (4.0 * m2 * (m2 * (u + v) * (u + v) - s * u * v))) / sqrt(als) / u);
        aj46 = m2 * (-pi * v * (s * u * (t + v) - 2.0 * m2 * v * (u + v)) / u / (m2 * (u + v) * (u + v) - s * u * v) +
               (s * (u - v) + 2.0 * m2 * v) * aj45) / u / als;
        aj47 = m2 * (2.0 * aj45 * (s * s * (u - v) * (u - v) + 6.0 * m2 * m2 * v * v - 2.0 * m2 * s * (t * u - 2.0 * v * (v - v))) -
                pi * (v * (-(s * s * s * u * u * (s * s - 3.0 * u * u + 10.0 * u * v - 6.0 * v * v + 2.0 * s * (v - u))) +
                4.0 * m2 * m2 * (u + v) * (8.0 * s * s * u * u - (u + v) * (4.0 * s * u - 3.0 * v * v)) +
                2.0 * m2 * s * s * u * (2.0 * u * u * u - u * u * u - 13.0 * u * v * v - 7.0 * u * v * v - 6.0 * v * v + s 
                * s * u * (u + 5.0 * v) - s * v * (8.0 * u * v - 5.0 * u * u + v * v)) - 2.0 * m2 * m2 * s * (8.0 * s * u 
                * u * v * (s + 2.0 * v) - (u + v) * (8.0 * s * s * u * u + 2.0 * u * u * (s + u) + 6.0 * u * u * v 
                + 6.0 * s * u * v * v - 2.0 * u * u * v * v - 15.0 * u * v * v - 3.0 * v * v)))) / (u * (s * u * v 
                - m2 * (u + v) * (u + v)) * (s * u * v - m2 * (u + v) * (u + v))) / (2.0 * als * als * u * u);
        aj48 = pi * (t * v * (4.0 * m2 - t) * (t * (u - v) + 2.0 * m2 * v)
                - 2.0 * m2 * sqrt(t * (t - 4.0 * m2)) * (t * u * (s + u) - 2.0 * m2 * v * v)
                * log((2.0 * m2 - t + sqrt(t * (t - 4.0 * m2))) / (2.0 * m2))
                / (pow(t * (t - 4.0 * m2), 2) / pow(u, 3)));
        aj49 = pi * ((t * (t - 4.0 * m2) * (-2.0 * m2 * (t - 6.0 * m2) * pow(v, 4)
                + t * ((-64.0 * pow(m2, 3) - t * ((s + t) * (s + t) - (s + 3.0 * t) * u)
                + 2.0 * m2 * (2.0 * pow(s, 2) + 5.0 * t * (2.0 * u - v))) * (-pow(v, 2))
                + (t * u - 8.0 * pow(m2, 2) - 2.0 * m2 * (3.0 * s - 2.0 * t + u)) * pow(v, 3))
                - pow(t, 2) * (t - 4.0 * m2) * (2.0 * (2.0 * m2 - u) * pow(v, 2) + pow(u, 2) * (s + u + v))) / v
                - 4.0 * m2 * sqrt(t * (t - 4.0 * m2)) * (-3.0 * m2 * (t - 2.0 * m2) * pow(v, 3)
                - t * u * (2.0 * s * (t - m2) * v + (t - 4.0 * m2) * (t * u + 2.0 * m2 * v)))
                * log((2.0 * m2 - t + sqrt(t * (t - 4.0 * m2))) / (-2.0 * m2)) 
                / (pow(t, 3) * pow(t - 4.0 * m2, 3) * pow(u, 4)));
        aj50 = -pi / pow(u, 3) / als * (v * (s * (v - u) - 2.0 * m2 * v) + m2 * (s * u * (t + v) - 2.0 * m2 * v * (v + u)) 
                / sqrt(als) * log(pow(v * sqrt(als) + (s * v - 2.0 * m2 * (u + v)), 2) / (4.0 * m2 * (m2 * pow(u + v, 2) 
                - s * u * v)));
        aj51 = pi * ((v * (-(pow(s, 3) * u * pow(u - v, 2) * v) - 4.0 * pow(m2, 2) * s
                - ((2.0 * s - u) * pow(u, 2) * (s + u) - pow(u, 2) * (11.0 * s + 3.0 * u) * v
                + u * (u - 2.0 * s) * pow(v, 2) - 6.0 * u * pow(v, 3) + pow(v, 4))
                + m2 * pow(s, 2) * (pow(u, 2) * (pow(u, 2) + pow(s + u, 2) - 8.0 * (s + u) * v)
                - pow(v, 2) * (2.0 * pow(u, 2) + 4.0 * u * v + pow(v, 2)))
                + 4.0 * pow(m2, 3) * ((u + v) * (u + v) * (3.0 * pow(v, 2) - 4.0 * s * u)
                + 4.0 * s * u * (s * u - v * (u + v)))) / (pow(als, 2) * (m2 * pow(u + v, 2) - s * u * v))
                + (2.0 * (-3.0 * (pow(m2, 2) + pow(s - m2, 2)) * u * pow(v, 2) + 3.0 * m2 * (s - 2.0 * m2)
                - pow(v, 3) - s * (2.0 * m2 + s) * pow(u, 2) * (v - t) + 2.0 * s * u * v * (s * u + (s - 2.0 * m2) * (v - t)))
                * log(pow((sqrt(als) * v + s * v - 2.0 * m2 * (u + v)), 2) / (4.0 * m2 * (-(m2 * pow(u + v, 2) - s * u * v)))) 
                * m2) / pow(als, 2.5)) / pow(u, 4);
    }

    sr1  = 4.0 * (2.0 * (2.0 * (2.0 * t - 3.0 * v + 5.0 * s - 2.0 * m2) * m2 * pls - ((3.0 * (t - v) + 2.0 * s) 
         * pls * s + 1.0)) * aj7 - (6.0 * (2.0 * m2 - s) * m2 * pls + pls * pow(s, 2.0) + 1.0) * aj6 + 2.0 * (pls * 
         pow(s, 2.0) - 1.0) * aj5 + 4.0 * ((4.0 * (t - v + s) + pls * s * t * t) * m2 - (2.0 * (pls * pow(t, 2.0) + 2.0) 
         * pow(m2, 2.0) + (s + t - v) * (s + t - v)) * aj40 - (pls * pow(s, 2.0) - 1.0) * aj4 - 2.0 * (2.0 * (((t - v) 
         * (t - 2.0 * v) + 6.0 * pow(s, 2.0) + 2.0 * (3.0 * t - 2.0 * v) * s) * pls + 4.0 - 4.0 * (2.0 * s + t) * 
         m2 * pls) * m2 - ((2.0 * (pow(t - v, 2.0) + pow(s, 2.0)) + (3.0 * t - 2.0 * v) * s) * pls * s + 2.0 * 
         (t - v + s))) * aj39 - 2.0 * (2.0 * m2 - s) * aj37 * pls - 2.0 * (2.0 * (4.0 * (tt * v - 1.0 - 2.0 * s * tt) 
         - (2.0 * t - v) * pls * t + 8.0 * m2 * tt) * m2 + 4.0 * pow(s, 2.0) * tt - t - 4.0 * (tt * v - 1.0) * s 
         + pls * s * pow(t, 2.0)) * aj35 - (2.0 * (((2.0 * tt * v - 9.0 - 4.0 * s * tt) * s + (2.0 * tt * v - 3.0) * v) 
         * pls * s - (4.0 * (tt * v - 3.0 - s * tt) * s * tt + 6.0 * tt * v - 5.0) - 16.0 * (2.0 * pls * s - tt) 
         * pow(m2, 2.0) * tt - 2.0 * (((2.0 * tt * v - 5.0) * v + t + 4.0 * (tt * v - 2.0 - 3.0 * s * tt) * s) * pls 
         - 4.0 * (tt * v - 1.0 - 2.0 * s * tt) * tt) * m2) * m2 + 2.0 * (2.0 * tt * v - 3.0 - 2.0 * s * tt) * s - (2.0 
         * (tt * v - 2.0) * v + 3.0 * t) + (2.0 * s + t) * pls * pow(s, 2.0)) * aj34 - 4.0 * pow(2.0 * m2 - s, 2.0) 
         * aj29 + 4.0 * (2.0 * (4.0 * m2 - 3.0 * s) * m2 * pls + pls * pow(s, 2.0) - 1.0) * aj28 * s + 2.0 * ((4.0 
         * (((6.0 * pow(s, 2.0) - pow(v, 2.0)) * tt - 2.0 * (tt * v - 2.0) * s) * pls + 2.0 * (tt * v - 1.0 - 2.0 * s 
         * tt) * tt - 4.0 * (2.0 * pls * s - tt) * pow(m2, 2.0) * tt) * m2 - (2.0 * (2.0 * (tt * v - 3.0 - s * tt) * s 
         * tt + 3.0 * tt * v - 1.0) - (2.0 * (tt * v - 7.0 - 2.0 * s * tt) * s + (2.0 * tt * v + 3.0) * v) * pls * s)) 
         * aj23 + 2.0 * (8.0 * (tt * v - 1.0 - 2.0 * s * tt + 2.0 * pow(m2, 2.0) * tt) * m2 + 4.0 * pow(s, 2.0) * tt 
         - t - 4.0 * (tt * v - 1.0) * s) * aj22 + (2.0 * (4.0 * (((6.0 * pow(s, 2.0) - pow(v, 2.0)) * tt - 2.0 * (tt * v 
         - 2.0) * s) * pls + 2.0 * (tt * v - 1.0 - 2.0 * s * tt) * tt - 4.0 * (2.0 * pls * s - tt) * pow(m2, 2.0) * tt) 
         * m2 + ((2.0 * tt * v - 9.0 - 4.0 * s * tt) * s + (2.0 * tt * v + 1.0) * v) * pls * s - (4.0 * (tt * v - 3.0 
         - s * tt) * s * tt + 6.0 * tt * v - 1.0)) * m2 + 2.0 * (2.0 * tt * v - 1.0 - 2.0 * s * tt) * s - (2.0 
         * (tt * v - 1.0) * v + t) + (t - 2.0 * v + 2.0 * s) * pls * pow(s, 2.0)) * aj21 + 4.0 * (pls * s * t 
         - 1.0 - 2.0 * m2 * pls * t) * aj12 - 2.0 * (4.0 * pow(m2, 2.0) * pls - 2.0 * m2 * pls * s + 1.0) * aj11;

    sr2  = -4.0 * ((2.0 * (tt * v - 1.0 + t * vt - 2.0 * (tt - vt - 6.0 * s * tt * vt) * s - 2.0 * (t * vt + tt 
         * v - 4.0 * (tt - vt) * s) * pls * pow(s, 2.0) - 4.0 * (5.0 * (tt - vt) * pls * s + 6.0 * tt * vt) * m2 * s) 
         * m2 - (tt * v - 1.0 + t * vt - 2.0 * (tt - vt - 2.0 * s * tt * vt) * s - (tt * v - 1.0 + t * vt - 2.0 * (tt 
         - vt) * s) * pls * pow(s, 2.0)) * s + 16.0 * ((tt * v + 1.0 + t * vt + 2.0 * (tt - vt) * s) * pls + 2.0 * tt 
         * vt) * pow(m2, 3.0)) * aj36 + (2.0 * (((2.0 * tt - 11.0 * vv - 8.0 * s * tt * vv) * s - (3.0 * t * vv - 4.0 
         * tt * v - 2.0)) * pls * s - (tt + 6.0 * vv + 6.0 * s * tt * vv)) * m2 - (((tt - 4.0 * vv - 2.0 * s * tt * vv) 
         * s - 2.0 * (t * vv - 1.0)) * (pls * pow(s, 2.0) + 1.0) + 32.0 * (s * tt + 1.0) * pow(m2, 3.0) * pls * vv) 
         - 8.0 * ((2.0 * tt * v + 1.0 - 2.0 * t * vv - (5.0 * s * tt + 4.0) * s * vv) * pls - 2.0 * tt * vv) 
         * pow(m2, 2.0)) * aj6 + (8.0 * (2.0 * (tt + vv) - 5.0 * s * tt * vv + 4.0 * pow(m2, 2.0) * tt * vv) 
         * pow(m2, 2.0) * pls * s - ((2.0 * s * vv - 3.0) * s * tt + tt * v - 1.0) * (pls * pow(s, 2.0) - 1.0) 
         + 2.0 * ((8.0 * pow(s, 2.0) * vv - 7.0 * s + 2.0 * v) * pls * s - (2.0 * s * vv + 1.0)) * m2 * tt) * aj4 
         + (2.0 * (4.0 * (((2.0 * (t * vv + 1.0) - 5.0 * (vt - vv) * s) * s - (pow(t, 2.0) * vv - 2.0 * t + v)) 
         * pls - 2.0 * (tt + 2.0 * vv - 3.0 * (vt - vv) * s * tt) - 2.0 * (((vt - vv) * t + 1.0 - 2.0 * (vt - vv) 
         * s) * pls + 2.0 * (vt - vv) * tt) * m2) * m2 - (t * vt - 7.0 * t * vv + 4.0 * tt * v + 2.0 + 12.0 * (vt - vv) 
         * pow(s, 2.0) * tt + 2.0 * (vt - 9.0 * vv - 2.0 * tt) * s + (3.0 * t - 2.0 * v - 8.0 * (vt - vv) * s * s 
         - ((2.0 * vt - 7.0 * vv) * t - 2.0) * s) * pls * s)) * m2 + 2.0 * (2.0 * (vt - vv) * s * tt + vt - 5.0 * vv) 
         * pow(s, 2.0) - (3.0 * t * t * vv - 4.0 * t + 2.0 * v) + ((vt - 9.0 * vv) * t + 6.0) * s - (s * vt - s * vv 
         - t * vv) * (2.0 * s + t) * pls * pow(s, 2.0)) * aj34 - (((pls * pow(s, 2.0) + 1.0) * (s * tt + 1.0) - 16.0 
         * pow(m2, 3.0) * pls * tt) * vv + 8.0 * (tt + vv + 2.0 * s * tt * vv) * pow(m2, 2.0) * pls - 2.0 * ((2.0 * tt 
         + 3.0 * vv + 4.0 * s * tt * vv) * pls * s + tt * vv) * m2) * aj32 - (2.0 * (4.0 * (2.0 * (tt * v + 1.0 + t 
         * vt + 2.0 * (tt - vt) * s) * pls + 2.0 * tt * vt) * m2 + (2.0 * (tt * v - 3.0 - t * vt) * s - (5.0 
         * (tt - vt) * pow(s, 2.0) - tt * v * v)) * pls + 2.0 * (2.0 * tt + vt - 3.0 * s * tt * vt)) * m2 + 7.0 
         * tt * v - 1.0 + t * vt - 6.0 * (3.0 * tt + vt - 2.0 * s * tt * vt) * s - ((7.0 * tt * v - 18.0 
         - 4.0 * t * vt) * s - (8.0 * (tt - vt) * pow(s, 2.0) - 3.0 * v)) * pls * s) * m2 + (3.0 * tt * v - 4.0) 
         * v + 2.0 * t + 2.0 * (5.0 * tt + vt - 2.0 * s * tt * vt) * pow(s, 2.0) - (9.0 * tt * v - 5.0 + t * vt) 
         * s - ((tt * v - 4.0) * v + 2.0 * t + 2.0 * (tt - vt) * s * s - (3.0 * tt * v - 7.0 - t * vt) * s) 
         * pls * pow(s, 2.0)) * aj24 + ((2.0 * (vt + vv - 2.0 * (vt - vv) * s * tt) * s - ((vt - vv) * t + 2.0) 
         - ((vt + vv) * t - 2.0 * (vt - vv) * s) * pls * pow(s, 2.0)) * s + 16.0 * ((t * vt + 1.0 - 2.0 * (vt - vv) 
         * s) * pls + 2.0 * (vt - vv) * tt) * pow(m2, 3.0) - 8.0 * (((2.0 * vt + vv) * t + 2.0 - 5.0 * (vt - vv) 
         * s) * pls * s - 2.0 * (tt + vt - 3.0 * (vt - vv) * s * tt)) * pow(m2, 2.0) + 2.0 * ((vt - vv) * t 
         + 2.0 + 12.0 * (vt - vv) * pow(s, 2.0) * tt - 2.0 * (3.0 * vt + vv + 2.0 * tt) * s + ((4.0 * vt + vv) 
         * t + 2.0 - 8.0 * (vt - vv) * s) * pls * pow(s, 2.0)) * m2) * aj21 - ((s * vv - 1.0 - 2.0 * m2 * vv) 
         * (pls * pow(s, 2.0) - 1.0) + 8.0 * pow(m2, 2.0) * pls * s * vv) * aj20 * tt - 2.0 * (4.0 * (tt 
         - 2.0 * vt) * pow(m2, 2.0) * pls + tt - 2.0 * ((tt - 2.0 * vt) * pls * s + 2.0 * tt * vt) * m2) 
         * aj16 - 2.0 * (4.0 * (2.0 * vt + vv) * pow(m2, 2.0) * pls + vv - 2.0 * ((2.0 * vt + vv) * pls * s 
         - 2.0 * tt * vt) * m2) * aj15 - 2.0 * (2.0 * tt * v - 1.0 - (3.0 * tt + vt) * s + pls * pow(s, 2.0) 
         + 8.0 * (tt + vt) * pow(m2, 3.0) * pls + 4.0 * ((2.0 * tt * v - 1.0 - t * vt - 2.0 * s * tt) * pls 
         - 2.0 * tt * vt) * pow(m2, 2.0) - ((3.0 * tt * v - 1.0 - 2.0 * t * vt - 2.0 * s * tt) * pls * s 
         - 2.0 * (2.0 * s * vt + 1.0) * tt) * m2) * aj13 - 2.0 * (2.0 * (2.0 * (((vt - vv) * t + 3.0) * pls 
         + 2.0 * (vt - vv) * tt - 2.0 * (vt - vv) * m2 * pls) * m2 - (2.0 * ((vt - vv) * s * tt - 2.0 * vv) + (t * vt 
         + 2.0) * pls * s)) * m2 - (t * vv - 1.0 - (vt - 3.0 * vv) * s - pls * pow(s, 2.0))) * aj11 
         - 2.0 * (4.0 * (tt + vv + 2.0 * s * tt * vv - 2.0 * pow(m2, 2.0) * tt * vv) * pow(m2, 2.0) 
         * pls + (s * vv - 1.0) * tt + ((tt + 3.0 * vv + 2.0 * s * tt * vv) * pls * s + 2.0 * tt * vv) * m2) * aj1;

    sr3  = -4.0 * (2.0 * ((((2.0 * m2 * tt - s * vv) * m2 * pls + (pls * std::pow(s, 2) - 1.0) * vv) * aj2
         - 2.0 * (pls * s - tt - 2.0 * m2 * pls) * aj19 + (4.0 * std::pow(m2, 2) * pls - 2.0 * m2 * pls * s +
         1.0) * aj18 * vv) * tt - ((4.0 * tt * v - 1.0 - 2.0 * s * tt) * pls * s - 2.0 * (2.0 * tt * v -
         1.0 - 2.0 * s * tt) * tt + 2.0 * ((4.0 * s - 3.0 * v) * pls - 4.0 * tt) * m2 * tt) * aj17 + (
         2.0 * tt - vv - 2.0 * s * tt * vv + pls * std::pow(s, 2) * vv - (2.0 * (2.0 * tt - 3.0 * vv - 2.0 * s * tt
         * vv) * tt - (tt - 3.0 * vv - 2.0 * s * tt * vv) * pls * s) * m2 - 2.0 * ((tt * v - 3.0 - 4.0 * s
         * vv) * pls + 4.0 * tt * vv) * std::pow(m2, 2) * tt) * aj16 + (2.0 * ((6.0 * std::pow(s, 2) + std::pow(v, 2) - (tt * v +
         4.0) * s * v) * pls - 4.0 * (s - v) * tt - 2.0 * ((4.0 * s - tt * std::pow(v, 2)) * pls - 2.0 * tt) * m2
         ) * m2 * tt - (((2.0 * tt * v - 1.0) * v + 2.0 * std::pow(s, 2) * tt - (3.0 * tt * v - 1.0) * s) * pls * s
         + 2.0 * (2.0 * tt * v - 1.0 - s * tt) * s * tt - (2.0 * std::pow(tt, 2) * std::pow(v, 2) - 2.0 * tt * v + 1.0))) *
         aj14) - (2.0 * ((2.0 * (6.0 * tt - 7.0 * vv - 2.0 * s * tt * vv) * s - (tt * v + 3.0) * (tt *
         v - 1.0)) * pls * s - (4.0 * (2.0 * tt - 3.0 * vv - s * tt * vv) * s * tt - (3.0 * std::pow(tt, 2) * v -
         6.0 * tt + 2.0 * vv)) - 8.0 * ((4.0 * s * vv - tt * v) * pls - 2.0 * tt * vv) * std::pow(m2, 2) * tt +
         2.0 * ((((2.0 * tt * v - 3.0) * v + 12.0 * std::pow(s, 2) * vv) * tt - 2.0 * (std::pow(tt, 2) * v + 6.0 * tt -
         4.0 * vv) * s) * pls + 4.0 * (2.0 * tt - vv - 2.0 * s * tt * vv) * tt) * m2) * m2 - (3.0 * tt *
         v - 4.0 + 2.0 * t * vv - 2.0 * (3.0 * tt - 2.0 * vv - 2.0 * s * tt * vv) * s - (tt * v - 4.0 + 2.0
         * t * vv - 2.0 * (tt - 2.0 * vv) * s) * pls * std::pow(s, 2))) * aj13 - ((pls * std::pow(s, 2) - 1.0) * vv +
         2.0 * (pls * s - vv) * m2 * tt) * aj10 * tt + (4.0 * ((2.0 * (3.0 * std::pow(s, 2) * std::pow(vv, 2) + tt * v)
         - (tt + 8.0 * vv) * s) * pls - 4.0 * (s * vv - 1.0) * tt * vv) * std::pow(m2, 2) * tt - (8.0 * ((4.0 *
         s * std::pow(vv, 2) - tt) * pls - 2.0 * tt * std::pow(vv, 2)) * std::pow(m2, 3) * tt - (tt - vv - 2.0 * s * tt * vv) * (
         pls * std::pow(s, 2) - 1.0)) + 2.0 * (((9.0 * tt - vv - 2.0 * s * tt * vv) * s * vv - (tt * v + 3.0) *
         tt) * pls * s - (2.0 * (2.0 * tt - vv - s * tt * vv) * s * tt * vv - (3.0 * std::pow(tt, 2) - tt * vv + std::pow(vv, 2)
         ))) * m2) * aj1);

    sr4  = -4.0 * (((2.0 * (2.0 * t * vv - 3.0 + s * vv) * s + 2.0 * t * t * vv - 5.0 * t + 4.0 * v) *
         ((pls * s * s + 1.0) - 16.0 * (2.0 * t * vv - 3.0 + 6.0 * s * vv) * m2 * m2 * m2 * pls - 2.0 * (((
         11.0 * s * vv + 18.0 * t * vv - 21.0) * s + 9.0 * t * t * vv - 20.0 * t + 12.0 * v) * pls * s +
         9.0 * t * vv - 10.0 + 7.0 * s * vv) * m2 + 4.0 * (((22.0 * t * vv - 25.0 + 20.0 * s * vv) * s +
         4.0 * t * t * vv - 18.0 * t + 13.0 * v) * pls + 8.0 * vv) * m2 * m2) * aj6 - 2.0 * (((11.0 *
         t - 8.0 * v + 5.0 * s) * s + 2.0 * (4.0 * t * t - 6.0 * t * v + 3.0 * v * v)) * pls * s + 4.0 * (t
         - v + s) + 16.0 * (2.0 * t - v + 4.0 * s) * m2 * m2 * pls - 2.0 * (((26.0 * t - 17.0 * v + 18.0 *
         s) * s + 2.0 * (4.0 * t * t - 6.0 * t * v + 3.0 * v * v)) * pls + 8.0) * m2) * aj7 - 2.0 * ((
         2.0 * (t - v) - 3.0 * s * s * vv + 2.0 * (t * vv + 1.0) * s) * pls * s - (2.0 * (t * vv - 1.0) -
         s * vv) + 4.0 * (5.0 * (s * vv - 1.0) * pls * s + vv - 8.0 * m2 * pls * s * vv) * m2) * aj4 * m2
         - 2.0 * (2.0 * ((2.0 * (7.0 * t - 6.0 * v + 4.0 * s) * s * s + (3.0 * t * t - 3.0 * t * v + 2.0 *
         v * v) * (t - v) + (13.0 * t * t - 17.0 * t * v + 8.0 * v * v) * s) * pls + 8.0 * (t - v + s))
         * m2 - (8.0 * (((4.0 * (t - v) + 5.0 * s) * s + 2.0 * t * t - 2.0 * t * v + v * v) * pls + 3.0)
         * m2 * m2 + ((2.0 * s * s + 3.0 * s * t - 2.0 * s * v + 4.0 * t * t - 4.0 * t * v + 2.0 * v * v) *
         pls * s + 2.0 * (s + t - v)) * (s + t - v) - 32.0 * m2 * m2 * m2 * pls * s) * aj39 + 4.0 * (2.0 *
         m2 - s) * aj38 * pls - 4.0 * (((7.0 * t - 6.0 * v + 9.0 * s - 8.0 * m2) * m2 * pls - ((3.0 * (t -
         v) + 2.0 * s) * pls * s + 1.0))) * aj37 - (((2.0 * (2.0 * (t * vv - 1.0) + s * vv) * s + 3.0 *
         (t * t * vv - 2.0 * t + v)) * s + t * t * t * vv - 3.0 * t * t + 3.0 * t * v - v * v) * (pls * s *
         s + 1.0) + 64.0 * (tt + vv) * m2 * m2 * m2 * m2 * pls * s + 16.0 * (((2.0 * (tt * v + 1.0 - 2.0 *
         t * vv) - (5.0 * tt + 7.0 * vv) * s) * s - ((t * vv - 3.0) * t + (tt * v + 1.0) * v)) * pls - 3.0 *
         (tt + vv)) * m2 * m2 * m2 - 4.0 * ((6.0 * (tt * v + 2.0 - 3.0 * t * vv) * s * s - (8.0 * s * s *
         s * tt + 18.0 * s * s * s * vv + 2.0 * t * t * t * vv - 10.0 * t * t + 13.0 * t * v - 5.0 * v * v) -
         (3.0 * (2.0 * t * vv - 7.0) * t + 2.0 * (2.0 * tt * v + 3.0) * v) * s) * pls + 2.0 * (tt * v +
         3.0 - 7.0 * t * vv - (4.0 * tt + 7.0 * vv) * s)) * m2 * m2 - 2.0 * ((2.0 * t * t * t * vv - 10.0 *
         t * t + 11.0 * t * v - 3.0 * v * v + 2.0 * (tt + 5.0 * vv) * s * s * s - (2.0 * tt * v + 13.0 -
         15.0 * t * vv) * s * s + ((7.0 * t * vv - 18.0) * t + 2.0 * (tt * v + 3.0) * v) * s) * pls * s +
         7.0 * t * t * vv - 10.0 * t + 3.0 * v + 2.0 * (tt + 5.0 * vv) * s * s - (2.0 * (tt * v + 5.0) -
         15.0 * t * vv) * s) * m2) * aj34 - (8.0 * (t * vv - 4.0 + 2.0 * s * vv - 2.0 * m2 * vv) * m2 *
         m2 * pls + (t * vv - 2.0 + s * vv) * (pls * s * s + 1.0) - 2.0 * (3.0 * (t * vv - 2.0 + s * vv) *
         pls * s + 2.0 * vv) * m2) * aj32 + (8.0 * ((4.0 * s * s * s * tt + v * v - 3.0 * (tt * v - 1.0) *
         s * s + ((2.0 * tt * v - 1.0) * v - 6.0 * t) * s) * pls - (tt * v - 3.0 - 4.0 * s * tt)) * m2 * m2
         - (16.0 * (((5.0 * s * s + v * v) * tt - 2.0 * (tt * v + 1.0) * s) * pls + 3.0 * tt) * m2 * m2 *
         m2 - ((pls * s * s + 1.0) * (2.0 * s * s - 2.0 * s * v + v * v) + 64.0 * m2 * m2 * m2 * m2 * pls *
         s * tt)) + 2.0 * ((2.0 * (tt * v - 4.0 - s * tt) * s * s + 4.0 * t * t - 3.0 * t * v - 3.0 * v * v
         - 2.0 * ((tt * v - 3.0) * v - 2.0 * t) * s) * pls * s - (2.0 * s * s * tt - 5.0 * v - 2.0 * (tt *
         v - 4.0) * s)) * m2) * aj23 - (16.0 * (((((5.0 * tt + 7.0 * vv) * s - 2.0 * tt * v) * s + (tt *
         v + 1.0) * v) * pls + 3.0 * (tt + vv)) * m2 * m2 * m2 - ((pls * s * s + 1.0) * (2.0 * s * s * vv
         - 2.0 * s + v) + 64.0 * (tt + vv) * m2 * m2 * m2 * m2 * pls) * s - 2.0 * ((3.0 * t * t * vv - 5.0 *
         t - 2.0 * tt * v * v - 2.0 * (tt + 5.0 * vv) * s * s + (2.0 * tt * v + 5.0 + t * vv) * s) * pls *
         s * s - (t * t * vv + v + 2.0 * (tt + 5.0 * vv) * s * s - (2.0 * (tt * v + 1.0) + t * vv) * s)) *
         m2 + 8.0 * (((3.0 * tt * v + 1.0 + t * vv - (4.0 * tt + 9.0 * vv) * s) * s + (t * vv - 1.0) * t -
         (2.0 * tt * v + 1.0) * v) * pls * s - (t * vv - tt * v + (4.0 * tt + 7.0 * vv) * s)) * m2 * m2) *
         aj21 + 2.0 * ((3.0 * s * vv - 2.0) * pls * s - vv - 4.0 * m2 * pls * s * vv) * aj20 * m2 - (4.0 *
         (2.0 * m2 - s) * (t * vv - 2.0) * m2 * pls + (pls * s * s + 1.0) * (t * vv - 1.0)) * aj15 + (4.0 *
         ((2.0 * (4.0 * (t * vv - 1.0) + 3.0 * s * vv) * s + 2.0 * t * t * vv - 7.0 * t + 7.0 * v) * pls +
         6.0 * vv) * m2 * m2 + (t * t * vv - 3.0 * t + 2.0 * v + (3.0 * t * vv - 2.0) * s) * (pls * s * s +
         1.0) - 16.0 * (2.0 * s + t) * m2 * m2 * m2 * pls * vv - 2.0 * (((9.0 * t * vv - 8.0 + 2.0 * s * vv)
         * s + 4.0 * t * t * vv - 11.0 * t + 9.0 * v) * pls * s + 7.0 * t * vv - 8.0 + 2.0 * s * vv) * m2) *
         aj11 - (8.0 * (t * vv - 3.0 + 3.0 * s * vv - 2.0 * m2 * vv) * m2 * m2 * pls + (t * vv - 2.0 + s *
         vv) * (pls * s * s + 1.0) - 2.0 * ((5.0 * t * vv - 7.0 + 4.0 * s * vv) * pls * s + 4.0 * vv) * m2)
         * aj1) * uu;

    sr5  = -4.0 * (((16.0 * (((2.0 * (t * uu + 3.0) - 7.0 * s * uu) * s - 2.0 * (t * t * uu + t + v))
         * pls - 3.0 * uu) * m2 * m2 * m2 - ((2.0 * (t * uu + 1.0 - s * uu) * s - (t * t * uu + t + v)) *
         (pls * s * s + 1.0) - 64.0 * m2 * m2 * m2 * m2 * pls * uu) * s - 8.0 * (((5.0 * t * uu + 11.0 -
         9.0 * s * uu) * s - (3.0 * t * t * uu + 5.0 * t + 4.0 * v)) * pls * s - 7.0 * (s * uu - 1.0)) *
         m2 * m2 + 2.0 * ((5.0 * t * uu + 12.0 - 10.0 * s * uu) * s - (2.0 * t * t * uu + 2.0 * t + v) -
         (2.0 * (5.0 * s * uu - 4.0 * t * uu - 6.0) * s + 4.0 * t * t * uu + 6.0 * t + 5.0 * v) * pls * s *
         s) * m2 - 4.0 * (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) *
         (2.0 * m2 - s + v) * dz1 * m2) * ds + 2.0 * (2.0 * ((t * uu + 3.0) * t + 3.0 * s * s * uu - (3.0 *
         t * uu + 2.0) * s) * pls * s - 3.0 * (t * uu + 1.0 - 2.0 * s * uu)) * m2 + (2.0 * (t * uu + 1.0 -
         s * uu) * s - (t * t * uu + t + v)) * (pls * s * s + 1.0) - 8.0 * (2.0 * (s - t) * s * uu + t * t *
         uu + t + v) * m2 * m2 * pls + ((pls * s * s + 1.0) * (2.0 * s * s - 2.0 * s * v + v * v) - 32.0 *
         m2 * m2 * m2 * pls * s - 2.0 * (2.0 * (4.0 * s * s - 3.0 * s * v + v * v) * pls * s + 8.0 * s - 3.0 *
         v) * m2 + 8.0 * ((5.0 * s * s - 2.0 * s * v + v * v) * pls + 3.0) * m2 * m2) * dz1 + 4.0 * (8.0 *
         m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) * (2.0 * m2 - s) * ds *
         ds * m2) * aj45 + (((pls * s * s + 1.0) * (s + t) - 48.0 * m2 * m2 * m2 * pls + 4.0 * (7.0 * t -
         4.0 * v + 21.0 * s) * m2 * m2 * pls - 2.0 * ((11.0 * t - 4.0 * v + 11.0 * s) * pls * s + 4.0)) *
         dz1 - 2.0 * (2.0 * (8.0 * m2 - 5.0 * s) * m2 * pls + pls * s * s + 1.0)) * aj8 - ((pls * s * s +
         1.0) * (s + t) - 48.0 * m2 * m2 * m2 * pls + 4.0 * (7.0 * t - 4.0 * v + 13.0 * s) * m2 * m2 * pls -
         2.0 * ((7.0 * t - 4.0 * v + 7.0 * s) * pls * s + 4.0)) * aj6 * dz1 - 16.0 * aj5 * m2 * pls * s + 4.0 *
         (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) * (2.0 * m2 - s) *
         aj46 * ds - (2.0 * (4.0 * ((3.0 * (3.0 * s * vt + 2.0) * s + t * t * vt + 3.0 * t + v) * pls * s +
         t * vt + 1.0 + 7.0 * s * vt) * m2 * m2 + ((pls * s * s + 1.0) * s * s + 32.0 * m2 * m2 * m2 * m2 *
         pls) * s * vt - 8.0 * (((7.0 * s * vt + 4.0) * s + t * t * vt + t + v) * pls + 3.0 * vt) * m2 * m2 *
         m2 + ((t * vt + 1.0 - 10.0 * s * vt) * s - (t * t * vt + t + v) - 2.0 * (5.0 * s * s * vt + 2.0 * s +
         2.0 * t) * pls * s * s) * m2) * dz1 + 16.0 * (6.0 * t * vt - 1.0 + 6.0 * s * vt - 8.0 * m2 * vt) *
         m2 * m2 * m2 * pls - ((2.0 * s * s + t * t) * uu - (2.0 * uu - vt) * s * t) * (pls * s * s + 1.0) +
         4.0 * ((t - 2.0 * v - 2.0 * (uu + vt) * t * t + ((4.0 * uu - 15.0 * vt) * t - 4.0 * (uu + 2.0 * vt) *
         s) * s) * pls - 2.0 * vt) * m2 * m2 - 2.0 * ((2.0 * (3.0 * (uu - vt) * t - 1.0 - (3.0 * uu + vt) * s) *
         s + 2.0 * t - v - 2.0 * (uu + vt) * t * t) * pls * s + (3.0 * uu + vt) * t - 1.0 - 6.0 * s * uu) * m2
         + (((pls * s * s + 1.0) * (2.0 * s * s - 2.0 * s * t + t * t) + 64.0 * m2 * m2 * m2 * m2 * pls) * s -
         16.0 * ((7.0 * s * s - 2.0 * s * t + 2.0 * t * t) * pls + 3.0) * m2 * m2 * m2 + 8.0 * ((9.0 * s * s -
         5.0 * s * t + 3.0 * t * t) * pls + 7.0) * m2 * m2 * s - 2.0 * (10.0 * s * s - 5.0 * s * t + 2.0 * t *
         t + 2.0 * (5.0 * s * s - 4.0 * s * t + 2.0 * t * t) * pls * s * s) * m2) * (uu - vt) * ds) * aj43 - ((
         (t - v + 3.0 * s) * (t - v) - 2.0 * (s * vt - 2.0) * s * s) * (pls * s * s + 1.0) - 64.0 * m2 * m2 *
         m2 * m2 * pls * s * vt + 16.0 * (((7.0 * s * vt - 8.0) * s + t * t * vt - 2.0 * t + 2.0 * v) * pls +
         3.0 * vt) * m2 * m2 * m2 + 2.0 * (((10.0 * (s * vt - 2.0) * s - (14.0 * t - 13.0 * v)) * s - (6.0 *
         t * t - 9.0 * t * v + 4.0 * v * v)) * pls * s - ((t * vt + 20.0 - 10.0 * s * vt) * s - (t * t * vt -
         8.0 * t + 9.0 * v))) * m2 - 4.0 * ((2.0 * (9.0 * s * vt - 16.0) * s * s - (5.0 * t * t - 7.0 * t * v
         + 4.0 * v * v) + (2.0 * t * t * vt - 17.0 * t + 16.0 * v) * s) * pls + 2.0 * (t * vt - 10.0 + 7.0 * s *
         vt)) * m2 * m2) * dz1 - (2.0 * (((6.0 * t * vt - 1.0 + 2.0 * s * vt) * s + 2.0 * t * t * vt - 6.0 * t +
         9.0 * v) * pls * s - t * vt) * m2 - (((t * vt - 1.0) * s - (t - 2.0 * v)) * (pls * s * s + 1.0) +
         128.0 * m2 * m2 * m2 * m2 * pls * vt - 16.0 * (6.0 * t * vt + 5.0 + 6.0 * s * vt) * m2 * m2 * m2 *
         pls) - 4.0 * (((15.0 * t * vt + 7.0 + 8.0 * s * vt) * s + 2.0 * t * t * vt - 3.0 * t + 7.0 * v) * pls
         + 2.0 * vt) * m2 * m2) + 4.0 * (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls *
         s * s * s + s) * (2.0 * m2 - s) * ds * ds * m2 + (16.0 * (((2.0 * (t * vt + 3.0) - 7.0 * s * vt) * s -
         2.0 * (t * t * vt + t + v)) * pls - 3.0 * vt) * m2 * m2 * m2 - ((2.0 * (t * vt + 1.0 - s * vt) * s -
         (t * t * vt + t + v)) * (pls * s * s + 1.0) - 64.0 * m2 * m2 * m2 * m2 * pls * vt) * s - 8.0 * (((
         5.0 * t * vt + 11.0 - 9.0 * s * vt) * s - (3.0 * t * t * vt + 5.0 * t + 4.0 * v)) * pls * s - 7.0 *
         (s * vt - 1.0)) * m2 * m2 + 2.0 * ((5.0 * t * vt + 12.0 - 10.0 * s * vt) * s - (2.0 * t * t * vt +
         2.0 * t + v) - (2.0 * (5.0 * s * vt - 4.0 * t * vt - 6.0) * s + 4.0 * t * t * vt + 6.0 * t + 5.0 * v) *
         pls * s * s) * m2 - 4.0 * (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s *
         s + s) * (2.0 * m2 - s + v) * dz1 * m2) * ds) * aj42 - 4.0 * (2.0 * m2 - s) * aj41 * dz1 * m2 * pls +
         8.0 * (s + t - 4.0 * m2) * aj4 * dz1 * m2 * pls * s + 2.0 * (4.0 * ((3.0 * (3.0 * s * vt + 2.0) * s +
         t * t * vt + 3.0 * t + v) * pls * s + t * vt + 1.0 + 7.0 * s * vt) * m2 * m2 + ((pls * s * s + 1.0) *
         s * s + 32.0 * m2 * m2 * m2 * m2 * pls) * s * vt - 8.0 * (((7.0 * s * vt + 4.0) * s + t * t * vt + t +
         v) * pls + 3.0 * vt) * m2 * m2 * m2 + ((t * vt + 1.0 - 10.0 * s * vt) * s - (t * t * vt + t + v) - 2.0 *
         (5.0 * s * s * vt + 2.0 * s + 2.0 * t) * pls * s * s) * m2) * aj36 * dz1 + (((t - v + 3.0 * s) * (t -
         v) - 2.0 * (s * vt - 2.0) * s * s) * (pls * s * s + 1.0) + 64.0 * (tt - vt) * m2 * m2 * m2 * m2 * pls *
         s + 16.0 * (((2.0 * (tt * v - 4.0) - (5.0 * tt - 7.0 * vt) * s) * s + t * t * vt - 2.0 * t + 2.0 * v) *
         pls - 3.0 * (tt - vt)) * m2 * m2 * m2 + 2.0 * (((2.0 * (tt * v - 10.0 - (tt - 5.0 * vt) * s) * s - (
         14.0 * t - 13.0 * v)) * s - (6.0 * t * t - 9.0 * t * v + 4.0 * v * v)) * pls * s + t * t * vt - 8.0 * t
         + 9.0 * v - 2.0 * (tt - 5.0 * vt) * s * s + (2.0 * (tt * v - 10.0) - t * vt) * s) * m2 - 4.0 * ((2.0 *
         (3.0 * tt * v - 16.0 - (4.0 * tt - 9.0 * vt) * s) * s * s - (5.0 * t * t - 7.0 * t * v + 4.0 * v * v) +
         (2.0 * t * t * vt - 17.0 * t + 16.0 * v) * s) * pls + 2.0 * (3.0 * tt * v - 10.0 + t * vt - (4.0 * tt -
         7.0 * vt) * s)) * m2 * m2) * aj34 * dz1 + 4.0 * (2.0 * m2 - s) * aj32 * dz1 * m2 * pls - 4.0 * (8.0 *
         m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) * (2.0 * m2 - s) * aj30 *
         ds + 4.0 * (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) * (2.0 *
         m2 - s) * aj28 * ds - 4.0 * (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s *
         s + s) * (2.0 * m2 - s) * aj26 * ds - ((16.0 * (((2.0 * (t * uu + 3.0) - 7.0 * s * uu) * s - 2.0 * (t *
         t * uu + t + v)) * pls - 3.0 * uu) * m2 * m2 * m2 - ((2.0 * (t * uu + 1.0 - s * uu) * s - (t * t * uu + t
         + v)) * (pls * s * s + 1.0) - 64.0 * m2 * m2 * m2 * m2 * pls * uu) * s - 8.0 * (((5.0 * t * uu + 11.0 -
         9.0 * s * uu) * s - (3.0 * t * t * uu + 5.0 * t + 4.0 * v)) * pls * s - 7.0 * (s * uu - 1.0)) * m2 * m2 +
         2.0 * ((5.0 * t * uu + 12.0 - 10.0 * s * uu) * s - (2.0 * t * t * uu + 2.0 * t + v) - (2.0 * (5.0 * s *
         uu - 4.0 * t * uu - 6.0) * s + 4.0 * t * t * uu + 6.0 * t + 5.0 * v) * pls * s * s) * m2 - 4.0 * (8.0 *
         m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) * (2.0 * m2 - s + v) * dz1 *
         m2) * ds + 2.0 * (2.0 * ((t * uu + 3.0) * t + 3.0 * s * s * uu - (3.0 * t * uu + 2.0) * s) * pls * s -
         3.0 * (t * uu + 1.0 - 2.0 * s * uu)) * m2 + (2.0 * (t * uu + 1.0 - s * uu) * s - (t * t * uu + t + v)) *
         (pls * s * s + 1.0) - 8.0 * (2.0 * (s - t) * s * uu + t * t * uu + t + v) * m2 * m2 * pls + ((pls * s *
         s + 1.0) * (2.0 * s * s - 2.0 * s * v + v * v) - 32.0 * m2 * m2 * m2 * pls * s - 2.0 * (2.0 * (4.0 * s *
         s - 3.0 * s * v + v * v) * pls * s + 8.0 * s - 3.0 * v) * m2 + 8.0 * ((5.0 * s * s - 2.0 * s * v + v * v)
         * pls + 3.0) * m2 * m2) * dz1 + 4.0 * (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls *
         s * s * s + s) * (2.0 * m2 - s) * ds * ds * m2) * aj25 + (((pls * s * s + 1.0) * (2.0 * s * s - 2.0 * s *
         t + t * t) + 64.0 * m2 * m2 * m2 * m2 * pls) * s - 16.0 * ((7.0 * s * s - 2.0 * s * t + 2.0 * t * t) *
         pls + 3.0) * m2 * m2 * m2 + 8.0 * ((9.0 * s * s - 5.0 * s * t + 3.0 * t * t) * pls + 7.0) * m2 * m2 * s -
         2.0 * (10.0 * s * s - 5.0 * s * t + 2.0 * t * t + 2.0 * (5.0 * s * s - 4.0 * s * t + 2.0 * t * t) * pls *
         s * s) * m2) * (uu - vt) * ds - (2.0 * ((2.0 * (3.0 * t * uu + 1.0 - 3.0 * s * uu) * s - (2.0 * t * t *
         uu + 6.0 * t - v)) * pls * s + 3.0 * (t * uu + 1.0 - 2.0 * s * uu)) * m2 - ((2.0 * (t * uu + 1.0 - s * uu)
         * s - (t * t * uu + t + v)) * (pls * s * s + 1.0) - 8.0 * (2.0 * (s - t) * s * uu + t * t * uu + t + v) *
         m2 * m2 * pls)) * aj24 - (8.0 * ((4.0 * s * s * s * tt + v * v - 3.0 * (tt * v - 1.0) * s * s + ((2.0 *
         tt * v - 1.0) * v - 6.0 * t) * s) * pls - (3.0 * (tt * v - 1.0) - 4.0 * s * tt)) * m2 * m2 - (16.0 * ((
         2.0 * (tt * v + 1.0) - 5.0 * s * tt) * pls * s - 3.0 * tt) * m2 * m2 * m2 + (pls * s * s + 1.0) * (2.0 *
         s * s - 2.0 * s * v + v * v) + 64.0 * m2 * m2 * m2 * m2 * pls * s * tt + 2.0 * (((2.0 * (tt * v - 4.0 -
         s * tt) * s + 4.0 * t + 5.0 * v) * s + 4.0 * t * t - t * v - 2.0 * v * v) * pls * s - (2.0 * s * s * tt -
         3.0 * v - 2.0 * (tt * v - 4.0) * s)) * m2) * aj23 * dz1 + ((16.0 * (((2.0 * (t * vt + 3.0) - 7.0 * s *
         vt) * s - 2.0 * (t * t * vt + t + v)) * pls - 3.0 * vt) * m2 * m2 * m2 - ((2.0 * (t * vt + 1.0 - s * vt) *
         s - (t * t * vt + t + v)) * (pls * s * s + 1.0) - 64.0 * m2 * m2 * m2 * m2 * pls * vt) * s - 8.0 * (((
         5.0 * t * vt + 11.0 - 9.0 * s * vt) * s - (3.0 * t * t * vt + 5.0 * t + 4.0 * v)) * pls * s - 7.0 * (s *
         vt - 1.0)) * m2 * m2 + 2.0 * ((5.0 * t * vt + 12.0 - 10.0 * s * vt) * s - (2.0 * t * t * vt + 2.0 * t + v)
         - (2.0 * (5.0 * s * vt - 4.0 * t * vt - 6.0) * s + 4.0 * t * t * vt + 6.0 * t + 5.0 * v) * pls * s * s) * m2
         - 4.0 * (8.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) * (2.0 * m2 - s +
         v) * dz1 * m2) * ds - 4.0 * ((2.0 * m2 - s + v) * dz1 * tt - (2.0 * m2 - s) * ds * ds) * (8.0 * m2 * m2 *
         pls * s - 6.0 * m2 * pls * s * s - 6.0 * m2 + pls * s * s * s + s) * m2) * aj21 + 4.0 * (2.0 * m2 - s) *
         aj16 * m2 * pls * vt - 4.0 * (2.0 * m2 - s) * aj15 * m2 * pls * vt - (32.0 * m2 * m2 * m2 * pls * vt - 8.0 *
         m2 * m2 * pls * s * vt - 16.0 * m2 * m2 * pls * t * vt + 4.0 * m2 * m2 * pls + 8.0 * m2 * pls * s * t * vt -
         2.0 * m2 * pls * s - 4.0 * m2 * vt + pls * s * s + 1.0) * aj13 + (32.0 * m2 * m2 * m2 * pls * vt - 8.0 * m2 *
         m2 * pls * s * vt - 16.0 * m2 * m2 * pls * t * vt - 20.0 * m2 * m2 * pls + 8.0 * m2 * pls * s * t * vt + 6.0 *
         m2 * pls * s - 4.0 * m2 * vt - pls * s * s - 1.0) * aj11;

    sr6  = 4.0 * (((((2.0 * (s * tt + 2.0) * s * vv + 3.0 * (t * vv - 1.0)) * s + t * t * vv - 2.0 * t + v) *
         (pls * s * s + 1.0) + 64.0 * m2 * m2 * m2 * m2 * pls * s * tt * vv - 16.0 * ((tt * v - 1.0 + 2.0 * t * vv +
         (7.0 * s * vv + 2.0) * s * tt) * pls + 3.0 * tt * vv) * m2 * m2 * m2 + 4.0 * (((2.0 * (3.0 * tt + 4.0 * vv
         + 9.0 * s * tt * vv) * s - (2.0 * tt * v - 1.0)) * s + (2.0 * t * vv - 1.0) * t + (2.0 * tt * v - 1.0) * v)
         * pls + 2.0 * (2.0 * tt + 3.0 * vv + 7.0 * s * tt * vv)) * m2 * m2 + 2.0 * (((2.0 * tt * v + 3.0 - 3.0 * t *
         vv - 2.0 * (tt + 6.0 * vv + 5.0 * s * tt * vv) * s) * s + (t * vv - 1.0) * t - (2.0 * tt * v - 1.0) * v) *
         pls * s + tt * v + 3.0 - 5.0 * t * vv - (3.0 * tt + 13.0 * vv + 10.0 * s * tt * vv) * s) * m2) * aj6 - 2.0 *
         (((pls * s * s + 1.0) * s * s + 32.0 * m2 * m2 * m2 * m2 * pls) * s * tt - 8.0 * ((7.0 * s * s * tt - 4.0 *
         s + tt * v * v) * pls + 3.0 * tt) * m2 * m2 * m2 - (2.0 * ((5.0 * s * tt - 2.0) * s - 2.0 * (t - v)) * pls *
         s * s + (10.0 * s * s + s * v + v * v) * tt) * m2 + 4.0 * (((tt * v + 2.0) * v - 2.0 * t + 3.0 * (3.0 * s *
         tt - 2.0) * s) * pls * s + (7.0 * s - v) * tt) * m2 * m2) * aj36 + (((pls * s * s + 1.0) * (s * tt + 1.0) -
         32.0 * m2 * m2 * m2 * pls * tt) * vv + 8.0 * (tt + vv + 3.0 * s * tt * vv) * m2 * m2 * pls - 4.0 * ((tt +
         vv + 2.0 * s * tt * vv) * pls * s + tt * vv) * m2) * aj33 + (((tt - 4.0 * vv - 2.0 * s * tt * vv) * s - 2.0 *
         (t * vv - 1.0)) * (pls * s * s + 1.0) - 16.0 * (tt - 4.0 * vv - 2.0 * s * tt * vv) * m2 * m2 * m2 * pls +
         4.0 * (4.0 * tt * v - 1.0 - 4.0 * t * vv + (tt - 18.0 * vv - 10.0 * s * tt * vv) * s) * m2 * m2 * pls - 2.0 *
         (((2.0 * tt - 15.0 * vv - 8.0 * s * tt * vv) * s - (5.0 * t * vv - 4.0 * tt * v - 2.0)) * pls * s - (2.0 *
         tt + 7.0 * vv + 5.0 * s * tt * vv)) * m2) * aj32 + (8.0 * (4.0 * s * tt + 1.0 - 4.0 * m2 * tt) * m2 * m2 *
         pls + (pls * s * s + 1.0) * (s * tt + 1.0) - 2.0 * (5.0 * pls * s * s + 1.0) * m2 * tt) * aj2 * vv + (4.0 *
         (2.0 * m2 - s) * (2.0 * tt - vv) * m2 * pls + (pls * s * s + 1.0) * (tt - vv)) * aj18 - (4.0 * ((2.0 * (
         8.0 * tt - vv - 3.0 * s * tt * vv) * s - (7.0 * tt * v - 3.0)) * pls - 6.0 * tt * vv) * m2 * m2 - ((2.0 *
         (tt * v - 1.0) - (3.0 * tt - 2.0 * vv) * s) * (pls * s * s + 1.0) - 32.0 * (s * vv - 1.0) * m2 * m2 * m2 *
         pls * tt) + 2.0 * ((9.0 * tt * v - 4.0 - t * vv - (13.0 * tt - 5.0 * vv - 2.0 * s * tt * vv) * s) * pls * s
         - (3.0 * (3.0 * tt - vv) - 2.0 * s * tt * vv)) * m2) * aj16 + (((tt * v - 1.0) * v + 2.0 * (2.0 * tt - vv) *
         s * s - (3.0 * tt * v - 2.0) * s) * (pls * s * s + 1.0) - 64.0 * m2 * m2 * m2 * m2 * pls * s * tt * vv + 16.0 *
         (((5.0 * s * s * vv + 2.0 * v) * tt - 2.0 * (4.0 * tt + vv) * s) * pls + 3.0 * tt * vv) * m2 * m2 * m2 + 4.0 *
         ((2.0 * (16.0 * tt - 3.0 * vv - 4.0 * s * tt * vv) * s * s + (5.0 * tt * v - 3.0) * v - (17.0 * tt * v + 6.0 -
         12.0 * t * vv) * s) * pls + 2.0 * (10.0 * tt - 3.0 * vv - 4.0 * s * tt * vv)) * m2 * m2 + 2.0 * (((14.0 * tt *
         v - 3.0 - 4.0 * t * vv + 2.0 * (s * tt * vv - 10.0 * tt + 4.0 * vv) * s) * s - ((4.0 * t * vv - 7.0) * t + (6.0 *
         tt * v - 1.0) * v)) * pls * s + 2.0 * (s * tt * vv - 10.0 * tt + 4.0 * vv) * s + 8.0 * tt * v - 5.0) * m2) *
         aj13 + (((pls * s * s + 1.0) * (s * tt + 1.0) - 32.0 * m2 * m2 * m2 * pls * tt) * vv + 8.0 * (tt + vv + 5.0 *
         s * tt * vv) * m2 * m2 * pls - 4.0 * (3.0 * (s * tt + 1.0) * pls * s + tt) * m2 * vv) * aj10 - (16.0 * ((
         2.0 * (tt - vv) + 5.0 * s * s * tt * vv * vv - 2.0 * (3.0 * tt + 2.0 * vv) * s * vv) * pls + 3.0 * tt * vv *
         vv) * m2 * m2 * m2 - (((tt - 3.0 * vv - 2.0 * s * tt * vv) * s - (t * vv - 2.0)) * (pls * s * s + 1.0) + 64.0 *
         m2 * m2 * m2 * m2 * pls * s * tt * vv * vv) + 2.0 * (((2.0 * (s * vv - 6.0) * s * tt * vv + 8.0 * tt - 5.0 *
         vv) * s + t * vv - 3.0 * tt * v) * pls * s + 2.0 * (s * vv - 6.0) * s * tt * vv + 6.0 * tt - 5.0 * vv) * m2 +
         4.0 * ((2.0 * (11.0 * tt + 2.0 * vv - 4.0 * s * tt * vv) * s * s * vv + 2.0 * t * vv - 1.0 - (13.0 * tt + 4.0 *
         vv - 4.0 * t * vv * vv) * s) * pls - 2.0 * (4.0 * s * vv - 7.0) * tt * vv) * m2 * m2) * aj1) * uu;

    sr7  = -4.0 * ((2.0 * ((2.0 * uu - 3.0 * vv + 7.0 * tt - 12.0 * s * tt * vv) * s - (2.0 * (tt * v - 1.0) +
         t * vv) - ((2.0 * uu - 5.0 * vv - 5.0 * tt + 2.0 * (uu + 4.0 * vv) * s * tt) * s + 4.0 * tt * v + 3.0 + (2.0 *
         uu - 3.0 * vv) * t) * pls * s * s) * m2 + ((pls * s * s + 1.0) * (2.0 * s * s * vv - 2.0 * s + v) * s - 256.0 *
         m2 * m2 * m2 * m2 * m2 * pls * uu) * tt + 64.0 * (tt + 2.0 * uu + (5.0 * uu - vv) * s * tt) * m2 * m2 * m2 *
         m2 * pls - 16.0 * (((4.0 * (2.0 * uu - vv) + 5.0 * tt + (10.0 * uu - vv) * s * tt) * s + 2.0 * tt * v + 1.0 +
         t * uu) * pls + (uu + 3.0 * vv) * tt) * m2 * m2 * m2 + 8.0 * (((5.0 * uu - 7.0 * vv + tt + (5.0 * uu + 4.0 *
         vv) * s * tt) * s + 3.0 * (tt * v + 1.0) + (2.0 * uu - vv) * t) * pls * s - (uu - vv + 2.0 * tt - (uu + 9.0 *
         vv) * s * tt)) * m2 * m2) * aj45 - (16.0 * (6.0 * tt + 5.0 * vv + 13.0 * s * tt * vv - 12.0 * m2 * tt * vv) *
         m2 * m2 * m2 * pls + ((s * vv - 2.0) * s * tt - (t * vv + 1.0)) * (pls * s * s + 1.0) - 4.0 * (2.0 * tt * v +
         3.0 + 2.0 * t * vv + (16.0 * s * tt * vv + 15.0 * tt + 18.0 * vv) * s) * m2 * m2 * pls + 2.0 * ((2.0 * tt * v +
         3.0 + 4.0 * t * vv + (8.0 * tt + 7.0 * vv + s * tt * vv) * s) * pls * s + tt - vv - 3.0 * s * tt * vv) * m2) *
         aj8 - 2.0 * (2.0 * (((2.0 * tt * v - 17.0 - 2.0 * s * tt) * s * s - (3.0 * t * t - 3.0 * t * v + 2.0 * v * v) -
         ((2.0 * tt * v - 13.0) * v + 20.0 * t) * s + 16.0 * (tt * v - 5.0 - 3.0 * s * tt + 4.0 * m2 * tt) * m2 * m2) *
         pls - 4.0 * (((3.0 * tt * v - 16.0 - 4.0 * s * tt) * s - ((tt * v - 4.0) * v + 7.0 * t)) * pls - tt)) * m2 +
         ((5.0 * t - 3.0 * v + 3.0 * s) * s + 2.0 * (2.0 * t * t - 2.0 * t * v + v * v)) * pls * s - 2.0 * t) * aj44 - (
         2.0 * ((2.0 * t * t * vv + 4.0 * t - 3.0 * v - 2.0 * (uu - vv) * s * s * s * tt - (2.0 * uu - 7.0 * vv + 2.0 *
         tt) * s * s - ((2.0 * uu - 5.0 * vv) * t - 2.0) * s) * pls * s + (2.0 * uu + vv + tt) * s - (t * vv - 4.0)) * m2
         - (16.0 * ((tt * v - 6.0 + (uu - 5.0 * vv) * t - (tt - 8.0 * uu + 14.0 * vv - 10.0 * (uu - vv) * s * tt) * s) *
         pls + (uu - vv) * tt) * m2 * m2 * m2 - (64.0 * (2.0 * (uu - 2.0 * vv - tt) + 5.0 * (uu - vv) * s * tt) * m2 *
         m2 * m2 * m2 * pls - ((pls * s * s + 1.0) * (s * s + t * t) * vv + 256.0 * (uu - vv) * m2 * m2 * m2 * m2 * m2 *
         pls * tt))) - 4.0 * ((2.0 * t * t * vv + 4.0 * t - v - 10.0 * (uu - vv) * s * s * s * tt - 5.0 * (2.0 * (uu -
         2.0 * vv) + tt) * s * s - (2.0 * tt * v - 11.0 + 4.0 * (uu - 3.0 * vv) * t) * s) * pls + 2.0 * (uu - vv + tt -
         (uu - vv) * s * tt)) * m2 * m2) * aj43 + (((pls * s * s + 1.0) * (s * tt + 1.0) - 32.0 * m2 * m2 * m2 * pls *
         tt) * vv + 8.0 * (tt + vv + 3.0 * s * tt * vv) * m2 * m2 * pls - 4.0 * ((tt + vv + 2.0 * s * tt * vv) * pls * s
         + tt * vv) * m2) * aj41 + 4.0 * ((s * s * tt * vv + s * vv + 1.0) * pls * s - (s * tt + 1.0) * vv - 2.0 * (3.0 *
         s * vv - 1.0 - 4.0 * m2 * vv) * m2 * pls * s * tt) * aj4 * m2 - 2.0 * ((2.0 * uu - 3.0 * vv + 7.0 * tt - 12.0 *
         s * tt * vv) * s - (2.0 * (tt * v - 1.0) + t * vv) - ((2.0 * uu - 5.0 * vv - 5.0 * tt + 2.0 * (uu + 4.0 * vv) *
         s * tt) * s + 4.0 * tt * v + 3.0 + (2.0 * uu - 3.0 * vv) * t) * pls * s * s) * m2 + ((pls * s * s + 1.0) * (2.0 *
         s * s * vv - 2.0 * s + v) * s - 256.0 * m2 * m2 * m2 * m2 * m2 * pls * uu) * tt + 64.0 * (tt + 2.0 * uu + (5.0 *
         uu - vv) * s * tt) * m2 * m2 * m2 * m2 * pls - 16.0 * (((4.0 * (2.0 * uu - vv) + 5.0 * tt + (10.0 * uu - vv) *
         s * tt) * s + 2.0 * tt * v + 1.0 + t * uu) * pls + (uu + 3.0 * vv) * tt) * m2 * m2 * m2 + 8.0 * (((5.0 * uu -
         7.0 * vv + tt + (5.0 * uu + 4.0 * vv) * s * tt) * s + 3.0 * (tt * v + 1.0) + (2.0 * uu - vv) * t) * pls * s -
         (uu - vv + 2.0 * tt - (uu + 9.0 * vv) * s * tt)) * m2 * m2) * aj25 + (((pls * s * s + 1.0) * (2.0 * s * s - 2.0 *
         s * v + v * v) - 256.0 * m2 * m2 * m2 * m2 * m2 * pls * uu) * tt + 64.0 * (tt + 2.0 * uu + 5.0 * s * tt * uu) *
         m2 * m2 * m2 * m2 * pls + 8.0 * ((5.0 * s * s * s * tt * uu + 2.0 * s * t * uu + tt * v * v + (6.0 * tt + 5.0 *
         uu) * s * s) * pls + 2.0 * tt - uu + s * tt * uu) * m2 * m2 - 16.0 * ((tt * v + 1.0 + t * uu + (5.0 * tt + 8.0 *
         uu + 10.0 * s * tt * uu) * s) * pls + tt * uu) * m2 * m2 * m2 - 2.0 * (((2.0 * tt * v + 5.0) * v - 4.0 * t + 2.0 *
         (4.0 * tt + uu + s * tt * uu) * s * s - (3.0 * tt * v + 2.0 - 2.0 * t * uu) * s) * pls * s + 2.0 * (3.0 * tt -
         uu) * s - 5.0 * tt * v) * m2) * aj24 - 2.0 * ((3.0 * s * vv - 2.0) * pls * s - vv - 4.0 * m2 * pls * s * vv) *
         aj20 * m2 * tt + 4.0 * (2.0 * m2 - s) * aj19 * pls * tt - 4.0 * ((2.0 * tt * v + 1.0 - s * tt) * pls * s - tt -
         (3.0 * (tt * v + 1.0) - 2.0 * s * tt - 4.0 * m2 * tt) * m2 * pls) * aj17 + (((pls * s * s + 1.0) * (tt + vv) -
         16.0 * m2 * m2 * m2 * pls * tt * vv + 2.0 * (tt - 3.0 * vv - 2.0 * s * tt * vv) * m2 * pls * s + 8.0 * (tt + vv +
         2.0 * s * tt * vv) * m2 * m2 * pls) * aj16 + 2.0 * (((3.0 * tt * v - 1.0 - 2.0 * s * tt) * s - 2.0 * ((2.0 * tt *
         v - 1.0) * v + 2.0 * t)) * pls * s + 2.0 * (tt * v + 1.0 - s * tt) - 8.0 * (3.0 * s * tt + 4.0 - 4.0 * m2 * tt) *
         m2 * m2 * pls - 2.0 * (((4.0 * tt * v - 5.0 - 6.0 * s * tt) * s - 3.0 * (t + tt * v * v)) * pls - 4.0 * tt) * m2)
         * aj14 - (2.0 * (8.0 * (3.0 * (tt + vv + s * tt * vv) - 4.0 * m2 * tt * vv) * m2 * m2 * pls + tt - 3.0 * vv - 2.0 *
         s * tt * vv + (s + 6.0 * t) * pls * s * vv) * m2 - (2.0 * tt * v + 1.0 + t * vv - (2.0 * tt + vv) * s) * (pls *
         s * s + 1.0) - 4.0 * ((2.0 * tt * v + 3.0 + 2.0 * t * vv + (5.0 * (tt + 2.0 * vv) + 2.0 * s * tt * vv) * s) *
         pls - 2.0 * tt * vv) * m2 * m2) * aj13 + (((pls * s * s + 1.0) * (s * tt + 1.0) - 32.0 * m2 * m2 * m2 * pls *
         tt) * vv + 8.0 * (tt + vv + 3.0 * s * tt * vv) * m2 * m2 * pls - 2.0 * ((5.0 * tt + 4.0 * vv + 2.0 * s * tt * vv)
         * pls * s + tt * vv) * m2) * aj1;

    sr8  = -4.0 * ((2.0 * (2.0 * (4.0 * ((2.0 * (7.0 * t * vv - 4.0 + 10.0 * s * vv) * s + 2.0 * t * t * vv - 7.0 * t 
         + 5.0 * v) * pls + 6.0 * vv - 16.0 * m2 * pls * s * vv) * m2 - (((4.0 * (9.0 * s * vv + 12.0 * t * vv - 8.0) * s 
         + 22.0 * t * t * vv - 37.0 * t + 16.0 * v) * s + 2.0 * t * t * t * vv - 9.0 * t * t + 12.0 * t * v - 4.0 * v * v) 
         * pls + 4.0 * (5.0 * t * vv - 2.0 + 5.0 * s * vv))) * m2 + (((14.0 * s * vv + 25.0 * t * vv - 20.0) * s + 16.0 
         * t * t * vv - 29.0 * t + 12.0 * v) * s + 5.0 * t * t * t * vv - 13.0 * t * t + 12.0 * t * v - 3.0 * v * v) 
         * pls * s + 2.0 * (5.0 * s * vv + 9.0 * t * vv - 6.0) * s + 12.0 * t * t * vv - 14.0 * t + 5.0 * v) * m2 
         - ((2.0 * (2.0 * (t * vv - 1.0) + s * vv) * s + 4.0 * t * t * vv - 6.0 * t + 3.0 * v) * s + 2.0 * t * t 
         * t * vv - 4.0 * t * t + 3.0 * t * v - v * v - ((2.0 * t - v) * (t - v) - 2.0 * s * s * s * vv - 4.0 
         * (t * vv - 1.0) * s * s - (2.0 * t * t * vv - 6.0 * t + 3.0 * v) * s) * pls * s * s)) * aj6 
         - 2.0 * (2.0 * (((2.0 * t - v) * (t - v) * (t - v) + 6.0 * s * s * s + 2.0 * (3.0 * t - 4.0 * v) 
         * s * s + (6.0 * t * t - 11.0 * t * v + 6.0 * v * v) * s) * pls + 4.0 * (t - v + s)) * m2 
         - (4.0 * ((4.0 * (2.0 * (t - v) + 3.0 * s) * s + 4.0 * t * t - 6.0 * t * v + 3.0 * v * v) * pls + 6.0) 
         * m2 * m2 + ((s * s - s * v + t * t - 2.0 * t * v + v * v) * pls * s + s + t - v) * (s + t - v) 
         - 64.0 * m2 * m2 * m2 * pls * s)) * aj7 - 4.0 * (2.0 * (t - v + 3.0 * s - 4.0 * m2) * m2 * pls - ((t - v + s) 
         * pls * s + 1.0)) * aj38 + 2.0 * (2.0 * ((4.0 * t - 3.0 * v) * (t - v) + 4.0 * s * s + (12.0 * t - 11.0 * v) * s) 
         * m2 * pls - (8.0 * (4.0 * t - 3.0 * v + 2.0 * s) * m2 * m2 * pls + (pls * s * s + 3.0 * pls * s * t 
         - 3.0 * pls * s * v + 2.0) * (s + t - v))) * aj37 + 2.0 * (4.0 * m2 * m2 * pls - 2.0 * m2 * pls * s 
         + 1.0) * (4.0 * m2 * vv - s * vv - t * vv + 1.0) * aj33 + 2.0 * (((2.0 * (t * vv - 1.0) + s * vv) * s 
         + t * t * vv - 2.0 * t + v) * (pls * s * s + 1.0) - 8.0 * (4.0 * t * vv - 7.0 + 2.0 * s * vv) * m2 * m2 * m2 
         * pls + 2.0 * ((20.0 * t * vv - 21.0 + 12.0 * s * vv) * s + 4.0 * t * t * vv - 11.0 * t + 6.0 * v) * m2 * m2 
         * pls - (((16.0 * t * vv - 15.0 + 9.0 * s * vv) * s + 7.0 * t * t * vv - 13.0 * t + 6.0 * v) * pls * s 
         + 2.0 * (3.0 * (t * vv - 1.0) + s * vv)) * m2) * aj32 - ((t * vv - 1.0 + s * vv) * (pls * s * s + 1.0) 
         - 32.0 * m2 * m2 * m2 * pls * vv + 4.0 * (2.0 * t * vv - 3.0 + 8.0 * s * vv) * m2 * m2 * pls 
         - 2.0 * ((3.0 * t * vv - 4.0 + 5.0 * s * vv) * pls * s + 3.0 * vv) * m2) * aj2 - ((t * vv - 1.0 + s * vv) 
         * (pls * s * s + 1.0) - 32.0 * m2 * m2 * m2 * pls * vv + 8.0 * (t * vv - 1.0 + 3.0 * s * vv) * m2 * m2 * pls 
         - 2.0 * ((2.0 * t * vv - 1.0 + 4.0 * s * vv) * pls * s + vv) * m2) * aj10 + (((3.0 * t * vv 
         - 4.0 + 2.0 * s * vv) * s + t * t * vv - 3.0 * t + 2.0 * v) * (pls * s * s + 1.0) - 128.0 * m2 * m2 * m2 
         * m2 * pls * s * vv * vv + 8.0 * ((4.0 * (t * vv - 6.0 + 3.0 * s * vv) * s * vv - (4.0 * t * vv - 11.0)) 
         * pls + 6.0 * vv * vv) * m2 * m2 * m2 + 2.0 * (((s * vv - 14.0) * s * s * vv - (5.0 * t * t * vv 
         - 12.0 * t + 8.0 * v) - (t * t * vv * vv + 15.0 * t * vv - 21.0) * s) * pls * s + (s * vv - 10.0) * s * vv 
         + t * t * vv * t * t - 9.0 * t * vv + 11.0) * m2 - 4.0 * ((2.0 * ((t * vv - 16.0 + 3.0 * s * vv) * s * vv 
         - (11.0 * t * vv - 16.0)) * s - (2.0 * t * t * vv - 7.0 * t + 6.0 * v)) * pls + 4.0 * (t * vv - 4.0 + s * vv) 
         * vv) * m2 * m2) * aj1) * uu * uu;

    sr9  = 4.0 * (2.0 * (((s * vv - 2.0) * s - (t * t * vv - v) - (t - 2.0 * v + 3.0 * s) * pls * s * s 
         + 8.0 * (5.0 * t * vv + 4.0 + 11.0 * s * vv - 12.0 * m2 * vv) * m2 * m2 * m2 * pls + ((9.0 * t 
         - 13.0 * v + 2.0 * s * s * vv + (2.0 * t * vv + 21.0) * s) * pls * s + 4.0 * (t * vv + 2.0)) * m2 
         - 4.0 * (((5.0 * t * vv + 12.0 + 6.0 * s * vv) * s + t * t * vv + 3.0 * t - 5.0 * v) * pls 
         + 2.0 * vv) * m2 * m2 + (4.0 * (((17.0 * t - 12.0 * v + 20.0 * s) * s + 3.0 * t * t - 5.0 * t * v + v * v) * pls 
         + 6.0) * m2 * m2 + (t - v + s) * t + (2.0 * s - v) * (s + t) * pls * s * s - 8.0 * (6.0 * t - 5.0 * v + 14.0 * s) 
         * m2 * m2 * m2 * pls - (((20.0 * t - 13.0 * v + 21.0 * s) * s + 7.0 * t * t - 13.0 * t * v + 2.0 * v * v) 
         * pls * s + 4.0 * (3.0 * t - v + 2.0 * s)) * m2) * dz1) * aj8 - (s * s + s * t - t * v + (s + t - v) * (s + t) 
         * pls * s * s - 8.0 * (6.0 * t - 5.0 * v + 6.0 * s) * m2 * m2 * m2 * pls + 4.0 * ((13.0 * t - 9.0 * v 
         + 10.0 * s) * s + 3.0 * t * t - 5.0 * t * v + v * v) * m2 * m2 * pls - (((11.0 * s + 18.0 * t - 10.0 * v) * s 
         + 7.0 * t * t - 10.0 * t * v + 2.0 * v * v) * pls * s + 2.0 * (t - 2.0 * v + 3.0 * s)) * m2) * aj6 * dz1 
         - (((3.0 * s * vv - 1.0) * s + t * t * vv - t + v + ((s * vv - 3.0) * s - (t * t * vv - t - v)) 
         * pls * s * s) * s - 16.0 * ((6.0 * s * s * vv - 2.0 * s * t * vv + v) * pls + 4.0 * vv) * m2 * m2 * m2 
         + 8.0 * (((t * vv - 8.0 + 8.0 * s * vv) * s - (t * t * vv - 2.0 * t - v)) * pls * s + 2.0 * t * vv - 1.0 
         + 8.0 * s * vv) * m2 * m2 - 2.0 * ((4.0 * t * vv - 3.0 + 11.0 * s * vv) * s + t * t * vv - t + v 
         + ((2.0 * (t * vv - 7.0) + 7.0 * s * vv) * s - (t * t * vv - 5.0 * v)) * pls * s * s) * m2 
         + (8.0 * ((5.0 * s - 2.0 * v) * pls * s * v - 2.0 * (3.0 * s - v) - 2.0 * ((2.0 * s - v) * pls * v 
         - 2.0) * m2) * m2 * m2 - (4.0 * s * s - 2.0 * s * v + v * v - (2.0 * s - v) * pls * s * s * v) * s 
         + 2.0 * (12.0 * s * s - 6.0 * s * v + v * v - 4.0 * (2.0 * s - v) * pls * s * s * v) * m2) * dz1) * aj45 
         + (16.0 * (((11.0 * s * vv + 16.0 * t * vv - 8.0) * s + 5.0 * t * t * vv - 2.0 * t - 2.0 * v) * pls 
         + 4.0 * vv - 8.0 * (3.0 * s + 2.0 * t - 2.0 * m2) * m2 * pls * vv) * m2 * m2 * m2 - (((t * vv - 2.0) * s + t) 
         * s + 3.0 * t * t * t * vv - 5.0 * t * t + 3.0 * t * v - v * v + (((t * vv - 2.0) * s + t) 
         * s - (t * t * vv - 3.0 * t * t + 3.0 * t * v - v * v)) * pls * s * s) + 2.0 * (((3.0 * t - v) 
         * (t - v) + s * s * vv + (6.0 * t * vv - 11.0) * s * s + (t * t * vv + 2.0 * t + v) * s) * pls 
         * s + (4.0 * t * vv - 7.0 + s * vv) * s + 11.0 * t * t * vv - 9.0 * t + 5.0 * v) * m2 - 8.0 * (((t * vv - 1.0) 
         * t * t + 4.0 * s * s * s * vv + (5.0 * t * vv - 2.0) * s * t + (10.0 * t * vv - 11.0) * s * s) * pls 
         + 8.0 * t * vv - 5.0 + 2.0 * s * vv) * m2 * m2 - (8.0 * (3.0 * (2.0 * t - v + 6.0 * s) * pls * s * s 
         + 2.0 * (2.0 * t - v + 5.0 * s) + 32.0 * m2 * m2 * pls * s) * m2 * m2 + (2.0 * t * t - 2.0 * t * v + v * v 
         + 2.0 * s * s - (2.0 * t * t - 2.0 * t * v + v * v - 2.0 * s * s) * pls * s * s) * s 
         - 16.0 * ((20.0 * s * s + v * v + 2.0 * (2.0 * t - v) * s) * pls + 6.0) * m2 * m2 * m2 
         - 2.0 * (2.0 * (2.0 * (2.0 * t - v) + 5.0 * s) * s + 2.0 * t * t - 2.0 * t * v + v * v 
         + 2.0 * ((2.0 * t - v + 7.0 * s) * s - (t * t - t * v + v * v)) * pls * s * s) * m2) * dz1) * aj43 
         - 2.0 * (2.0 * (2.0 * (t * vv - 2.0 + 3.0 * s * vv - 4.0 * m2 * vv) * m2 * pls - ((t * vv - 2.0 + s * vv) 
         * pls * s + 2.0 * vv)) * m2 + t * vv - 1.0 + s * vv - pls * s * s + (((s + t) * s * s - 48.0 * m2 * m2 * m2 
         + 4.0 * (3.0 * t - 2.0 * v + 9.0 * s) * m2 * m2) * pls - 2.0 * ((3.0 * t - 2.0 * v + 5.0 * s) * pls * s 
         + 2.0) * m2) * dz1) * aj41 + ((pls * s * s - 1.0) * (s * s * vv + s - t) - 64.0 * m2 * m2 * m2 * pls * s * vv 
         + 8.0 * ((5.0 * s + 2.0 * t) * pls * s + 1.0) * m2 * m2 * vv - 2.0 * ((5.0 * s * s * vv + 3.0 * s - v) 
         * pls * s + s * vv + 1.0) * m2 + 2.0 * (16.0 * m2 * m2 * pls * s - 6.0 * m2 * pls * s * s + 3.0 * m2 
         * pls * s * v - 6.0 * m2 + pls * s * s * s - pls * s * s * t - s + t) * (4.0 * m2 - s - t) * dz1) * aj4 
         + (8.0 * (3.0 * (2.0 * t - v + 6.0 * s) * pls * s * s + 2.0 * (2.0 * t - v + 5.0 * s) + 32.0 * m2 * m2 
         * pls * s) * m2 * m2 + (2.0 * t * t - 2.0 * t * v + v * v + 2.0 * s * s - (2.0 * t * t - 2.0 * t * v + v * v 
         - 2.0 * s * s) * pls * s * s) * s - 16.0 * ((20.0 * s * s + v * v + 2.0 * (2.0 * t - v) * s) * pls 
         + 6.0) * m2 * m2 * m2 - 2.0 * (2.0 * (2.0 * (2.0 * t - v) + 5.0 * s) * s + 2.0 * t * t - 2.0 * t * v + v * v 
         + 2.0 * ((2.0 * t - v + 7.0 * s) * s - (t * t - t * v + v * v)) * pls * s * s) * m2) * aj36 * dz1 
         - 2.0 * (2.0 * ((3.0 * t - 2.0 * v + 3.0 * s) * pls * s + 4.0) * m2 + 48.0 * m2 * m2 * m2 * pls 
         - 36.0 * m2 * m2 * pls * s - 12.0 * m2 * m2 * pls * t + 8.0 * m2 * m2 * pls * v - s - t) * aj32 * dz1 
         + (((3.0 * s * vv - 1.0) * s + t * t * vv - t + v + ((s * vv - 3.0) * s - (t * t * vv - t - v)) 
         * pls * s * s) * s - 16.0 * ((6.0 * s * s * vv - 2.0 * s * t * vv + v) * pls + 4.0 * vv) * m2 * m2 * m2 
         + 8.0 * (((t * vv - 8.0 + 8.0 * s * vv) * s - (t * t * vv - 2.0 * t - v)) * pls * s + 2.0 * t * vv - 1.0 
         + 8.0 * s * vv) * m2 * m2 - 2.0 * ((4.0 * t * vv - 3.0 + 11.0 * s * vv) * s + t * t * vv - t + v 
         + ((2.0 * (t * vv - 7.0) + 7.0 * s * vv) * s - (t * t * vv - 5.0 * v)) * pls * s * s) * m2 
         + (8.0 * ((5.0 * s - 2.0 * v) * pls * s * v - 2.0 * (3.0 * s - v) - 2.0 * ((2.0 * s - v) * pls * v - 2.0) * m2) 
         * m2 * m2 - (4.0 * s * s - 2.0 * s * v + v * v - (2.0 * s - v) * pls * s * s * v) * s + 2.0 * (12.0 * s * s 
         - 6.0 * s * v + v * v - 4.0 * (2.0 * s - v) * pls * s * s * v) * m2) * dz1) * aj25 + ((2.0 * (s - v) * s 
         + 2.0 * t * t + v * v + (2.0 * (s - v) * s - (2.0 * t * t - v * v)) * pls * s * s) * t 
         + 256.0 * m2 * m2 * m2 * pls * s - 16.0 * ((12.0 * s * s + v * v + 4.0 * (3.0 * t - v) * s) * pls + 6.0) 
         * m2 * m2 * m2 + 8.0 * ((6.0 * s * s * s + t * v * v + 2.0 * (2.0 * t - v) * (t - v) * s 
         + 2.0 * (7.0 * t - 2.0 * v) * s * s) * pls + 2.0 * (5.0 * t - v + 2.0 * s)) * m2 * m2 
         + 2.0 * ((3.0 * (t - v) * t * v - 2.0 * s * s * s - 2.0 * (6.0 * t - v) * s * s - (2.0 * t * t 
         - 9.0 * t * v + v * v) * s) * pls * s - (2.0 * (4.0 * t - v + s) * s + 10.0 * t * t - 4.0 * t * v + v * v)) * m2) 
         * aj23 * dz1 - ((s * vv + 1.0 - 2.0 * m2 * vv) * (pls * s * s - 1.0) + 8.0 * m2 * m2 * pls * s * vv 
         - 2.0 * (s + t - 4.0 * m2) * (pls * s * s - 1.0) * dz1) * aj20 - ((pls * s * s + 1.0) * (t * vv - 1.0) 
         - 16.0 * m2 * m2 * m2 * pls * vv + 8.0 * (t * vv - 2.0 + 2.0 * s * vv) * m2 * m2 * pls 
         - 2.0 * ((3.0 * t * vv - 5.0 + s * vv) * pls * s + vv) * m2) * aj16 + (16.0 * (3.0 * t * vv + 1.0 + s * vv 
         - 4.0 * m2 * vv) * m2 * m2 * m2 * pls - (t * t * vv - v - (t * vv - 2.0) * s) * (pls * s * s + 1.0) 
         + 8.0 * (((s * vv - 6.0) * s - (t * t * vv - 3.0 * v)) * pls + vv) * m2 * m2 + 2.0 * ((4.0 * (t - 2.0 * v) 
         - s * s * vv - 3.0 * (t * vv - 3.0) * s) * pls * s - (t * vv - 6.0 + s * vv)) * m2) * aj13 
         - ((t * vv - 2.0 + s * vv) * (pls * s * s + 1.0) - 32.0 * m2 * m2 * m2 * pls * vv 
         + 8.0 * (t * vv - 2.0 + 2.0 * s * vv) * m2 * m2 * pls - 2.0 * ((2.0 * t * vv - 3.0 + 5.0 * s * vv) 
         * pls * s + 2.0 * vv) * m2) * aj1) * uu;

    sr10 = 4.0 * (2.0 * (2.0 * (2.0 * (t - v + 3.0 * s - 4.0 * m2) * m2 * pls - ((t - v + s) * pls * s + 1.0)) 
         * aj9 - (4.0 * m2 * m2 * pls - 2.0 * m2 * pls * s + 1.0) * aj8 - 2.0 * (2.0 * m2 - s) * (2.0 * m2 - s) 
         * aj51 + (4.0 * (2.0 * t * uu - 1.0 + 2.0 * s * uu - 4.0 * m2 * uu) * m2 - ((4.0 * t * uu - 1.0) * s - (t - v))) 
         * aj50 + (pls * s * s - 1.0) * aj5 + 2.0 * (16.0 * (t - v + 2.0 * s) * m2 * m2 * m2 * pls - (32.0 * m2 * m2 * m2 
         * pls + t * t) + ((s + t - v) * (s + t - v) * pls * s + 4.0 * t) * m2 - 2.0 * ((5.0 * s + t - v) * (s + t - v) 
         * pls + 2.0) * m2 * m2) * aj49 - (2.0 * (4.0 * ((4.0 * t - 3.0 * v + 6.0 * s) * pls - 2.0 * uu - 8.0 * m2 
         * pls) * m2 + 2.0 * (2.0 * t * uu - 1.0 + 2.0 * s * uu) - (6.0 * s + 2.0 * t - v) * (s + t - v) * pls) * m2 
         - ((4.0 * t * uu - 1.0) * s - (t - v) - (s + t - v) * (s + t - v) * pls * s)) * aj48 + 4.0 * (2.0 * m2 - s) 
         * (2.0 * m2 - s) * aj47 - (4.0 * (2.0 * t * uu - 1.0 + 2.0 * s * uu + 3.0 * pls * s * s - 4.0 
         * (pls * s + uu) * m2) * m2 - ((4.0 * t * uu - 3.0) * s - (t - v) + 2.0 * pls * s * s * s)) * aj46 
         - (2.0 * (2.0 * (2.0 * t * uu - 1.0) * s * uu - (2.0 * t * uu - 3.0) - (6.0 * s * t * uu - 3.0 * s 
         + 2.0 * t * t * uu + 2.0 * t + 3.0 * v) * pls * s + 64.0 * m2 * m2 * m2 * pls * uu - 16.0 * ((2.0 * t 
         * uu + 1.0 + 4.0 * s * uu) * pls - uu * uu) * m2 * m2) * m2 + (2.0 * s * uu - 1.0) * s + 2.0 * t * t 
         * uu + t + v + (t + v - s) * pls * s * s + 8.0 * (((8.0 * t * uu + 1.0 + 3.0 * s * uu) * s + t * t 
         * uu + t + v) * pls - 2.0 * (t * uu + 1.0 + s * uu) * uu) * m2 * m2) * aj45 + 2.0 * ((s * s 
         + 2.0 * t * t + (t + v) * s) * pls * s - 2.0 * t + 2.0 * (4.0 * (t + 2.0 * v + s) * m2 - (3.0 * s + t) 
         * (s + t + v)) * m2 * pls) * aj44 + (2.0 * (2.0 * (2.0 * t * uu - 1.0) * s * uu - (2.0 * t * uu - 3.0) 
         - (6.0 * s * t * uu - 3.0 * s + 2.0 * t * t * uu + 2.0 * t - v) * pls * s + 64.0 * m2 * m2 * m2 * pls * uu 
         - 8.0 * ((4.0 * t * uu + 1.0 + 8.0 * s * uu) * pls - 2.0 * uu * uu) * m2 * m2) * m2 + (2.0 * s * uu + 1.0) * s
         + 2.0 * t * t * uu - t + v + (t - v - s) * pls * s * s + 4.0 * (((16.0 * t * uu + 1.0 + 6.0 * s * uu) * s 
         + 2.0 * t * t * uu + t - 2.0 * v) * pls - 4.0 * (t * uu + 1.0 + s * uu) * uu) * m2 * m2) 
         * aj43 - (pls * s * s - 1.0) * aj4 - 4.0 * (2.0 * m2 - s) * (2.0 * m2 - s) * aj31 - 4.0 * (2.0 * (4.0 * m2 
         - 3.0 * s) * m2 * pls + pls * s * s - 1.0) * aj30 * s - 8.0 * (2.0 * m2 - s) * (2.0 * m2 - s) * aj27 
         + 2.0 * (4.0 * (2.0 * t * uu - 1.0 + 2.0 * s * uu + 3.0 * pls * s * s - 4.0 * (pls * s + uu) * m2) * m2 
         - ((4.0 * t * uu - 3.0) * s - (t - v) + 2.0 * pls * s * s * s)) * aj26 + (2.0 * (2.0 * (2.0 * t * uu - 1.0) 
         * s * uu - (2.0 * t * uu - 3.0) - (6.0 * s * t * uu - 3.0 * s + 2.0 * t * t * uu + 2.0 * t + 3.0 * v) * pls 
         * s + 64.0 * m2 * m2 * m2 * pls * uu - 16.0 * ((2.0 * t * uu + 1.0 + 4.0 * s * uu) * pls - uu * uu) * m2 * m2) 
         * m2 + (2.0 * s * uu - 1.0) * s + 2.0 * t * t * uu + t + v + (t + v - s) * pls * s * s 
         + 8.0 * (((8.0 * t * uu + 1.0 + 3.0 * s * uu) * s + t * t * uu + t + v) * pls - 2.0 * (t * uu + 1.0 + s * uu) 
         * uu) * m2 * m2) * aj25 + 2.0 * (16.0 * ((2.0 * t * uu + 1.0 + 4.0 * s * uu) * pls - uu * uu - 4.0 * m2 * pls 
         * uu) * m2 * m2 * m2 - (t * t * uu + v + s * s * uu - (s - t) * pls * s * s) + ((2.0 * t * t * uu + 2.0 * t 
         + 5.0 * v + 6.0 * (t * uu - 1.0) * s) * pls * s - 2.0 * ((2.0 * t * uu - 1.0) * s * uu - (t * uu - 2.0)) * m2 
         - 4.0 * (((8.0 * t * uu + 1.0 + 3.0 * s * uu) * s + t * t * uu + t + v) * pls - 2.0 * (t * uu + 1.0 + s * uu) 
         * uu) * m2 * m2) * aj24 - 2.0 * (2.0 * m2 - s) * aj17 * pls + 2.0 * ((2.0 * t + v) * pls * s - 1.0 - 2.0 
         * (t + 2.0 * v - 2.0 * m2) * m2 * pls) * aj14 - (4.0 * (2.0 * m2 - s) * m2 * pls + pls * s * s + 1.0) * aj13);

    // fsir = aj6
    // std::cout << "fsir: " << fsir << std::endl;
    if (fsir < 0) nn += 1;
    if (ikey == -1) {
        fir = -4.0 * ((aj[4] + aj[6] + m2 * aj[0] * vv * vv + aj[13]) + (t - 2.0 * m2) * (aj[22] + aj[12] * vv) + (s - 2.0 * m2) * (aj[37] + vv * aj[3]) - (s + t - 2.0 * m2) * (vv * aj[6] + aj[24]));
        fsir = fsir - alfa / (4.0 * pi * pi) * fir * sig(t, pl, 0);
    }
    //TODO:/FIXME: DO WE WANT THIS TO RETURN FSIR???
    return fsir;
}
