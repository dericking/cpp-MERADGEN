#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

// Constants and includes
const double pi = 4.0 * atan(1.0);
const double alfa = 0.729735e-2;
const double m = 0.511000e-3;
const double m2 = 0.261112e-6;
const double barn = 0.389379e6;

// Function prototypes
double sig(double t, double pl, int i);
double vacpol(double t);
double L1f(double xt, double xs, double xm2);
double xsBt(double pl, double xs, double xt, double xu);
double dgg1(double xs, double xt, double xu);
double dgg2(double xs, double xt, double xu);
double dcanc(double xvmin, double xs, double xt, double xu);
double fspen(double x);
void simps(double a1, double b1, double h1, double reps1, double aeps1, double (*funct)(double), double &ai, double &aih, double &aiabs);
void simpsx(double a, double b, int np, double ep, double (*func)(double), double &res);
void vectrec(std::vector<double> &vpgen);
void merad_init(double elab);
void grid_init();
double fsirv(double v);
double fsirv1(double v);
void zd(double t, double t1, double v);
double urand();

void meradgen(double ppl, std::vector<double> &vpgen, std::vector<double> &r) {
    double y, ppl, u0, u, vmax, vmin;
    double fsirv, sig, xsvt, xsbt, dcanc;
    double distsiv[61], distarv[61];
    double distsit1[121], distart1[121];
    double distsiz[61], distarz[61];
    double t1p[5];
    double vvn, vvo, sin, sio;
    double tt1n, tt1o, sirad1;
    double zmin, zmax, zzn, zzo, det;
    double sitot, xs0, xsvr, xsV, xsB, xsF, xsR;
    double xsadd, sirand, urand;
    double fs, v, t1min, t1max;
    double t1z, t1z1, t1z2, stept1, t1p;
    double sigmar;
    int nv = 60, nt1 = 30, nz = 60;
    int nn, iv, ikey = 0;
    int nt1, it1, itt, iz, iz1, iz2;
    int nz, i;

    double pl = ppl;
    double t = vpgen[4] * vpgen[4] - vpgen[1] * vpgen[1] - vpgen[2] * vpgen[2] - vpgen[3] * vpgen[3];
    vmax = (s + t) / 2.0;
    vmin = 2.0 * Egmin * m * 1e-1;

    if (itest == 2) goto label22;
    if (itest == 3) goto label44;
    if (ikey == 0) {
        if (itest == 1) {
            double step = (vmax - vmin) / nbin;
            for (int i = 0; i <= nbin; ++i) {
                bin[i] = vmin + i * step;
            }
        }
        vvn = vmin;
        sin = 0.0;
        nn = 0;
        distsiv[0] = 0.0;
        distarv[0] = vvn;
        xs0 = sig(t, pl, 0);
        for (iv = 1; iv <= nv; ++iv) {
            vvo = vvn;
            sio = sin;
            vvn = vmin + (vmax - vmin) * grv(iv);
            sin = fsir(t, 0.0, vvn, 0.0, pl, nn, 2);
            distsiv[iv] = distsiv[iv - 1] + (sin + sio) * (vvn - vvo) / 2.0;
            distarv[iv] = vvn;
        }
        if (itest == 1) ikey = 1;
        u0 = -s - t;
        xs0 = sig(t, pl, 0);
        xsvr = sig(t, pl, 1);
        xsB = xsBt(pl, s, t, u0) + xsBt(pl, s, u0, t);
        xsF = xs0 * dcanc(vmin, s, t, u0);
        double xsadd;
        simpsx(1e-22, vmin, 10000, 1e-3, fsirv, xsadd);
        sinonr = xs0 + xsvr + xsB + xsF + xsadd;
    }
    sirad = distsiv[nv];
    sitot = sirad + sinonr;
    if (itest == 0) weight = sitot / xs0;
    if (itest == 1) weight = sirad;
    sigmar = sirad;

    sirand = r[0] * sitot;
    if (sirand <= sinonr) {
        vgen = 0.0;
        t1gen = t;
        zgen = 0.0;
        ich = 0;
        for (i = 0; i < 4; ++i) {
            vprad[i] = vpgen[i];
            phirad[i] = 0.0;
        }
        return;
    }
    ich = 1;
    for (iv = 1; iv <= nv; ++iv) {
        if (distsiv[iv] > sirand - sinonr) {
            vgen = distarv[iv - 1] + (distarv[iv] - distarv[iv - 1]) *
                   (sirand - sinonr - distsiv[iv - 1]) /
                   (distsiv[iv] - distsiv[iv - 1]);
            if (itest == 1) return;
            goto label22;
        }
    }

label22:
    if (ikey == 0) {
        u = 4.0 * m2 + vgen - s - t;
        t1min = (2.0 * m2 * t + vgen * (t - vgen - sqrt((t - vgen) * (t - vgen) - 4.0 * m2 * t))) /
                (2.0 * (m2 + vgen));
        t1max = m2 * t * t / (m2 + vgen) / t1min;
        t1z = (-s * t * (t + u) + 2.0 * m2 * vgen * vgen) / ((s - vgen) * (s - vgen) - 4.0 * m2 * s);
        t1z1 = (-u * t * (t + s) + 2.0 * m2 * vgen * vgen) / ((u - vgen) * (u - vgen) - 4.0 * m2 * u);
        t1z2 = (-s * vgen * (t + s) - 2.0 * m2 * (2.0 * t * u + (u - 2.0 * (s + vgen)) * vgen)) / ((u - vgen) * (u - vgen) - 4.0 * m2 * u);
        t1p[0] = t1min;
        t1p[1] = t1z;
        t1p[2] = std::min(t1z1, t1z2);
        t1p[3] = std::max(t1z1, t1z2);
        t1p[4] = t1max;
        tt1n = t1min;
        sin = 0.0;
        distsit1[0] = 0.0;
        distart1[0] = tt1n;
        for (i = 1; i <= 4; ++i) {
            for (it1 = 1; it1 <= nt1; ++it1) {
                tt1o = tt1n;
                sio = sin;
                tt1n = t1p[i - 1] + (t1p[i] - t1p[i - 1]) * grt1(it1);
                zd(t, tt1n, vgen);
                sin = fsir(t, tt1n, vgen, 0.0, pl, nn, 1);
                distsit1[(i - 1) * nt1 + it1] = distsit1[(i - 1) * nt1 + it1 - 1] + (sin + sio) * (tt1n - tt1o) / 2.0;
                distart1[(i - 1) * nt1 + it1] = tt1n;
            }
        }
        sirad = distsit1[4 * nt1];
        if (itest == 2) weight = sirad;
    }
    sirand = r[1] * sirad;
    ich = 1;
    for (it1 = 1; it1 <= 4 * nt1; ++it1) {
        if (distsit1[it1] > sirand) {
            t1gen = distart1[it1 - 1] + (distart1[it1] - distart1[it1 - 1]) *
                    (sirand - distsit1[it1 - 1]) /
                    (distsit1[it1] - distsit1[it1 - 1]);
            if (itest == 2) return;
            goto label44;
        }
    }

label44:
    if (ikey == 0) {
        ich = 1;
        zd(t, t1gen, vgen);
        det = (bz * bz - az * cz);
        zmax = (-bz + sqrt(det)) / az;
        zmin = cz / az / zmax;
        if (itest == 3) ikey = 1;
        double step = (zmax - zmin) / nbin;
        for (i = 0; i <= nbin; ++i) {
            bin[i] = zmin + i * step;
        }
        zzn = zmin;
        sin = 0.0;
        nn = 0;
        distsiz[0] = 0.0;
        distarz[0] = zzn;
        for (iz = 1; iz <= nz; ++iz) {
            zzo = zzn;
            sio = sin;
            zzn = zmin + (zmax - zmin) * grz(iz);
            sin = fsir(t, t1gen, vgen, zzn, pl, nn, 0);
            distsiz[iz] = distsiz[iz - 1] + (sin + sio) * (zzn - zzo) / 2.0;
            distarz[iz] = zzn;
        }
        sirad = distsiz[nz];
        if (itest == 3) weight = sirad;
    }
    sirand = r[2] * sirad;
    for (iz = 1; iz <= nz; ++iz) {
        if (distsiz[iz] > sirand) {
            zgen = distarz[iz - 1] + (distarz[iz] - distarz[iz - 1]) *
                   (sirand - distsiz[iz - 1]) /
                   (distsiz[iz] - distsiz[iz - 1]);
            if (itest == 3) return;
            goto label55;
        }
    }

label55:
    if (itest != 0 && itest != 4) return;
    vectrec(vpgen);
    return;
}

void vectrec(std::vector<double> &vpgen) {
    double al1, al2, al3, al4, al5, al6, al7, al8;
    double sls, sl1, sl3, sl8, phi, sp, cp;
    double t = vpgen[4] * vpgen[4] - vpgen[1] * vpgen[1] - vpgen[2] * vpgen[2] - vpgen[3] * vpgen[3];
    phi = atan2(vpgen[2], vpgen[1]);
    al1 = (s - vgen) * (s - vgen) - 4.0 * m2 * s;
    al2 = s + 2.0 * t - vgen - 4.0 * m2;
    al3 = -s * t * (s + t - vgen - 4.0 * m2) - m2 * vgen * vgen;
    al4 = s * (s - vgen - 4.0 * m2) - (s + vgen) * zgen;
    al5 = vgen * zgen * (s - vgen - zgen) - m2 * (vgen + zgen) * (vgen + zgen);
    al6 = s * (vgen - zgen) - vgen * (vgen + zgen);
    al7 = (s + 2.0 * t1gen - zgen - 4.0 * m2) * al1 - al2 * al4;
    al8 = 16.0 * al3 * al5 - al7 * al7;
    sls = sqrt(als);
    sl1 = sqrt(al1);
    sl3 = sqrt(al3);
    sl8 = sqrt(al8);
    if (urand() > 0.5) sl8 = -sl8;
    sp = sin(phi);
    cp = cos(phi);
    vprad[0] = -(sls * sl1 * sl8 * sp + (4.0 * al3 * al4 - s * al2 * al7) * cp) / (4.0 * al1 * sls * sl3);
    vprad[1] = (sls * sl1 * sl8 * cp - (4.0 * al3 * al4 - s * al2 * al7) * sp) / (4.0 * al1 * sls * sl3);
    vprad[2] = (als * al1 - s * (al7 + al2 * al4)) / (2.0 * sqrt(s) * al1 * sls);
    vprad[3] = -zgen / 2.0 / sqrt(s);
    phirad[0] = (sls * sl1 * sl8 * sp + (4.0 * al3 * al6 - s * al2 * al7) * cp) / (4.0 * al1 * sls * sl3);
    phirad[1] = (-sls * sl1 * sl8 * cp + (4.0 * al3 * al6 - s * al2 * al7) * sp) / (4.0 * al1 * sls * sl3);
    phirad[2] = sqrt(s) * (al7 + al2 * al6) / (2.0 * al1 * sls);
    phirad[3] = (vgen + zgen) / 2.0 / sqrt(s);
    return;
}

void merad_init(double elab) {
    En = elab;
    s = 2.0 * (En * m + m2);
    als = s * (s - 4.0 * m2);
    coeb = 4.0 * pi * alfa * alfa * barn * (s - 2.0 * m2) / als;
    coer = alfa * alfa * alfa * barn * (s - 2.0 * m2) / als / pi / 4.0;
    Egmin = En * 1e-2;
    grid_init();
}

void grid_init() {
    int nv = 60, nt1 = 30, nz = 60;
    for (int i = 1; i <= 30; ++i) {
        grv[i - 1] = (i - 1) / 29.0 / 4.0;
    }
    for (int i = 31; i <= 45; ++i) {
        grv[i - 1] = 0.25 + (i - 30) / 15.0 / 4.0;
    }
    for (int i = 46; i <= nv; ++i) {
        grv[i - 1] = 0.5 + (i - 45) / 15.0 / 2.0;
    }
    for (int i = 1; i <= 7; ++i) {
        grt1[i - 1] = 0.1 * i * i / 49.0 / 2.0;
        grt1[30 - i] = 1.0 + grt1[0] - grt1[i - 1];
    }
    for (int i = 8; i <= 15; ++i) {
        grt1[i - 1] = (0.1 + 0.9 * (i - 7) / 8.0) / 2.0;
        grt1[30 - i] = 1.0 + grt1[0] - grt1[i - 1];
    }
    for (int i = 1; i <= 30; ++i) {
        grz[i - 1] = 0.5 * i * i / 900.0;
    }
    for (int i = 1; i <= 30; ++i) {
        grz[60 - i] = 1.0 - 0.49 * i * i / 900.0;
    }
}

double sig(double t, double pl, int i) {
    double u = 4.0 * m2 - s - t;
    double ss = s - 2.0 * m2;
    double u1 = (ss * ss + u * u) / 2.0 + 2.0 * m2 * (s + 2.0 * t - 3.0 * m2);
    double u2 = (ss * ss + t * t) / 2.0 + 2.0 * m2 * (s + 2.0 * u - 3.0 * m2);
    double u3 = (s - 2.0 * m2) * (s - 6.0 * m2);
    double pl1 = -t * (-t * s * s / 2.0 / als - ss);
    double pl2 = -t * (-t * s * s / 2.0 / als - 2.0 * m2) + als / 2.0 - ss * s;
    double pl3 = -(ss * ss * ss * ss + 4.0 * m2 * (ss * ss * t + m2 * (-ss * ss + 2.0 * t * t - 4.0 * m2 * t))) / als;
    if (i == 0) {
        return coeb * (u1 / t / t + u2 / u / u + u3 / u / t + pl * (pl1 / t / t + pl2 / u / u + pl3 / u / t));
    } else {
        double dsvt = vacpol(-t) + L1f(t, s, m2);
        double dsvu = vacpol(-u) + L1f(u, s, m2);
        double du1 = dsvt;
        double du2 = dsvu;
        double du3 = dsvt + dsvu;
        double dp1 = dsvt;
        double dp2 = dsvu;
        double dp3 = dsvt + dsvu;
        return alfa / pi * coeb * (du1 * u1 / t / t + du2 * u2 / u / u + du3 * u3 / u / t / 2.0 + pl * (dp1 * pl1 / t / t + dp2 * pl2 / u / u + dp3 * pl3 / u / t / 2.0));
    }
}

double vacpol(double t) {
    double am2[3] = {0.26110e-6, 0.111637e-1, 3.18301};
    double suml = 0.0;
    for (int i = 0; i < 3; ++i) {
        double a2 = 2.0 * am2[i];
        double sqlmi = sqrt(t * t + 2.0 * a2 * t);
        double allmi = log((sqlmi + t) / (sqlmi - t)) / sqlmi;
        suml += 2.0 * (t + a2) * allmi / 3.0 - 10.0 / 9.0 + 4.0 * a2 * (1.0 - a2 * allmi) / 3.0 / t;
    }
    double aaa, bbb, ccc, sumh;
    if (t < 1.0) {
        aaa = -1.345e-9;
        bbb = -2.302e-3;
        ccc = 4.091;
    } else if (t < 64.0) {
        aaa = -1.512e-3;
        bbb = -2.822e-3;
        ccc = 1.218;
    } else {
        aaa = -1.1344e-3;
        bbb = -3.0680e-3;
        ccc = 9.9992e-1;
    }
    sumh = -(aaa + bbb * log(1.0 + ccc * t)) * 2.0 * pi / alfa;
    return suml + sumh;
}

double L1f(double xt, double xs, double xm2) {
    return -2.0 * log(fabs(xt) / xs) * (log(fabs(xt) / xm2) - 1.0) + log(fabs(xt) / xm2) + pow(log(fabs(xt) / xm2), 2) + 4.0 * (pow(pi, 2) / 12.0 - 1.0);
}

double xsBt(double pl, double xs, double xt, double xu) {
    return 2.0 * pow(alfa, 3) / pow(xt, 2) * barn * ((1.0 + pl) * pow(xu, 2) / xs * dgg1(xs, xt, xu) - (1.0 - pl) * pow(xs, 2) / xu * dgg2(xs, xt, xu));
}

double dgg1(double xs, double xt, double xu) {
    double LS = log(xs / fabs(xt));
    double LX = log(xu / xt);
    double dgg = pow(LS, 2) * (pow(xs, 2) + pow(xu, 2)) / 2.0 / xt - LS * xu - (pow(LX, 2) + pow(pi, 2)) * pow(xu, 2) / xt;
    return 2.0 * log(xs / fabs(xu)) * log(sqrt(fabs(xu / xs))) - xt / pow(xu, 2) * dgg;
}

double dgg2(double xs, double xt, double xu) {
    double LS = log(xs / fabs(xt));
    double LX = log(xu / xt);
    double dgg = pow(LS, 2) * pow(xs, 2) / xt + LX * xs - (pow(LX, 2) + pow(pi, 2)) * (pow(xs, 2) + pow(xu, 2)) / 2.0 / xt;
    return 2.0 * log(xs / fabs(xu)) * log(sqrt(fabs(xu / xs))) - xt / pow(xs, 2) * dgg;
}

double dcanc(double xvmin, double xs, double xt, double xu) {
    double lm = log(-xt / m2);
    double lr = log(-xu / xs);
    double del1s = -2.5 * pow(lm, 2) + (3.0 - 2.0 * lr) * lm - pow(lr, 2) / 2.0 - (lm - 1.0) * log(xs * (xs + xt) / pow(xt, 2)) - pow(pi, 2) / 3.0 + 1.0;
    double del1h = -pow(lm, 2) / 2.0 + (log(pow(xt, 2) * pow(xs + xt, 2) * (xs - xvmin) / xs / (xvmin - xt) / xvmin / pow(xs + xt - xvmin, 2)) + 1.0) * lm - pow(log(-xvmin / xt), 2) / 2.0 - pow(log(1.0 - xvmin / xt), 2) + log((xs + xt) / (xs + xt - xvmin)) * log((xs + xt) * (xs + xt - xvmin) / pow(xt, 2)) + log((xvmin - xs) / xt) * log(1.0 - xvmin / xs) + log(-xvmin / xt) + fspen((xs - xvmin) / xs) - fspen((xt - xvmin) / xt) + 2.0 * (fspen(xvmin / xs) - fspen(xvmin / xt) - fspen(xvmin / (xt + xs))) - pow(pi, 2) / 6.0;
    return alfa / pi * (4.0 * log(xvmin / sqrt(m2 * xs)) * (log(xt * xu / m2 / xs) - 1.0) + del1s + del1h);
}

double fspen(double x) {
    const double f1 = 1.644934;
    if (x == 0.0) return fspens(x);
    if (x == 0.5) return fspens(x);
    if (x == 1.0) return f1 - log(x) * log(1.0 - x + 1e-10) - fspens(1.0 - x);
    if (x == 2.0) return f1 - 0.5 * log(x) * log(pow(x - 1.0, 2) / x) + fspens(1.0 - 1.0 / x);
    return 2.0 * f1 - 0.5 * pow(log(x), 2) - fspens(1.0 / x);
}

void simps(double a1, double b1, double h1, double reps1, double aeps1, double (*funct)(double), double &ai, double &aih, double &aiabs) {
    double h = (b1 > a1) ? h1 : -h1;
    double s = (h1 > 0) ? 1.0 : -1.0;
    double a = a1, b = b1;
    ai = 0.0;
    aih = 0.0;
    aiabs = 0.0;
    double p[5] = {4.0, 4.0, 2.0, 4.0, 1.0};
    double reps = fabs(reps1);
    double aeps = fabs(aeps1);
    double f[7] = {1e16, 1e16, 1e16, 1e16, 1e16, 1e16, 1e16};
    double x = a;
    double c = 0.0;
    f[0] = funct(x) / 3.0;
    while (true) {
        double x0 = x;
        if ((x0 + 4.0 * h - b) * s > 0.0) {
            h = (b - x0) / 4.0;
            if (h == 0.0) break;
            for (int k = 1; k < 7; ++k) {
                f[k] = 1e16;
            }
            c = 1.0;
        }
        double di2 = f[0];
        double di3 = fabs(f[0]);
        for (int k = 1; k < 5; ++k) {
            x += h;
            if ((x - b) * s > 0.0) x = b;
            if (f[k] == 1e16) f[k] = funct(x) / 3.0;
            di2 += p[k] * f[k];
            di3 += p[k] * fabs(f[k]);
        }
        double di1 = (f[0] + 4.0 * f[2] + f[4]) * 2.0 * h;
        di2 *= h;
        di3 *= h;
        double eps = fabs((aiabs + di3) * reps);
        if (eps < aeps) eps = aeps;
        double delta = fabs(di2 - di1);
        if (delta < eps) {
            if (delta < eps / 8.0) {
                h *= 2.0;
                f[0] = f[4];
                f[1] = f[5];
                f[2] = f[6];
                for (int k = 3; k < 7; ++k) {
                    f[k] = 1e16;
                }
            } else {
                f[0] = f[4];
                f[2] = f[5];
                f[4] = f[6];
                f[1] = 1e16;
                f[3] = 1e16;
                f[5] = 1e16;
                f[6] = 1e16;
                di1 = di2 + (di2 - di1) / 15.0;
                ai += di1;
                aih += di2;
                aiabs += di3;
            }
        } else {
            h /= 2.0;
            f[6] = f[4];
            f[5] = f[3];
            f[4] = f[2];
            f[3] = f[1];
            f[2] = 1e16;
            f[1] = 1e16;
            x = x0;
            c = 0.0;
        }
        if (c != 0.0) break;
    }
}

void simpsx(double a, double b, int np, double ep, double (*func)(double), double &res) {
    double step = (b - a) / np;
    double ra, r2, r3;
    simps(a, b, step, ep, 1e-18, func, ra, res, r2, r3);
}

double fsirv(double v) {
    return fsir(t, 0.0, v, 0.0, pl, nn, -1);
}

double fsirv1(double v) {
    return fsir(t, 0.0, v, 0.0, pl, nn, 2);
}

void zd(double t, double t1, double v) {
    double u = v - s - t + 4.0 * m2;
    az = (v - t) * (v - t) - 4.0 * m2 * t;
    bz = -(v * (2.0 * m2 * (t + t1) + t1 * (v - t))) + s * (-t * t + t1 * v + t * (t1 + v));
    cz = pow(s * (t - t1) + t1 * v, 2) - 4.0 * m2 * (s * pow(t - t1, 2) + t1 * v * v);
    az1 = az;
    bz1 = -(t * (s + t - 4.0 * m2) * (t - t1)) + (t * (2.0 * t - t1) - 2.0 * m2 * (t + t1) + s * (t + t1)) * v - t * v * v;
    cz1 = pow((s + t) * (t - t1) - t * v, 2) + 4.0 * m2 * (-pow((s + t) * (t - t1), 2) + (t - t1) * (t + t1) * v - t1 * v * v);
    az2 = az;
    bz2 = (4.0 * m2 - s - t) * t * (4.0 * m2 - t1) + (6.0 * m2 * t - s * t - 2.0 * m2 * t1 + s * t1 - t * t1) * v + (-4.0 * m2 + s) * v * v;
    cz2 = u * (-4.0 * m2 * pow(s + t1 - 4.0 * m2, 2) + (16.0 * m2 * m2 + t1 * t1 - 4.0 * m2 * (s + 2.0 * t1)) * u) - 2.0 * (2.0 * m2 - t1) * (4.0 * m2 - s - t1) * u * v + pow(s + t1 - 4.0 * m2, 2) * v * v;
}

double urand() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

