#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <algorithm>
#ifdef POWER
#include <builtins.h>
#endif

//----------------------------------------------------------------------
const int N = 20000;
const double L = 10.0;
const int D = 4;
const int X = 0;
const int Y = 1;
const int Z = 2;
double qx[N], qy[N], qz[N];
double px[N], py[N], pz[N];
const double dt = 0.001;
const double C2 = 0.1;
const double CUTOFF2 = L * L * 0.9 * 0.9;
//----------------------------------------------------------------------
double
myrand(void) {
  return static_cast<double>(rand()) / (static_cast<double>(RAND_MAX));
}
//----------------------------------------------------------------------
double
myclock(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}
//----------------------------------------------------------------------
void
calcforce(void) {
  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qx[i];
      double dy = qy[j] - qy[i];
      double dz = qz[j] - qz[i];
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      px[i] += df * dx;
      py[i] += df * dy;
      pz[i] += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
  }
}
//----------------------------------------------------------------------
void
calcforce_hand(void) {
  for (int i = 0; i < N - 1; i++) {
    const double qix = qx[i];
    const double qiy = qy[i];
    const double qiz = qz[i];
    double pix = px[i];
    double piy = py[i];
    double piz = pz[i];
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qix;
      double dy = qy[j] - qiy;
      double dz = qz[j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      pix += df * dx;
      piy += df * dy;
      piz += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
    px[i] = pix;
    py[i] = piy;
    pz[i] = piz;
  }
}
//----------------------------------------------------------------------
void
calcforce_hand_if(void) {
  for (int i = 0; i < N - 1; i++) {
    const double qix = qx[i];
    const double qiy = qy[i];
    const double qiz = qz[i];
    double pix = px[i];
    double piy = py[i];
    double piz = pz[i];
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qix;
      double dy = qy[j] - qiy;
      double dz = qz[j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CUTOFF2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      pix += df * dx;
      piy += df * dy;
      piz += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
    px[i] = pix;
    py[i] = piy;
    pz[i] = piz;
  }
}
//----------------------------------------------------------------------
double r2_flag[N];
//----------------------------------------------------------------------
void
calcforce_hand_if_pre(void) {
  for (int i = 0; i < N - 1; i++) {
    const double qix = qx[i];
    const double qiy = qy[i];
    const double qiz = qz[i];
    double pix = px[i];
    double piy = py[i];
    double piz = pz[i];
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qix;
      double dy = qy[j] - qiy;
      double dz = qz[j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      r2_flag[j] = (r2 > CUTOFF2) ? 0.0 : 1.0;
    }
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qix;
      double dy = qy[j] - qiy;
      double dz = qz[j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt * r2_flag[j];
      pix += df * dx;
      piy += df * dy;
      piz += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
    px[i] = pix;
    py[i] = piy;
    pz[i] = piz;
  }
}
//----------------------------------------------------------------------
void
calcforce_hand_if_min(void) {
  for (int i = 0; i < N - 1; i++) {
    const double qix = qx[i];
    const double qiy = qy[i];
    const double qiz = qz[i];
    double pix = px[i];
    double piy = py[i];
    double piz = pz[i];
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qix;
      double dy = qy[j] - qiy;
      double dz = qz[j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      r2 = std::min(r2, CUTOFF2);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      pix += df * dx;
      piy += df * dy;
      piz += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
    px[i] = pix;
    py[i] = piy;
    pz[i] = piz;
  }
}

//----------------------------------------------------------------------
void
calcforce_hand_if_condop(void) {
  for (int i = 0; i < N - 1; i++) {
    const double qix = qx[i];
    const double qiy = qy[i];
    const double qiz = qz[i];
    double pix = px[i];
    double piy = py[i];
    double piz = pz[i];
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qix;
      double dy = qy[j] - qiy;
      double dz = qz[j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      r2 = (r2 < CUTOFF2) ? r2 : CUTOFF2;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      pix += df * dx;
      piy += df * dy;
      piz += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
    px[i] = pix;
    py[i] = piy;
    pz[i] = piz;
  }
}

//----------------------------------------------------------------------
void
calcforce_hand_if2(void) {
  for (int i = 0; i < N - 1; i++) {
    const double qix = qx[i];
    const double qiy = qy[i];
    const double qiz = qz[i];
    double pix = px[i];
    double piy = py[i];
    double piz = pz[i];
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qix;
      double dy = qy[j] - qiy;
      double dz = qz[j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
#ifdef POWER
      df = __fsel(r2 - CUTOFF2, df, 0.0);
#else
      df = (r2 > CUTOFF2) ? 0.0 : df;
#endif
      pix += df * dx;
      piy += df * dy;
      piz += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
    px[i] = pix;
    py[i] = piy;
    pz[i] = piz;
  }
}
//----------------------------------------------------------------------
void
calcforce_if(void) {
  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qx[i];
      double dy = qy[j] - qy[i];
      double dz = qz[j] - qz[i];
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CUTOFF2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      px[i] += df * dx;
      py[i] += df * dy;
      pz[i] += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
  }
}
//----------------------------------------------------------------------
void
calcforce_if2(void) {
  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      double dx = qx[j] - qx[i];
      double dy = qy[j] - qy[i];
      double dz = qz[j] - qz[i];
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df =  ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
#ifdef POWER
      df = __fsel(r2 - CUTOFF2, df, 0.0);
#else
      df = (r2 > CUTOFF2) ? 0.0 : df;
#endif
      px[i] += df * dx;
      py[i] += df * dy;
      pz[i] += df * dz;
      px[j] -= df * dx;
      py[j] -= df * dy;
      pz[j] -= df * dz;
    }
  }
}

//----------------------------------------------------------------------
void
measure(void(*pfunc)(), const char *name) {
  double st = myclock();
  pfunc();
  double t = myclock() - st;
  printf("N=%d, %s %f [sec]\n", N, name, t);
}
//----------------------------------------------------------------------
int
main(void) {
  for (int i = 0; i < N; i++) {
    px[i] = 0.0;
    py[i] = 0.0;
    pz[i] = 0.0;
    qx[i] = L * myrand();
    qy[i] = L * myrand();
    qz[i] = L * myrand();
  }
  measure(&calcforce, "calcforce");
  measure(&calcforce_hand, "calcforce_hand");
  measure(&calcforce_if, "calcforce_if");
  measure(&calcforce_hand_if, "calcforce_hand_if");
  measure(&calcforce_hand_if_min, "calcforce_hand_if_min");
  measure(&calcforce_hand_if_min, "calcforce_hand_if_condop");
}
//----------------------------------------------------------------------
