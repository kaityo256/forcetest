#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <algorithm>
#ifdef POWER
#include <builtins.h>
#endif

//----------------------------------------------------------------------
const int N = 20000;
const double L = 10.0;
const int D = 3;
const int X = 0;
const int Y = 1;
const int Z = 2;
double q[N][D],p[N][D];
const double dt = 0.001;
const double C2 = 0.1;
const double CUTOFF2 = L*L*0.9*0.9; 
//----------------------------------------------------------------------
double
myrand(void){
  return static_cast<double>(rand())/(static_cast<double>(RAND_MAX));
}
//----------------------------------------------------------------------
double
myclock(void){
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec*1e-6;
}
//----------------------------------------------------------------------
void
calcforce(void){
  for(int i=0;i<N-1;i++){
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - q[i][X];
      double dy = q[j][Y] - q[i][Y];
      double dz = q[j][Z] - q[i][Z];
      double r2 = (dx*dx + dy*dy + dz*dz);
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
      p[i][X] += df*dx;
      p[i][Y] += df*dy;
      p[i][Z] += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
  }
}
//----------------------------------------------------------------------
void
calcforce_hand(void){
  for(int i=0;i<N-1;i++){
      const double qx = q[i][X];
      const double qy = q[i][Y];
      const double qz = q[i][Z];
      double px = p[i][X];
      double py = p[i][Y];
      double pz = p[i][Z];
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - qx;
      double dy = q[j][Y] - qy;
      double dz = q[j][Z] - qz;
      double r2 = (dx*dx + dy*dy + dz*dz);
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
      px += df*dx;
      py += df*dy;
      pz += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
      p[i][X] = px;
      p[i][Y] = py;
      p[i][Z] = pz;
  }
}
//----------------------------------------------------------------------
void
calcforce_hand_if(void){
  for(int i=0;i<N-1;i++){
      const double qx = q[i][X];
      const double qy = q[i][Y];
      const double qz = q[i][Z];
      double px = p[i][X];
      double py = p[i][Y];
      double pz = p[i][Z];
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - qx;
      double dy = q[j][Y] - qy;
      double dz = q[j][Z] - qz;
      double r2 = (dx*dx + dy*dy + dz*dz);
      if(r2 > CUTOFF2) continue;
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
      px += df*dx;
      py += df*dy;
      pz += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
      p[i][X] = px;
      p[i][Y] = py;
      p[i][Z] = pz;
  }
}
//----------------------------------------------------------------------
double r2_flag[N];
//----------------------------------------------------------------------
void
calcforce_hand_if_pre(void){
  for(int i=0;i<N-1;i++){
      const double qx = q[i][X];
      const double qy = q[i][Y];
      const double qz = q[i][Z];
      double px = p[i][X];
      double py = p[i][Y];
      double pz = p[i][Z];
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - qx;
      double dy = q[j][Y] - qy;
      double dz = q[j][Z] - qz;
      double r2 = (dx*dx + dy*dy + dz*dz);
      r2_flag[j] = (r2 > CUTOFF2)? 0.0: 1.0;
    }
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - qx;
      double dy = q[j][Y] - qy;
      double dz = q[j][Z] - qz;
      double r2 = (dx*dx + dy*dy + dz*dz);
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt*r2_flag[j];
      px += df*dx;
      py += df*dy;
      pz += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
      p[i][X] = px;
      p[i][Y] = py;
      p[i][Z] = pz;
  }
}
//----------------------------------------------------------------------
void
calcforce_hand_if_min(void){
  for(int i=0;i<N-1;i++){
      const double qx = q[i][X];
      const double qy = q[i][Y];
      const double qz = q[i][Z];
      double px = p[i][X];
      double py = p[i][Y];
      double pz = p[i][Z];
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - qx;
      double dy = q[j][Y] - qy;
      double dz = q[j][Z] - qz;
      double r2 = (dx*dx + dy*dy + dz*dz);
      r2 = std::min(r2,CUTOFF2);
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
      px += df*dx;
      py += df*dy;
      pz += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
      p[i][X] = px;
      p[i][Y] = py;
      p[i][Z] = pz;
  }
}

//----------------------------------------------------------------------
void
calcforce_hand_if_condop(void){
  for(int i=0;i<N-1;i++){
      const double qx = q[i][X];
      const double qy = q[i][Y];
      const double qz = q[i][Z];
      double px = p[i][X];
      double py = p[i][Y];
      double pz = p[i][Z];
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - qx;
      double dy = q[j][Y] - qy;
      double dz = q[j][Z] - qz;
      double r2 = (dx*dx + dy*dy + dz*dz);
      r2 = (r2<CUTOFF2)? r2:CUTOFF2;
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
      px += df*dx;
      py += df*dy;
      pz += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
      p[i][X] = px;
      p[i][Y] = py;
      p[i][Z] = pz;
  }
}

//----------------------------------------------------------------------
void
calcforce_hand_if2(void){
  for(int i=0;i<N-1;i++){
      const double qx = q[i][X];
      const double qy = q[i][Y];
      const double qz = q[i][Z];
      double px = p[i][X];
      double py = p[i][Y];
      double pz = p[i][Z];
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - qx;
      double dy = q[j][Y] - qy;
      double dz = q[j][Z] - qz;
      double r2 = (dx*dx + dy*dy + dz*dz);
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
#ifdef POWER
      df = __fsel(r2 - CUTOFF2,df,0.0);
#else
      df = (r2 > CUTOFF2)? 0.0: df;
#endif
      px += df*dx;
      py += df*dy;
      pz += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
      p[i][X] = px;
      p[i][Y] = py;
      p[i][Z] = pz;
  }
}
//----------------------------------------------------------------------
void
calcforce_if(void){
  for(int i=0;i<N-1;i++){
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - q[i][X];
      double dy = q[j][Y] - q[i][Y];
      double dz = q[j][Z] - q[i][Z];
      double r2 = (dx*dx + dy*dy + dz*dz);
      if(r2 > CUTOFF2) continue;
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
      p[i][X] += df*dx;
      p[i][Y] += df*dy;
      p[i][Z] += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
  }
}
//----------------------------------------------------------------------
void
calcforce_if2(void){
  for(int i=0;i<N-1;i++){
    for(int j=i+1;j<N;j++){
      double dx = q[j][X] - q[i][X];
      double dy = q[j][Y] - q[i][Y];
      double dz = q[j][Z] - q[i][Z];
      double r2 = (dx*dx + dy*dy + dz*dz);
      double r6 = r2*r2*r2;
      double df =  ((24.0*r6-48.0)/(r6*r6*r2)+C2*8.0)*dt;
#ifdef POWER
      df = __fsel(r2 - CUTOFF2,df,0.0);
#else
      df = (r2 > CUTOFF2)? 0.0: df;
#endif
      p[i][X] += df*dx;
      p[i][Y] += df*dy;
      p[i][Z] += df*dz;
      p[j][X] -= df*dx;
      p[j][Y] -= df*dy;
      p[j][Z] -= df*dz;
    }
  }
}

//----------------------------------------------------------------------
void
measure(void(*pfunc)(),const char *name){
  double st = myclock();
  pfunc();
  double t = myclock()-st;
  printf("N=%d, %s %f [sec]\n",N,name,t);
}
//----------------------------------------------------------------------
int
main(void){
  for(int i=0;i<N;i++){
    p[i][X] = 0.0;
    p[i][Y] = 0.0;
    p[i][Z] = 0.0;
    q[i][X] = L*myrand();
    q[i][Y] = L*myrand();
    q[i][Z] = L*myrand();
  } 
  measure(&calcforce,"calcforce");
  measure(&calcforce_hand,"calcforce_hand");
  measure(&calcforce_if,"calcforce_if");
  measure(&calcforce_hand_if,"calcforce_hand_if");
  measure(&calcforce_hand_if_min,"calcforce_hand_if_min");
  measure(&calcforce_hand_if_min,"calcforce_hand_if_condop");
/*
  measure(&calcforce_if2,"calcforce_if2");
  measure(&calcforce_hand,"calcforce_hand");
  measure(&calcforce_hand_if,"calcforce_hand_if");
  measure(&calcforce_hand_if2,"calcforce_hand_if2");
  measure(&calcforce_hand_if_min,"calcforce_hand_if_min");
  measure(&calcforce_hand_if_pre,"calcforce_hand_if_pre");
*/
}
//----------------------------------------------------------------------
