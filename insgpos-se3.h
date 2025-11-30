/*------------------------------------------------------------------------------
 * insgpos-se3.h : INS constants, types and function prototypes
 *
 * author   : sujinglan
 *-----------------------------------------------------------------------------*/
#ifndef INSLIB_H
#define INSLIB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* constants ------------------------------------------------------------------*/
#define INSS_NONE        0  /* ins updates status: no update */
#define INSS_TIME        1  /* ins updates status: INS mechanization and propagate states and covariance */
#define INSS_GNSS        2  /* ins updates status: update GNSS position solution */
#define INSS_INIT        3  /* ins updates status: INS initial done */
#define INSS_SPPC        4  /* ins updates status: INS/GNSS-SPP */
#define INSS_RTDC        5  /* ins updates status: INS/GNSS-RTD */
#define INSS_FLTC        6  /* ins updates status: INS/GNSS-RTK float */
#define INSS_FIXC        7  /* ins updates status: INS/GNSS-RTK fix */

#define SOLQ_NONE   0       /* solution status: no solution */
#define SOLQ_FIX    1       /* solution status: fix */
#define SOLQ_FLOAT  2       /* solution status: float */
#define SOLQ_DGPS   4       /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE 5       /* solution status: single */
#define SOLQ_DR     7       /* solution status: dead reconing */

#define INS_RANDOM_WALK  1  /* stochastic process settings: random walk */
#define INS_RANDOM_CONS  2  /* stochastic process settings: random const */
#define INS_GAUSS_MARKOV 3  /* stochastic process settings: gauss-markov */

#define UDINS_IMUS       1
#define UDINS_GNSS       2

#define DATA_TYPE_IMUS   1
#define DATA_TYPE_GNSS   2

#define Mg2M             9.80665E-6
#define SCULL_CORR       1
#define MIN_CNT          5
#define VAR_POS          SQR(30.0)
#define VAR_VEL          SQR(30.0)
#define VAR_ATT          SQR(5.0*D2R)
#define VAR_BAS          SQR(1000.0*Mg2M)
#define VAR_BGS          SQR(0.5*D2R/3600.0)
#define VAR_RLM          SQR(0.1)
#define VAR_SA           SQR(1E-6)
#define VAR_SG           SQR(1E-8)

#define PI               3.1415926535897932  /* pi */
#define D2R              (PI/180.0)          /* deg to rad */
#define R2D              (180.0/PI)          /* rad to deg */
#define OMGE             7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */

#define RE_WGS84         6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84         (1.0/298.257223563) /* earth flattening (WGS84) */

#define RE               6378137.0           /* earth semi-major axis (WGS84) (m) */
#define WGS84_FE         (1.0/298.257223563) /* earth flattening (WGS84) */
#define MU               3.986004418E14      /* Earth's gravitational constant and its WGS-84 value */
#define J2               1.082627E-3         /* Earth's second gravitational constant */

#define INP(opt)         (3)
#define INV(opt)         (3)
#define INA(opt)         (3)
#define INBA(opt)        (3)
#define INBG(opt)        (3)
#define INRL(opt)        ((opt)->estrl?3:0)
#define INSA(opt)        ((opt)->estsa?3:0)
#define INSG(opt)        ((opt)->estsg?3:0)
#define INOL(opt)        ((opt)->estol?3:0)
#define INX(opt)         (INP(opt)+INV(opt)+INA(opt)+INBA(opt)+INBG(opt)+INRL(opt)+INSA(opt)+INSG(opt)+INOL(opt))

#define IIA(opt)         (0)
#define IIV(opt)         (INA(opt))
#define IIP(opt)         (INA(opt)+INV(opt))
#define IIBA(opt)        (INA(opt)+INV(opt)+INP(opt))
#define IIBG(opt)        (INA(opt)+INV(opt)+INP(opt)+INBA(opt))

#define IIRL(opt)        (INP(opt)+INV(opt)+INA(opt)+INBA(opt)+INBG(opt))
#define IISA(opt)        (INP(opt)+INV(opt)+INA(opt)+INBA(opt)+INBG(opt)+INRL(opt))
#define IISG(opt)        (INP(opt)+INV(opt)+INA(opt)+INBA(opt)+INBG(opt)+INRL(opt)+INSA(opt))
#define IIOL(opt)        (INP(opt)+INV(opt)+INA(opt)+INBA(opt)+INBG(opt)+INRL(opt)+INSA(opt)+INSG(opt))

#define SQR(x)           ((x)*(x))
#define SQRT(x)          ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

/* type definitions ----------------------------------------------------------*/
typedef struct {            /* time struct */
    time_t time;            /* time (s) expressed by standard time_t */
    double sec;             /* fraction of second under 1 s */
} gtime_t;

typedef struct {            /* solution type */
    gtime_t time;           /* time (GPST) */
    double rr[6];           /* position/velocity (m|m/s) {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    double qr[6];           /* position variance/covariance (m^2) */
                            /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
                            /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    double qv[6];           /* velocity variance/covariance (m^2/s^2) */
    double qa[6];           /* attitude variance/covariance (rad^2)*/
    double att[3];          /* attitude (rad) {roll,pitch,yaw} */
    uint8_t type;           /* type (0:xyz-ecef,1:enu-baseline) */
    uint8_t stat;           /* solution status (SOLQ_???) */
    uint8_t ns;
    double age;             /* age of differential (s) */
    double accl[3];         /* IMU accl measurement */
    double gyro[3];         /* IMU gyro measurement */
} sol_t;

typedef struct {            /* IMU data record type */
    gtime_t time;           /* IMU time */
    double gyro[3];         /* angular rate measurements in b-frame (rad/s,FRD) */
    double accl[3];         /* force measurements in b-frame (m/s^2,FRD) */
    double temp;            /* temperature (C) */
    int stat;               /* sensor status (0:ok) */
} imud_t;

typedef struct {            /* INS options struct */
    int estrl;              /* estimate lever arm terms of IMU body-frame to GNSS antenna */
    int estsa;              /* estimate scale factor of accelerometers */
    int estsg;              /* estimate scale factor of gyroscopes */
    int estol;              /* estimate  lever arm terms of IMU body-frame to odometer */
    double psd[7];          /* gyro/accl/accelerometer-bias/gyro-bias/lever-arm(IMU-GNSS)/sa/sg PSD:
                             * {rad^2/s,m^2 s^-3,m^2 s^-5,rad^2 s^-3}
                             * */
    double ba[3],bg[3];     /* IMU accelerometer-bias and gyro-bias */
    double Ma[9];           /* non-orthogonal between sensor axes and body frame for accelerometers
                             *    | sax    ra_xy  ra_xz |      | sgx    rg_xy  rg_xz |
                             * Ma=| ra_yx  say    ra_yz |   Mg=| rg_yx  sgy    rg_yz |
                             *    | ra_zx  ra_zy  saz   |      | rg_zx  rg_zy  sgz   |
                             * */
    double Mg[9];           /* non-orthogonal between sensor axes and body frame for gyroscopes */
    double Gg[9];           /* g-dependent bias for a gyro triad */
    double rl[3];           /* lever arm terms of IMU-frame to GNSS antenna */
    double ol[3];           /* lever arm terms of odometer to IMU body-frame */
    double ia[3];           /* IMU mounting angle related to vehicle body-frame ([roll,pitch,yaw]/rad) */
    double oa[3];           /* odometer mounting angle related to vehicle body-frame ([roll,pitch,yaw]/rad) */
} ins_opt_t;

typedef struct {
    gtime_t time;           /* INS filter time (GPST) */
    int stat;               /* INS filter status */
    int nx;                 /* number of estimated states */
    double tt;              /* time difference between current and previous (s) */
    double *x,*P;           /* EKF estimated states/covariance matrix */
    double *Q,*F;           /* EKF filter propagate system noise/transition matrix */
    double *xp,*Pp;         /* EKF filter propagate states/covariance matrix */
    double *x0,*P0;         /* EKF filter states/covariance matrix before propagate */
    double *dx;
} ins_flt_t;

typedef struct ins {        /* INS control/result type */
    sol_t sol;              /* INS solution data to GNSS antenna */
    gtime_t time;           /* INS current update solution time */
    int stat;               /* INS solution status */
    double dt;              /* INS time difference between current and previous (s) */
    double re[3];           /* INS position (e-frame) */
    double ve[3];           /* INS velocity (e-frame) */
    double ae[3];           /* INS acceleration (e-frame) */
    double Cbe[9];          /* INS position/velocity/acceleration (ECEF)/b-frame to e-frame trans matrix */
    double ba[3],bg[3];     /* IMU accelerometer-bias and gyro-bias */
    double Ma[9];           /* non-orthogonal between sensor axes and body frame for accelerometers
                             *    | sax    ra_xy  ra_xz |      | sgx    rg_xy  rg_xz |
                             * Ma=| ra_yx  say    ra_yz |   Mg=| rg_yx  sgy    rg_yz |
                             *    | ra_zx  ra_zy  saz   |      | rg_zx  rg_zy  sgz   |
                             * */
    double Mg[9];           /* non-orthogonal between sensor axes and body frame for gyroscopes */
    double Gg[9];           /* g-dependent bias for a gyro triad */
    double Ca[9];           /* DCM of IMU mounting angle related to vehicle body-frame */
    double Co[9];           /* DCM of odometer mounting angle related to vehicle body-frame */
    double fb0[3],omgb0[3]; /* uncorrected specific-force (b-frame)/angular rate (b-frame) */
    double fb[3],omgb[3];   /* corrected specific-force (b-frame)/angular rate (b-frame) */
    double rl[3];           /* lever arm terms of IMU body-frame to GNSS antenna */
    double ol[3];           /* lever arm terms of odometer to IMU body-frame */
    ins_flt_t flt;          /* INS EKF filter */
    ins_opt_t opt;          /* INS process options */
    imud_t  imu[2];         /* IMU measurement data */
    struct ins *pins;       /* precious epoch INS states */
} ins_t;

typedef struct {             /* data record type */
    gtime_t time;            /* data record time (GPST) */
    int type;                /* data type (see: DATA_TYPE_xxxx) */
    void *data;              /* pointer of input measurement data */
} data_t;

typedef struct {             /* igpos data type */
    int n,nmax;              /* number of observation data/allocated */
    data_t *data;            /* igpos data records */
} datas_t;

/* function definitions-------------------------------------------------------*/
extern int udins_init(ins_t *ins, const sol_t *sols, const imud_t *imu);
extern int udins_statse3(ins_t *ins, const imud_t *imu);
extern void udins_imui(ins_t *ins, gtime_t time, imud_t *imu);
extern int udins_gnssse3(ins_t *ins, const sol_t *sol);
extern void udins_uapvse3(ins_t *ins, const double *x);
extern int insposse3(int udflg, ins_t *ins, const imud_t *imu, const sol_t *sol);
extern void insinit(ins_t *ins, const ins_opt_t *opt);
extern void insfree(ins_t *ins);

extern void insemdl(const double *accl, const double *gyro, double *cor_accl, double *cor_gyro, const ins_t * ins);
extern void insgrty(const double *re, double *ge);
extern void insgrtp(const double *pos, double *g);
extern void insscullcorr(ins_t *ins, const ins_opt_t *opt, double dt, double *dv, double *da);
extern double vel2head(const double *vel);
extern void insp2gant(const double *ri, const double *rl, const double *Cbe, double *rg);

extern int is2gs(int stat);
extern int gs2is(int stat);
extern void ins2x(const ins_t *ins, double *x);

extern void quat2dcm(const double *qba, double *Cba);
extern void dcm2quat(const double *dcm, double* quat);
extern void rve2quat(const double *rvec, double *q);
extern void quatmul(const double *qab, const double *qca, double *qcb);
extern void normquat(double *q);
extern void normrpy(double *rpy);
extern void skewsym3(const double *ang, double *C);
extern void ned2xyz(const double *pos, double *Cne);
extern void rpy2dcm(const double *rpy, double *Cnb);
extern void dcm2rpy(const double *Cnb, double *rpy);

extern void so3_exp(const double *omega, double *C);
extern void so3_log(const double *C, double *omega, double *theta);

extern void matmul3(const char *tr, const double *A, const double *B, double *C);
extern void matmul3v(const char *tr, const double *A, const double *b, double *c);
extern void matt(const double *A, int n, int m, double *At);
extern void seteye(double* A, int n);
extern void setzero(double *A, int n, int m);
extern void matsb(double *A, int m, int n, const double *B, int p, int q, int isr, int isc);

extern void outisol(FILE *fp, const sol_t *sol);

extern void matcpy(double *A, const double *B, int n, int m);
extern void ecef2pos(const double *r, double *pos);
extern double dot(const double *a, const double *b, int n);
extern double timediff(gtime_t t1, gtime_t t2);
extern double *mat(int n, int m);
extern int *imat(int n, int m);
extern double *zeros(int n, int m);
extern double *eye(int n);
extern void matmul33(const char *tr,const double *A,const double *B,const double *C,
                     int n,int p,int q,int m,double *D);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
extern int matinv(double *A, int n);
extern double norm(const double *a, int n);
extern int filter(double *x, double *P, const double *H, const double *v,
                  const double *R, int n, int m, const int *ix_, int nx_);
extern void time2str(gtime_t t, char *s, int n);
extern void time2epoch(gtime_t t, double *ep);
extern void cross3(const double *a, const double *b, double *c);
extern void covenu(const double *pos, const double *P, double *Q);
extern void xyz2enu(const double *pos, double *EE);
extern gtime_t gpst2time(int week, double sec);

#ifdef __cplusplus
}
#endif
#endif