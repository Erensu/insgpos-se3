/*------------------------------------------------------------------------------
 * insgpos-se3.c: INS positioning functions
 *
 * author  : sujinglan
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2021/06/17 1.0 new
 *-----------------------------------------------------------------------------*/
#include "insgpos-se3.h"

/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1, xp=x+K*v, Pp=(I-K*H')*P
*
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
static int filter_(const double *x, const double *P, const double *H,
                   const double *v, const double *R, int n, int m,
                   double *xp, double *Pp, double *Ko)
{
    double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
    double *P1=mat(n,n),*P2=mat(n,n);
    int i,j,info;

    matcpy(Q,R,m,m);
    matcpy(xp,x,n,1);
    matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */

        matmul33("NNT",I,P,I,n,n,n,n,P1);
        matmul33("NNT",K,R,K,n,m,m,n,P2);
        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) Pp[i+j*n]=P1[i+j*n]+P2[i+j*n];
        }
    }
    free(F); free(Q); free(K); free(I);
    free(P1); free(P2);
    return info;
}
extern int filter(double *x, double *P, const double *H, const double *v,
                  const double *R, int n, int m, const int *ix_, int nx_)
{
    double *x_,*xp_,*P_,*Pp_,*H_;
    int i,j,k,info,*ix;

    ix=imat(n,1);
    if (ix_&&nx_) {
        memcpy(ix,ix_,sizeof(int)*nx_);
        k=nx_;
    }
    else {
        for (i=k=0;i<n;i++) if (x[i]!=0.0&&P[i+i*n]>0.0) ix[k++]=i;
    }
    x_=mat(k,1); xp_=mat(k,1); P_=mat(k,k); Pp_=mat(k,k); H_=mat(k,m);
    for (i=0;i<k;i++) {
        x_[i]=x[ix[i]];
        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
    }
    info=filter_(x_,P_,H_,v,R,k,m,xp_,Pp_,NULL);
    for (i=0;i<k;i++) {
        x[ix[i]]=xp_[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
    }
    free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
    return info;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a, int n)
{
    return sqrt(dot(a,a,n));
}
/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors
* args   : double *a,*b     I   vector a,b (3 x 1)
*          double *c        O   outer product (a x b) (3 x 1)
* return : none
*-----------------------------------------------------------------------------*/
extern void cross3(const double *a, const double *b, double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}
/* multiply matrix -----------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

    for (i=0;i<n;i++) for (j=0;j<k;j++) {
            d=0.0;
            switch (f) {
                case 1: for (x=0;x<m;x++) d+=A[i+x*n]==0.0||B[x+j*m]==0.0?0.0:A[i+x*n]*B[x+j*m]; break;
                case 2: for (x=0;x<m;x++) d+=A[i+x*n]==0.0||B[j+x*k]==0.0?0.0:A[i+x*n]*B[j+x*k]; break;
                case 3: for (x=0;x<m;x++) d+=A[x+i*m]==0.0||B[x+j*m]==0.0?0.0:A[x+i*m]*B[x+j*m]; break;
                case 4: for (x=0;x<m;x++) d+=A[x+i*m]==0.0||B[j+x*k]==0.0?0.0:A[x+i*m]*B[j+x*k]; break;
            }
            if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
        }
}
/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
{
    double big,s,tmp,*vv=mat(n,1);
    int i,imax=0,j,k;

    *d=1.0;
    for (i=0;i<n;i++) {
        big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
        if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
            if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
        }
        if (j!=imax) {
            for (k=0;k<n;k++) {
                tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
            }
            *d=-(*d); vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j+j*n]==0.0) {free(vv); return -1;}
        if (j!=n-1) {
            tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
        }
    }
    free(vv);
    return 0;
}
/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double *A, int n, const int *indx, double *b)
{
    double s;
    int i,ii=-1,ip,j;

    for (i=0;i<n;i++) {
        ip=indx[i]; s=b[ip]; b[ip]=b[i];
        if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
        b[i]=s;
    }
    for (i=n-1;i>=0;i--) {
        s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
    }
}
/* inverse of matrix ---------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double d,*B;
    int i,j,*indx;

    indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
    if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) A[i+j*n]=0.0;
        A[j+j*n]=1.0;
        lubksb(B,n,indx,A+j*n);
    }
    free(indx); free(B);
    return 0;
}
/* D=A*B*C---------------------------------------------------------------------*/
extern void matmul33(const char *tr,const double *A,const double *B,const double *C,
                     int n,int p,int q,int m,double *D)
{
    char tr_[8];
    double *T=mat(n,q);
    matmul(tr,n,q,p,1.0,A,B,0.0,T);
    sprintf(tr_,"N%c",tr[2]);
    matmul(tr_,n,m,q,1.0,T,C,0.0,D);
    free(T);
}
/* new matrix ------------------------------------------------------------------
* allocate memory of matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *mat(int n, int m)
{
    double *p;

    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)calloc(sizeof(double),n*m))) {
        return NULL;
    }
    return p;
}
/* new integer matrix ----------------------------------------------------------
* allocate memory of integer matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern int *imat(int n, int m)
{
    int *p;

    if (n<=0||m<=0) return NULL;
    if (!(p=(int *)calloc(sizeof(int),n*m))) {
        return NULL;
    }
    return p;
}
/* zero matrix -----------------------------------------------------------------
* generate new zero matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *zeros(int n, int m)
{
    double *p;

#if NOCALLOC
    if ((p=mat(n,m))) for (n=n*m-1;n>=0;n--) p[n]=0.0;
#else
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)calloc(sizeof(double),n*m))) {
        return NULL;
    }
#endif
    return p;
}
/* identity matrix -------------------------------------------------------------
* generate new identity matrix
* args   : int    n         I   number of rows and columns of matrix
* return : matrix pointer (if n<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *eye(int n)
{
    double *p;
    int i;

    if ((p=zeros(n,n))) for (i=0;i<n;i++) p[i+i*n]=1.0;
    return p;
}
/* time difference -------------------------------------------------------------
* difference between gtime_t structs
* args   : gtime_t t1,t2    I   gtime_t structs
* return : time difference (t1-t2) (s)
*-----------------------------------------------------------------------------*/
extern double timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
}
/* convert calendar day/time to time -------------------------------------------
* convert calendar day/time to gtime_t struct
* args   : double *ep       I   day/time {year,month,day,hour,min,sec}
* return : gtime_t struct
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern gtime_t epoch2time(const double *ep)
{
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    gtime_t time={0};
    int days,sec,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];

    if (year<1970||2099<year||mon<1||12<mon) return time;

    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec=(int)floor(ep[5]);
    time.time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
    time.sec=ep[5]-sec;
    return time;
}
/* gps time to time ------------------------------------------------------------
* convert week and tow in gps time to gtime_t struct
* args   : int    week      I   week number in gps time
*          double sec       I   time of week in gps time (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2time(int week, double sec)
{
    static const double gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
    gtime_t t=epoch2time(gpst0);

    if (sec<-1E9||1E9<sec) sec=0.0;
    t.time+=(time_t)86400*7*week+(int)sec;
    t.sec=sec-(int)sec;
    return t;
}
/* copy matrix -----------------------------------------------------------------
* copy matrix
* args   : double *A        O   destination matrix A (n x m)
*          double *B        I   source matrix B (n x m)
*          int    n,m       I   number of rows and columns of matrix
* return : none
*-----------------------------------------------------------------------------*/
extern void matcpy(double *A, const double *B, int n, int m)
{
    memcpy(A,B,sizeof(double)*n*m);
}
/* set unit matrix-----------------------------------------------------------*/
extern void seteye(double* A, int n)
{
    int i,j;
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
            A[i+j*n]=(i==j?1.0:0.0);
        }
}
/* multiply 3d matrix -------------------------------------------------------*/
extern void matmul3(const char *tr, const double *A, const double *B, double *C)
{
    matmul(tr,3,3,3,1.0,A,B,0.0,C);
}
/* multiply 3d matrix and vector ---------------------------------------------*/
extern void matmul3v(const char *tr, const double *A, const double *b, double *c)
{
    char t[]="NN";
    t[0]=tr[0];
    matmul(t,3,1,3,1.0,A,b,0.0,c);
}
/* transpose matrix-----------------------------------------------------------
 * args   : double  *A      I   matrix
 *          int      n      I   rows of transpose matrix
 *          int      m      I   cols of transpose matrix
 *          double  *At     O   transpose matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void matt(const double *A, int n, int m, double *At)
{
    int i,j;
    for (i=0;i<n;i++) for (j=0;j<m;j++) At[i+j*n]=A[j+i*m];
}
/* asign a block matrix to another matrix by giving index -----------------
 * args    :  double *A   IO  input matrix for asign in elements
 *            int m,n     I   rows and cols of input matrix
 *            double *B   I   asign matrix
 *            int p,q     I   rows and cols of asign matrix
 *            int isr,isc I   start row and col to asign matrix
 * return  : none
 * ----------------------------------------------------------------------*/
extern void matsb(double *A, int m, int n, const double *B, int p, int q, int isr, int isc)
{
    int i,j;
    for (i=isr;i<isr+p;i++) {
        for (j=isc;j<isc+q;j++) A[i+j*m]=B[(i-isr)+(j-isc)*p];
    }
}
/* set all matrix elements to zero ------------------------------------------*/
extern void setzero(double *A, int n, int m)
{
    memset(A,0,sizeof(double)*n*m);
}
/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;

    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos)
{
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;

    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}
/* transform covariance to local tangental coordinate --------------------------
* transform ecef covariance to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *P        I   covariance in ecef coordinate
*          double *Q        O   covariance in local tangental coordinate
* return : none
*-----------------------------------------------------------------------------*/
extern void covenu(const double *pos, const double *P, double *Q)
{
    double EE[9],EP[9];

    xyz2enu(pos,EE);
    matmul("NN",3,3,3,1.0,EE,P,0.0,EP);
    matmul("NT",3,3,3,1.0,EP,EE,0.0,Q);
}
/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void xyz2enu(const double *pos, double *EE)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);

    EE[0]=-sinl;      EE[3]=cosl;       EE[6]=0.0;
    EE[1]=-sinp*cosl; EE[4]=-sinp*sinl; EE[7]=cosp;
    EE[2]=cosp*cosl;  EE[5]=cosp*sinl;  EE[8]=sinp;
}
/* time to string --------------------------------------------------------------
* convert gtime_t struct to string
* args   : gtime_t t        I   gtime_t struct
*          char   *s        O   string ("yyyy/mm/dd hh:mm:ss.ssss")
*          int    n         I   number of decimals
* return : none
*-----------------------------------------------------------------------------*/
extern void time2str(gtime_t t, char *s, int n)
{
    double ep[6];

    if (n<0) n=0; else if (n>12) n=12;
    if (1.0-t.sec<0.5/pow(10.0,n)) {t.time++; t.sec=0.0;};
    time2epoch(t,ep);
    sprintf(s,"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f",ep[0],ep[1],ep[2],
            ep[3],ep[4],n<=0?2:n+3,n<=0?0:n,ep[5]);
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern void time2epoch(gtime_t t, double *ep)
{
    const int mday[]={ /* # of days in a month */
            31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
            31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days,sec,mon,day;

    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(t.time/86400);
    sec=(int)(t.time-(time_t)days*86400);
    for (day=days%1461,mon=0;mon<48;mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
    ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t.sec;
}
/* INS states convert to vector states------------------------------------------*/
extern void ins2x(const ins_t *ins, double *x)
{
    so3_log(ins->Cbe,x+IIA(opt),NULL);
    matcpy(x+ IIV(&ins->opt),ins->ve,1,3);
    matcpy(x+ IIP(&ins->opt),ins->re,1,3);
    matcpy(x+IIBA(&ins->opt),ins->ba,1,3);
    matcpy(x+IIBG(&ins->opt),ins->bg,1,3);

    if (ins->opt.estrl) matcpy(x+IIRL(&ins->opt),ins->rl,1,3);
    if (ins->opt.estol) matcpy(x+IIOL(&ins->opt),ins->ol,1,3);
    if (ins->opt.estsa) {
        x[IISA(&ins->opt)+0]=ins->Ma[0];
        x[IISA(&ins->opt)+1]=ins->Ma[4];
        x[IISA(&ins->opt)+2]=ins->Ma[8];
    }
    if (ins->opt.estsg) {
        x[IISA(&ins->opt)+0]=ins->Mg[0];
        x[IISA(&ins->opt)+1]=ins->Mg[4];
        x[IISA(&ins->opt)+2]=ins->Mg[8];
    }
}
/* SE3 convert to INS states ---------------------------------------------------*/
extern int se32x(const double *T, double *x)
{
    double C[9];
    int i,j;
    for (i=0;i<3;i++) for (j=0;j<3;j++) C[i+j*3]=T[i+j*5];
    so3_log(C,x+IIA(opt),NULL);
    for (i=0;i<3;i++) {
        x[IIV(opt)+i]=T[i+3*5];
        x[IIP(opt)+i]=T[i+4*5];
    }
}
/* ins error model -----------------------------------------------------------*/
extern void insemdl(const double *accl, const double *gyro, double *cor_accl, double *cor_gyro, const ins_t *ins)
{
    int i,j;
    double Mai[9],Mgi[9],*I=eye(3),T[9]={0},Gf[3]={0};

    for (i=0;i<3;i++) for (j=0;j<3;j++) Mai[i+j*3]=I[i+j*3]+ins->Ma[i+j*3];
    for (i=0;i<3;i++) for (j=0;j<3;j++) Mgi[i+j*3]=I[i+j*3]+ins->Mg[i+j*3];

    if (!matinv(Mai,3)&&!matinv(Mgi,3)) {
        matmul("NN",3,1,3,1.0,Mai,accl,0.0,T);
        matmul("NN",3,1,3,1.0,Mgi,gyro,0.0,T+3);
    }
    for (i=0;i<3;i++) cor_accl[i]=T[i]-ins->ba[i];
    matmul("NN",3,1,3,1.0,ins->Gg,accl,0.0,Gf);
    for (i=0;i<3;i++) cor_gyro[i]=T[3+i]-ins->bg[i]-Gf[i];
    free(I);
}
/* gravity at earth surface --------------------------------------------------
 * gravity at earth surface
 * args   : double *pos      I   position {lat,lon,height} (rad/m)
 * return : gravity acceleration at earth surface (m/s^2)
 *----------------------------------------------------------------------------*/
static double gravity0(const double *pos)
{
    double sinp,e2=WGS84_FE*(2.0-WGS84_FE);
    sinp=sin(pos[0]); /* (2.85) */
    return 9.7803253359*(1.0+0.001931853*sinp*sinp)/sqrt(1-e2*sinp*sinp);
}
/* gravity model ---------------------------------------------------------------
 * gravity model in e-frame
 * args   : double *re       I   position (ecef) (m)
 *          double *ge       O   gravity acceleration (ecef) (m/s^2)
 * return : none
 *-----------------------------------------------------------------------------*/
extern void insgrty(const double *re, double *ge)
{
    double pos[3],gn[3]={0},Cne[9];
    ecef2pos(re,pos);
    gn[2]=gravity0(pos)*(1.0-2.0*pos[2]/RE);
    ned2xyz(pos,Cne);
    matmul3v("N",Cne,gn,ge);
}
/* calculates  acceleration due to gravity resolved about ecef-frame ---------
 * args  : double *pos     I   cartesian position of body frame w.r.t. ecef frame,
 *                             resolved about ecef-frame axes (m)
 *         double *g       O   Acceleration due to gravity (m/s^2) in ecef-frame
 * return: none
 * --------------------------------------------------------------------------*/
extern void insgrtp(const double *pos, double *g)
{
    double r,z,gamma;

    /* calculate distance from center of the earth */
    if ((r=norm(pos,3))<RE_WGS84/2.0) {
        g[0]=g[1]=0.0; g[2]=-9.81;
        return;
    }
    z=-MU/(r*r*r);
    gamma=1.5*J2*SQR(RE)/SQR(r);
    g[0]=z*(pos[0]+gamma*(1.0-5.0*SQR(pos[2]/r))*pos[0]);
    g[1]=z*(pos[1]+gamma*(1.0-5.0*SQR(pos[2]/r))*pos[1]);
    g[2]=z*(pos[2]+gamma*(3.0-5.0*SQR(pos[2]/r))*pos[2]);
    g[0]+=SQR(OMGE)*pos[0];
    g[1]+=SQR(OMGE)*pos[1];
}
/* estimate heading using the navigation frame velocity-----------------------
 * args   :  double *vel  I  velocity in navigation frame (ned-frame)
 * return :  heading (rad)
 * --------------------------------------------------------------------------*/
extern double vel2head(const double *vel)
{
    return atan2(vel[1],fabs(vel[0])<1E-4?1E-4:vel[0]);
}
/* 3d skew symmetric matrix --------------------------------------------------*/
extern void skewsym3(const double *ang, double *C)
{
    C[0]=0.0;     C[3]=-ang[2]; C[6]=ang[1];
    C[1]=ang[2];  C[4]=0.0;     C[7]=-ang[0];
    C[2]=-ang[1]; C[5]=ang[0];  C[8]=0.0;
}
/* quaternion convert to dcm--------------------------------------------------*/
extern void quat2dcm(const double *qba, double *Cba)
{
    double qn[4];
    matcpy(qn,qba,1,4);
    normquat(qn);
    Cba[0]=SQR(qn[0])+SQR(qn[1])-SQR(qn[2])-SQR(qn[3]);
    Cba[1]=2.0*(qn[1]*qn[2]+qn[0]*qn[3]);
    Cba[2]=2.0*(qn[1]*qn[3]-qn[0]*qn[2]);
    Cba[3]=2.0*(qn[1]*qn[2]-qn[0]*qn[3]);
    Cba[4]=SQR(qn[0])-SQR(qn[1])+SQR(qn[2])-SQR(qn[3]);
    Cba[5]=2.0*(qn[2]*qn[3]+qn[0]*qn[1]);
    Cba[6]=2.0*(qn[1]*qn[3]+qn[0]*qn[2]);
    Cba[7]=2.0*(qn[2]*qn[3]-qn[0]*qn[1]);
    Cba[8]=SQR(qn[0])-SQR(qn[1])-SQR(qn[2])+SQR(qn[3]);
}
/* dcm convert to quaternion-------------------------------------------------*/
extern void dcm2quat(const double *C, double* q)
{
    double tr=C[0]+C[4]+C[8];
    double P1=1.0+tr;
    double P2=1.0+2.0*C[0]-tr;
    double P3=1.0+2.0*C[4]-tr;
    double P4=1.0+2.0*C[8]-tr;
    int i;

    if (P1>=P2&&P1>=P3&&P1>=P4) {
        q[0]=0.5*sqrt(P1);
        q[1]=(C[5]-C[7])/4.0/q[0];
        q[2]=(C[6]-C[2])/4.0/q[0];
        q[3]=(C[1]-C[3])/4.0/q[0];
    }
    else if (P2>=P3&&P2>=P4) {
        q[1]=0.5*sqrt(P2);
        q[2]=(C[1]+C[3])/4.0/q[1];
        q[3]=(C[6]+C[2])/4.0/q[1];
        q[0]=(C[5]-C[7])/4.0/q[1];
    }
    else if (P3>=P4) {
        q[2]=0.5*sqrt(P3);
        q[3]=(C[5]+C[7])/4.0/q[2];
        q[0]=(C[6]-C[2])/4.0/q[2];
        q[1]=(C[1]+C[3])/4.0/q[2];
    }
    else {
        q[3]=0.5*sqrt(P4);
        q[0]=(C[1]-C[3])/4.0/q[3];
        q[1]=(C[6]+C[2])/4.0/q[3];
        q[2]=(C[5]+C[7])/4.0/q[3];
    }
    if (q[0]<0) {
        for (i=0;i<4;i++) q[i]=-q[i];
    }
    normquat(q);
}
/* quaternion multiply-------------------------------------------------------*/
extern void quatmul(const double *qab, const double *qca, double *qcb)
{
    setzero(qcb,4,1);
    qcb[0]=(qab[0])*qca[0]+(-qab[1])*qca[1]+(-qab[2])*qca[2]+(-qab[3])*qca[3];
    qcb[1]=(qab[1])*qca[0]+( qab[0])*qca[1]+(-qab[3])*qca[2]+( qab[2])*qca[3];
    qcb[2]=(qab[2])*qca[0]+( qab[3])*qca[1]+( qab[0])*qca[2]+(-qab[1])*qca[3];
    qcb[3]=(qab[3])*qca[0]+(-qab[2])*qca[1]+( qab[1])*qca[2]+( qab[0])*qca[3];
}
/* normlization quaternion---------------------------------------------------*/
extern void normquat(double *q)
{
    int i;
    double e=norm(q,4);
    for (i=0;i<4;i++) q[i]/=e;
}
/* rotate vector convert to quaternion---------------------------------------*/
extern void rve2quat(const double *rv, double *q)
{
    double a=norm(rv,3);
    if (a<1E-12) {
        q[0]=1.0;
        q[1]=q[2]=q[3]=0.0;
        return;
    }
    q[0]=cos(0.5*a);
    q[1]=sin(0.5*a)/(0.5*a)*0.5*rv[0];
    q[2]=sin(0.5*a)/(0.5*a)*0.5*rv[1];
    q[3]=sin(0.5*a)/(0.5*a)*0.5*rv[2];
    normquat(q);
}
/* conversion matrix of ned frame to ecef frame--------------------------------
 * conversion matrix of ned to ecef
 * args   : double *pos     I   position {lat,lon,height} (rad/m)
 *          double *Cne     O   conversion matrix between frame
 * return : none
 * --------------------------------------------------------------------------*/
extern void ned2xyz(const double *pos, double *Cne)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);

    Cne[0]=-sinp*cosl; Cne[3]=-sinl; Cne[6]=-cosp*cosl;
    Cne[1]=-sinp*sinl; Cne[4]= cosl; Cne[7]=-cosp*sinl;
    Cne[2]=      cosp; Cne[5]= 0.0;  Cne[8]=-sinp;
}
/* ned frame to body frame---------------------------------------------------
 * args    : double rpy      I  attitude {roll,picth,yaw} (rad)
 *           double Cnb      O  matrix of ned to body frame
 * return  : none
 * -------------------------------------------------------------------------*/
extern void rpy2dcm(const double *rpy, double *Cnb)
{
    double sin_phi=sin(rpy[0]),cos_phi=cos(rpy[0]),
            sin_theta=sin(rpy[1]),cos_theta=cos(rpy[1]),
            sin_psi=sin(rpy[2]),cos_psi=cos(rpy[2]);
    Cnb[0]= cos_theta*cos_psi;
    Cnb[3]= cos_theta*sin_psi;
    Cnb[6]=-sin_theta;
    Cnb[1]=-cos_phi*sin_psi+sin_phi*sin_theta*cos_psi;
    Cnb[4]= cos_phi*cos_psi+sin_phi*sin_theta*sin_psi;
    Cnb[7]= sin_phi*cos_theta;
    Cnb[2]= sin_phi*sin_psi+cos_phi*sin_theta*cos_psi;
    Cnb[5]=-sin_phi*cos_psi+cos_phi*sin_theta*sin_psi;
    Cnb[8]= cos_phi*cos_theta;
}
/* converts a coordinate transformation matrix to the corresponding set of
 * Euler angles -------------------------------------------------------------
 * args    : double *Cnb      I matrix of ned-frame to body-frame
 *           double *rpy      O attitude {roll,pitch,yaw}
 * return  : none
 * --------------------------------------------------------------------------*/
extern void dcm2rpy(const double *Cnb, double *rpy)
{
    rpy[0]= atan2(Cnb[7],Cnb[8]);  /* roll */
    rpy[1]=-asin (Cnb[6]);         /* pitch */
    rpy[2]= atan2(Cnb[3],Cnb[0]);  /* yaw */
}
/* normalization RPY----------------------------------------------------------*/
extern void normrpy(double *rpy)
{
    rpy[0]=fmod(rpy[0],2.0*PI);
    rpy[1]=fmod(rpy[1],2.0*PI);
    rpy[2]=fmod(rpy[2],2.0*PI);

    if (fabs(rpy[0])>fabs(2.0*PI-fabs(rpy[0]))) {
        if      (rpy[0]> 0.0) rpy[0]-=2.0*PI;
        else if (rpy[0]<=0.0) rpy[0]+=2.0*PI;
    }
    if (fabs(rpy[1])>fabs(2.0*PI-fabs(rpy[1]))) {
        if      (rpy[1]> 0.0) rpy[1]-=2.0*PI;
        else if (rpy[1]<=0.0) rpy[1]+=2.0*PI;
    }
    if (fabs(rpy[2])>fabs(2.0*PI-fabs(rpy[2]))) {
        if      (rpy[2]> 0.0) rpy[2]-=2.0*PI;
        else if (rpy[2]<=0.0) rpy[2]+=2.0*PI;
    }
}
/* so3 exponential map--------------------------------------------------------*/
extern void so3_exp(const double *omega, double *C)
{
    double I[9]={1,0,0,0,1,0,0,0,1};
    double a,n[3],S[9],S2[9];
    int i;

    if ((a=norm(omega,3))<1E-8) {
        skewsym3(omega,S);
        for (i=0;i<9;i++) C[i]=I[i]+S[i];
        return;
    }
    for (i=0;i<3;i++) n[i]=omega[i]/a;
    skewsym3(n,S);
    matmul("NN",3,3,1,1.0,n,n,0.0,S2);

    for (i=0;i<9;i++) {
        C[i]=cos(a)*I[i]+(1.0-cos(a))*S2[i]+sin(a)*S[i];
    }
}
/* so3 logarithm map----------------------------------------------------------*/
extern void so3_log(const double *C, double *omega, double *theta)
{
    double tr=C[0]+C[4]+C[8];
    double a=acos((tr-1.0)/2.0);
    double c=fabs(a)>1E-6?a/(2.0*sin(a)):0.5,S[9];
    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) S[i+j*3]=c*(C[i+j*3]-C[j+i*3]);
    }
    omega[0]=S[5];
    omega[1]=S[6];
    omega[2]=S[1];
    if (theta) *theta=norm(omega,3);
}
/* INS position convert to GNSS antenna---------------------------------------*/
extern void insp2gant(const double *ri, const double *rl, const double *Cbe, double *rg)
{
    double T[3];
    matmul3v("N",Cbe,rl,T);
    rg[0]=ri[0]+T[0];
    rg[1]=ri[1]+T[1];
    rg[2]=ri[2]+T[2];
}
/* rotational and sculling motion correction --------------------------------*/
extern void insscullcorr(ins_t *ins, const ins_opt_t *opt, double dt, double *dv, double *da)
{
    double dap[3]={0},dvp[3]={0},dak[3],dvk[3],dv1[3],dv2[3],dv3[3];
    double domg,a1,a2,dv4[3];
    int i;
    for (i=0;i<3;i++) {
        if (ins->pins) dap[i]=ins->pins->omgb[i]*ins->dt;
        if (ins->pins) dvp[i]=ins->pins->fb[i]*ins->dt;
    }
    for (i=0;i<3;i++) {
        dak[i]=ins->omgb[i]*dt; dvk[i]=ins->fb[i]*dt;
    }
    cross3(dak,dvk,dv1);
    cross3(dap,dvk,dv2);
    cross3(dvp,dak,dv3);

    domg=norm(dak,3);
    if (fabs(domg)>1E-6) {
        a1=(1.0-cos(domg))/SQR(domg);
        a2=1.0/SQR(domg)*(1.0-sin(domg)/domg);
    }
    else {
        a1=0.5-SQR(domg)/24.0+SQR(SQR(domg))/720.0;
        a2=1.0/6.0-SQR(domg)/120.0+SQR(SQR(domg))/5040.0;
    }
    cross3(dak,dv1,dv4);

    for (i=0;i<3;i++) {
        dv[i]=a1*dv1[i]+a2*dv4[i]+1.0/12.0*(dv2[i]+dv3[i]);
    }
    cross3(dap,dak,da);
    for (i=0;i<3;i++) da[i]=1.0/12.0*da[i];
}
/* INS solution convert to GNSS solution status-------------------------------*/
extern int is2gs(int stat)
{
    if      (stat==INSS_TIME) return SOLQ_DR;
    else if (stat==INSS_SPPC) return SOLQ_SINGLE;
    else if (stat==INSS_RTDC) return SOLQ_DGPS;
    else if (stat==INSS_FLTC) return SOLQ_FLOAT;
    else if (stat==INSS_FIXC) return SOLQ_FIX;
    else return SOLQ_NONE;
}
/* GNSS solution convert to INS solution status-------------------------------*/
extern int gs2is(int stat)
{
    if      (stat==SOLQ_DR    ) return INSS_TIME;
    else if (stat==SOLQ_SINGLE) return INSS_SPPC;
    else if (stat==SOLQ_DGPS  ) return INSS_RTDC;
    else if (stat==SOLQ_FLOAT ) return INSS_FLTC;
    else if (stat==SOLQ_FIX   ) return INSS_FIXC;
    else return INSS_NONE;
}
/* solution to covariance ----------------------------------------------------*/
static void soltocov(const sol_t *sol, double *P)
{
    P[0]     =sol->qr[0]; /* xx or ee */
    P[4]     =sol->qr[1]; /* yy or nn */
    P[8]     =sol->qr[2]; /* zz or uu */
    P[1]=P[3]=sol->qr[3]; /* xy or en */
    P[5]=P[7]=sol->qr[4]; /* yz or nu */
    P[2]=P[6]=sol->qr[5]; /* zx or ue */
}
/* sqrt of covariance --------------------------------------------------------*/
static double sqvar(double covar)
{
    return covar<0.0?-sqrt(-covar):sqrt(covar);
}
/* output INS solution to file------------------------------------------------*/
extern void outisol(FILE *fp, const sol_t *sol)
{
    double pos[3],P[9]={0},Q[9];
    const char *sep="  ";
    char buff[1024],s[64],*p=(char *)buff;

    time2str(sol->time,s,3);
    ecef2pos(sol->rr,pos);
    soltocov(sol,P); covenu(pos,P,Q);

    p+=sprintf(p,"%s%s%14.9f%s%14.9f",s,sep,pos[0]*R2D,sep,pos[1]*R2D);
    p+=sprintf(p,"%s%10.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s",
               sep,pos[2],sep,sol->stat,sep,sol->ns,sep,SQRT(Q[4]),sep,
               SQRT(Q[0]),sep,SQRT(Q[8]),sep,sqvar(Q[1]),sep,sqvar(Q[2]),
               sep,sqvar(Q[5]),sep,sol->age,sep);

    p+=sprintf(p,"%s%8.3f%s%8.3f%s%8.3f",sep,sol->att[0]*R2D,sep,sol->att[1]*R2D,sep,sol->att[2]*R2D);
    p+=sprintf(p,"%s%8.3f%s%8.3f%s%8.3f",sep,sol->rr[3],sep,sol->rr[4],sep,sol->rr[5]);
    sprintf(p,"\n");
    fprintf(fp,"%s",buff);
    fflush(fp);
}
/* initial INS struct---------------------------------------------------------*/
extern void insinit(ins_t *ins, const ins_opt_t *opt)
{
    ins->flt.nx=INX(opt);
    ins->flt.x =mat(ins->flt.nx,1);
    ins->flt.xp=mat(ins->flt.nx,1);
    ins->flt.x0=mat(ins->flt.nx,1);
    ins->flt.dx=mat(ins->flt.nx,1);
    ins->flt.P =mat(ins->flt.nx,ins->flt.nx);
    ins->flt.F =mat(ins->flt.nx,ins->flt.nx);
    ins->flt.Q =mat(ins->flt.nx,ins->flt.nx);
    ins->flt.Pp=mat(ins->flt.nx,ins->flt.nx);
    ins->flt.P0=mat(ins->flt.nx,ins->flt.nx);
    ins->opt=*opt;
    matcpy(ins->rl,opt->rl,1,3);
    matcpy(ins->ol,opt->ol,1,3);
    matcpy(ins->ba,opt->ba,1,3);
    matcpy(ins->bg,opt->bg,1,3);
    matcpy(ins->Ma,opt->Ma,3,3);
    matcpy(ins->Mg,opt->Mg,3,3);
    matcpy(ins->Gg,opt->Gg,1,3);
    ins->pins=calloc(1,sizeof(ins_t));
    rpy2dcm(opt->ia,ins->Ca);
    rpy2dcm(opt->oa,ins->Co);
}
/* INS states convert to solution----------------------------------------------*/
static void ins2sol(int type, const ins_t *ins, sol_t *solo)
{
    const ins_flt_t *flt=&ins->flt;
    int i,nx=ins->flt.nx;
    double pos[3],Cne[9],Cnb[9],*P=flt->P;

    /* update INS solutions */
    solo->time=flt->stat?flt->time:ins->time;
    solo->stat=is2gs(ins->stat);

    matcpy(solo->rr+0,ins->re,1,3);
    matcpy(solo->rr+3,ins->ve,1,3);

    ecef2pos(ins->re,pos);
    ned2xyz(pos,Cne);
    matmul3("TN",ins->Cbe,Cne,Cnb);
    dcm2rpy(Cnb,solo->att);

    /* update INS covariance */
    solo->qr[3]=P[IIP(opt)+0+(IIP(opt)+1)*nx];
    solo->qr[4]=P[IIP(opt)+1+(IIP(opt)+2)*nx];
    solo->qr[5]=P[IIP(opt)+2+(IIP(opt)+0)*nx];
    solo->qv[3]=P[IIV(opt)+0+(IIV(opt)+1)*nx];
    solo->qv[4]=P[IIV(opt)+1+(IIV(opt)+2)*nx];
    solo->qv[5]=P[IIV(opt)+2+(IIV(opt)+0)*nx];
    solo->qa[3]=P[IIA(opt)+0+(IIA(opt)+1)*nx];
    solo->qa[4]=P[IIA(opt)+1+(IIA(opt)+2)*nx];
    solo->qa[5]=P[IIA(opt)+2+(IIA(opt)+0)*nx];
    for (i=0;i<3;i++) {
        solo->qr[i]=P[IIP(opt)+i+(IIP(opt)+i)*nx];
        solo->qv[i]=P[IIV(opt)+i+(IIV(opt)+i)*nx];
        solo->qa[i]=P[IIA(opt)+i+(IIA(opt)+i)*nx];
    }
}
/* free INS struct-----------------------------------------------------------*/
extern void insfree(ins_t *ins)
{
    free(ins->flt.x ); ins->flt.x =NULL;
    free(ins->flt.P ); ins->flt.P =NULL;
    free(ins->flt.F ); ins->flt.F =NULL;
    free(ins->flt.Q ); ins->flt.Q =NULL;
    free(ins->flt.xp); ins->flt.xp=NULL;
    free(ins->flt.x0); ins->flt.x0=NULL;
    free(ins->flt.P0); ins->flt.P0=NULL;
    free(ins->flt.Pp); ins->flt.Pp=NULL;
    free(ins->flt.dx); ins->flt.dx=NULL;
    free(ins->pins); ins->pins=NULL;
    memset(ins,0,sizeof(ins_t));
}
/* INS states convert to SE3---------------------------------------------------*/
extern int x2se3(const double *x, double *T)
{
    double Cbe[9];
    int i,j;

    so3_exp(x,Cbe); seteye(T,5);
    for (i=0;i<3;i++) for (j=0;j<3;j++) T[i+j*5]=Cbe[i+j*3];
    for (i=0;i<3;i++) T[i+3*5]=x[IIV(opt)+i];
    for (i=0;i<3;i++) T[i+4*5]=x[IIP(opt)+i];
}
/* update INS filter system noise matrix---------------------------------------*/
static void udins_Qs(const ins_t *ins, const ins_opt_t *opt, double dt, const double *F, double *Q)
{
    int i,j,nx=ins->flt.nx;
    double *Qn=zeros(12,12),*Gn=mat(12,nx),*Qp=mat(nx,nx),I[9]={1,0,0,0,1,0,0,0,1};

    for (i=0;i< 3;i++) Qn[i+i*12]=opt->psd[0]*fabs(dt);
    for (i=3;i< 6;i++) Qn[i+i*12]=opt->psd[1]*fabs(dt);
    for (i=6;i< 9;i++) Qn[i+i*12]=opt->psd[2]*fabs(dt);
    for (i=9;i<12;i++) Qn[i+i*12]=opt->psd[3]*fabs(dt);
    matsb(Gn,nx,12,I,3,3 ,0,0);
    matsb(Gn,nx,12,I,3,3 ,3,3);
    matsb(Gn,nx,12,I,3,3 ,9,6);
    matsb(Gn,nx,12,I,3,3,12,9);
    matmul33("NNT",Gn,Qn,Gn,nx,12,12,nx,Qp);

    matmul33("NNT",F,Qp,F,nx,nx,nx,nx,Q);
    for(i=0;i<nx;i++) for (j=0;j<nx;j++) Q[i+j*nx]=(Q[i+j*nx]+Qp[i+j*nx])/2.0;
    free(Qn); free(Gn); free(Qp);

    for (i=IIRL(opt);i<IIRL(opt)+INRL(opt);i++) Q[i+i*nx]=opt->psd[4]*fabs(dt);
    for (i=IISA(opt);i<IISA(opt)+INSA(opt);i++) Q[i+i*nx]=opt->psd[5]*fabs(dt);
    for (i=IISG(opt);i<IISG(opt)+INSG(opt);i++) Q[i+i*nx]=opt->psd[6]*fabs(dt);
}
/* propagate matrix for stochastic parameters----------------------------------*/
static void stochasticPhi(int opt, int ix, int nix, int nx, double dt, double *F)
{
    double corr_time=3600.0;
    int i;
    for (i=ix;i<ix+nix;i++) {
        if      (opt==INS_RANDOM_CONS ) F[i+i*nx]=1.0;
        else if (opt==INS_RANDOM_WALK ) F[i+i*nx]=1.0;
        else if (opt==INS_GAUSS_MARKOV) F[i+i*nx]=exp(-dt/corr_time);
    }
}
/* update INS filter transmit matrix--------------------------------------------*/
static void udins_Fp(const ins_t *ins, const ins_opt_t *opt, double dt, double *F)
{
    int i,j,nx=ins->flt.nx;
    double we[3]={0,0,OMGE},Omge[9],Ww[9],Wf[9],We[9],I[9]={1,0,0,0,1,0,0,0,1};
    double omega[3]={0},T[9];

    seteye(F,nx);
    skewsym3(ins->omgb,Ww);

    /* attitude transmit matrix */
    for (i=IIA(opt);i<IIA(opt)+INA(opt);i++) {
        for (j= IIA(opt);j< IIA(opt)+ INA(opt);j++) F[i+j*nx]-=Ww[i-IIA(opt)+(j-IIA(opt))*3]*dt;
        for (j=IIBG(opt);j<IIBG(opt)+INBG(opt);j++) F[i+j*nx] =I[i-IIA(opt)+(j-IIBG(opt))*3]*dt;
        for (j=IISG(opt);j<IISG(opt)+INSG(opt);j++) F[i+j*nx] =I[i-IIA(opt)+(j-IISG(opt))*3]*ins->omgb[j-IISG(opt)]*dt;
    }
    /* velocity transmit matrix */
    matmul("TN",3,1,3,1.0,ins->Cbe,we,0.0,omega);
    omega[0]+=ins->omgb[0];
    omega[1]+=ins->omgb[1];
    omega[2]+=ins->omgb[2];
    skewsym3(omega,Ww);
    skewsym3(ins->fb,Wf);
    skewsym3(we,Omge);
    matmul33("TNN",ins->Cbe,Omge,Omge,3,3,3,3,T);
    matmul("NN",3,3,3,1.0,T,ins->Cbe,0.0,We);

    for (i=IIV(opt);i<IIV(opt)+INV(opt);i++) {
        for (j= IIV(opt);j< IIV(opt)+ INV(opt);j++) F[i+j*nx]-= Ww[i-IIV(opt)+(j-IIV(opt))*3]*dt;
        for (j= IIA(opt);j< IIA(opt)+ INA(opt);j++) F[i+j*nx] =-Wf[i-IIV(opt)+(j-IIA(opt))*3]*dt;
        for (j= IIP(opt);j< IIP(opt)+ INP(opt);j++) F[i+j*nx] =-We[i-IIV(opt)+(j-IIP(opt))*3]*dt;
        for (j=IIBA(opt);j<IIBA(opt)+INBA(opt);j++) F[i+j*nx] =I[i-IIV(opt)+(j-IIBA(opt))*3]*dt;
        for (j=IISA(opt);j<IISA(opt)+INSA(opt);j++) F[i+j*nx] =I[i-IIV(opt)+(j-IISA(opt))*3]*ins->fb[j-IISA(opt)]*dt;
    }
    /* position transmit matrix */
    matmul("TN",3,1,3,1.0,ins->Cbe,we,0.0,omega);
    omega[0]=ins->omgb[0]-omega[0];
    omega[1]=ins->omgb[1]-omega[1];
    omega[2]=ins->omgb[2]-omega[2];
    skewsym3(omega,Ww);

    for (i=IIP(opt);i<IIP(opt)+INP(opt);i++) {
        for (j=IIP(opt);j<IIP(opt)+INP(opt);j++) F[i+j*nx]-=Ww[i-IIP(opt)+(j-IIP(opt))*3]*dt;
        for (j=IIV(opt);j<IIV(opt)+INV(opt);j++) F[i+j*nx] =(i-IIP(opt))==(j-IIV(opt))?dt:0.0;
    }
    /* for stochastic parameters */
    stochasticPhi(INS_GAUSS_MARKOV,IIBA(opt),3,nx,dt,F);
    stochasticPhi(INS_GAUSS_MARKOV,IIBG(opt),3,nx,dt,F);
    if (opt->estrl) {
        stochasticPhi(INS_RANDOM_CONS,IIRL(opt),3,nx,dt,F);
    }
    if (opt->estsa) {
        stochasticPhi(INS_GAUSS_MARKOV,IISA(opt),3,nx,dt,F);
    }
    if (opt->estsg) {
        stochasticPhi(INS_GAUSS_MARKOV,IISG(opt),3,nx,dt,F);
    }
}
/* initial INS filter covariance----------------------------------------------*/
static void udins_P0(ins_t *ins, gtime_t time, const ins_opt_t *opt, const sol_t *sol)
{
    ins_flt_t *flt=&ins->flt;
    int i,nx=flt->nx;

    setzero(flt->x, 1,nx);
    setzero(flt->P,nx,nx);
    flt->time=time;
    for (i= IIP(opt);i< IIP(opt)+ INP(opt);i++) flt->P[i+i*nx]=sol->qr[i-IIP(opt)]?sol->qr[i-IIP(opt)]:VAR_POS;
    for (i= IIV(opt);i< IIV(opt)+ INV(opt);i++) flt->P[i+i*nx]=sol->qv[i-IIV(opt)]?sol->qv[i-IIV(opt)]:VAR_VEL;
    for (i= IIA(opt);i< IIA(opt)+ INA(opt);i++) flt->P[i+i*nx]=VAR_ATT;
    for (i=IIBA(opt);i<IIBA(opt)+INBA(opt);i++) flt->P[i+i*nx]=VAR_BAS;
    for (i=IIBG(opt);i<IIBG(opt)+INBG(opt);i++) flt->P[i+i*nx]=VAR_BGS;
    for (i=IIRL(opt);i<IIRL(opt)+INRL(opt);i++) flt->P[i+i*nx]=VAR_RLM;
    for (i=IISA(opt);i<IISA(opt)+INSA(opt);i++) flt->P[i+i*nx]=VAR_SA;
    for (i=IISG(opt);i<IISG(opt)+INSG(opt);i++) flt->P[i+i*nx]=VAR_SG;
}
/* update INS filter covariance-----------------------------------------------*/
static void udins_Pp(ins_t *ins, const ins_opt_t *opt, double dt)
{
    int i,j,nx=ins->flt.nx;
    double *Q=mat(nx,nx),*F=mat(nx,nx),*T=mat(nx,nx);

    udins_Fp(ins,opt,dt,F);
    udins_Qs(ins,opt,dt,F,Q);

    matmul("NN",nx,nx,nx,1.0,F,ins->flt.P,0.0,T);
    matmul("NT",nx,nx,nx,1.0,T,F,0.0,ins->flt.P);

    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) ins->flt.P[i+j*nx]=ins->flt.P[i+j*nx]+Q[i+j*nx];
    }
    matcpy(ins->flt.F,F,nx,nx);
    matcpy(ins->flt.Q,Q,nx,nx);
    free(Q); free(F); free(T);
}
/* update INS filter states---------------------------------------------------*/
static void udins_xp(ins_t *ins, const ins_opt_t *opt, double dt)
{
    int i;
    for (i=0;i<ins->flt.nx;i++) ins->flt.x[i]=1E-30;
}
/* update INS filter---------------------------------------------------------*/
static int udins_ufltse3(ins_t *ins, const imud_t *imu)
{
    if (ins->stat==INSS_NONE) return 0;
    ins->flt.tt=timediff(imu->time,ins->flt.time);
    ins->flt.time=imu->time;

    int nx=ins->flt.nx;
    matcpy(ins->flt.P0,ins->flt.P,nx,nx);

    udins_xp(ins,&ins->opt,ins->flt.tt);
    udins_Pp(ins,&ins->opt,ins->flt.tt);

    ins2x(ins,ins->flt.xp);
    matcpy(ins->flt.Pp,ins->flt.P,nx,nx);
    return 1;
}
/* jacobian of position measurement by attitude error term ------------------*/
static void jacobian_p_att(const double *Cbe, const double *lever, double *dpdatt)
{
    double wl[9];
    skewsym3(lever,wl);
    matmul("NN",3,3,3,-1.0,Cbe,wl,0.0,dpdatt);
}
/* GNSS solution measurement residuals----------------------------------------*/
static int gnssres(ins_t *ins, const sol_t *sol, const double *x, const double *P,
                   double *v, double *H, double *R)
{
    double rg[3],r1[9],Rv[9]={0},Rp[9]={0},var[6],Rpr[9],Rvr[9];
    int i,j,nv=0,nx=ins->flt.nx;

    /* INS position to GNSS antenna position */
    insp2gant(ins->re,ins->rl,ins->Cbe,rg);
    jacobian_p_att(ins->Cbe,ins->rl,r1);

    /* initial GNSS solution measurement residuals and variance */
    for (i=0;i<6;i++) var[i]=1E8;
    for (i=0;i<6;i++) v[i]=0.0;
    for (i=0;i<6;i++) for (j=0;j< 6;j++) R[i+j*6]=0.0;
    for (i=0;i<6;i++) for (j=0;j<nx;j++) H[j+i*nx]=0.0;

    /* check GNSS solution validation */
    if (norm(sol->rr+0,3)<=0.0) {
        return 0;
    }
    /* GNSS position measurement residuals */
    for (i=0;i<3;i++) {
        if (!sol->rr[i]) continue;
        var[i]=!sol->qr[i]?VAR_POS:sol->qr[i];
        v[nv]=sol->rr[i]-rg[i];
        for (j=IIA(opt);j<IIA(opt)+INA(opt);j++) H[j+nv*nx]=r1[i+(j-IIA(opt))*3];
        for (j=IIP(opt);j<IIP(opt)+INP(opt);j++) H[j+nv*nx]=ins->Cbe[i+(j-IIP(opt))*3];
        nv++;
    }
    /* GNSS position/velocity measurement variance */
    for (i=0;i<3;i++) {
        Rp[i+i*3]=var[i+0];
        Rv[i+i*3]=var[i+3];
    }
    matmul33("TNN",ins->Cbe,Rp,ins->Cbe,3,3,3,3,Rpr);
    matmul33("TNN",ins->Cbe,Rv,ins->Cbe,3,3,3,3,Rvr);
    for (i=0;i<3;i++) for (j=0;j<3;j++) R[i+j*6]=Rpr[i-0+(j-0)*3];
    for (i=3;i<6;i++) for (j=3;j<6;j++) R[i+j*6]=Rvr[i-3+(j-3)*3];
    return 6;
}
/* INS filter close loop correction-------------------------------------------*/
extern void udins_uapvse3(ins_t *ins, const double *x)
{
    double I[9]={1,0,0,0,1,0,0,0,1},omgbc[3],fbc[3],Cbec[9];
    double fe[3],ge[3],Omge[9]={0,OMGE,0,-OMGE,0,0,0,0,0},cori[3],W[9];
    double dp[3],dv[3];
    int i,j;

    double dT[25],Tp[25],Tc[25];

    x2se3(x,dT); seteye(Tp,5);
    for (i=0;i<3;i++) for (j=0;j<3;j++) Tp[i+j*5]=ins->Cbe[i+j*3];
    for (i=0;i<3;i++) Tp[i+3*5]=ins->ve[i];
    for (i=0;i<3;i++) Tp[i+4*5]=ins->re[i];
    matmul("NN",5,5,5,1.0,Tp,dT,0.0,Tc);

    for (i=0;i<3;i++) for (j=0;j<3;j++) ins->Cbe[i+j*3]=Tc[i+j*5];
    for (i=0;i<3;i++) {
        ins->ve[i]=Tc[i+3*5];
        ins->re[i]=Tc[i+4*5];
    }
    for (i=0;i<3;i++) ins->ba[i]=ins->ba[i]-x[IIBA(opt)+i];
    for (i=0;i<3;i++) ins->bg[i]=ins->bg[i]-x[IIBG(opt)+i];
    insemdl(ins->fb,ins->omgb,fbc,omgbc,ins);
    matcpy(ins->fb,fbc,1,3);
    matcpy(ins->omgb,omgbc,1,3);

    matmul3v("N",ins->Cbe,fbc,fe);
    insgrtp(ins->re,ge);
    matmul3v("N",Omge,ins->ve,cori);
    for (i=0;i<3;i++) ins->ae[i]=fe[i]+ge[i]-2.0*cori[i];

    if (INRL(&ins->opt)) {
        for (i=0;i<3;i++) ins->rl[i]-=x[IIRL(&ins->opt)+i];
    }
    if (INOL(&ins->opt)) {
        for (i=0;i<3;i++) ins->ol[i]-=x[IIOL(&ins->opt)+i];
    }
    if (INSG(&ins->opt)) {
        for (i=IISG(&ins->opt);i<IISG(&ins->opt)+INSG(&ins->opt);i++) ins->Mg[i-IISG(&ins->opt)+(i-IISG(&ins->opt))*3]+=x[i];
    }
    if (INSA(&ins->opt)) {
        for (i=IISA(&ins->opt);i<IISA(&ins->opt)+INSA(&ins->opt);i++) ins->Ma[i-IISA(&ins->opt)+(i-IISA(&ins->opt))*3]+=x[i];
    }
}
/* INS filter fusion with GNSS solution---------------------------------------*/
static int udflt_gnss(ins_t *ins, const sol_t *sol)
{
    int nx=ins->flt.nx,nv,flag=0;
    double *v,*H,*R,*xp,*Pp,dt=timediff(sol->time,ins->time),outage=0.5;
    ins_flt_t *flt=&ins->flt;
    sol_t sol_=*sol;

    flt->stat=0;
    if (fabs(dt)>outage) return -1;
    xp=mat(nx,nx); Pp=mat(nx,nx); v=zeros(1,6); H=mat(6,nx); R=zeros(6,6);

    /* predict INS states */
    matcpy(xp,flt->x,1,nx);
    matcpy(Pp,flt->P,nx,nx);

    /* GNSS solution measurement residuals */
    if (!(nv=gnssres(ins,&sol_,flt->x,flt->P,v,H,R))) {
        free(v); free(H);
        free(R); free(xp); free(Pp);
        return 0;
    }
    /* Kalman filter update */
    if (filter(xp,Pp,H,v,R,nx,nv,NULL,0)) {
        free(v); free(H);
        free(R); free(xp); free(Pp);
        return 0;
    }
    /* INS filter close loop correction */
    udins_uapvse3(ins,xp);

    ins->stat=gs2is(sol->stat);
    ins->flt.time=sol->time;
    ins->flt.stat=1;
    matcpy(ins->flt.x,xp, 1,nx);
    matcpy(ins->flt.P,Pp,nx,nx);

    free(v); free(H); free(R); free(xp); free(Pp);
    return ins->flt.stat;
}
/* update INS using GNSS solution---------------------------------------------*/
extern int udins_gnssse3(ins_t *ins, const sol_t *sol)
{
    if (ins->stat==INSS_NONE) return -1;
    if (udflt_gnss(ins,sol)<0) ins->stat=INSS_NONE;
    return ins->stat;
}
/* update INS solution--------------------------------------------------------*/
static void udins_usol(int udflg, int type, ins_t *ins, const imud_t *imu, const sol_t *solg, sol_t *solo)
{
    if (udflg==UDINS_IMUS) {
        if (ins->stat!=INSS_NONE) ins2sol(type,ins,solo);
    }
    else if (udflg==UDINS_GNSS) {
        if (ins->stat==INSS_NONE) *solo=*solg;
        else {
            *solo=*solg;
            ins2sol(type,ins,solo);
        }
    }
}
/* update INS states forward-----------------------------------------------------
 * args:  ins_t *ins  IO  INS struct data
 *        imud_t *imu I   IMU measurement data
 * return: status (1:ok,0:fail)
 * -----------------------------------------------------------------------------*/
extern int udins_statse3(ins_t *ins, const imud_t *imu)
{
    double dt,dqb[4],qk_1[4],Ck_1[9],qk[4],dqe[4],qtmp[4],da[3]={0},dv[3]={0};
    double domgb[3]={0},domge[3]={0},dvfk[3],dvbk[4],dCe[9];
    double wv[3],Omge[3],ge[3],vek_1[3],rm[3],vm[3];
    int i;

    if (ins->stat==INSS_NONE) return 0;
    if (ins->pins) *ins->pins=*ins;

    if ((dt=timediff(imu->time,ins->time))>5.0||fabs(dt)<1E-6) {
        ins->dt=timediff(imu->time,ins->time);
        ins->time=imu->time;
        return 0;
    }
    matcpy(ins->omgb0,imu->gyro,1,3);
    matcpy(ins->fb0,imu->accl,1,3);
    ins->imu[0]=ins->imu[1];
    ins->imu[1]=*imu;
    insemdl(imu->accl,imu->gyro,ins->fb,ins->omgb,ins);

    /* update attitude */
#if SCULL_CORR
    insscullcorr(ins,&ins->opt,dt,dv,da);
#endif
    for (i=0;i<3;i++) domgb[i]=ins->omgb[i]*dt+da[i];
    rve2quat(domgb,dqb);

    matcpy(Ck_1,ins->Cbe,3,3);
    dcm2quat(ins->Cbe,qk_1);
    domge[2]=-OMGE*dt;
    rve2quat(domge,dqe);
    quatmul(qk_1,dqb,qtmp);
    quatmul(dqe,qtmp,qk);
    normquat(qk);
    quat2dcm(qk,ins->Cbe);

    /* update velocity */
    dvbk[0]=ins->fb[0]*dt+dv[0];
    dvbk[1]=ins->fb[1]*dt+dv[1];
    dvbk[2]=ins->fb[2]*dt+dv[2];
    domge[2]=-OMGE*dt/2.0;
    rve2quat(domge,dqe);
    quat2dcm(dqe,dCe);
    matmul33("NNN",dCe,Ck_1,dvbk,3,3,3,1,dvfk);

    matcpy(vek_1,ins->ve,1,3);

    Omge[0]=0.0; Omge[1]=0.0; Omge[2]=OMGE;
    cross3(Omge,vek_1,wv);

    /* middle point velocity and gravity */
    for (i=0;i<3;i++) rm[i]=ins->re[i]+vek_1[i]*dt/2.0;
    insgrty(rm,ge);
    for (i=0;i<3;i++) {
        vm[i]=vek_1[i]+1.0/2.0*(dvfk[i]+(ge[i]-2.0*wv[i])*dt);
    }
    for (i=0;i<3;i++) rm[i]=ins->re[i]+vm[i]*dt/2.0;
    insgrty(rm,ge);
    cross3(Omge,vm,wv);
    for (i=0;i<3;i++) {
        ins->ae[i]=dvfk[i]/dt+(ge[i]-2.0*wv[i]);
        ins->ve[i]+=ins->ae[i]*dt;
    }
    /* update position */
    for (i=0;i<3;i++) ins->re[i]+=0.5*(vek_1[i]+ins->ve[i])*dt;
    ins->dt=timediff(imu->time,ins->time);
    ins->time=imu->time;
    ins->stat=INSS_TIME;
    return 1;
}
/* interpolate IMU measurement data---------------------------------------------*/
static int imuinterpolate(const imud_t *imu1, const imud_t *imu2, gtime_t time, imud_t *imu)
{
    double thresdt=1.0,v[3],w[3];
    double t0=timediff(imu2->time,imu1->time);
    double t1=timediff(time,imu1->time),t2=timediff(time,imu2->time);
    int i;

    memset(imu,0,sizeof(*imu));
    if (!imu1->time.time) return 0;
    if (!imu2->time.time) return 0;
    if (fabs(t0)>thresdt||fabs(t1)>thresdt||fabs(t2)>thresdt) return 0;

    /* rate of IMU accl/gyro */
    for (i=0;i<3;i++) {
        v[i]=(imu2->accl[i]-imu1->accl[i])/t0;
        w[i]=(imu2->gyro[i]-imu1->gyro[i])/t0;
    }
    /* linear interpolate of IMU accl/gyro measurement */
    for (i=0;i<3;i++) {
        imu->accl[i]=imu2->accl[i]+v[i]*t2;
        imu->gyro[i]=imu2->gyro[i]+w[i]*t2;
    }
    imu->time=time;
    return 1;
}
extern void udins_imui(ins_t *ins, gtime_t time, imud_t *imu)
{
    if (ins->stat==INSS_NONE) return;
    imuinterpolate(&ins->imu[0],&ins->imu[1],time,imu);
}
/* update GNSS antenna velocity-----------------------------------------------*/
static int sol2vel(const sol_t *sol, double *vr)
{
    static sol_t sols[MIN_CNT];
    double dt,vt[MIN_CNT][3]={0},vthres=5.0,pqthres=30.0,vqthres=3.0;
    double thresdt=3.0;
    int i,j;

    /* check solutions status */
    if (sol->stat==SOLQ_NONE) return 0;
    if (fabs(timediff(sol->time,sols[MIN_CNT-1].time))<1E-2) return 0;
    if (sqrt(sol->qr[0]+sol->qr[1]+sol->qr[2])>pqthres||
        sqrt(sol->qv[0]+sol->qv[1]+sol->qv[2])>vqthres) {
        return 0;
    }
    for (i=0;i<MIN_CNT-1;i++) sols[i]=sols[i+1];
    sols[i]=*sol;

    /* update velocity by position difference */
    for (i=0;i<MIN_CNT-1;i++) {
        if ((dt=timediff(sols[i+1].time,sols[i].time))>thresdt||fabs(dt)<1E-2) continue;
        for (j=0;j<3;j++) {
            vt[i][j]=(sols[i+1].rr[j]-sols[i].rr[j])/dt;
        }
    }
    /* initial INS only velocity greater than vthres */
    for (i=0;i<MIN_CNT-1;i++) {
        if (norm(vt[i],3)<=0.0||norm(vt[i],3)<vthres) return 0;
    }
    matcpy(vr,vt[MIN_CNT-2],1,3);
    return 1;
}
/* check INS initial----------------------------------------------------------*/
static int isinitins(ins_t *ins, sol_t *sols, const imud_t *imu)
{
    double vr[3],thresdt=0.1;

    if (ins->stat!=INSS_NONE) return 0;
    if (fabs(timediff(imu->time,sols->time))>thresdt) return 0;
    if (!sol2vel(sols,vr)) return 0;
    if (norm(sols->rr+3,3)<=0) matcpy(sols->rr+3,vr,1,3);
    return 1;
}
/* initial INS states---------------------------------------------------------*/
static void udins_is(ins_t *ins, const sol_t *sols, const imud_t *imu, const double dt)
{
    double llh[3],C[9],rpy[3]={0},vn[3],Cbn[9],T[9],TT[9],rr[6];
    double Omge[9]={0,OMGE,0,-OMGE,0,0,0,0,0},gyro[3]={0};
    int i;

    for (i=0;i<3;i++) rr[i]=sols->rr[i]+sols->rr[i+3]*dt;
    for (i=3;i<6;i++) rr[i]=sols->rr[i];

    ecef2pos(rr,llh); ned2xyz(llh,C);
    matmul("TN",3,1,3,1.0,C,rr+3,0.0,vn);

    rpy[1]=atan2(-vn[2],sqrt(SQR(vn[0]))+SQR(vn[1]));
    rpy[2]=vel2head(vn);
    rpy2dcm(rpy,C); matt(C,3,3,Cbn);
    ned2xyz(llh,C);
    matmul("NN",3,3,3,1.0,C,Cbn,0.0,ins->Cbe);

    matcpy(ins->re,rr,3,1);
    matmul("NN",3,1,3,-1.0,ins->Cbe,ins->rl,1.0,ins->re);

    skewsym3(gyro,T);
    matmul("NN",3,1,3,1.0,T,ins->rl,0.0,TT);
    matmul("NN",3,1,3,1.0,ins->Cbe,TT,0.0,T);

    matmul33("NNN",Omge,ins->Cbe,ins->rl,3,3,3,1,TT);
    ins->ve[0]=rr[3]-T[0]+TT[0];
    ins->ve[1]=rr[4]-T[1]+TT[1];
    ins->ve[2]=rr[5]-T[2]+TT[2];
}
/* initial INS states/filter--------------------------------------------------*/
extern int udins_init(ins_t *ins, const sol_t *sols, const imud_t *data)
{
    sol_t sol_=*sols;

    /* check whether initial INS */
    if (!isinitins(ins,&sol_,data)) {
        return 0;
    }
    /* initial INS states and covariance */
    udins_is(ins,&sol_,data,timediff(data->time,sol_.time));
    udins_P0(ins,data->time,&ins->opt,&sol_);
    ins->stat=INSS_INIT;
    ins->time=sol_.time;
    ins->flt.time=sol_.time;
    return 1;
}
/* INS positioning in SE(2)3---------------------------------------------------
 * args:    int udflg    I   INS update flag (see UDINS_xxxx)
 *          ins_t *ins   IO  INS control/result struct
 *          imud_d *imu  I   IMU measurement data
 *          sol_t *sol   I   GNSS/attitude solution data
 *          nav_t *nav   I   navigation data
 *          obsd_t *obs  I   GNSS observation data
 *          int n        I   number of GNSS observation data
 *          odom_t *odo  I   odometer measurement data
 * return : status (0:no solution,1:valid solution)
 * ---------------------------------------------------------------------------*/
extern int insposse3(int udflg, ins_t *ins, const imud_t *imu, const sol_t *sol)
{
    imud_t imut=*imu;
    switch (udflg) {
        case UDINS_IMUS:
            udins_init(ins,&ins->sol,&imut);
            udins_statse3(ins,&imut);
            udins_ufltse3(ins,&imut);
            break;

        case UDINS_GNSS:
            udins_imui(ins,sol->time,&imut);
            udins_statse3(ins,&imut);
            udins_ufltse3(ins,&imut);
            udins_gnssse3(ins,sol);
            break;
    }
    udins_usol(udflg,1,ins,imu,sol,&ins->sol);
    return ins->stat;
}


