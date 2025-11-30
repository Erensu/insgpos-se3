#include "insgpos-se3.h"

int main(int argc, const char *argv[])
{
    const char *rawfile="..\\test_data";
    FILE *fp;

    if (!(fp=fopen(rawfile,"r"))) return 0;

    FILE *fp_sol=fopen(".\\sol.out","w");
    ins_t ins={0};
    ins_opt_t iopt={0};

    iopt.rl[0]=-0.6;
    iopt.rl[1]=0.001;
    iopt.rl[2]=-1.35;
    iopt.psd[0]=SQR(0.25*D2R/sqrt(3600.0));
    iopt.psd[1]=SQR(2.5/sqrt(3600.0));
    iopt.psd[2]=1E-9;
    iopt.psd[3]=1E-13;
    iopt.psd[4]=1E-8;
    iopt.psd[5]=1E-8;
    iopt.psd[6]=1E-8;
    iopt.estrl=0;
    iopt.estsa=0;
    iopt.estsg=0;

    insinit(&ins,&iopt);

    while (!feof(fp)) {
        char buff[1024];
        if (!fgets(buff,sizeof(buff),fp)) break;

        double sow,week,val[6];
        int type;
        sol_t sol={0};
        imud_t imu={0};
        if (sscanf(buff,"%lf %lf %d %lf %lf %lf %lf %lf %lf\n",&week,&sow,&type,&val[0],&val[1],&val[2],&val[3],&val[4],&val[5])<9) continue;

        if (type==1) {
            imu.time=gpst2time(week,sow);
            matcpy(imu.accl,val,1,3);
            matcpy(imu.gyro,val+3,1,3);
            insposse3(UDINS_IMUS,&ins,&imu,&sol);
        }
        else if (type==2) {
            sol.time=gpst2time(week,sow);
            sol.stat=SOLQ_FIX;
            matcpy(sol.rr,val,1,3);
            sol.qr[0]=SQR(val[3]);
            sol.qr[1]=SQR(val[4]);
            sol.qr[2]=SQR(val[5]);
            insposse3(UDINS_GNSS,&ins,&imu,&sol);
        }
        if (ins.stat!=INSS_NONE) {
            outisol(fp_sol,&ins.sol);
            fflush(fp_sol);
        }
        else if (ins.stat==INSS_NONE) {
            if (type==2) {
                outisol(fp_sol,&sol);
                fflush(fp_sol);
            }
        }
    }
    insfree(&ins);
    fclose(fp);
    fclose(fp_sol);
    return 1;
}