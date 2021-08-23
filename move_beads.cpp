# include <math.h>
# include "slink.h"
# include <map>

void move_beads(Bead *bead, SlipLink *sl, double h, double ksi, double kBT, long int *in_ran1, int n_ch, int n_m, int istep)
{

int i,j,l,txjj;
double w1,w2,lw1,du1,u1,Fsum;
double cte2,sigma1;
int bead_index;
FILE *beadpos;
sigma1=sqrt(2.0*h*kBT/ksi);
cte2=1.0/ksi;
int bead_size = n_m*n_ch;
//beadpos=fopen("beadpos.dat","w");
//fprintf(stdout,"istep=%d\n",istep);
for (i=0;i<bead_size;i++) {
 	bead_index = i; //g_bead_index[i];
        for (l=0;l<3;l++) {
        w1=ran1(in_ran1);
        w2=ran1(in_ran1);
        lw1=log(w1);
        du1=sqrt(-2.0*lw1);
        u1=du1*cos(2.0*M_PI*w2);


        /* sum of forces multiplied by 1/ksi=cte2 */
        Fsum=(bead[bead_index].F_SP[l]+bead[bead_index].F_SL[l])*cte2;
        /* Temporary update of the beads coordinates */
        bead[bead_index].XYZ[l]=bead[bead_index].XYZold[l]+Fsum*h+sigma1*u1;
        }
}
}
