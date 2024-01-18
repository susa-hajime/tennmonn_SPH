/*----------------------------------------*/
/*　光子消失問題  (非保存)          　　      */
/*----------------------------------------*/
#include<stdio.h>
#include <stdlib.h>
#include<math.h>

#define num 100000
#define length 100.
#define tend 3.
#define nout 2
#define kedge 50.

/* alpha = 100 or 1 */
#define alpha 100.
//#define alpha 1.

int main(void)
{
  int i,ngrid,count;
  double y[num],k[num],x[num],tau[num],dy[num],ya[num],dx,t,dt,dtout,tout;
  FILE *output_file;
  output_file=fopen("output_simple_100.dat","w"); /* 出力ファイルオープン*/
//  output_file=fopen("output_simple_1.dat","w"); /* 出力ファイルオープン*/

  dtout = tend/(double)nout;
  tout = 0.;

  /* ----------------  初期条件 ----------------- */
  dx = 0.1; 
  ngrid = (int)(length/dx + 0.1);
  
  t=0; dt=1e-4;
  for(i=0; i < num ; i++){
    y[i] = 1.;
    k[i] = 0.;
    x[i] = dx*(double)i;
  }
  k[0]=kedge;

  /* -------------- 初期条件終わり ------------- */

  /* -------- 時間積分 --------- */

  count = 0.;
  while(t < tend){

  if( count == 1 ){
    tout =tout + dtout;
    count=0;
    }
 
  tau[0]=0.;
  for(i=0; i < ngrid ; i++){
    tau[i+1]=tau[i] + 0.5*(y[i]+y[i+1])*(x[i+1]-x[i])*alpha;
    k[i]=kedge*exp(-tau[i]);
  }
  for(i=0; i < ngrid ; i++){
    dy[i] = dt * ( -k[i]*y[i]*alpha + pow(1.-y[i],2) );
  }

  for(i=0; i < ngrid ; i++){
    y[i] = y[i]+dy[i];
    if(x[i] > kedge*(1.-exp(-t))){
      ya[i]=1.;
      }
      else{
      ya[i]=0.;
    }
  }

    t=t+dt;

    if( t > tout ){
      count = 1;
      fprintf(output_file,"# t= %f\n",t);            
      for(i=0; i < ngrid ; i++){
        fprintf(output_file,"%f %f %f %f %f\n",(float)x[i],(float)y[i],(float)k[i],(float)tau[i],(float)ya[i]);
      }
    }
  }

}
