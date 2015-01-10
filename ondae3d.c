#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dissipation.h"
/** pi  = dphi/dt*/
/** psi = dphi/dr*/
/**compilar con "gcc -o archivo archivo.c -lm"*/
int main()
{   
    /**Definimos variables*/
    float t,dx,dt, v,a0,x0,s0,idx,diss;
    int nx,nt,z,n;
    FILE *fp,*fp1,*fp2;
    /**ancho de la malla*/
    printf("Introduzca el ancho de la malla\n");
    scanf("%d",&nx);
    /**nx=30*/
    /**Matrices que daran informacion de las ecuaciones diferenciales*/
    float X[nx+2],PI[nx+2], PHI[nx+2],PSI[nx+2],DPSI[nx+2],DPI[nx+2],SPSI[nx+2],SPI[nx+2],SPHI[nx+2],OPHI[nx+2],OPI[nx+2],OPSI[nx+2];

    float PI_K1[nx+2], PHI_K1[nx+2],PSI_K1[nx+2];

    /**Introducimos valores de avance espacial y temporal, la velocidad de la onda, y los parametros iniciales a0,x0,s0 de la gaussiana*/
    printf("Introduzca la division espacial\n");
    scanf("%f",&dx);
    /**dx=0.05;*/
    printf("Introduzca la division temporal\n");
    scanf("%f",&dt);
    /**dt=0.05;*/
    printf("Introduzca la velocidad\n");
    scanf("%f",&v);
    /**v=0.1;*/
    printf("Introduzca el numero total de pasos temporales\n");
    scanf("%d",&nt);
    /**nt=15;*/
    a0 = 1.0;
    x0 = (nx/2.0)*dx;
    s0 = 1.0;
    t=0.;
    idx=1.0/(2*dx);
    diss = 0.01;
    /**Fijamos valores iniciales*/
    for(z=2;z<nx+2;z++){
        X[z]=(z)*dx;
    }
    /*SIMETRIA*/
    for(z=0;z<2;z++){
        X[z]=(z-2)*dx;
    }

    for(z=2;z<nx+2;z++){
        PHI[z]=a0*exp((-pow((X[z]-x0),2))/pow(s0,2));
        PSI[z]=-2.*a0*((X[z]-x0)/pow(s0,2))*exp(-pow(X[z]-x0,2)/pow(s0,2));
        PI[z]=0;
    }
    /**Simetria**/
    for(z=0;z<2;z++){
        PHI[z]=PHI[z+1];
        PSI[z]=-PSI[z+1];
        PI[z]=PI[z+1];
    }
    /*Creamos documentos de escritura*/
    fp=fopen("salidaphi3.txt","wt");
    fp1=fopen("salidapsi3.txt","wt");
    fp2=fopen("salidapi3.txt","wt");
    /*Se escribe el resultado en el documento correspondiente*/
    for (z=0;z<nx+2;z++){
      fprintf(fp, "%f %f",X[z],PHI[z]);
      fprintf(fp,"\n");
    }
      fprintf(fp,"\n\n");
    for (z=0;z<nx+2;z++){
      fprintf(fp1, "%f %f",X[z],PSI[z]);
      fprintf(fp1,"\n");
    }
      fprintf(fp1,"\n\n");
    for (z=0;z<nx+2;z++){
      fprintf(fp2, "%f %f",X[z],PI[z]);
      fprintf(fp2,"\n");
    }
      fprintf(fp2,"\n\n");
    

/**Loop principal*/
    for(z=0;z<nt;z++){
//Avance temporal
      t=t+dt;
//Salvamos el paso temporal anterior
      for(n=0;n<nx+2;n++){ 
	OPHI[n]=PHI[n];
	OPSI[n]=PSI[n];
	OPI[n]=PI[n];
      }
      //Calculando puntos nuevos para un tiempo t 
      //ICN de tres pasos
      //Primer paso (medio intervalo)
      for(n=2;n<nx;n++){
	//Primero calculamos las derivadas
	DPSI[n]=(PSI[n+1]-PSI[n-1])*idx;
	DPI[n]=(PI[n+1]-PI[n-1])*idx;
      }
      for(n=2;n<nx+1;n++){
	//Evaluamos las fuentes
	SPHI[n]=PI[n];
	SPSI[n]=DPI[n];
	SPI[n]=(pow(v,2))*(DPSI[n]+(2./X[n])*PSI[n]);
      }
      //frontera derecha
      SPI[nx-1]=-(v*idx)*(3*PI[nx-1]-4*PI[nx-2]+PI[nx-3])-v*PI[nx-1]/X[nx-1];
      SPHI[nx-1]=-(v*idx)*(3*PHI[nx-1]-4*PHI[nx-2]+PHI[nx-3])-v*PHI[nx-1]/X[nx-1];
      SPSI[nx-1]=-(v*idx)*(3*PSI[nx-1]-4*PSI[nx-2]+PSI[nx-3])-v*PSI[nx-1]/X[nx-1];
      //Disipacion
      //dissipation(nx,dx,diss,PHI,SPHI,+1,PSI,SPSI,-1,PI,SPI,+1);
      //Actualizamos
      for(n=2;n<nx;n++){
	PHI[n]=OPHI[n]+0.5*dt*SPHI[n];
	PSI[n]=OPSI[n]+0.5*dt*SPSI[n];
	PI[n]=OPI[n]+0.5*dt*SPI[n];
      }
      //Frontera izquierda, aplicamos simetria para PI y PHI y antisimetria para PSI
      for(n=0;n<2;n++){
          PHI[n]=PHI[n+1];
          PSI[n]=-PSI[n+1];
          PI[n]=PI[n+1];
      }
      //segundo paso (medio intervalo)
      for(n=2;n<nx;n++){
	DPSI[n]=(PSI[n+1]-PSI[n-1])*idx;
	DPI[n]=(PI[n+1]-PI[n-1])*idx;
      }
      for(n=2;n<nx;n++){
	SPHI[n]=PI[n];
	SPSI[n]=DPI[n];
	SPI[n]=(pow(v,2))*(DPSI[n] +(2./X[n])*PSI[n]);
      }
      //frontera derecha
      SPI[nx-1]=-(v*idx)*(3*PI[nx-1]-4*PI[nx-2]+PI[nx-3])-v*PI[nx-1]/X[nx-1];
      SPHI[nx-1]=-(v*idx)*(3*PHI[nx-1]-4*PHI[nx-2]+PHI[nx-3])-v*PHI[nx-1]/X[nx-1];
      SPSI[nx-1]=-(v*idx)*(3*PSI[nx-1]-4*PSI[nx-2]+PSI[nx-3])-v*PSI[nx-1]/X[nx-1];
      //dissipation(nx,dx,diss,PHI,SPHI,+1,PSI,SPSI,-1,PI,SPI,+1);
      for(n=2;n<nx;n++){
	PHI[n]=OPHI[n]+0.5*dt*SPHI[n];
	PSI[n]=OPSI[n]+0.5*dt*SPSI[n];
	PI[n]=OPI[n]+0.5*dt*SPI[n];
      }
      //Frontera izquierda, simetria
      for(n=0;n<2;n++){
          PHI[n]=PHI[n+1];
          PSI[n]=-PSI[n+1];
          PI[n]=PI[n+1];
      }
      //Tercer paso (completo)
      for(n=2;n<nx;n++){
	DPSI[n]=(PSI[n+1]-PSI[n-1])*idx;
	DPI[n]=(PI[n+1]-PI[n-1])*idx;
      }
      for(n=2;n<nx;n++){
	SPHI[n]=PI[n];
	SPSI[n]=DPI[n];
	//SPI[n]=(pow(v,2))*DPSI[n];
	SPI[n]=(pow(v,2))*(DPSI[n] +(2./X[n])*PSI[n]);
      }
      //frontera derecha
      SPI[nx-1]=-(v*idx)*(3*PI[nx-1]-4*PI[nx-2]+PI[nx-3])-v*PI[nx-1]/X[nx-1];
      SPHI[nx-1]=-(v*idx)*(3*PHI[nx-1]-4*PHI[nx-2]+PHI[nx-3])-v*PHI[nx-1]/X[nx-1];
      SPSI[nx-1]=-(v*idx)*(3*PSI[nx-1]-4*PSI[nx-2]+PSI[nx-3])-v*PSI[nx-1]/X[nx-1];
      //dissipation(nx,dx,diss,PHI,SPHI,+1,PSI,SPSI,-1,PI,SPI,+1);
      for(n=2;n<nx;n++){
	PHI[n]=OPHI[n]+0.5*dt*SPHI[n];
	PSI[n]=OPSI[n]+0.5*dt*SPSI[n];
	PI[n]=OPI[n]+0.5*dt*SPI[n];
      }

      //Frontera izquierda
      for(n=0;n<2;n++){
          PHI[n]=PHI[n+1];
          PSI[n]=-PSI[n+1];
          PI[n]=PI[n+1];
      }
    
    //Se escribe el resultado en el documento correspondiente
    for(n=0;n<nx+1;n++){
      fprintf(fp, "%f %f",X[n],PHI[n]);
      fprintf(fp,"\n");
    }
      fprintf(fp,"\n\n");
    for(n=0;n<nx+1;n++){
      fprintf(fp1, "%f %f",X[n],PSI[n]);
      fprintf(fp1,"\n");
    }
      fprintf(fp1,"\n\n");
    for(n=0;n<nx+1;n++){
      fprintf(fp2, "%f %f",X[n],PI[n]);
      fprintf(fp2,"\n");
    }
      fprintf(fp2,"\n\n");
    }
    //Clean exit
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
}
