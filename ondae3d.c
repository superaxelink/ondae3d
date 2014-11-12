#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/** pi  = dphi/dt*/
/** psi = dphi/dx*/
/**compilar con "gcc -o archivo archivo.c -lm"*/
int main()
{   
    /**Definimos variables*/
    float t,dx,dt, v,a0,x0,s0;
    int nx,nt,z,n;
    FILE *fp,*fp1,*fp2;
    /**ancho de la malla*/
    printf("Introduzca el ancho de la malla\n");
    scanf("%d",&nx);
    /**nx=30*/
    /**Matrices que daran informacion de las ecuaciones diferenciales*/
    float R[nx],PI[nx], PHI[nx],PSI[nx],DPSI[nx],DPI[nx],SPSI[nx],SPI[nx],SPHI[nx],OPHI[nx],OPI[nx],OPSI[nx];

    float PI_K1[nx], PHI_K1[nx],PSI_K1[nx];

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
    t=0;
/**Fijamos valores iniciales*/
    for(z=0;z<nx;z++){
        R[z]=(z+0.5)*dx;
    }
    for(z=0;z<nx;z++){
        PHI[z] = a0*exp((-pow((R[z]-x0),2))/pow(s0,2));
        PI[z]=0.0;
    }
    for(z=0;z<nx;z++){
        PSI[z]=-2.*a0*((R[z]-x0)/pow(s0,2))*exp(-pow(R[z]-x0,2)/pow(s0,2));
    }
    /**Creamos documentos de escritura*/
    fp=fopen("salidaphi3.txt","wt");
    fp1=fopen("salidapsi3.txt","wt");
    fp2=fopen("salidapi3.txt","wt");
    /**Se escribe el resultado en el documento correspondiente*/
    for (z=0;z<nx;z++){
      fprintf(fp, "%f %f",R[z],PHI[z]);
      fprintf(fp,"\n");
    }
      fprintf(fp,"\n\n");
    for (z=0;z<nx;z++){
      fprintf(fp1, "%f %f",R[z],PSI[z]);
      fprintf(fp1,"\n");
    }
      fprintf(fp1,"\n\n");
    for (z=0;z<nx;z++){
      fprintf(fp2, "%f %f",R[z],PI[z]);
      fprintf(fp2,"\n");
    }
      fprintf(fp2,"\n\n");
    

/**Loop principal*/
    for(z=0;z<nt;z++){
/**Avance temporal*/
      t=t+dt;
/**Salvamos el paso temporal anterior*/
      for(n=0;n<nx;n++){ 
	OPHI[n]=PHI[n];
	OPSI[n]=PSI[n];
	OPI[n]=PI[n];
      }

      /**Calculando puntos nuevos para un tiempo t*/ 

      //ICN de tres pasos
      
      //Primero
      for(n=1;n<nx-1;n++){
	/*Primero calculamos las derivadas*/
	DPSI[n]=(8*PSI[n+1]+PSI[n-2]-PSI[n+2]-8*PSI[n-1])/(12.*dx);
	DPI[n]=(8*PI[n+1]+PI[n-2]-PI[n+2]-8*PI[n-1])/(12.*dx);
      }

      for(n=1;n<nx-1;n++){
	/**Evaluamos las fuentes*/
	SPHI[n]=PI[n];
	SPSI[n]=DPI[n];
	SPI[n]=(pow(v,2))*DPSI[n];
        /*SPI[n]=(pow(v,2))*(DPSI[n]+(2/R[n])*PSI[n]);*/

	//Primer paso (medio intervalo)
	PHI[n]=OPHI[n]+0.5*dt*SPHI[n];
	PSI[n]=OPSI[n]+0.5*dt*SPSI[n];
	PI[n]=OPI[n]+0.5*dt*SPI[n];
      }

      /**frontera derecha*/
      SPI[nx-1]=-v/dx*(PI[nx-1]-PI[nx-2]);//0.;
      SPHI[nx-1]=-v/dx*(PHI[nx-1]-PHI[nx-2]);
      SPSI[nx-1]=-v/dx*(PSI[nx-1]-PSI[nx-2]);

      PHI[nx-1]=OPHI[nx-1]+0.5*dt*SPHI[nx-1];
      PSI[nx-1]=OPSI[nx-1]+0.5*dt*SPSI[nx-1];
      PI[nx-1]=OPI[nx-1]+0.5*dt*SPI[nx-1];

      /**Frontera izquierda*/
      SPI[0]=v/dx*(PI[1]-PI[0]);//0.;
      SPHI[0]=v/dx*(PHI[1]-PHI[0]);
      SPSI[0]=v/dx*(PSI[1]-PSI[0]);

      PHI[0]=OPHI[0]+0.5*dt*SPHI[0];
      PSI[0]=OPSI[0]+0.5*dt*SPSI[0];
      PI[0]=OPI[0]+0.5*dt*SPI[0];



      for(n=1;n<nx-1;n++){
	/*Primero calculamos las derivadas*/
	DPSI[n]=(8*PSI[n+1]+PSI[n-2]-PSI[n+2]-8*PSI[n-1])/(12.*dx);
	DPI[n]=(8*PI[n+1]+PI[n-2]-PI[n+2]-8*PI[n-1])/(12.*dx);
      }

      for(n=1;n<nx-1;n++){
	/**Evaluamos las fuentes*/
	SPHI[n]=PI[n];
	SPSI[n]=DPI[n];
        SPI[n]=(pow(v,2))*DPSI[n];
        /*SPI[n]=(pow(v,2))*(DPSI[n]+(2/R[n])*PSI[n]);*/
	//Agrego pasos intermedios para hacerlo estable
	// y aumentar el orden

	//Segundo paso (medio intervalo)
	PHI[n]=OPHI[n]+0.5*dt*SPHI[n];
	PSI[n]=OPSI[n]+0.5*dt*SPSI[n];
	PI[n]=OPI[n]+0.5*dt*SPI[n];
      }

        /**frontera derecha*/
      SPI[nx-1]=-v/dx*(PI[nx-1]-PI[nx-2]);//0.;
      SPHI[nx-1]=-v/dx*(PHI[nx-1]-PHI[nx-2]);
      SPSI[nx-1]=-v/dx*(PSI[nx-1]-PSI[nx-2]);

      PHI[nx-1]=OPHI[nx-1]+0.5*dt*SPHI[nx-1];
      PSI[nx-1]=OPSI[nx-1]+0.5*dt*SPSI[nx-1];
      PI[nx-1]=OPI[nx-1]+0.5*dt*SPI[nx-1];

      /**Frontera izquierda*/
      SPI[0]=v/dx*(PI[1]-PI[0]);//0.;
      SPHI[0]=v/dx*(PHI[1]-PHI[0]);
      SPSI[0]=v/dx*(PSI[1]-PSI[0]);

      PHI[0]=OPHI[0]+0.5*dt*SPHI[0];
      PSI[0]=OPSI[0]+0.5*dt*SPSI[0];
      PI[0]=OPI[0]+0.5*dt*SPI[0];

      for(n=1;n<nx-1;n++){
	/*Primero calculamos las derivadas*/
	DPSI[n]=(8*PSI[n+1]+PSI[n-2]-PSI[n+2]-8*PSI[n-1])/(12.*dx);
	DPI[n]=(8*PI[n+1]+PI[n-2]-PI[n+2]-8*PI[n-1])/(12.*dx);
      }

      for(n=1;n<nx-1;n++){
	/**Evaluamos las fuentes*/
	SPHI[n]=PI[n];
	SPSI[n]=DPI[n];
	SPI[n]=(pow(v,2))*DPSI[n];
        /*SPI[n]=(pow(v,2))*(DPSI[n]+(2/R[n])*PSI[n]);*/

	//Agrego pasos intermedios para hacerlo estable
	// y aumentar el orden

	//Tercer paso (intervalo completo)
	PHI[n]=OPHI[n]+dt*SPHI[n];
	PSI[n]=OPSI[n]+dt*SPSI[n];
	PI[n]=OPI[n]+dt*SPI[n];
      }

        /**frontera derecha*/
      SPI[nx-1]=-v/dx*(PI[nx-1]-PI[nx-2]);//0.;
      SPHI[nx-1]=-v/dx*(PHI[nx-1]-PHI[nx-2]);
      SPSI[nx-1]=-v/dx*(PSI[nx-1]-PSI[nx-2]);

      PHI[nx-1]=OPHI[nx-1]+dt*SPHI[nx-1];
      PSI[nx-1]=OPSI[nx-1]+dt*SPSI[nx-1];
      PI[nx-1]=OPI[nx-1]+dt*SPI[nx-1];

        /**Frontera izquierda*/
      SPI[0]=v/dx*(PI[1]-PI[0]);//0.;
      SPHI[0]=v/dx*(PHI[1]-PHI[0]);
      SPSI[0]=v/dx*(PSI[1]-PSI[0]);
      PHI[0]=OPHI[0]+dt*SPHI[0];
      PSI[0]=OPSI[0]+dt*SPSI[0];
      PI[0]=OPI[0]+dt*SPI[0];

    /**Se escribe el resultado en el documento correspondiente*/
    for (n=0;n<nx;n++){
      fprintf(fp, "%f %f",R[n],PHI[n]);
      fprintf(fp,"\n");
    }
      fprintf(fp,"\n\n");
    for (n=0;n<nx;n++){
      fprintf(fp1, "%f %f",R[n],PSI[n]);
      fprintf(fp1,"\n");
    }
      fprintf(fp1,"\n\n");
    for (n=0;n<nx;n++){
      fprintf(fp2, "%f %f",R[n],PI[n]);
      fprintf(fp2,"\n");
    }
      fprintf(fp2,"\n\n");


    }

    //Clean exit

    fclose(fp);
    fclose(fp1);
    fclose(fp2);

}
