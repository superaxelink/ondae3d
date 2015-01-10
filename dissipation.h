dissipation(int Nx, float dx, float diss, float var[], float svar[], int sym,float var1[], float svar1[], int sym1,float var2[], float svar2[], int sym2)
{

    int i;
    float idx;
    idx= 1.0/dx;
    
//! ***   SECOND ORDER   ***

//!    Interior points:  Second order evolution
//!    requires fourth order dissipation.

    for(i=1;i<Nx-1;i++){
        svar[i] = svar[i]-diss*idx*(var[i+2]-4.0*var[i+1]+6.0*var[i]-4.0*var[i-1]+var[i-2]);
        svar1[i] = svar1[i]-diss*idx*(var1[i+2]-4.0*var1[i+1]+6.0*var1[i]-4.0*var1[i-1]+var1[i-2]);
        svar2[i] = svar2[i]-diss*idx*(var2[i+2]-4.0*var2[i+1]+6.0*var2[i]-4.0*var2[i-1]+var[i-2]);
    }
// Ghost points using symmetries.
     for(i=0;i<2;i++){
        svar[i-1] = sym*svar[i];
        svar1[i-1] = sym1*svar1[i];
        svar2[i-1] = sym2*svar2[i];
     }
}
