using namespace std;
typedef  double Real;
#include <iostream>
#include <string.h>
#include <math.h>
#include "mat_vet.h"
#include "Runge_Kutta.h"
void azzera(int ns, Real *A, Real *b, Real *c, Real *bc=0, int *Ind = 0);
void calcola_c(int ns, Real *c, Real *A);

// vuole in ingresso:
// primo vettore che sommo, passo h, vettorer dei coeff della comb lineare,matrice K,numnero degli stadi, lunghezza di ogni vettore, vettore di uscita(v e y possono essere lo stesso valore)
void step(Real *v, Real h, Real *coef, Real *K, int stadi, int d, Real *y){
    Real *aux = new Real [d];
    matmat(coef,K,aux,1,stadi,d);

    // se il vettore Un è nullo non sommare niente
    if (v==0){
        for(int ii=0;ii<d;ii++) y[ii]=h*aux[ii];
    }
    else for(int ii=0;ii<d;ii++) y[ii] =v[ii]+ h*aux[ii];

    delete [] aux;
}

// azzera i coefficenti della mat di butcher
void azzera(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind){
    for(int i=0;i<ns*ns;i++)A[i]=0.0;
    for(int i=0;i<ns;i++)b[i]=0.0;
    for(int i=0;i<ns;i++)c[i]=0.0;
    if(bc!=0)for(int i=0;i<ns;i++)bc[i]=0.0;
    if(Ind!=0)for(int i=0;i<ns;i++)Ind[i]=1;
}

void calcola_c(int ns, Real *c, Real *A){
// calcolo c come somma sulle righe di A
    for(int k=0; k<ns; k++)
    {
        Real sum=*(A+k*ns);
        for (int j=1; j<ns; j++)
        {
            sum+=*(A+ns*k+j);
        }
        c[k]=sum;
        //      cout << c[k]<<endl;
    }
}

// bc in realta serve solo per i metodi embedded e Ind serve per i metodi di
// dirk


// METODI RK ESPLICITI

void EE(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Eulero Esplicito  p =1, ns = 1;
b[0]=1.0;
calcola_c(ns,c,A);
}

void RK4(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Runge-Kutta 4  p =4, ns = 4;
b[0]=b[3]=1.0/6.0;
b[1]=b[2]=1.0/3.0;
*(A+4)=*(A+9)=0.5;
*(A+14)=1.0;
calcola_c(ns,c,A);
}

void HEUN(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Heun p =2, ns = 2;
b[0]=b[1]=0.5;
*(A+2)=1.0;
calcola_c(ns,c,A);
}

void EMOD(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Eulero Modificato  p =2, ns = 2;
b[1]=1.0;
*(A+2)=0.5;
calcola_c(ns,c,A);
}

void RKFEHL5(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind){
    azzera(ns,A,b,c,bc,Ind);
    // RK5  p =5, ns = 6;
    b[0]=16.0/135.0;
    b[1]=0.0;
    b[2]=6656.0/12825.0;
    b[3]=28561.0/56430.0;
    b[4]=-9.0/50.0;
    b[5]=2.0/55.0;
    //riga 2
    *(A+6)=0.25;
    //riga 3
    *(A+12)=3.0/32.0;
    *(A+13)=9.0/32.0;
    //riga 4
    *(A+18)=1932.0/2197.0;
    *(A+19)=-7200.0/2197.0;
    *(A+20)=7296.0/2197.0;
    //riga 5
    *(A+24)=439.0/216.0;
    *(A+25)=-8.0;
    *(A+26)=3680.0/513.0;
    *(A+27)=-845.0/4104.0;
    //riga 6
    *(A+30)=-8.0/27.0;
    *(A+31)=2.0;
    *(A+32)=-3544.0/2565.0;
    *(A+33)=1859.0/4104.0;
    *(A+34)=-11.0/40.0;
calcola_c(ns,c,A);

}

// METODI RK IMMERSI

void RKFEHL54(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// RK Fehlberg 54  p = 5, ns = 6;
// esplicito immerso
        b[0]=16.0/135.0;
        b[1]=0.0;
        b[2]=6656.0/12825.0;
        b[3]=28561.0/56430.0;
        b[4]=-9.0/50.0;
        b[5]=2.0/55.0;
        bc[0]=25.0/216.0;
        bc[1]=0.0;
        bc[2]=1408.0/2565.0;
        bc[3]=2197.0/4104.0;
        bc[4]=-1.0/5.0;
        bc[5]=0.0;
        //riga 2
        *(A+6)=0.25;
        //riga 3
        *(A+12)=3.0/32.0;
        *(A+13)=9.0/32.0;
        //riga 4
        *(A+18)=1932.0/2197.0;
        *(A+19)=-7200.0/2197.0;
        *(A+20)=7296.0/2197.0;
        //riga 5
        *(A+24)=439.0/216.0;
        *(A+25)=-8.0;
        *(A+26)=3680.0/513.0;
        *(A+27)=-845.0/4104.0;
        //riga 6
        *(A+30)=-8.0/27.0;
        *(A+31)=2.0;
        *(A+32)=-3544.0/2565.0;
        *(A+33)=1859.0/4104.0;
        *(A+34)=-11.0/40.0;
calcola_c(ns,c,A);
}

void EEHEUN(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Eulero Esplicito - Heun  p = 2, ns = 2;
// esplicito immerso
        b[0]=0.5;
        b[1]=0.5;
        bc[0]=1.0;
        bc[1]=0.0;
        //riga 2
        *(A+2)=1.0;
calcola_c(ns,c,A);
}

void DP87(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Dormand Prince 87  p = 8, ns = 13;
// esplicito immerso
        b[0]=14005451.0/335480064.0;
        b[5]=-59238493.0/1068277825.0;
        b[6]=181606767.0/758867731.0;
        b[7]=561292985.0/797845732.0;
        b[8]=-1041891430.0/1371343529.0;
        b[9]=760417239.0/1151165299.0;
        b[10]=118820643.0/751138087.0;
        b[11]=-528747749.0/2220607170.0;
        b[12]=0.25;
        bc[0]=13451932.0/455176623.0;
        bc[5]=-808719846.0/976000145.0;
        bc[6]=1757004468.0/5645159321.0;
        bc[7]=656045339.0/265891186.0;
        bc[8]=-3867574721.0/1518517206.0;
        bc[9]=465885868.0/322736535.0;
        bc[10]=53011238.0/667516719.0;
        bc[11]=2.0/45.0;
        *(A+13)=1.0/18.0;
        //riga 3
        *(A+26)=1.0/48.0;
        *(A+27)=1.0/16.0;
        //riga 4
        *(A+39)=1.0/32.0;
        *(A+41)=3.0/32.0;
        //riga 5
        *(A+52)=5.0/16.0;
        *(A+54)=-75.0/64.0;
        *(A+55)=75.0/64.0;
        //riga 6
        *(A+65)=3.0/80.0;
        *(A+68)=3.0/16.0;
        *(A+69)=3.0/20.0;
        //riga 7
        *(A+78)=29443841.0/614563906.0;
        *(A+81)=77736538.0/692538347.0;
        *(A+82)=-28693883.0/1125000000.0;
        *(A+83)=23124283.0/1800000000.0;
        //riga 8
        *(A+91)=16016141.0/946692911.0;
        *(A+94)=61564180.0/158732637.0;
        *(A+95)=22789713.0/633445777.0;
        *(A+96)=545815736.0/2771057229.0;
        *(A+97)=-180193667.0/1043307555.0;
        //riga 9
        *(A+104)=39632708.0/573591083.0;
        *(A+107)=-433636366.0/683701615.0;
        *(A+108)=-421739975.0/2616292301.0;
        *(A+109)=100302831.0/723423059.0;
        *(A+110)=790204164.0/839813087.0;
        *(A+111)=800635310.0/3783071287.0;
        //riga 10
        *(A+117)=246121993.0/1340847787.0;
        *(A+120)=-37695042795.0/15268766246.0;
        *(A+121)=-309121744.0/1061227803.0;
        *(A+122)=-12992083.0/490766935.0;
        *(A+123)=6005943493.0/2108947869.0;
        *(A+124)=393006217.0/1396673457.0;
        *(A+125)=123872331.0/1001029789.0;
        //riga 11
        *(A+130)=-1028468189.0/846180014.0;
        *(A+133)=8478235783.0/508512852.0;
        *(A+134)=1311729495.0/1432422823.0;
        *(A+135)=-10304129995.0/1701304382.0;
        *(A+136)=-48777925059.0/3047939560.0;
        *(A+137)=15336726248.0/1032824649.0;
        *(A+138)=-45442868181.0/3398467696.0;
        *(A+139)=3065993473.0/597172653.0;
        //riga 12
        *(A+143)=185892177.0/718116043.0;
        *(A+146)=-3185094517.0/667107341.0;
        *(A+147)=-477755414.0/1098053517.0;
        *(A+148)=-703635378.0/230739211.0;
        *(A+149)=5731566787.0/1027545527.0;
        *(A+150)=5232866602.0/850066563.0;
        *(A+151)=-4093664535.0/808688257.0;
        *(A+152)=3962137247.0/1805957418.0;
        *(A+153)=65686358.0/487910083.0;
        //riga 13
        *(A+156)=403863854.0/491063109.0;
        *(A+159)=-5068492393.0/434740067.0;
        *(A+160)=-411421997.0/543043805.0;
        *(A+161)=652783627.0/914296604.0;
        *(A+162)=11173962825.0/925320556.0;
        *(A+163)=-13158990841.0/6184727034.0;
        *(A+164)=3936647629.0/1978049680.0;
        *(A+165)=-160528059.0/685178525.0;
        *(A+166)=248638103.0/1413531060.0;
calcola_c(ns,c,A);
}

// METODI IMPLICITI

void EI(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind){
azzera(ns,A,b,c,bc,Ind);
// EULERO IMPLICITO / Radau1, p = 1, ns = 1;
// implicito DIRK
b[0]=1.0;
A[0]=1.0;
calcola_c(ns,c,A);
}

void CN(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Crank - Nicolson p = 2, ns = 2;
// implicito DIRK
b[0]=b[1]=0.5;
*(A+2)=*(A+3)=0.5;
if(Ind!=0)Ind[0]=0;
calcola_c(ns,c,A);
}

void GAUSS1(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Gauss 1   p = 2, ns = 1;
// implicito DIRK
        b[0]=1;
        *A=0.5;
calcola_c(ns,c,A);
}

void GAUSS2(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Gauss 2   p = 4, ns = 2;
// implicito
        b[0]=0.5;
        b[1]=0.5;
        *A=1.0/4.0;
        *(A+1)=(3.0-2.0*sqrt(3))/12.0;
        *(A+2)=(3.0+2.0*sqrt(3))/12.0;
        *(A+3)=1.0/4.0;
calcola_c(ns,c,A);
}

void GAUSS3(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Gauss 3   p = 6, ns = 3;
// implicito
        b[0]=5.0/18.0;
        b[1]=4.0/9.0;
        b[2]=5.0/18.0;
        *A=5.0/36.0;
        *(A+1)=2.0/9.0-sqrt(15)/15.0;
        *(A+2)=5.0/36.0-sqrt(15)/30.0;
        *(A+3)=5.0/36.0+sqrt(15)/24.0;
        *(A+4)=2.0/9.0;
        *(A+5)=5.0/36.0-sqrt(15)/24.0;
        *(A+6)=5.0/36.0+sqrt(15)/30.0;
        *(A+7)=2.0/9.0+sqrt(15)/15.0;
        *(A+8)=5.0/36.0;
calcola_c(ns,c,A);
}

void RADAU2(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Radau 2   p = 3, ns = 2;
// implicito
        b[0]=3.0/4.0;
        b[1]=1.0/4.0;
        *A=5.0/12.0;
        *(A+1)=-1.0/12.0;
        *(A+2)=3.0/4.0;
        *(A+3)=1.0/4.0;
calcola_c(ns,c,A);
}

void RADAU3(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind){
azzera(ns,A,b,c,bc,Ind);
// Radau 3   p = 5, ns = 3;
// implicito
        Real s6=sqrt(6);
        b[0]=(16.0-s6)/36.0;
        b[1]=(16.0+s6)/36.0;
        b[2]=1.0/9.0;
        *A=(88.0-7.0*s6)/360.0;
        *(A+1)=(296.0-169.0*s6)/1800.0;
        *(A+2)=(-2.0+3.0*s6)/225.0;
        *(A+3)=(296.0+169.0*s6)/1800.0;
        *(A+4)=(88.0+7.0*s6)/360.0;;
        *(A+5)=(-2.0-3.0*s6)/225.0;
        *(A+6)=(16.0-s6)/36.0;
        *(A+7)=(16.0+s6)/36.0;
        *(A+8)=1.0/9.0;
calcola_c(ns,c,A);
}

void LOBATTO_3C_ORD2(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind)
{

    azzera(ns,A,b,c,bc,Ind);

    bc[0]=1.0;
    bc[1]=0;

    b[0]=0.5;
    b[1]=0.5;

    *A=0.5;
    *(A+1)=-0.5;
    *(A+2)=0.5;
    *(A+3)=0.5;

    calcola_c(ns,c,A);

}


void LOBATTO_3C_ORD4(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind)
{

    azzera(ns,A,b,c,bc,Ind);

    bc[0]=-0.5;
    bc[1]=2.0;
    bc[2]=-0.5;

    b[0]=1.0/6.0;
    b[1]=2.0/3.0;
    b[2]=1.0/6.0;



    *A=1.0/6.0;
    *(A+1)=-1.0/3.0;
    *(A+2)=1.0/6.0;

    *(A+3)=1.0/6.0;
    *(A+4)=5.0/12.0;
    *(A+5)=-1.0/12.0;

    *(A+6)=1.0/6.0;
    *(A+7)=2.0/3.0;
    *(A+8)=1.0/6.0;

    calcola_c(ns,c,A);
}







