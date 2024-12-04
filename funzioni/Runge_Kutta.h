using namespace std;
typedef  double Real;
void step(Real *v, Real h, Real *a, Real *K, int stadi, int d, Real *y);
// metodi espliciti
void EE(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);
void RK4(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);
void HEUN(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);
void EMOD(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);
void RKFEHL5(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);

// metodi DIRK: diagonalmente impliciti
// con elementi diagonali se non nulli tutti uguali a stesso valore

// metodi impliciti NON diagonali


// metodi espliciti immersi
void RKFEHL54(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);
void EEHEUN(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);
void DP87(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);

// metodi impliciti DIRK immersi(adattivi)
void EI(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);
void CN(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);
void GAUSS1(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);

// metodi impliciti NON DIAGONALI immersi (adattivi)
void GAUSS2(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);
void GAUSS3(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);
void RADAU2(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);
void RADAU3(int ns, Real *A, Real *b, Real *c, Real *bc,int *Ind);

// metodi impliciti NON DIAGONALI immersi (adattivi)
void LOBATTO_3C_ORD2(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);
void LOBATTO_3C_ORD4(int ns, Real *A, Real *b, Real *c, Real *bc, int *Ind);
