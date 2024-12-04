
typedef double Real;
void gauss(Real *A, Real x[], Real b[], int n);
void lu(Real *A, int P[], int n);
void risist(Real *A, int P[],Real x[], Real b[], int n);
void matmat(Real *A, Real *B, Real *Ris, int ra, int ca, int cb);
Real norm_2(Real v[],int n);
void copia(Real *a, Real *b,int d);
