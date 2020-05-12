#include <cmath>
#include <iomanip>
#include <iostream>

inline void mvm(double A[],double x[],double b[], int size)
{
  for(int i = 0; i < size; i++)
  {
    b[i] = 0.;
    for(int j = 0; j < size; j++)
    {
      b[i] += A[i*size+j] * x[j];
    }
  }
}

inline double ip(double a[], double b[], int size)
{
  double value=0.;
  for(int i = 0; i < size; i++)
  {
    value += a[i]*b[i];
  }
  return value;
}

int main()
{
  int size = 3;
  double alpha, beta;
  double r0, rk, rk1;
  double eps = pow(2,-50);

  double *r = new double[size];
  double *p = new double[size];
  double *b = new double[size];
  double *x = new double[size];
  double *Ax = new double[size];
  double *Ap = new double[size];
  double *A = new double[size*size];


  A[0*size+0]=1.;A[0*size+1]=2.;A[0*size+2]=-3.;
  A[1*size+0]=2.;A[1*size+1]=5.;A[1*size+2]=-4.;
  A[2*size+0]=-3.;A[2*size+1]=-4.;A[2*size+2]=8.;


  b[0] = -4;b[1] = -3;b[2] = 12;

  for(int i = 0;i < size; i++) x[i]=0.;

  mvm(A, x, Ax, size);
  for(int i = 0; i < size; i++)
  {
    r[i] = b[i] - Ax[i];
  }
  for(int i = 0; i < size; i++)
  {
    p[i] = r[i];
  }

  r0 = sqrt(ip(r,r,size));

  for(int k = 0; k < 100000; k++)
  {
    mvm(A, p, Ap, size);
    alpha = ip(r, r, size) / ip(p, Ap, size);
    rk = ip(r, r, size);

    for(int i = 0; i<size; i++)
    {
      x[i] = x[i] + alpha*p[i];
    }
    for(int i = 0; i<size; i++)
    {
      r[i] = r[i]-alpha*Ap[i];
    }

    rk1 = sqrt(ip(r,r,size));
    std::cout << "step: " << k << "  rk: " << rk1 << std::endl;
    if(rk1 / r0 <= eps)
      break;

    beta = ip(r,r,size) / rk;
    for(int i = 0; i < size; i++)
    {
      p[i] = r[i] + beta*p[i];
    }
  }

  std::cout << eps << std::endl;

  for(int i = 0; i < size; i++)
  {
    std::cout << x[i] << std::endl;
  }

  delete [] r,p,b,x,Ax,Ap,A;

  return 0;
}
