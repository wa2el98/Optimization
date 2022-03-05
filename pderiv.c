// pderiv.c

#include <math.h>
#include "proj.h"


void find_minimum(double xstar[], double x0[],  int n, int m, FUN_PTR original_f, FUN_PTR h[], double epsilon)
 {
	 int min= 2147483647;
	 int tmp;
	 for(int i=0; i<m;i++)
	 {
		 tmp=(*h[i])(x0);
		 if(tmp<min)
			 min=tmp;
	 }
	 
 
 }



void copy_vector(double dest[], double source[], int n)
{
  int i;
  
  for(i=0; i < n; i++)
   dest[i] = source [i];  

} // copy_vector



double approx_partial_derivative(double (*obj_f)(double x[]),
     int i, double x[])
{
   double temp1,temp2, xi_orig, result, h, prev;
   double initial_h_const = 1.0/1024.0;
   double eps_const = 0.0000001;
   double local_x[NMAX];
  
  copy_vector(local_x, x, n+m);
   

   xi_orig = local_x[i];
   h = fabs(x[i]*initial_h_const);

   if(h < initial_h_const)
     h = initial_h_const; 

   local_x[i] =  xi_orig + h;

  temp1 = (*obj_f)(local_x); 


   local_x[i] =  xi_orig - h;

  temp2 = (*obj_f)(local_x); 

 
  result =  (temp1 - temp2)/(2*h);

do{
  h = h / 2.0;

  prev = result;

  local_x[i] =  xi_orig + h;

  temp1 = (*obj_f)(local_x); 

   local_x[i] =  xi_orig - h;

  temp2 = (*obj_f)(local_x); 
 
  result =  (temp1 - temp2)/(2*h);

} while ( fabs(prev - result) > ( fabs(prev) *  eps_const));

   local_x[i] =  xi_orig;

  return result;

} // approx_partial_derivative


double approx_twice_partial_derivative(double (*obj_f)(double x[]),
     int i, int j, double x[])
{
   double temp1,temp2,  xi_orig, xj_orig,result, hi, hj;
   double eps_const = 1.0/1048576.0;
   double local_x[NMAX];
  
  copy_vector(local_x, x, n);


   xj_orig = local_x[j];

   hj = fabs(x[j]*eps_const);

   if(hj < eps_const)
     hj = eps_const; 

 
   local_x[j] =  xj_orig + hj;

   temp1 =  approx_partial_derivative(obj_f, i, local_x);

   local_x[j] =  xj_orig - hj;

   temp2 =  approx_partial_derivative(obj_f, i, local_x);

  result =  (temp1 - temp2)/(2*hj);

   local_x[j] =  xj_orig;

  return result;

} // approx_partial_derivative
