#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "time.h"
#include "problem_size.hpp"
 
 int main(int argc, char const *argv[])
{

    clock_t start, end;

    float *ddivdx       = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *ddivdx_pl    = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *ddivdy       = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *ddivdy_pl    = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *ddivdz       = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *ddivdz_pl    = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *rhogvx       = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *rhogvx_pl    = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *rhogvy       = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *rhogvy_pl    = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *rhogvz       = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *rhogvz_pl    = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *rhogw        = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *rhogw_pl     = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *coef_intp    = (float*) malloc( (ADM_iall*ADM_jall*3        *ADM_nxyz*TJ*ADM_lall   ) * sizeof(float) );
    float *coef_intp_pl = (float*) malloc( (ADM_gall_pl      *3        *ADM_nxyz*      ADM_lall_pl) * sizeof(float) );
    float *coef_diff    = (float*) malloc( (ADM_iall*ADM_jall*6        *ADM_nxyz*      ADM_lall   ) * sizeof(float) );
    float *coef_diff_pl = (float*) malloc( (                  ADM_vlink*ADM_nxyz*      ADM_lall_pl) * sizeof(float) );
    float *RGSQRTH      = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *RGSQRTH_pl   = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *RGAM         = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *RGAM_pl      = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *RGAMH        = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall   ) * sizeof(float) );
    float *RGAMH_pl     = (float*) malloc( (ADM_gall_pl      *ADM_kall*ADM_lall_pl) * sizeof(float) );
    float *C2WfactGz    = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*6*ADM_lall   ) * sizeof(float) );
    float *C2WfactGz_pl = (float*) malloc( (ADM_gall_pl      *ADM_kall*6*ADM_lall_pl) * sizeof(float) );
    float *sclt         = (float*) malloc( (ADM_iall*ADM_jall*TJ) * sizeof(float) );
    float *sclt_pl      = (float*) malloc( (ADM_gall_pl      ) * sizeof(float) );
    float *rhogvx_vm    = (float*) malloc( (ADM_iall*ADM_jall) * sizeof(float) );
    float *rhogvx_vm_pl = (float*) malloc( (ADM_gall_pl      ) * sizeof(float) );
    float *rhogvy_vm    = (float*) malloc( (ADM_iall*ADM_jall) * sizeof(float) );
    float *rhogvy_vm_pl = (float*) malloc( (ADM_gall_pl      ) * sizeof(float) );
    float *rhogvz_vm    = (float*) malloc( (ADM_iall*ADM_jall) * sizeof(float) );
    float *rhogvz_vm_pl = (float*) malloc( (ADM_gall_pl      ) * sizeof(float) );
    float *rhogw_vm     = (float*) malloc( (ADM_iall*ADM_jall*ADM_kall) * sizeof(float) );
    float *rhogw_vm_pl  = (float*) malloc( (ADM_gall_pl      *ADM_kall) * sizeof(float) );
    float *GRD_rdgz     = (float*) malloc( (ADM_kall) * sizeof(float) );
    float sclt_rhogw;
    float sclt_rhogw_pl;
    int n;

    int imin = ADM_imin;
    int imax = ADM_imax;
    int jmin = ADM_jmin;
    int jmax = ADM_jmax;
    int kall = ADM_kall;
    int kmin = ADM_kmin;
    int kmax = ADM_kmax;
    int lall = ADM_lall;

    // Use following To init the values from dumpio.c
    /*
    // ########< read input data >######## //
    // set intial values using dumpio
  void *ORG_ddivdx      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_ddivdy      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_ddivdz      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_rhogvx      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_rhogvy      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_rhogvz      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_rhogw       =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_coef_intp   =  malloc( (ADM_iall*ADM_jall*3*ADM_nxyz*TJ*ADM_lall) * sizeof(float_type));
  void *ORG_coef_diff   =  malloc( (ADM_iall*ADM_jall*6*ADM_nxyz*ADM_lall) * sizeof(float_type));
  void *ORG_RGSQRTH     =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_RGAM        =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_RGAMH       =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  void *ORG_C2WfactGz   =  malloc( (ADM_iall*ADM_jall*ADM_kall*6*ADM_lall) * sizeof(float_type));

  int32_t EX_fid;  
  char    *EX_fname = (char*) malloc(1024 * sizeof(char));
  // IO_FREAD = fread mode
  int32_t IO_FREAD = 0;
  dumpio_syscheck();
  dumpio_mk_fname(EX_fname,"snapshot.dc_divdamp3d","pe",SET_prc_me-1,6);
  dumpio_fopen(EX_fname,EX_fid);
  // MUST READ ALL VALUES (INCLUDING _p ARRAYS)
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_ddivdx);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_ddivdy);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_ddivdz);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_rhogvx);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_rhogvy);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_rhogvz);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_rhogw);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*3*ADM_nxyz*2*ADM_lall   , ORG_coef_intp);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*6*ADM_nxyz*  ADM_lall   , ORG_coef_diff);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall       , ORG_RGSQRTH);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall       , ORG_RGAM);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall       , ORG_RGAMH);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*6*ADM_lall     , ORG_C2WfactGz);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     );
  dumpio_fclose(EX_fid);
  
  // can a storage type be assigned a value from an array??
  // example: scl = ORG_scl

  for(int l=0; l<ADM_lall; ++l) // how to retrieve the l dimension??
  {
   for(int i=0; i<ADM_iall; ++i)
    for(int j=0; j<ADM_jall; ++j)
      for(int k=0; k<ADM_kall; ++k)
      {
        ddivdx[index3D(i,j,k)]= *(float_type*)((ORG_ddivdx + index4D(i,j,k,l)));
        ddivdy[index3D(i,j,k)]= *(float_type*)((ORG_ddivdy + index4D(i,j,k,l)));
        ddivdz[index3D(i,j,k)]= *(float_type*)((ORG_ddivdz + index4D(i,j,k,l)));
        rhogvx[index3D(i,j,k)]= *(float_type*)((ORG_rhogvx + index4D(i,j,k,l)));
        rhogvy[index3D(i,j,k)]= *(float_type*)((ORG_rhogvy + index4D(i,j,k,l)));
        rhogvz[index3D(i,j,k)]= *(float_type*)((ORG_rhogvz + index4D(i,j,k,l)));
        RGSQRTH[index3D(i,j,k)]=*(float_type*)((ORG_RGSQRTH+ index4D(i,j,k,l)));
        RGAM[index3D(i,j,k)]   =*(float_type*)((ORG_RGAM   + index4D(i,j,k,l)));
        RGAMH[index3D(i,j,k)]  =*(float_type*)((ORG_RGAMH  + index4D(i,j,k,l)));
        for(int d=0; d<6; ++d)      
            C2WfactGz[index4D(i,j,k,d)]=* (float_type*)((ORG_C2WfactGz + index5D(i,j,k,d,l)));
        
      }
   for(int i=0; i<ADM_iall; ++i)
    for(int j=0; j<ADM_jall; ++j)
    {
      // the k dimension is killed
      for(int a=0; a<3; ++a)
         for(int d=0; d<ADM_nxyz; ++d)      
            for(int t=0; t<TJ; ++t)      
              coef_intp[index6D(i,j,1,a,d,t)]=* (float_type*)((ORG_coef_intp + index6D(i,j,a,d,t,l)));
      // the k dimension is killed
      for(int a=0; a<6; ++a)
         for(int d=0; d<ADM_nxyz; ++d)      
            coef_diff[index5D(i,j,1,a,d)]=* (float_type*)((ORG_coef_diff + index5D(i,j,a,d,l)));
    }
  }  // l loop
  for(int i=0; i<ADM_iall; ++i)
    for(int j=0; j<ADM_jall; ++j) {
      rhogvx_vm[index2D(i,j)]=(i+j)/2;
      rhogvy_vm[index2D(i,j)]=(i+j)/4;
      rhogvz_vm[index2D(i,j)]=(i+j)/8;
      for(int k=0; k<ADM_kall; ++k)
        rhogw_vm[index3D(i,j,k)]=(i+j+k)/2;
      for(int t=0; t<TJ; ++t)  
        sclt[index3D(i,j,t)]=(i+j+t)/2;
    }
  for(int k=0; k<ADM_kall; ++k) 
    GRD_rdgz[k]=(k)/2;   

  free(ORG_ddivdx);     
  //free();
  free(ORG_ddivdy);
  //free();
  free(ORG_ddivdz);
  //free();
  free(ORG_rhogvx);
  //free();
  free(ORG_rhogvy);
  //free();
  free(ORG_rhogvz);
  //free();
  free(ORG_rhogw);
  //free();
  free(ORG_coef_intp);
  //free();
  free(ORG_coef_diff);
  //free();
  free(ORG_RGSQRTH);
  //free();
  free(ORG_RGAM);
  //free();
  free(ORG_RGAMH);
  //free();
  free(ORG_C2WfactGz);
  //free();
  
  //###############################################################################
  */
  //scl.print();
  for(int l=0; l<ADM_lall; ++l) // how to retrieve the l dimension??
  {
   for(int i=0; i<ADM_iall; ++i)
    for(int j=0; j<ADM_jall; ++j)
      for(int k=0; k<ADM_kall; ++k)
      {
        ddivdx[index3D(i,j,k)]= (i+j+k)/2;
        ddivdy[index3D(i,j,k)]= (i+j+k)/2;
        ddivdz[index3D(i,j,k)]= (i+j+k)/4;
        rhogvx[index3D(i,j,k)]= (i+j+k)/4;
        rhogvy[index3D(i,j,k)]= (i+j+k)/8;
        rhogvz[index3D(i,j,k)]= (i+j+k)/8;
        RGSQRTH[index3D(i,j,k)]=(i+j+k)/2;
        RGAM[index3D(i,j,k)]   =(i+j+k)/2;
        RGAMH[index3D(i,j,k)]  =(i+j+k)/4;
        for(int d=0; d<6; ++d)      
            C2WfactGz[index4D(i,j,k,d)]=(i+j+k+d)/2;
        
      }
   for(int i=0; i<ADM_iall; ++i)
    for(int j=0; j<ADM_jall; ++j)
    {
      // the k dimension is killed
      for(int a=0; a<3; ++a)
         for(int d=0; d<ADM_nxyz; ++d)      
            for(int t=0; t<TJ; ++t)      
              coef_intp[index6D(i,j,1,a,d,t)]=(i+j+a+d)/2;
      // the k dimension is killed
      for(int a=0; a<6; ++a)
         for(int d=0; d<ADM_nxyz; ++d)      
            coef_diff[index5D(i,j,1,a,d)]=(i+j+d)/2;
    }
  }  // l loop
  for(int i=0; i<ADM_iall; ++i)
    for(int j=0; j<ADM_jall; ++j) {
      rhogvx_vm[index2D(i,j)]=(i+j)/2;
      rhogvy_vm[index2D(i,j)]=(i+j)/4;
      rhogvz_vm[index2D(i,j)]=(i+j)/8;
      for(int k=0; k<ADM_kall; ++k)
        rhogw_vm[index3D(i,j,k)]=(i+j+k)/2;
      for(int t=0; t<TJ; ++t)  
        sclt[index3D(i,j,t)]=(i+j+t)/2;
    }
  for(int k=0; k<ADM_kall; ++k) 
    GRD_rdgz[k]=k/2.0;    

  start = clock();
  // -1 to account for zero-index
  imin = ADM_imin -1;
  imax = ADM_imax -1;
  jmin = ADM_jmin -1;
  jmax = ADM_jmax -1;
  kall = ADM_kall -1;
  kmin = ADM_kmin -1;
  kmax = ADM_kmax -1;
  lall = ADM_lall -1; 

  #pragma omp parallel default(none),private(i,j,k,l,sclt_rhogw), \
    shared(imin,imax,jmin,jmax,kall,kmin,kmax,lall,ADM_have_sgp,GRD_rdgz, \
    ddivdx,ddivdy,ddivdz,rhogvx,rhogvy,rhogvz,rhogw,sclt,coef_intp,coef_diff, \
    rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm,C2WfactGz,RGAMH,RGSQRTH,RGAM)
  for(int l = 0; l<ADM_lall; ++l) {  

       #pragma omp for schedule(static)
       for (int i = imin-1; i < imax+1; ++i)
        for (int j = jmin-1; j < jmax+1; ++j)
          for (int k = kmin+1; k < kmax; ++k)
            rhogw_vm[index3D(i,j,k)] = ( C2WfactGz[index5D(i,j,k,1,l)] * rhogvx[index4D(i,j,k,l)] 
                                       + C2WfactGz[index5D(i,j,k,2,l)] * rhogvx[index4D(i,j,k-1,l)] 
                                       + C2WfactGz[index5D(i,j,k,3,l)] * rhogvy[index4D(i,j,k  ,l)] 
                                       + C2WfactGz[index5D(i,j,k,4,l)] * rhogvy[index4D(i,j,k-1,l)] 
                                       + C2WfactGz[index5D(i,j,k,5,l)] * rhogvz[index4D(i,j,k  ,l)] 
                                       + C2WfactGz[index5D(i,j,k,6,l)] * rhogvz[index4D(i,j,k-1,l)] )
                                       * RGAMH[index4D(i,j,k,l)]                          
                                       + rhogw[index4D(i,j,k,l)] * RGSQRTH[index4D(i,j,k,l)];            

       #pragma omp for schedule(static)
       for (int i = imin-1; i < imax+1; ++i)
        for (int j = jmin-1; j < jmax+1; ++j) {
          rhogw_vm[index3D(i,j,kmin)] = 0.0;
          rhogw_vm[index3D(i,j,kmax+1)] = 0.0;
        }

       for (int k = kmin; k < kmax; ++k) {

          #pragma omp for schedule(static)
          for (int i = imin-1; i<imax+1; ++i)
            for (int j = jmin-1; j<jmax+1; ++j) {
             rhogvx_vm[index2D(i,j)] = rhogvx[index4D(i,j,k,l)] * RGAM[index4D(i,j,k,l)];
             rhogvy_vm[index2D(i,j)] = rhogvy[index4D(i,j,k,l)] * RGAM[index4D(i,j,k,l)];
             rhogvz_vm[index2D(i,j)] = rhogvz[index4D(i,j,k,l)] * RGAM[index4D(i,j,k,l)];
          }
          // affirm nowait
          #pragma omp for nowait schedule(static)
          for (int i = imin-1; i<imax; ++i)
            for (int j = jmin-1; j<jmax; ++j) {
             sclt_rhogw = ( ( rhogw_vm[index3D(i,j,k+1)] + rhogw_vm[index3D(i+1,j,k+1)] + rhogw_vm[index3D(i+1,j+1,k+1)] ) 
                          - ( rhogw_vm[index3D(i,j,k  )] + rhogw_vm[index3D(i+1,j,k  )] + rhogw_vm[index3D(i+1,j+1,k  )] ) 
                          ) / 3.0 * GRD_rdgz[k];

             sclt[index3D(i,j,TI)] = coef_intp[index6D(i,j,1,XDIR,TI,l)] * rhogvx_vm[index2D(i,j)] 
                          + coef_intp[index6D(i,j,2,XDIR,TI,l)] * rhogvx_vm[index2D(i+1,j  )] 
                          + coef_intp[index6D(i,j,3,XDIR,TI,l)] * rhogvx_vm[index2D(i+1,j+1)] 
                          + coef_intp[index6D(i,j,1,YDIR,TI,l)] * rhogvy_vm[index2D(i  ,j  )] 
                          + coef_intp[index6D(i,j,2,YDIR,TI,l)] * rhogvy_vm[index2D(i+1,j  )] 
                          + coef_intp[index6D(i,j,3,YDIR,TI,l)] * rhogvy_vm[index2D(i+1,j+1)] 
                          + coef_intp[index6D(i,j,1,ZDIR,TI,l)] * rhogvz_vm[index2D(i  ,j  )] 
                          + coef_intp[index6D(i,j,2,ZDIR,TI,l)] * rhogvz_vm[index2D(i+1,j  )] 
                          + coef_intp[index6D(i,j,3,ZDIR,TI,l)] * rhogvz_vm[index2D(i+1,j+1)] 
                          + sclt_rhogw;
          }

          #pragma omp for schedule(static)
          for (int i = imin-1; i<imax; ++i)
            for (int j = jmin-1; j<jmax; ++j) {
             sclt_rhogw = ( ( rhogw_vm[index3D(i,j,k+1)] + rhogw_vm[index3D(i+1,j+1,k+1)] + rhogw_vm[index3D(i,j+1,k+1)] )
                          - ( rhogw_vm[index3D(i,j,k  )] + rhogw_vm[index3D(i+1,j+1,k  )] + rhogw_vm[index3D(i,j+1,k  )] )
                          ) / 3.0 * GRD_rdgz[k];

             sclt[index3D(i,j,TJ)] = coef_intp[index6D(i,j,1,XDIR,TJ,l)] * rhogvx_vm[index2D(i  ,j  )]
                          + coef_intp[index6D(i,j,2,XDIR,TJ,l)] * rhogvx_vm[index2D(i+1,j+1)]
                          + coef_intp[index6D(i,j,3,XDIR,TJ,l)] * rhogvx_vm[index2D(i  ,j+1)]
                          + coef_intp[index6D(i,j,1,YDIR,TJ,l)] * rhogvy_vm[index2D(i  ,j  )]
                          + coef_intp[index6D(i,j,2,YDIR,TJ,l)] * rhogvy_vm[index2D(i+1,j+1)]
                          + coef_intp[index6D(i,j,3,YDIR,TJ,l)] * rhogvy_vm[index2D(i  ,j+1)]
                          + coef_intp[index6D(i,j,1,ZDIR,TJ,l)] * rhogvz_vm[index2D(i  ,j  )]
                          + coef_intp[index6D(i,j,2,ZDIR,TJ,l)] * rhogvz_vm[index2D(i+1,j+1)]
                          + coef_intp[index6D(i,j,3,ZDIR,TJ,l)] * rhogvz_vm[index2D(i  ,j+1)]
                          + sclt_rhogw;
          }

          if ( ADM_have_sgp[l] ) //pentagon
             sclt[index3D(imin-1,jmin-1,TI)] = sclt[index3D(imin,jmin-1,TJ)];
          
          #pragma omp for nowait schedule(static)
          for (int i = imin; i<imax; ++i)
            for (int j = jmin; j<jmax; ++j) {
             ddivdx[index4D(i,j,k,l)] = coef_diff[index5D(i,j,1,XDIR,l)] * ( sclt[index3D(i  ,j  ,TI)] + sclt[index3D(i  ,j  ,TJ)] )
                                      + coef_diff[index5D(i,j,2,XDIR,l)] * ( sclt[index3D(i  ,j  ,TJ)] + sclt[index3D(i-1,j  ,TI)] )
                                      + coef_diff[index5D(i,j,3,XDIR,l)] * ( sclt[index3D(i-1,j  ,TI)] + sclt[index3D(i-1,j-1,TJ)] )
                                      + coef_diff[index5D(i,j,4,XDIR,l)] * ( sclt[index3D(i-1,j-1,TJ)] + sclt[index3D(i-1,j-1,TI)] )
                                      + coef_diff[index5D(i,j,5,XDIR,l)] * ( sclt[index3D(i-1,j-1,TI)] + sclt[index3D(i  ,j-1,TJ)] )
                                      + coef_diff[index5D(i,j,6,XDIR,l)] * ( sclt[index3D(i  ,j-1,TJ)] + sclt[index3D(i  ,j  ,TI)] );
          }

          #pragma omp for nowait schedule(static)
          for (int i = imin; i<imax; ++i)
            for (int j = jmin; j<jmax; ++j) {
             ddivdy[index4D(i,j,k,l)] = coef_diff[index5D(i,j,1,YDIR,l)] * ( sclt[index3D(i  ,j  ,TI)] + sclt[index3D(i  ,j  ,TJ)] )
                                      + coef_diff[index5D(i,j,2,YDIR,l)] * ( sclt[index3D(i  ,j  ,TJ)] + sclt[index3D(i-1,j  ,TI)] )
                                      + coef_diff[index5D(i,j,3,YDIR,l)] * ( sclt[index3D(i-1,j  ,TI)] + sclt[index3D(i-1,j-1,TJ)] )
                                      + coef_diff[index5D(i,j,4,YDIR,l)] * ( sclt[index3D(i-1,j-1,TJ)] + sclt[index3D(i-1,j-1,TI)] )
                                      + coef_diff[index5D(i,j,5,YDIR,l)] * ( sclt[index3D(i-1,j-1,TI)] + sclt[index3D(i  ,j-1,TJ)] )
                                      + coef_diff[index5D(i,j,6,YDIR,l)] * ( sclt[index3D(i  ,j-1,TJ)] + sclt[index3D(i  ,j  ,TI)] );
          }

          #pragma omp for nowait schedule(static)
          for (int i = imin; i<imax; ++i)
            for (int j = jmin; j<jmax; ++j) {
             ddivdz[index4D(i,j,k,l)] = coef_diff[index5D(i,j,1,ZDIR,l)] * ( sclt[index3D(i  ,j  ,TI)] + sclt[index3D(i  ,j  ,TJ)] )
                                      + coef_diff[index5D(i,j,2,ZDIR,l)] * ( sclt[index3D(i  ,j  ,TJ)] + sclt[index3D(i-1,j  ,TI)] )
                                      + coef_diff[index5D(i,j,3,ZDIR,l)] * ( sclt[index3D(i-1,j  ,TI)] + sclt[index3D(i-1,j-1,TJ)] )
                                      + coef_diff[index5D(i,j,4,ZDIR,l)] * ( sclt[index3D(i-1,j-1,TJ)] + sclt[index3D(i-1,j-1,TI)] )
                                      + coef_diff[index5D(i,j,5,ZDIR,l)] * ( sclt[index3D(i-1,j-1,TI)] + sclt[index3D(i  ,j-1,TJ)] )
                                      + coef_diff[index5D(i,j,6,ZDIR,l)] * ( sclt[index3D(i  ,j-1,TJ)] + sclt[index3D(i  ,j  ,TI)] );
          }
           //workshare 
          #pragma omp workshare
          {
            for (int i = 0; i < ADM_iall; ++i) {
              ddivdx[index4D(i,jmin-1,k,l)] = 0.0;
              ddivdy[index4D(i,jmin-1,k,l)] = 0.0;
              ddivdz[index4D(i,jmin-1,k,l)] = 0.0;
              ddivdx[index4D(i,jmax+1,k,l)] = 0.0;
              ddivdy[index4D(i,jmax+1,k,l)] = 0.0;
              ddivdz[index4D(i,jmax+1,k,l)] = 0.0;
            }            
            for (int j = 0; j < ADM_jall; ++j) {
              ddivdx[index4D(imin-1,j,k,l)] = 0.0;
              ddivdy[index4D(imin-1,j,k,l)] = 0.0;
              ddivdz[index4D(imin-1,j,k,l)] = 0.0;
              ddivdx[index4D(imax+1,j,k,l)] = 0.0;
              ddivdy[index4D(imax+1,j,k,l)] = 0.0;
              ddivdz[index4D(imax+1,j,k,l)] = 0.0;
            }
          }
        } // end k loop
         #pragma omp workshare
         {
           for (int i = 0; i < ADM_iall; ++i)
            for (int j = 0; j < ADM_jall; ++j) {
               ddivdx[index4D(i,j,kmin-1,l)] = 0.0;
               ddivdx[index4D(i,j,kmax+1,l)] = 0.0;
               ddivdy[index4D(i,j,kmin-1,l)] = 0.0;
               ddivdy[index4D(i,j,kmax+1,l)] = 0.0;
               ddivdz[index4D(i,j,kmin-1,l)] = 0.0;
               ddivdz[index4D(i,j,kmax+1,l)] = 0.0;
           }
         }
    } // l loop
    
    if ( ADM_have_pl )  {
       for (int l = 0; l < ADM_lall_pl; ++l) {
        for (int ij = 0; ij < ADM_gall_pl; ++ij) {
          for (int k = ADM_kmin+1; k < ADM_kmax; ++k)
             rhogw_vm_pl[index2D(ij,k)] = ( C2WfactGz_pl[index4D(ij,k,1,l)] * rhogvx_pl[index3D(ij,k  ,l)]
                                  + C2WfactGz_pl[index4D(ij,k,2,l)] * rhogvx_pl[index3D(ij,k-1,l)]
                                  + C2WfactGz_pl[index4D(ij,k,3,l)] * rhogvy_pl[index3D(ij,k  ,l)]
                                  + C2WfactGz_pl[index4D(ij,k,4,l)] * rhogvy_pl[index3D(ij,k-1,l)]
                                  + C2WfactGz_pl[index4D(ij,k,5,l)] * rhogvz_pl[index3D(ij,k  ,l)]
                                  + C2WfactGz_pl[index4D(ij,k,6,l)] * rhogvz_pl[index3D(ij,k-1,l)]
                                  ) * RGAMH_pl[index3D(ij,k,l)]
                                  + rhogw_pl[index3D(ij,k,l)] * RGSQRTH_pl[index3D(ij,k,l)];
        }
          
        for (int ij = 0; ij < ADM_gall_pl; ++ij) {
             rhogw_vm_pl[index2D(ij,ADM_kmin  )] = 0.0;
             rhogw_vm_pl[index2D(ij,ADM_kmax+1)] = 0.0;
        }

        n = ADM_gslf_pl;

        for (int k = ADM_kmin; k < ADM_kmax; ++k) {
          for (int v = 0; v < ADM_gall_pl; ++v) {
              rhogvx_vm_pl[v] = rhogvx_pl[index3D(v,k,l)] * RGAM_pl[index3D(v,k,l)];
              rhogvy_vm_pl[v] = rhogvy_pl[index3D(v,k,l)] * RGAM_pl[index3D(v,k,l)];
              rhogvz_vm_pl[v] = rhogvz_pl[index3D(v,k,l)] * RGAM_pl[index3D(v,k,l)];
          }

          for (int v = ADM_gmin_pl; v < ADM_gmax_pl; ++v)  {
                int ij   = v;
                int ijp1 = v + 1;
                if( ijp1 == ADM_gmax_pl+1 )
                    ijp1 = ADM_gmin_pl;

                sclt_rhogw_pl = ( ( rhogw_vm_pl[index2D(n,k+1)] + rhogw_vm_pl[index2D(ij,k+1)] + rhogw_vm_pl[index2D(ijp1,k+1)] )
                                - ( rhogw_vm_pl[index2D(n,k  )] + rhogw_vm_pl[index2D(ij,k  )] + rhogw_vm_pl[index2D(ijp1,k  )] )
                                ) / 3.0 * GRD_rdgz[k];

                sclt_pl[ij] = coef_intp_pl[index4D(v,1,XDIR,l)] * rhogvx_vm_pl[n   ]
                            + coef_intp_pl[index4D(v,2,XDIR,l)] * rhogvx_vm_pl[ij  ]
                            + coef_intp_pl[index4D(v,3,XDIR,l)] * rhogvx_vm_pl[ijp1]
                            + coef_intp_pl[index4D(v,1,YDIR,l)] * rhogvy_vm_pl[n   ]
                            + coef_intp_pl[index4D(v,2,YDIR,l)] * rhogvy_vm_pl[ij  ]
                            + coef_intp_pl[index4D(v,3,YDIR,l)] * rhogvy_vm_pl[ijp1]
                            + coef_intp_pl[index4D(v,1,ZDIR,l)] * rhogvz_vm_pl[n   ]
                            + coef_intp_pl[index4D(v,2,ZDIR,l)] * rhogvz_vm_pl[ij  ]
                            + coef_intp_pl[index4D(v,3,ZDIR,l)] * rhogvz_vm_pl[ijp1]
                            + sclt_rhogw_pl;
             }

             for (int i = 0; i < ADM_iall; ++i) {
               ddivdx_pl[index3D(i,k,l)] = 0.0;
               ddivdy_pl[index3D(i,k,l)] = 0.0;
               ddivdz_pl[index3D(i,k,l)] = 0.0;
             }

             for(int v =ADM_gmin_pl; v < ADM_gmax_pl; ++v) {
                int ij   = v;
                int ijm1 = v - 1;
                if( ijm1 == ADM_gmin_pl-1 )
                  ijm1 = ADM_gmax_pl; // cyclic condition

                ddivdx_pl[index3D(n,k,l)] = ddivdx_pl[index3D(n,k,l)] + coef_diff_pl[index3D(v-1,XDIR,l)] * ( sclt_pl[ijm1] + sclt_pl[ij] );
                ddivdy_pl[index3D(n,k,l)] = ddivdy_pl[index3D(n,k,l)] + coef_diff_pl[index3D(v-1,YDIR,l)] * ( sclt_pl[ijm1] + sclt_pl[ij] );
                ddivdz_pl[index3D(n,k,l)] = ddivdz_pl[index3D(n,k,l)] + coef_diff_pl[index3D(v-1,ZDIR,l)] * ( sclt_pl[ijm1] + sclt_pl[ij] );
             }
          }

          for (int ij = 0; ij < ADM_gall_pl; ++ij) {
             ddivdx_pl[index3D(ij,ADM_kmin-1,l)] = 0.0;
             ddivdx_pl[index3D(ij,ADM_kmax+1,l)] = 0.0;
             ddivdy_pl[index3D(ij,ADM_kmin-1,l)] = 0.0;
             ddivdy_pl[index3D(ij,ADM_kmax+1,l)] = 0.0;
             ddivdz_pl[index3D(ij,ADM_kmin-1,l)] = 0.0;
             ddivdz_pl[index3D(ij,ADM_kmax+1,l)] = 0.0;
          }
         }  // end l loop
      }
    else {
      for (int i = 0; i < ADM_iall; ++i)
        for (int j = 0; j < ADM_jall; ++j)
          for (int k = 0; k < ADM_kall; ++k) {
             ddivdx_pl[index3D(i,j,k)] = 0.0;
             ddivdy_pl[index3D(i,j,k)] = 0.0;
             ddivdz_pl[index3D(i,j,k)] = 0.0;
     }
    }
  
  end = clock();
  float time1 = ((float)(end-start))/CLOCKS_PER_SEC;
  printf("Runtime: %f", time1);
  
  free(ddivdx);     
  free(ddivdx_pl);
  free(ddivdy);
  free(ddivdy_pl);
  free(ddivdz);
  free(ddivdz_pl);
  free(rhogvx);
  free(rhogvx_pl);
  free(rhogvy);
  free(rhogvy_pl);
  free(rhogvz);
  free(rhogvz_pl);
  free(rhogw);
  free(rhogw_pl);
  free(coef_intp);
  free(coef_intp_pl);
  free(coef_diff);
  free(coef_diff_pl);
  free(RGSQRTH);
  free(RGSQRTH_pl);
  free(RGAM);
  free(RGAM_pl);
  free(RGAMH);
  free(RGAMH_pl);
  free(C2WfactGz);
  free(C2WfactGz_pl);
  free(GRD_rdgz);

    



    // a CPU side storage to compare results
    // can pass pointer to compare to another data
    //storage_type_6d ref(metadata_vt, 0.0, "ref");
    /*
// compare ref 

    verifier verif(1e-13);
    array<array<uint_t, 2>, 3> halos{{ {halo_size,halo_size}, {halo_size,halo_size}, {halo_size,halo_size} }};
    bool result = verif.verify(grid, vt, ref, halos);

// report timing
#ifdef BENCHMARK
        std::cout << divdamp3d->print_meter() << std::endl;
#endif

    ASSERT_TRUE(result);
}
*/

    return 0;
}


/**
@}
*/
