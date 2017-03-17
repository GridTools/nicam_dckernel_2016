#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "time.h"
#include "problem_size.hpp"
 
 int main(int argc, char const *argv[])
{

    clock_t start_c, end_c;

    float *dscl         =  (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float));
    float *scl          =  (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float));
    float *kh           =  (float*) malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float));
    float *coef_intp    =  (float*) malloc( (ADM_iall*ADM_jall*3*ADM_nxyz*TJ*ADM_lall) * sizeof(float));
    float *coef_diff    =  (float*) malloc( (ADM_iall*ADM_jall*6*ADM_nxyz*ADM_lall) * sizeof(float));
    float *vt           =  (float*) malloc( (ADM_iall*ADM_jall*ADM_nxyz*TJ) * sizeof(float));

    float *dscl_pl      =  (float*) malloc( (ADM_gall_pl*ADM_kall*ADM_lall_pl) * sizeof(float));
    float *scl_pl       =  (float*) malloc( (ADM_gall_pl*ADM_kall*ADM_lall_pl) * sizeof(float));
    float *kh_pl        =  (float*) malloc( (ADM_gall_pl*ADM_kall*ADM_lall_pl) * sizeof(float));
    float *coef_intp_pl =  (float*) malloc( (ADM_gall_pl*3*ADM_nxyz*ADM_lall_pl) * sizeof(float));
    float *coef_diff_pl =  (float*) malloc( (ADM_vlink*ADM_nxyz*ADM_lall_pl) * sizeof(float));
    float *vt_pl        =  (float*) malloc( (ADM_gall_pl*ADM_nxyz) * sizeof(float));

    bool read_from_dump = false;

    void *ORG_dscl;
    void *ORG_scl;
    void *ORG_kh;
    void *ORG_coef_intp;
    void *ORG_coef_diff;

    if(read_from_dump) {
    // ########< read input data >######## //
    // set intial values using dumpio
        ORG_dscl      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float));
        ORG_scl       =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float));
        ORG_kh        =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float));
        ORG_coef_intp =  malloc( (ADM_iall*ADM_jall*3*ADM_nxyz*TJ*ADM_lall) * sizeof(float));
        ORG_coef_diff =  malloc( (ADM_iall*ADM_jall*6*ADM_nxyz*ADM_lall) * sizeof(float));

        int32_t EX_fid;
        char    *EX_fname = (char*) malloc(1024 * sizeof(char));
        // IO_FREAD = fread mode
        int32_t IO_FREAD = 0;
        dumpio_syscheck();
        dumpio_mk_fname(EX_fname,"snapshot.dc_diffusion","pe",SET_prc_me-1,6);
        dumpio_fopen(EX_fname,EX_fid);
        dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_dscl);
        dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_scl);
        dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_kh);
        dumpio_read_data( EX_fid, ADM_iall*ADM_jall*3*ADM_nxyz*2*ADM_lall   , ORG_coef_intp);
        dumpio_read_data( EX_fid, ADM_iall*ADM_jall*6*ADM_nxyz*  ADM_lall   , ORG_coef_diff);
    }

    double dx = 1/(double)ADM_iall;
    double dy = 1/(double)ADM_jall;

    for(int l=0; l<ADM_lall; ++l)
    {
      for(int i=0; i<ADM_iall; ++i) {
        for(int j=0; j<ADM_jall; ++j) {
          double x = dx * (double)(i);
          double y = dy * (double)(j);

          for(int k=0; k<ADM_kall; ++k)
          {
              if(read_from_dump) {
                  scl[index3D(i,j,k)]  = *(float*)((ORG_scl + index4D(i,j,k,l)));
                  dscl[index3D(i,j,k)] = *(float*)((ORG_dscl + index4D(i,j,k,l)));
              }
              else {
                  scl[index3D(i,j,k)]= 2.4 + 7.6 * (k / (double)ADM_kall + cos(PI * (x + 2.5 * y)) + sin(4 * PI * (x + 1.5 * y))) / 4.;
              }
          }
        // note: for coef_diff dims<2> is killed (i.e. no K dimension)
          for(int c=0; c<6; ++c)
            for(int d=0; d<ADM_nxyz; ++d)
            {
              if(read_from_dump) {
                coef_diff[index5D(i,j,1,c,d)]=* (float*)((ORG_coef_diff + index5D(i,j,c,d,l)));
              }
              else{
                coef_diff[index5D(i,j,1,c,d)]=2.2 + 5.5 * (c / 6.0 + d / 3. + cos(PI * (x + 2.5 * y)) + sin(1.2 * PI * (x + 1.5 * y))) / 4.;
              }
            }
        // note: for coef_intp dims<2> is killed (i.e. no K dimension)
          for(int a=0; a<3; ++a)
            for(int d=0; d<ADM_nxyz; ++d)
              for(int t=0; t<TJ; ++t)
              {
                  if(read_from_dump) {
                      coef_intp[index6D(i,j,1,a,d,t)]=* (float*)((ORG_coef_intp + index6D(i,j,a,d,t,l)));
                  }
                  else {
                      coef_intp[index6D(i,j,1,a,d,t)]=3. + 2.5 * (a / 3.0 + d / 3. + t/2.0 + cos(PI * (x + 2.5 * y)) + sin(2 * PI * (x + 3.5 * y))) / 4.;
                  }
              }
          for(int k=0; k<ADM_kall; ++k) {
              if(read_from_dump) {
                  kh[index3D(i,j,k)]=* (float*)((ORG_kh + index4D(i,j,k,l)));
              }
              else {
                  kh[index3D(i,j,k)]=3. + 2.5 * (k / (double)ADM_kall + cos(PI * (x + 2.5 * y)) + sin(2 * PI * (x + 3.5 * y))) / 4.;
              }
          }
        }
      }
    }  // l loop  

  start_c = clock();
  // -1 to account for zero-index
  int imin = ADM_imin-1;
  int imax = ADM_imax-1;
  int jmin = ADM_jmin-1;
  int jmax = ADM_jmax-1;
  int kall = ADM_kall-1;
  int lall = ADM_lall-1;
  // all constant indeces lowered by minus one
  #pragma omp parallel default(none) private(i,j,k,l,d) shared(imin,imax,jmin,jmax,kall,lall,ADM_have_sgp,dscl,scl,kh,vt,coef_intp,coef_diff)
  {
  for(int l = 0; l<ADM_lall; ++l) {
    for (int k = 0; k < ADM_kall; ++k) {

      #pragma omp for nowait schedule(static) collapse(2)
      for (int d = XDIR; d < ZDIR; ++d)
        for (int i = imin-1; i < imax; ++i)
          for (int j = jmin-1; j < jmax; ++j)
            vt[index4D(i,j,d,TI)] = ( ( + 2.0 * coef_intp[index6D(i,j,0,d,TI,l)]
                                - 1.0 * coef_intp[index6D(i,j,1,d,TI,l)]
                                - 1.0 * coef_intp[index6D(i,j,2,d,TI,l)] ) * scl[index4D(i  ,j  ,k,l)]
                            + ( - 1.0 * coef_intp[index6D(i,j,0,d,TI,l)]
                                + 2.0 * coef_intp[index6D(i,j,1,d,TI,l)]
                                - 1.0 * coef_intp[index6D(i,j,2,d,TI,l)] ) * scl[index4D(i+1,j  ,k,l)]
                            + ( - 1.0 * coef_intp[index6D(i,j,0,d,TI,l)]
                                - 1.0 * coef_intp[index6D(i,j,1,d,TI,l)]
                                + 2.0 * coef_intp[index6D(i,j,2,d,TI,l)] ) * scl[index4D(i+1,j+1,k,l)]
                            ) / 3.0;

      #pragma omp for schedule(static) collapse(2)
      for (int d = XDIR; d < ZDIR; ++d)
        for (int i = imin-1; i < imax; ++i)
          for (int j = jmin-1; j < jmax; ++j)
            vt[index4D(i,j,d,TJ)] = ( ( + 2.0 * coef_intp[index6D(i,j,0,d,TJ,l)]
                                - 1.0 * coef_intp[index6D(i,j,1,d,TJ,l)]
                                - 1.0 * coef_intp[index6D(i,j,2,d,TJ,l)] ) * scl[index4D(i  ,j  ,k,l)]
                            + ( - 1.0 * coef_intp[index6D(i,j,0,d,TJ,l)]
                                + 2.0 * coef_intp[index6D(i,j,1,d,TJ,l)]
                                - 1.0 * coef_intp[index6D(i,j,2,d,TJ,l)] ) * scl[index4D(i+1,j+1,k,l)]
                            + ( - 1.0 * coef_intp[index6D(i,j,0,d,TJ,l)]
                                - 1.0 * coef_intp[index6D(i,j,1,d,TJ,l)]
                                + 2.0 * coef_intp[index6D(i,j,2,d,TJ,l)] ) * scl[index4D(i  ,j+1,k,l)]
                            ) / 3.0;

      if ( ADM_have_sgp[l] ) {// pentagon
        #pragma omp master
        for (int d = XDIR-1; d < ZDIR-1; ++d)
          vt[index4D(imin-1,jmin-1,d,TI)] = vt[index4D(imin,jmin-1,d,TJ)]; 
      }
      
      #pragma omp for schedule(static)
      for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j)
          dscl[index4D(i,j,k,l)] = ( coef_diff[index5D(i,j,0,XDIR,l)] * ( vt[index4D(i  ,j  ,XDIR,TI)] + vt[index4D(i  ,j  ,XDIR,TJ)] )
                          + coef_diff[index5D(i,j,0,YDIR,l)] * ( vt[index4D(i  ,j  ,YDIR,TI)] + vt[index4D(i  ,j  ,YDIR,TJ)] )
                          + coef_diff[index5D(i,j,0,ZDIR,l)] * ( vt[index4D(i  ,j  ,ZDIR,TI)] + vt[index4D(i  ,j  ,ZDIR,TJ)] )
                          ) * 0.5 * ( kh[index4D(i  ,j  ,k,l)] + kh[index4D(i+1,j+1,k,l)] );

      #pragma omp for schedule(static)
      for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j)
          dscl[index4D(i,j,k,l)] = dscl[index4D(i,j,k,l)]
                        + ( coef_diff[index5D(i,j,1,XDIR,l)] * ( vt[index4D(i  ,j  ,XDIR,TJ)] + vt[index4D(i-1,j  ,XDIR,TI)] )
                          + coef_diff[index5D(i,j,1,YDIR,l)] * ( vt[index4D(i  ,j  ,YDIR,TJ)] + vt[index4D(i-1,j  ,YDIR,TI)] )
                          + coef_diff[index5D(i,j,1,ZDIR,l)] * ( vt[index4D(i  ,j  ,ZDIR,TJ)] + vt[index4D(i-1,j  ,ZDIR,TI)] )
                          ) * 0.5 * ( kh[index4D(i  ,j  ,k,l)] + kh[index4D(i  ,j+1,k,l)] );

      #pragma omp for schedule(static)
      for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j)
          dscl[index4D(i,j,k,l)] = dscl[index4D(i,j,k,l)]
                        + ( coef_diff[index5D(i,j,2,XDIR,l)] * ( vt[index4D(i-1,j  ,XDIR,TI)] + vt[index4D(i-1,j-1,XDIR,TJ)] )
                          + coef_diff[index5D(i,j,2,YDIR,l)] * ( vt[index4D(i-1,j  ,YDIR,TI)] + vt[index4D(i-1,j-1,YDIR,TJ)] )
                          + coef_diff[index5D(i,j,2,ZDIR,l)] * ( vt[index4D(i-1,j  ,ZDIR,TI)] + vt[index4D(i-1,j-1,ZDIR,TJ)] )
                          ) * 0.5 * ( kh[index4D(i-1,j  ,k,l)] + kh[index4D(i  ,j  ,k,l)] );

      #pragma omp for schedule(static)
      for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j)
          dscl[index4D(i,j,k,l)] = dscl[index4D(i,j,k,l)]
                        + ( coef_diff[index5D(i,j,3,XDIR,l)]  * ( vt[index4D(i-1,j-1,XDIR,TJ)] + vt[index4D(i-1,j-1,XDIR,TI)] )
                          + coef_diff[index5D(i,j,3,YDIR,l)] * ( vt[index4D(i-1,j-1,YDIR,TJ)] + vt[index4D(i-1,j-1,YDIR,TI)] )
                          + coef_diff[index5D(i,j,3,ZDIR,l)] * ( vt[index4D(i-1,j-1,ZDIR,TJ)] + vt[index4D(i-1,j-1,ZDIR,TI)] )
                          ) * 0.5 * ( kh[index4D(i-1,j-1,k,l)] + kh[index4D(i  ,j  ,k,l)] );

      #pragma omp for schedule(static)
      for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j)
          dscl[index4D(i,j,k,l)] = dscl[index4D(i,j,k,l)]
                        + ( coef_diff[index5D(i,j,4,XDIR,l)] * ( vt[index4D(i-1,j-1,XDIR,TI)] + vt[index4D(i  ,j-1,XDIR,TJ)] )
                          + coef_diff[index5D(i,j,4,YDIR,l)] * ( vt[index4D(i-1,j-1,YDIR,TI)] + vt[index4D(i  ,j-1,YDIR,TJ)] )
                          + coef_diff[index5D(i,j,4,ZDIR,l)] * ( vt[index4D(i-1,j-1,ZDIR,TI)] + vt[index4D(i  ,j-1,ZDIR,TJ)] )
                          ) * 0.5 * ( kh[index4D(i  ,j-1,k,l)] + kh[index4D(i  ,j  ,k,l)] );

      #pragma omp for schedule(static)
      for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j)
          dscl[index4D(i,j,k,l)] = dscl[index4D(i,j,k,l)]
                        + ( coef_diff[index5D(i,j,5,XDIR,l)] * ( vt[index4D(i  ,j-1,XDIR,TJ)] + vt[index4D(i  ,j  ,XDIR,TI)] )
                          + coef_diff[index5D(i,j,5,YDIR,l)] * ( vt[index4D(i  ,j-1,YDIR,TJ)] + vt[index4D(i  ,j  ,YDIR,TI)] )
                          + coef_diff[index5D(i,j,5,ZDIR,l)] * ( vt[index4D(i  ,j-1,ZDIR,TJ)] + vt[index4D(i  ,j  ,ZDIR,TI)] )
                          ) * 0.5 * ( kh[index4D(i  ,j  ,k,l)] + kh[index4D(i+1,j  ,k,l)] );

      #pragma omp workshare
      {                        
        for (int i = 0; i < ADM_iall; ++i) {
          dscl[index4D(i,jmin-1,k,l)] = 0.0;
          dscl[index4D(i,jmax+1,k,l)] = 0.0;
        }
        for (int j = 0; j < ADM_jall; ++j) {
          dscl[index4D(imin-1,j,k,l)] = 0.0;
          dscl[index4D(imax+1,j,k,l)] = 0.0;
        }
      }
    } // k loop
   } // l loop
  } // end omp parallel
  
  if ( ADM_have_pl ) {
      int n = ADM_gslf_pl;
      for (int l = 0; l < ADM_lall_pl; ++l) {
        for (int k = 0; k < ADM_kall; ++k) {
          for (int v = ADM_gmin_pl; v < ADM_gmax_pl; ++v) {
            int ij   = v;
            int ijp1 = v + 1;
            if( ijp1 == ADM_gmax_pl+1 ) 
              ijp1 = ADM_gmin_pl;
            for (int d = 0; d < ADM_nxyz; ++d) {
              vt_pl[index2D(ij,d)] = ( ( + 2.0 * coef_intp_pl[index4D(v,0,d,l)]
                                - 1.0 * coef_intp_pl[index4D(v,1,d,l)]
                                - 1.0 * coef_intp_pl[index4D(v,2,d,l)] ) * scl_pl[index3D(n   ,k,l)]
                            + ( - 1.0 * coef_intp_pl[index4D(v,0,d,l)]
                                + 2.0 * coef_intp_pl[index4D(v,1,d,l)]
                                - 1.0 * coef_intp_pl[index4D(v,2,d,l)] ) * scl_pl[index3D(ij  ,k,l)]
                            + ( - 1.0 * coef_intp_pl[index4D(v,0,d,l)]
                                - 1.0 * coef_intp_pl[index4D(v,1,d,l)]
                                + 2.0 * coef_intp_pl[index4D(v,2,d,l)] ) * scl_pl[index3D(ijp1,k,l)]
                            ) / 3.0;
            }
        }
        for (int i = 0; i < ADM_iall; ++i)
          dscl_pl[index3D(i,k,l)] = 0.0;
        for (int v = ADM_gmin_pl; v < ADM_gmax_pl; ++v) {
            int ij   = v;
            int ijm1 = v - 1;
            if( ijm1 == ADM_gmin_pl-1 ) 
              ijm1 = ADM_gmax_pl; //cyclic condition
            dscl_pl[index3D(n,k,l)] = dscl_pl[index3D(n,k,l)]
                          + ( coef_diff_pl[index3D(v-1,XDIR,l)] * ( vt_pl[index2D(ijm1,XDIR)] + vt_pl[index2D(ij,XDIR)] )
                            + coef_diff_pl[index3D(v-1,YDIR,l)] * ( vt_pl[index2D(ijm1,YDIR)] + vt_pl[index2D(ij,YDIR)] )
                            + coef_diff_pl[index3D(v-1,ZDIR,l)] * ( vt_pl[index2D(ijm1,ZDIR)] + vt_pl[index2D(ij,ZDIR)] )
                            ) * 0.5 * ( kh_pl[index3D(n,k,l)] + kh_pl[index3D(ij,k,l)] );
          }
        } //end k
      } //end l
     } // end if 
  else
      memset(dscl_pl,0.0,sizeof(ADM_gall_pl*ADM_kall*ADM_lall_pl));
  
  
  end_c = clock();
  float time1 = ((float)(end_c-start_c))/CLOCKS_PER_SEC;
  printf("Runtime: %f", time1);
  
  if(read_from_dump) {
    free(ORG_dscl);
    free(ORG_scl);
    free(ORG_kh);
    free(ORG_coef_intp);
    free(ORG_coef_diff);
  }
  
  free(dscl);
  free(scl);
  free(kh);
  free(coef_intp);
  free(coef_diff);
  free(vt);
  free(dscl_pl);
  free(scl_pl);
  free(kh_pl);
  free(coef_intp_pl);
  free(coef_diff_pl);
  free(vt_pl);



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
