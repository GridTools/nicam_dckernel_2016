#define PEDANTIC_DISABLED
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "time.h"
#include "tools/verifier.hpp"
#include <stencil-composition/stencil-composition.hpp>
#include "problem_size.hpp"
#ifdef __CUDACC__
#include <boundary-conditions/apply_gpu.hpp>
#else
#include <boundary-conditions/apply.hpp>
#endif

using namespace gridtools;
using namespace expressions;
using gridtools::direction;
using gridtools::sign;
using gridtools::minus_;
using gridtools::zero_;
using gridtools::plus_;
  
 
 // single interval for vertical direction
 typedef gridtools::interval< level< 0, -1 >, level< 1, -1 > > x_interval;
 
 // only one axis, but different interval for different equations   
 typedef gridtools::interval< level< 0, -2 >, level< 1, 1 > > axis;


struct flux_1 {
  // Note: This operation runs only up to Kall-1    
  // Make sure to account for that
  // rhogw_vm   --> [I, J, K]            
  typedef accessor<0, enumtype::inout, extent<0,0,0,0>, 4 > rhogw_vm;
  // C2WfactGz  --> [I, J, K, A]+[L]            
  typedef accessor<1, enumtype::in, extent<0,0,0,0>, 6 > C2WfactGz;
  // rhogvx     --> [I, J, K]+[L]            
  typedef accessor<2, enumtype::in, extent<0,0,0,0>, 5 > rhogvx;
  // rhogvy     --> [I, J, K]+[L]            
  typedef accessor<3, enumtype::in, extent<0,0,0,0>, 5 > rhogvy;
  // rhogvz     --> [I, J, K]+[L]            
  typedef accessor<4, enumtype::in, extent<0,0,0,0>, 5 > rhogvz;
  // RGAMH      --> [I, J, K]+[L]            
  typedef accessor<5, enumtype::in, extent<0,0,0,0>, 5 > RGAMH;
  // rhogw      --> [I, J, K]+[L]            
  typedef accessor<6, enumtype::in, extent<0,0,0,0>, 5 > rhogw;
  // RGSQRTH    --> [I, J, K]+[L]            
  typedef accessor<7, enumtype::in, extent<0,0,0,0>, 5 > RGSQRTH;
  
  

  typedef boost::mpl::vector<rhogw_vm, C2WfactGz, rhogvx, rhogvy, rhogvz, RGAMH, rhogw, RGSQRTH> arg_list;

  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<4>::Index a;
        
        eval(rhogw_vm{}) =  ( eval(C2WfactGz{}) * eval(rhogvx{}) 
                            + eval(C2WfactGz{a+1}) * eval(rhogvx{k-1}) 
                            + eval(C2WfactGz{a+2}) * eval(rhogvy{}) 
                            + eval(C2WfactGz{a+3}) * eval(rhogvy{k-1}) 
                            + eval(C2WfactGz{a+4}) * eval(rhogvz{}) 
                            + eval(C2WfactGz{a+5}) * eval(rhogvz{k-1}) 
                            ) * eval(RGAMH{})                  // horizontal contribution
                          + eval(rhogw{}) * eval(RGSQRTH{});   // vertical   contribution

    }
}; // end flux_function

struct flux_2 {
    
  // Note: This operation runs only from 2 to Kall-1 (kmin to kmax)   
  // Make sure to account for that
  // rhogvx_vm  --> [I, J]            
  typedef accessor<0, enumtype::inout, extent<0,0,0,0>, 3 > rhogvx_vm;
  // rhogvy_vm  --> [I, J]            
  typedef accessor<1, enumtype::inout, extent<0,0,0,0>, 3 > rhogvy_vm;
  // rhogvz_vm  --> [I, J]
  typedef accessor<2, enumtype::inout, extent<0,0,0,0>, 3 > rhogvz_vm;
  // rhogvx     --> [I, J, K]+[L]            
  typedef accessor<3, enumtype::in, extent<0,0,0,0>, 5 > rhogvx;
  // rhogvy     --> [I, J, K]+[L]            
  typedef accessor<4, enumtype::in, extent<0,0,0,0>, 5 > rhogvy;
  // rhogvz     --> [I, J, K]+[L]            
  typedef accessor<5, enumtype::in, extent<0,0,0,0>, 5 > rhogvz;
  // RGAM      --> [I, J, K]+[L]            
  typedef accessor<6, enumtype::in, extent<0,0,0,0>, 5 > RGAM;
    
  typedef boost::mpl::vector<rhogvx_vm, rhogvy_vm, rhogvz_vm, rhogvx, rhogvy, rhogvz, RGAM> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        
        eval(rhogvx_vm{}) = eval(rhogvx{}) * eval(RGAM{});
        eval(rhogvy_vm{}) = eval(rhogvy{}) * eval(RGAM{});
        eval(rhogvz_vm{}) = eval(rhogvz{}) * eval(RGAM{});        
        
    }
}; // end flux_function_2

struct flux_3 {
    
  // Note: This operation runs only from 2 to Kall-1 (kmin to kmax)   
  // Make sure to account for that
  // sclt       --> [I, J, T]            
  typedef accessor<0, enumtype::inout, extent<0,-1,0,-1>, 4 > sclt;
  // rhogw_vm   --> [I, J, K]            
  typedef accessor<1, enumtype::in, extent<0,-1,0,-1>, 4 > rhogw_vm;
  // GRD_rdgz   --> [K]            
  typedef accessor<2, enumtype::in, extent<0,-1,0,-1>, 2 > GRD_rdgz;
  // coef_intp   --> [I, J, A, D, T] + [L]           
  typedef accessor<3, enumtype::in, extent<0,-1,0,-1>, 7 > coef_intp;
  // rhogvx_vm  --> [I, J]            
  typedef accessor<4, enumtype::inout, extent<0,-1,0,-1>, 3 > rhogvx_vm;
  // rhogvy_vm  --> [I, J]            
  typedef accessor<5, enumtype::inout, extent<0,-1,0,-1>, 3 > rhogvy_vm;
  // rhogvz_vm  --> [I, J]
  typedef accessor<6, enumtype::inout, extent<0,-1,0,-1>, 3 > rhogvz_vm;
  
  typedef boost::mpl::vector<sclt, rhogw_vm, GRD_rdgz, coef_intp, rhogvx_vm, rhogvy_vm, rhogvz_vm> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<1>::Index k_GRD_rdgz;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index a;
        dimension<3>::Index t;
        dimension<4>::Index d;
        dimension<5>::Index t_coef_intp;
        
        eval(sclt{}) = eval(coef_intp{})        * eval(rhogvx_vm{})
                     + eval(coef_intp{a+1})     * eval(rhogvx_vm{i+1})
                     + eval(coef_intp{a+2})     * eval(rhogvx_vm{i+1,j+1})
                     + eval(coef_intp{d+1})     * eval(rhogvy_vm{})
                     + eval(coef_intp{a+1,d+1}) * eval(rhogvy_vm{i+1})
                     + eval(coef_intp{a+1,d+1}) * eval(rhogvy_vm{i+1,j+1})
                     + eval(coef_intp{d+2})     * eval(rhogvz_vm{})
                     + eval(coef_intp{a+1,d+2}) * eval(rhogvz_vm{i+1})
                     + eval(coef_intp{a+2,d+2}) * eval(rhogvz_vm{i+1,j+1})
                     + ( ( eval(rhogw_vm{k+1})  + eval(rhogw_vm{i+1,k+1}) + eval(rhogw_vm{i+1,j+1,k+1}))
                       - ( eval(rhogw_vm{})     + eval(rhogw_vm{i+1})     + eval(rhogw_vm{i+1,j+1}))
                          ) / 3.0 * eval(GRD_rdgz{});

        eval(sclt{t+1}) = eval(coef_intp{t+1})      * eval(rhogvx_vm{})
                     + eval(coef_intp{a+1,t+1})     * eval(rhogvx_vm{i+1,j+1})
                     + eval(coef_intp{a+2,t+1})     * eval(rhogvx_vm{j+1})
                     + eval(coef_intp{d+1,t+1})     * eval(rhogvy_vm{})
                     + eval(coef_intp{a+1,d+1,t+1}) * eval(rhogvy_vm{i+1,j+1})
                     + eval(coef_intp{a+2,d+1,t+1}) * eval(rhogvy_vm{j+1})
                     + eval(coef_intp{d+2,t+1})     * eval(rhogvz_vm{})
                     + eval(coef_intp{a+1,d+2,t+1}) * eval(rhogvz_vm{i+1,j+1})
                     + eval(coef_intp{a+2,d+2,t+1}) * eval(rhogvz_vm{j+1})
                     + ( ( eval(rhogw_vm{k+1})      + eval(rhogw_vm{i+1,j+1,k+1}) + eval(rhogw_vm{j+1,k+1}))
                       - ( eval(rhogw_vm{})         + eval(rhogw_vm{i+1,j+1})      + eval(rhogw_vm{j+1}))
                          ) / 3.0 * eval(GRD_rdgz{});             

        
    }
}; // end flux_3

struct flux_4 {
    
  // ddivdx        --> [I, J, K]+[L]            
  typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > ddivdx;
  // ddivdy        --> [I, J, K]+[L]                       
  typedef accessor<1, enumtype::inout, extent<1,-1,1,-1>, 5 > ddivdy;
  // ddivdz        --> [I, J, K]+[L]              
  typedef accessor<2, enumtype::inout, extent<1,-1,1,-1>, 5 > ddivdz;
  // coef_diff     --> [I, J, C, D]+[L]
  typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 6 > coef_diff;
  // sclt       --> [I, J, T]            
  typedef accessor<4, enumtype::in, extent<1,-1,1,-1>, 4 > sclt;
    
  typedef boost::mpl::vector<ddivdx, ddivdy, ddivdz, coef_diff, sclt> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index a;
        dimension<3>::Index t;
        dimension<4>::Index d;

        eval(ddivdx{}) = eval(coef_diff{})        * ( eval(sclt{})           + eval(sclt{t+1}) )
                       + eval(coef_diff{a+1})     * ( eval(sclt{t+1})        + eval(sclt{i-1,j}) ) 
                       + eval(coef_diff{a+2})     * ( eval(sclt{i-1,j})      + eval(sclt{i-1,j-1,t+1}) ) 
                       + eval(coef_diff{a+3})     * ( eval(sclt{i-1,j-1,t+1})+ eval(sclt{i-1,j-1}) ) 
                       + eval(coef_diff{a+4})     * ( eval(sclt{i-1,j-1})    + eval(sclt{j-1,t+1}) ) 
                       + eval(coef_diff{a+5})     * ( eval(sclt{j-1,t+1})    + eval(sclt{}) );
        
        eval(ddivdy{}) = eval(coef_diff{d+1})     * ( eval(sclt{})           + eval(sclt{t+1}) ) 
                       + eval(coef_diff{a+1,d+1}) * ( eval(sclt{t+1})        + eval(sclt{i-1}) ) 
                       + eval(coef_diff{a+2,d+1}) * ( eval(sclt{i-1})        + eval(sclt{i-1,j-1,t+1}) ) 
                       + eval(coef_diff{a+3,d+1}) * ( eval(sclt{i-1,j-1,t+1})+ eval(sclt{i-1,j-1}) ) 
                       + eval(coef_diff{a+4,d+1}) * ( eval(sclt{i-1,j-1})    + eval(sclt{j-1,t+1}) ) 
                       + eval(coef_diff{a+5,d+1}) * ( eval(sclt{j-1,t+1})    + eval(sclt{}) );

        eval(ddivdz{}) = eval(coef_diff{d+2})     * ( eval(sclt{})           + eval(sclt{t+1}) ) 
                       + eval(coef_diff{a+1,d+2}) * ( eval(sclt{t+1})        + eval(sclt{i-1}) ) 
                       + eval(coef_diff{a+2,d+2}) * ( eval(sclt{i-1})        + eval(sclt{i-1,j-1,t+1}) ) 
                       + eval(coef_diff{a+3,d+2}) * ( eval(sclt{i-1,j-1,t+1})+ eval(sclt{i-1,j-1}) ) 
                       + eval(coef_diff{a+4,d+2}) * ( eval(sclt{i-1,j-1})    + eval(sclt{j-1,t+1}) ) 
                       + eval(coef_diff{a+5,d+2}) * ( eval(sclt{j-1,t+1})    + eval(sclt{}) );
    }
}; // end flux_4

//// Boundary condition functor
//template <typename T>
//struct direction_bc_input {
//    T value;

//    GT_FUNCTION
//    direction_bc_input()
//        : value(1)
//    {}

//    GT_FUNCTION
//    direction_bc_input(T v)
//        : value(v)
//    {}

//    // relative coordinates
//  template <sign I, sign J, sign K, typename DataField0>
//    GT_FUNCTION
//  void operator()(direction<I,J,K>,
//                    DataField0 & data_field0,
//                    uint_t i, uint_t j, uint_t k) const {
//        data_field0(i,j,k) = value;
//    }
//}; // end bounrady condition functor

int main(int argc, char const *argv[])
{
  
    using namespace enumtype;
    clock_t start, end;

    // [layout_map]
    //CUDA case
#ifdef __CUDACC__    
    typedef gridtools::layout_map<1>            layout_1d;
    typedef gridtools::layout_map<3,0>          layout_2d;
    typedef gridtools::layout_map<3,0,1>        layout_3d;
    typedef gridtools::layout_map<3,0,1,2>      layout_3d_l;
    typedef gridtools::layout_map<3,0,-1,1,2>   layout_4d;
    typedef gridtools::layout_map<4,3,-1,2,1,0> layout_5d;
#else    
    typedef gridtools::layout_map<1>            layout_1d;
    typedef gridtools::layout_map<0,1>          layout_2d;
    typedef gridtools::layout_map<0,1,2,3>      layout_3d_l;
    typedef gridtools::layout_map<0,1,2>        layout_3d;
    typedef gridtools::layout_map<2,3,-1,0,1>   layout_4d;
    typedef gridtools::layout_map<3,4,-1,0,1,2> layout_5d;
#endif

    // CUDA backend: use block strategy
#ifdef __CUDACC__
    typedef backend<Cuda,structured,Block> backend_t;
#else
    typedef backend<Host,structured,Block> backend_t;
#endif
    
    typedef backend_t::storage_info<0, layout_1d> storage_info_1d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_1d_t>::type storage_type_1d;    

    typedef backend_t::storage_info<1, layout_2d> storage_info_2d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_2d_t>::type storage_type_2d;    
    
    typedef backend_t::storage_info<2, layout_3d> storage_info_3d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_t>::type storage_type_3d;

    typedef backend_t::storage_info<3, layout_3d_l> storage_info_3d_t_l;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_t_l>::type storage_type_3d_l;

    typedef backend_t::storage_info<4, layout_4d> storage_info_4d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_4d_t>::type storage_type_4d;

    typedef backend_t::storage_info<5, layout_4d> storage_info_4d_coef_diff_t;
    typedef typename backend_t::storage_type<float_type, storage_info_4d_coef_diff_t>::type storage_type_4d_coef_diff;
     
    typedef backend_t::storage_info<6, layout_5d> storage_info_5d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_5d_t>::type storage_type_5d;

    // GRD_rdgz
    storage_info_1d_t           metadata_1d(ADM_kall);
    // rhogvx_vm & rhogvy_vm & rhogvz_vm
    storage_info_2d_t           metadata_2d(ADM_iall,ADM_jall);
    // rhogw_vm & sclt
    storage_info_3d_t           metadata_3d(ADM_iall,ADM_jall,ADM_kall);
    // rhogvx & rhogvy & rhogvz & RGAMH & rhogw & RGSQRTH & RGAM & ddivdx & ddivdy & ddivdz
    storage_info_3d_t_l         metadata_3d_l(ADM_iall,ADM_jall,ADM_kall);
    // coef_diff
    storage_info_4d_coef_diff_t metadata_4d_coef_diff(ADM_iall,ADM_jall,6,ADM_nxyz);
    // C2WfactGz
    //storage_info_4d_t           metadata_4d(ADM_iall+1,ADM_jall+1,1,ADM_nxyz,TJ);
    storage_info_4d_t           metadata_4d(ADM_iall,ADM_jall,ADM_kall,6);
    // coef_intp
    storage_info_5d_t           metadata_5d(ADM_iall,ADM_jall,3,ADM_nxyz,TJ);
    
    // flux_1 functor
    storage_type_3d                                   rhogw_vm(metadata_3d, 0.0, "rhogw_vm");
    field <storage_info_4d_t,ADM_lall>::type          C2WfactGz(metadata_4d, 0.0, "C2WfactGz");
    field <storage_info_3d_t_l,ADM_lall>::type        rhogvx(metadata_3d_l, 0.0, "rhogvx");
    field <storage_info_3d_t_l,ADM_lall>::type        rhogvy(metadata_3d_l, 0.0, "rhogvy");
    field <storage_info_3d_t_l,ADM_lall>::type        rhogvz(metadata_3d_l, 0.0, "rhogvz");
    field <storage_info_3d_t_l,ADM_lall>::type        RGAMH(metadata_3d_l, 0.0, "RGAMH");
    field <storage_info_3d_t_l,ADM_lall>::type        rhogw(metadata_3d_l, 0.0, "rhogw");
    field <storage_info_3d_t_l,ADM_lall>::type        RGSQRTH(metadata_3d_l, 0.0, "RGSQRTH");
    
    // flux_2 functor
    storage_info_2d_t                                 rhogvx_vm(metadata_2d, 0.0, "rhogvx_vm");
    storage_info_2d_t                                 rhogvy_vm(metadata_2d, 0.0, "rhogvy_vm");
    storage_info_2d_t                                 rhogvz_vm(metadata_2d, 0.0, "rhogvz_vm");
    field <storage_info_3d_t_l,ADM_lall>::type        RGAM(metadata_3d_l, 0.0, "RGAM");
    
    // flux_3 functor
    storage_info_3d_t                                 sclt(metadata_3d, 0.0, "sclt");
    storage_info_3d_t                                 rhogw_vm(metadata_3d, 0.0, "rhogw_vm");
    storage_info_1d_t                                 GRD_rdgz(metadata_1d, 0.0, "GRD_rdgz");
    field <storage_info_5d_t,ADM_lall>::type          coef_intp(metadata_5d, 0.0, "coef_intp");

    // flux_4 functor
    field <storage_type_3d_l,ADM_lall>::type          ddivdx(metadata_3d_l, 0.0, "ddivdx");
    field <storage_type_3d_l,ADM_lall>::type          ddivdy(metadata_3d_l, 0.0, "ddivdy");
    field <storage_type_3d_l,ADM_lall>::type          ddivdz(metadata_3d_l, 0.0, "ddivdz");
    field <storage_type_4d_coef_diff,ADM_lall>::type  coef_diff(metadata_4d_coef_diff, 0.0, "coef_diff");
        
    // Use following To fill "scl" up the values from dumpio.c
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
   for(int i=0; i<metadata_3d.template dims<0>(); ++i)
    for(int j=0; j<metadata_3d.template dims<1>(); ++j)
      for(int k=0; k<metadata_3d.template dims<2>(); ++k)
      {
        ddivdx.get_value<0>(i,j,k)= *(float_type*)((ORG_ddivdx + index4D(i,j,k,l)));
        ddivdy.get_value<0>(i,j,k)= *(float_type*)((ORG_ddivdy + index4D(i,j,k,l)));
        ddivdz.get_value<0>(i,j,k)= *(float_type*)((ORG_ddivdz + index4D(i,j,k,l)));
        rhogvx.get_value<0>(i,j,k)= *(float_type*)((ORG_rhogvx + index4D(i,j,k,l)));
        rhogvy.get_value<0>(i,j,k)= *(float_type*)((ORG_rhogvy + index4D(i,j,k,l)));
        rhogvz.get_value<0>(i,j,k)= *(float_type*)((ORG_rhogvz + index4D(i,j,k,l)));
        RGSQRTH.get_value<0>(i,j,k)=*(float_type*)((ORG_RGSQRTH+ index4D(i,j,k,l)));
        RGAM.get_value<0>(i,j,k)   =*(float_type*)((ORG_RGAM   + index4D(i,j,k,l)));
        RGAMH.get_value<0>(i,j,k)  =*(float_type*)((ORG_RGAMH  + index4D(i,j,k,l)));
        for(int d=0; d<metadata_4d.template dims<3>(); ++d)      
            C2WfactGz.get_value<0>(i,j,k,d)=* (float_type*)((ORG_C2WfactGz + index5D(i,j,k,d,l)));
        
      }
   for(int i=0; i<metadata_5d.template dims<0>(); ++i)
    for(int j=0; j<metadata_5d.template dims<1>(); ++j)
    {
      // the k dimension is killed
      for(int a=0; a<metadata_5d.template dims<3>(); ++a)
         for(int d=0; d<metadata_5d.template dims<4>(); ++d)      
            for(int t=0; t<metadata_5d.template dims<5>(); ++t)      
              coef_intp.get_value<0>(i,j,1,a,d,t)=* (float_type*)((ORG_coef_intp + index6D(i,j,a,d,t,l)));
      // the k dimension is killed
      for(int a=0; a<metadata_4d_coef_diff.template dims<3>(); ++a)
         for(int d=0; d<metadata_4d_coef_diff.template dims<4>(); ++d)      
            coef_diff.get_value<0>(i,j,1,a,d)=* (float_type*)((ORG_coef_diff + index5D(i,j,a,d,l)));
    }
  }  // l loop
  for(int i=0; i<metadata_3d.template dims<0>(); ++i)
    for(int j=0; j<metadata_3d.template dims<1>(); ++j)
      for(int k=0; k<metadata_3d.template dims<2>(); ++k)
        rhogw.get_value<0>(i,j,k)=* (float_type*)((ORG_rhogw + index4D(i,j,k,l)));

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
   for(int i=0; i<metadata_3d_l.template dims<0>(); ++i)
    for(int j=0; j<metadata_3d_l.template dims<1>(); ++j)
      for(int k=0; k<metadata_3d_l.template dims<2>(); ++k)
      {
        ddivdx.get_value<0>(i,j,k) = (i+j+k)/2;
        ddivdy.get_value<0>(i,j,k) = (i+j+k)/2;
        ddivdz.get_value<0>(i,j,k) = (i+j+k)/4;
        rhogvx.get_value<0>(i,j,k) = (i+j+k)/4;
        rhogvy.get_value<0>(i,j,k) = (i+j+k)/8;
        rhogvz.get_value<0>(i,j,k) = (i+j+k)/8;
        RGSQRTH.get_value<0>(i,j,k)= (i+j+k)/2;
        RGAM.get_value<0>(i,j,k)   = (i+j+k)/2;
        RGAMH.get_value<0>(i,j,k)  = (i+j+k)/2;
        rhogw.get_value<0>(i,j,k)  = (i+j+k)/4;
        for(int d=0; d<metadata_4d.template dims<3>(); ++d)      
            C2WfactGz.get_value<0>(i,j,k,d)=(i+j+k+d)/2;
        
      }
   for(int i=0; i<metadata_5d.template dims<0>(); ++i)
    for(int j=0; j<metadata_5d.template dims<1>(); ++j)
    {
      // the k dimension is killed
      for(int a=0; a<metadata_5d.template dims<3>(); ++a)
         for(int d=0; d<metadata_5d.template dims<4>(); ++d)      
            for(int t=0; t<metadata_5d.template dims<5>(); ++t)      
              coef_intp.get_value<0>(i,j,1,a,d,t)=(i+j+a+d)/2;
      // the k dimension is killed
      for(int a=0; a<metadata_4d_coef_diff.template dims<3>(); ++a)
         for(int d=0; d<metadata_4d_coef_diff.template dims<4>(); ++d)      
            coef_diff.get_value<0>(i,j,1,a,d)=(i+j+d)/2;
    }
  }  // l loop
  for(int i=0; i<metadata_3d.template dims<0>(); ++i)
    for(int j=0; j<metadata_3d.template dims<1>(); ++j) {
      rhogvx_vm.get_value<0>(i,j)=(i+j)/2;
      rhogvy_vm.get_value<0>(i,j)=(i+j)/4;
      rhogvz_vm.get_value<0>(i,j)=(i+j)/8;
      for(int k=0; k<metadata_3d.template dims<2>(); ++k)
        rhogw_vm.get_value<0>(i,j,k)=(i+j+k)/2;
      for(int t=0; t<metadata_3d.template dims<2>(); ++t)  
        sclt.get_value<0>(i,j,t)=(i+j+t)/2;
    }
  for(int k=0; k<metadata_1d.template dims<0>(); ++k) 
    GRD_rdgz.get_value<0>(k)=(k)/2;   
      

    // functor_1
    typedef arg<0, storage_type_3d>                           p_rhogw_vm;
    typedef arg<1, field <storage_info_4d_t,ADM_lall>::type>  p_C2WfactGz;
    typedef arg<2, field <storage_info_3d_t,ADM_lall>::type>  p_rhogvx;
    typedef arg<3, field <storage_info_3d_t,ADM_lall>::type>  p_rhogvy;
    typedef arg<4, field <storage_info_3d_t,ADM_lall>::type>  p_rhogvz;
    typedef arg<5, field <storage_info_3d_t,ADM_lall>::type>  p_RGAMH;
    typedef arg<6, field <storage_info_3d_t,ADM_lall>::type>  p_rhogw;
    typedef arg<7, field <storage_info_3d_t,ADM_lall>::type>  p_RGSQRTH;
    
    // functor_2
    typedef arg<8,  storage_info_2d_t>                        p_rhogvx_vm;
    typedef arg<9,  storage_info_2d_t>                        p_rhogvy_vm;
    typedef arg<10, storage_info_2d_t>                        p_rhogvz_vm;
    typedef arg<11, field <storage_info_3d_t,ADM_lall>::type> p_rhogvx;
    typedef arg<12, field <storage_info_3d_t,ADM_lall>::type> p_rhogvy;
    typedef arg<13, field <storage_info_3d_t,ADM_lall>::type> p_rhogvz;
    typedef arg<14, field <storage_info_3d_t,ADM_lall>::type> p_RGAM;

    // functor_3
    typedef arg<15, storage_info_3d_t>                        p_sclt;
    typedef arg<16, storage_info_3d_t>                        p_rhogw_vm;
    typedef arg<17, storage_info_1d_t>                        p_GRD_rdgz;
    typedef arg<18, field <storage_info_5d_t,ADM_lall>::type> p_coef_intp;
    typedef arg<19, storage_info_2d_t>                        p_rhogvx_vm;
    typedef arg<20, storage_info_2d_t>                        p_rhogvy_vm;
    typedef arg<21, storage_info_2d_t>                        p_rhogvz_vm;

    // fucntor_4
    typedef arg<22, field <storage_type_3d,ADM_lall>::type>           p_ddivdx;
    typedef arg<23, field <storage_type_3d,ADM_lall>::type>           p_ddivdy;
    typedef arg<24, field <storage_type_3d,ADM_lall>::type>           p_ddivdz;
    typedef arg<25, field <storage_type_4d_coef_diff,ADM_lall>::type> p_coef_diff;
    typedef arg<26, storage_info_3d_t>                                p_sclt;
    
    // accessor list of all place holders used
    typedef boost::mpl::vector<p_rhogw_vm, p_C2WfactGz, p_rhogvx, p_rhogvy, p_rhogvz, p_RGAMH, p_rhogw, p_RGSQRTH, p_rhogvx_vm, p_rhogvy_vm, p_rhogvz_vm, p_RGAM, p_sclt, p_rhogw_vm, p_GRD_rdgz, p_coef_intp, p_ddivdx, p_ddivdy, p_ddivdz, p_coef_diff> accessor_list;

    gridtools::domain_type<accessor_list> domain ((p_rhogw_vm() = rhogw_vm), (p_C2WfactGz() = C2WfactGz), (p_rhogvx() = rhogvx), (p_rhogvy()  = rhogvy), (p_rhogvz() = rhogvz), (p_RGAMH() = RGAMH), (p_rhogw() = rhogw), (p_RGSQRTH() = RGSQRTH), (p_rhogvx_vm() = rhogvx_vm), (p_rhogvy_vm() = rhogvy_vm), (p_rhogvz_vm() = rhogvz_vm), (p_RGAM() = RGAM), (p_sclt() = sclt), (p_rhogw_vm() = rhogw_vm), (p_GRD_rdgz() = GRD_rdgz), (p_coef_intp() = coef_intp), (p_ddivdx() = ddivdx), (p_ddivdy() = ddivdy), (p_ddivdz() = ddivdz), (p_coef_diff() = coef_diff));

    // i-1,j-1 ---- i+1, j+1
    // extra lower halo added to account for operation dependecy
    uint_t halo_size = 1; 
    uint_t di[5] = {halo_size*2, halo_size, halo_size*2,(uint_t) ADM_iall-halo_size-1,(uint_t) ADM_iall};
    uint_t dj[5] = {halo_size*2, halo_size, halo_size*2,(uint_t) ADM_jall-halo_size-1,(uint_t) ADM_jall};
    gridtools::grid<axis> grid(di,dj);
    

    // single splitter
    grid.value_list[0] = 0;
    grid.value_list[1] = ADM_kall-1;

    auto divdamp3d = make_computation<backend_t>
        (
         domain, grid,
         make_mss
         (
          execute<forward>(), // parallel means no vertical dependency, forward from 0 to K-1, otherwise it is backward
          make_esf<flux_1>(p_rhogw_vm(), p_C2WfactGz(), p_rhogvx(), p_rhogvy(), p_rhogvz(), p_RGAMH(), p_rhogw(), p_RGSQRTH()),//define in same order as accessors unique identifier          
          make_esf<flux_2>(p_rhogvx_vm(), p_rhogvy_vm(), p_rhogvz_vm(), p_rhogvx(), p_rhogvy(), p_rhogvz(), p_RGAM()),//define in same order as accessors unique identifier          
          make_esf<flux_3>(p_sclt(), p_rhogw_vm(), p_GRD_rdgz(), p_coef_intp(), p_rhogvx_vm(), p_rhogvy_vm(), p_rhogvz_vm()),//define in same order as accessors unique identifier          
          make_esf<flux_4>(p_ddivdx(), p_ddivdy(), p_ddivdz(), p_coef_diff(), p_sclt())//define in same order as accessors unique identifier                 
          )
         );

    // allocate and init storages
    divdamp3d->ready();

    // copy to device
    divdamp3d->steady();

    start = clock();
    for(int iter=0; iter < SET_iteration; ++iter) {  
      divdamp3d->run(); 
      //ToDo: call on DEBUG_valuecheck
    }
    end = clock();
    float time1 = ((float)(end-start))/CLOCKS_PER_SEC;
    printf("Runtime: %f", time1);

    // run boundary condition
    //gridtools::array<gridtools::halo_descriptor, 3> halos;
    //halos[0] = gridtools::halo_descriptor(0,0,0,ADM_iall-1,ADM_iall);
    //halos[1] = gridtools::halo_descriptor(0,0,0,ADM_jall-1,ADM_jall);
    //halos[2] = gridtools::halo_descriptor(0,0,0,ADM_kall-1,ADM_kall);
//#ifdef __CUDACC__
//    gridtools::boundary_apply_gpu<direction_bc_input<float_type> >(halos, direction_bc_input<float_type>(0.0)).apply(dscl);
//#else
//    gridtools::boundary_apply<direction_bc_input<float_type> >(halos, direction_bc_input<float_type>(0.0)).apply(dscl);
//#endif

    // copy data back and deallocate
    divdamp3d->finalize();

    //dscl.print();

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
