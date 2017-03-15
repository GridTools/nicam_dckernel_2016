
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


struct flux_function {
      
  // vt --> [I, J, D, T]            
  typedef accessor<0, enumtype::inout, extent<0,-1,0,-1>, 5 > vt; //<-1,0,-1,0>
  //typedef accessor<0, enumtype::inout, extent<0,-1,0,-1>, 5 > vt;
  // scl --> [I, J, K]+[L] 
  typedef accessor<1, enumtype::in, extent<0, -1, 0, -1>, 5 > scl; //<-1, 1, -1, 1>
  //typedef accessor<1, enumtype::in, extent<0, -1, 0, -1>, 5 > scl; //<-1, 1, -1, 1>
  // coef_intp   [I, J, A, D, T]+[L]
  typedef accessor<2, enumtype::in, extent<0, -1, 0, -1>, 7 > coef_intp; //<-1, 0, -1, 0>
  //typedef accessor<2, enumtype::in, extent<0, -1, 0, -1>, 7 > coef_intp; //<-1, 0, -1, 0>
  

  typedef boost::mpl::vector<vt, scl, coef_intp> arg_list;

  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index d_vt;
        dimension<3>::Index a;
        dimension<4>::Index t;
        dimension<4>::Index d_coef_intp;
        dimension<5>::Index t_coef_intp;
        

        
        // TI
        eval(vt{}) =      (( + 2.0 * eval(coef_intp{}) 
                             - 1.0 * eval(coef_intp{a+1}) 
                             - 1.0 * eval(coef_intp{a+2})  ) * eval(scl{}) 
                         + ( - 1.0 * eval(coef_intp{})
                             + 2.0 * eval(coef_intp{a+1})
                             - 1.0 * eval(coef_intp{a+2}) ) * eval(scl{i+1}) //i+1}) 
                         + ( - 1.0 * eval(coef_intp{})
                             - 1.0 * eval(coef_intp{a+1})
                             + 2.0 * eval(coef_intp{a+2}) ) * eval(scl{i+1,j+1}) //i+1,j+1}) 
                         ) / 3.0;

        // TJ    
        eval(vt{t+1}) =   (( + 2.0 * eval(coef_intp{t+1}) 
                             - 1.0 * eval(coef_intp{a+1,t+1}) 
                             - 1.0 * eval(coef_intp{a+2,t+1})  ) * eval(scl{}) 
                         + ( - 1.0 * eval(coef_intp{t_coef_intp+1})
                             + 2.0 * eval(coef_intp{a+1,t_coef_intp+1})
                             - 1.0 * eval(coef_intp{a+2,t_coef_intp+1}) ) * eval(scl{i+1,j+1}) //i+1,j+1}) 
                         + ( - 1.0 * eval(coef_intp{t_coef_intp+1})
                             - 1.0 * eval(coef_intp{a+1,t_coef_intp+1})
                             + 2.0 * eval(coef_intp{a+2,t_coef_intp+1}) ) * eval(scl{j+1}) //j+1}) 
                         ) / 3.0;  

    }
}; // end flux_function

struct consume_function_1 {
    
  // dscl      --> [I, J, K]+[L]            
  typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;
  //typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;
  // vt        --> [I, J, D, T]           
  typedef accessor<1, enumtype::in, extent<1,-1,1,-1>, 5 > vt;
  //typedef accessor<1, enumtype::in, extent<1,-1,1,-1>, 5 > vt;
  // coef_diff -->  [I, J, C, D]+[L]              
  typedef accessor<2, enumtype::in, extent<1,-1,1,-1>, 6 > coef_diff;
  //typedef accessor<2, enumtype::in, extent<1,-1,1,-1>, 6 > coef_diff;
  // kh        --> [I, J, K]+[L]
  typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 5 > kh;
  //typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 5 > kh;
    
  typedef boost::mpl::vector<dscl, vt, coef_diff, kh> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index c; // would be set to 1
        dimension<3>::Index d;
        dimension<4>::Index t;
        dimension<4>::Index d_coef_diff;
        
        eval(dscl{}) = ( eval(coef_diff{})    * ( eval(vt{})    + eval(vt{t+1})) 
                       + eval(coef_diff{d_coef_diff+1}) * ( eval(vt{d+1}) + eval(vt{d+1,t+1})) 
                       + eval(coef_diff{d_coef_diff+2}) * ( eval(vt{d+2}) + eval(vt{d+2,t+1})) 
                      ) * 0.5 * ( eval(kh{})  +   eval(kh{i+1,j+1}));
        
    }
}; // end consume_function_1

struct consume_function_2 {
    
  // dscl      --> [I, J, K]+[L]            
  typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;
  //typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;
  // vt        --> [I, J, D, T]           
  //typedef accessor<1, enumtype::in, extent<-1, 0, 0, 0>, 5 > vt;
  typedef accessor<1, enumtype::in, extent<1,-1,1,-1>, 5 > vt;
  // coef_diff -->  [I, J, C, D]+[L]              
  //typedef accessor<2, enumtype::in, extent<>, 6 > coef_diff;
  typedef accessor<2, enumtype::in, extent<1,-1,1,-1>, 6 > coef_diff;
  // kh        --> [I, J, K]+[L]
  //typedef accessor<3, enumtype::in,    extent<0, 0, 0, 1>, 5 > kh;
  typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 5 > kh;
    
  typedef boost::mpl::vector<dscl, vt, coef_diff, kh> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index c; // would be set to 2
        dimension<3>::Index d;
        dimension<4>::Index t;
        dimension<4>::Index d_coef_diff;
        
        eval(dscl{}) = ( eval(coef_diff{})    * ( eval(vt{t+1})     + eval(vt{i-1})) 
                       + eval(coef_diff{d_coef_diff+1}) * ( eval(vt{d+1,t+1}) + eval(vt{i-1,d+1})) 
                       + eval(coef_diff{d_coef_diff+2}) * ( eval(vt{d+2,t+1}) + eval(vt{i-1,d+2})) 
                      ) * 0.5 * ( eval(kh{})  +   eval(kh{j+1}));
        
    }
}; // end consume_function_2

struct consume_function_3 {
    
  // dscl      --> [I, J, K]+[L]            
  //typedef accessor<0, enumtype::inout, extent<>, 5 > dscl;
  typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;
  // vt        --> [I, J, D, T]           
  //typedef accessor<1, enumtype::in, extent<-1, 0, -1, 0>, 5 > vt;
  typedef accessor<1, enumtype::in, extent<1,-1,1,-1>, 5 > vt;
  // coef_diff -->  [I, J, C, D]+[L]              
  //typedef accessor<2, enumtype::in, extent<>, 6 > coef_diff;
  typedef accessor<2, enumtype::in, extent<1,-1,1,-1>, 6 > coef_diff;
  // kh        --> [I, J, K]+[L]
  //typedef accessor<3, enumtype::in,    extent<-1, 0, 0, 0>, 5 > kh;
  typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 5 > kh;
    
  typedef boost::mpl::vector<dscl, vt, coef_diff, kh> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index c; // would be set to 3
        dimension<3>::Index d;
        dimension<4>::Index t;
        dimension<4>::Index d_coef_diff;
        
        eval(dscl{}) = eval(dscl{}) +
                       ( eval(coef_diff{})    * ( eval(vt{i-1})      + eval(vt{i-1,j-1,t+1})) 
                       + eval(coef_diff{d_coef_diff+1}) * ( eval(vt{i-1,d+1}) + eval(vt{i-1,j-1,d+1,t+1})) 
                       + eval(coef_diff{d_coef_diff+2}) * ( eval(vt{i-1,d+2}) + eval(vt{i-1,j-1,d+2,t+1})) 
                      ) * 0.5 * ( eval(kh{i-1})  +   eval(kh{}));
        
    }
}; // end consume_function_3

struct consume_function_4 {
    
  // dscl      --> [I, J, K]+[L]            
  //typedef accessor<0, enumtype::inout, extent<>, 5 > dscl;
  typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;  
  // vt        --> [I, J, D, T]           
  //typedef accessor<1, enumtype::in, extent<-1, 0, -1, 0>, 5 > vt;
  typedef accessor<1, enumtype::in, extent<1,-1,1,-1>, 5 > vt;
  // coef_diff -->  [I, J, C, D]+[L]              
  //typedef accessor<2, enumtype::in, extent<>, 6 > coef_diff;
  typedef accessor<2, enumtype::in, extent<1,-1,1,-1>, 6 > coef_diff;
  // kh        --> [I, J, K]+[L]
  //typedef accessor<3, enumtype::in,    extent<-1, 0, -1, 0>, 5 > kh;
  typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 5 > kh;
    
  typedef boost::mpl::vector<dscl, vt, coef_diff, kh> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index c; // would be set to 4
        dimension<3>::Index d;
        dimension<4>::Index t;
        dimension<4>::Index d_coef_diff;
        
        eval(dscl{}) = eval(dscl{}) +
                       ( eval(coef_diff{})    * ( eval(vt{i-1,j-1,t+1})      + eval(vt{i-1,j-1})) 
                       + eval(coef_diff{d_coef_diff+1}) * ( eval(vt{i-1,j-1,d+1,t+1}) + eval(vt{i-1,j-1,d+1})) 
                       + eval(coef_diff{d_coef_diff+2}) * ( eval(vt{i-1,j-1,d+2,t+1}) + eval(vt{i-1,j-1,d+2})) 
                      ) * 0.5 * ( eval(kh{i-1,j-1})  +   eval(kh{}));
        
    }
}; // end consume_function_4

struct consume_function_5 {
    
  // dscl      --> [I, J, K]+[L]            
  //typedef accessor<0, enumtype::inout, extent<>, 5 > dscl;
  typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;  
  // vt        --> [I, J, D, T]           
  //typedef accessor<1, enumtype::in, extent<-1, 0, -1, 0>, 5 > vt;
  typedef accessor<1, enumtype::in, extent<1,-1,1,-1>, 5 > vt;
  // coef_diff -->  [I, J, C, D]+[L]              
  //typedef accessor<2, enumtype::in, extent<>, 6 > coef_diff;
  typedef accessor<2, enumtype::in, extent<1,-1,1,-1>, 6 > coef_diff;
  // kh        --> [I, J, K]+[L]
  //typedef accessor<3, enumtype::in,    extent<0, 0, -1, 0>, 5 > kh;
  typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 5 > kh;
    
  typedef boost::mpl::vector<dscl, vt, coef_diff, kh> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index c; // would be set to 5
        dimension<3>::Index d;
        dimension<4>::Index t;
        dimension<4>::Index d_coef_diff;
        
        eval(dscl{}) = eval(dscl{}) +
                       ( eval(coef_diff{})    * ( eval(vt{i-1,j-1})      + eval(vt{j-1,t+1})) 
                       + eval(coef_diff{d_coef_diff+1}) * ( eval(vt{i-1,j-1,d+1}) + eval(vt{j-1,d+1,t+1})) 
                       + eval(coef_diff{d_coef_diff+2}) * ( eval(vt{i-1,j-1,d+2}) + eval(vt{j-1,d+2,t+1})) 
                      ) * 0.5 * ( eval(kh{j-1})  +   eval(kh{}));
        
    }
}; // end consume_function_5

struct consume_function_6 {
    
  // dscl      --> [I, J, K]+[L]            
  //typedef accessor<0, enumtype::inout, extent<>, 5 > dscl;
  typedef accessor<0, enumtype::inout, extent<1,-1,1,-1>, 5 > dscl;  
  // vt        --> [I, J, D, T]           
  //typedef accessor<1, enumtype::in, extent<0, 0, -1, 0>, 5 > vt;
  typedef accessor<1, enumtype::in, extent<1,-1,1,-1>, 5 > vt;
  // coef_diff -->  [I, J, C, D]+[L]              
  //typedef accessor<2, enumtype::in, extent<>, 6 > coef_diff;
  typedef accessor<2, enumtype::in, extent<1,-1,1,-1>, 6 > coef_diff;
  // kh        --> [I, J, K]+[L]
  //typedef accessor<3, enumtype::in,    extent<0, 1, 0, 0>, 5 > kh;
  typedef accessor<3, enumtype::in,    extent<1,-1,1,-1>, 5 > kh;
    
  typedef boost::mpl::vector<dscl, vt, coef_diff, kh> arg_list;
  
  template <typename evaluation>
  GT_FUNCTION
  static void Do(evaluation const & eval, x_interval) {

        dimension<1>::Index i;
        dimension<2>::Index j;
        dimension<3>::Index k;
        dimension<3>::Index c; // would be set to 6
        dimension<3>::Index d;
        dimension<4>::Index t;
        dimension<4>::Index d_coef_diff;
        
        eval(dscl{}) = eval(dscl{}) +
                       ( eval(coef_diff{})    * ( eval(vt{j-1,t+1})      + eval(vt{})) 
                       + eval(coef_diff{d_coef_diff+1}) * ( eval(vt{j-1,d+1,t+1}) + eval(vt{})) 
                       + eval(coef_diff{d_coef_diff+2}) * ( eval(vt{j-1,d+2,t+1}) + eval(vt{})) 
                      ) * 0.5 * ( eval(kh{})  +   eval(kh{i+1}));
        
    }
}; // end consume_function_6

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
    typedef gridtools::layout_map<3,0,1>        layout_3d_K;
    typedef gridtools::layout_map<3,0,-1,1,2>   layout_4d;
    typedef gridtools::layout_map<4,3,-1,2,1,0> layout_5d;
#else    
    typedef gridtools::layout_map<0,1,2>        layout_3d_K;
    typedef gridtools::layout_map<2,3,-1,0,1>   layout_4d;
    typedef gridtools::layout_map<3,4,-1,0,1,2> layout_5d;
#endif

    // CUDA backend: use block strategy
#ifdef __CUDACC__
    typedef backend<Cuda,structured,Block> backend_t;
#else
    typedef backend<Host,structured,Block> backend_t;
#endif
    
    
    typedef backend_t::storage_info<0, layout_3d_K> storage_info_3d_K_t;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_K_t>::type storage_type_3d_K;

    typedef backend_t::storage_info<1, layout_4d> storage_info_4d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_4d_t>::type storage_type_4d;

    typedef backend_t::storage_info<2, layout_4d> storage_info_4d_coef_diff_t;
    typedef typename backend_t::storage_type<float_type, storage_info_4d_coef_diff_t>::type storage_type_4d_coef_diff;
     
    typedef backend_t::storage_info<3, layout_5d> storage_info_5d_t;
    typedef typename backend_t::storage_type<float_type, storage_info_5d_t>::type storage_type_5d;

    typedef backend_t::storage_info<4, layout_3d_K> storage_info_3d_K_kh_t;
    typedef typename backend_t::storage_type<float_type, storage_info_3d_K_kh_t>::type storage_type_3d_K_kh;

    // scl & dscl
    //storage_info_3d_K_t         metadata_3d_k(ADM_iall+1,ADM_jall+1,ADM_kall);
    storage_info_3d_K_t         metadata_3d_k(ADM_iall,ADM_jall,ADM_kall);
    // coef_diff
    storage_info_4d_coef_diff_t metadata_4d_coef_diff(ADM_iall,ADM_jall,6,ADM_nxyz);
    //storage_info_4d_coef_diff_t metadata_4d_coef_diff(ADM_iall,ADM_jall,1,6,ADM_nxyz);
    // vt  
    //storage_info_4d_t           metadata_4d(ADM_iall+1,ADM_jall+1,1,ADM_nxyz,TJ);
    storage_info_4d_t           metadata_4d(ADM_iall,ADM_jall,1,ADM_nxyz,TJ);
    // coef_intp
    storage_info_5d_t           metadata_5d(ADM_iall,ADM_jall,3,ADM_nxyz,TJ);
    //storage_info_5d_t           metadata_5d(ADM_iall,ADM_jall,1,3,ADM_nxyz,TJ);
    // kh
    //storage_info_3d_K_kh_t      metadata_3d_k_kh(ADM_iall+2,ADM_jall+2,ADM_kall);
    storage_info_3d_K_kh_t      metadata_3d_k_kh(ADM_iall,ADM_jall,ADM_kall);
    
    // flux functor
    storage_type_4d                                   vt(metadata_4d, 0.0, "vt");
    field <storage_type_3d_K,ADM_lall>::type          scl(metadata_3d_k, 0.0, "scl");
    field <storage_type_5d,ADM_lall>::type            coef_intp(metadata_5d, 0.0, "coef_intp");

    // consume functor
    field <storage_type_3d_K_kh,ADM_lall>::type       kh(metadata_3d_k_kh, 0.0, "kh");
    field <storage_type_3d_K,ADM_lall>::type          dscl(metadata_3d_k, 0.0, "dscl");
    field <storage_type_4d_coef_diff,ADM_lall>::type  coef_diff(metadata_4d_coef_diff, 0.0, "coef_diff");
    
    // Use following To fill "scl" up the values from dumpio.c
    
    /*
    // ########< read input data >######## //
    // set intial values using dumpio
  void *ORG_dscl      =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  //float_type *ORG_dscl_pl      =  malloc( (ADM_gall_pl*ADM_kall*ADM_lall_pl) * sizeof(ORG_dscl_pl[0]));
  void *ORG_scl       =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));
  //float_type *ORG_scl_pl       =  malloc( (ADM_gall_pl*ADM_kall*ADM_lall_pl) * sizeof(ORG_scl_pl[0]));
  void *ORG_kh        =  malloc( (ADM_iall*ADM_jall*ADM_kall*ADM_lall) * sizeof(float_type));        
  //float_type *ORG_kh_pl        =  malloc( (ADM_gall_pl*ADM_kall*ADM_lall_pl) * sizeof(ORG_kh_pl[0]));        
  void *ORG_coef_intp =  malloc( (ADM_iall*ADM_jall*3*ADM_nxyz*TJ*ADM_lall) * sizeof(float_type));        
  //float_type *ORG_coef_intp_pl =  malloc( (ADM_gall_pl*3*ADM_nxyz*ADM_lall_pl) * sizeof(ORG_coef_intp_pl[0]));        
  void *ORG_coef_diff =  malloc( (ADM_iall*ADM_jall*6*ADM_nxyz*ADM_lall) * sizeof(float_type));        
  //float_type *ORG_coef_diff_pl =  malloc( (ADM_vlink*ADM_nxyz*ADM_lall_pl) * sizeof(ORG_coef_diff_pl[0]));        

  int32_t EX_fid;  
  char    *EX_fname = (char*) malloc(1024 * sizeof(char));
  // IO_FREAD = fread mode
  int32_t IO_FREAD = 0;
  dumpio_syscheck();
  dumpio_mk_fname(EX_fname,"snapshot.dc_diffusion","pe",SET_prc_me-1,6);
  dumpio_fopen(EX_fname,EX_fid);
  // MUST READ ALL VALUES (INCLUDING _p ARRAYS)
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_dscl);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     ORG_dscl_pl);
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_scl);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     ORG_scl_pl);
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_kh);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     ORG_kh_pl);
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*3*ADM_nxyz*2*ADM_lall   , ORG_coef_intp);
  //dumpio_read_data( &EX_fid, ADM_gall_pl      *3*ADM_nxyz*  ADM_lall_pl, ORG_coef_intp_pl);
  dumpio_read_data( EX_fid, ADM_iall*ADM_jall*6*ADM_nxyz*  ADM_lall   , ORG_coef_diff);
  //dumpio_read_data( &EX_fid,           ADM_vlink*ADM_nxyz*  ADM_lall_pl, ORG_coef_diff_pl);
  //dumpio_fclose(EX_fid);
  // can a storage type be assigned a value from an array??
  // example: scl = ORG_scl

  for(int l=0; l<ADM_lall; ++l) // how to retrieve the l dimension??
  {
   for(int i=0; i<metadata_3d_k.template dims<0>(); ++i)
    for(int j=0; j<metadata_3d_k.template dims<1>(); ++j)
      for(int k=0; k<metadata_3d_k.template dims<2>(); ++k)
      {
        scl.get_value<0>(i,j,k)= *(float_type*)((ORG_scl + index4D(i,j,k,l)));
        dscl.get_value<0>(i,j,k)=*(float_type*)((ORG_dscl + index4D(i,j,k,l)));
      }
   for(int i=0; i<metadata_4d_coef_diff.template dims<0>(); ++i)
    for(int j=0; j<metadata_4d_coef_diff.template dims<1>(); ++j)
    {
      // note: for coef_diff dims<2> is killed (i.e. no K dimension)
      for(int c=0; c<metadata_4d_coef_diff.template dims<3>(); ++c)
        for(int d=0; d<metadata_4d_coef_diff.template dims<4>(); ++d)      
          coef_diff.get_value<0>(i,j,1,c,d)=* (float_type*)((ORG_coef_diff + index5D(i,j,c,d,l)));
      // note: for coef_intp dims<2> is killed (i.e. no K dimension)
      for(int a=0; a<metadata_5d.template dims<3>(); ++a)
         for(int d=0; d<metadata_5d.template dims<4>(); ++d)      
            for(int t=0; t<metadata_5d.template dims<5>(); ++t)      
              coef_intp.get_value<0>(i,j,1,a,d,t)=* (float_type*)((ORG_coef_intp + index6D(i,j,a,d,t,l)));
    }  
    for(int i=0; i<metadata_3d_k_kh.template dims<0>(); ++i)
       for(int j=0; j<metadata_3d_k_kh.template dims<1>(); ++j)
        for(int k=0; k<metadata_3d_k_kh.template dims<2>(); ++k)
          kh.get_value<0>(i,j,k)=* (float_type*)((ORG_kh + index4D(i,j,k,l)));
  }  // l loop

  free(ORG_dscl);     
  //free(ORG_dscl_pl);
  free(ORG_scl);
  //free(ORG_scl_pl);
  free(ORG_kh);
  //free(ORG_kh_pl);
  free(ORG_coef_intp);
  //free(ORG_coef_intp_pl);
  free(ORG_coef_diff);
  //free(ORG_coef_diff_pl);
  //###############################################################################
  */
    //scl.print();
  for(int l=0; l<ADM_lall; ++l) // how to retrieve the l dimension??
  {
   for(int i=0; i<metadata_3d_k.template dims<0>(); ++i)
    for(int j=0; j<metadata_3d_k.template dims<1>(); ++j)
      for(int k=0; k<metadata_3d_k.template dims<2>(); ++k)
      {
        scl.get_value<0>(i,j,k)= (i+j+k)/2.;
        dscl.get_value<0>(i,j,k)=i+j+k;
      }
   for(int i=0; i<metadata_4d_coef_diff.template dims<0>(); ++i)
    for(int j=0; j<metadata_4d_coef_diff.template dims<1>(); ++j)
    {
      // note: for coef_diff dims<2> is killed (i.e. no K dimension)
      for(int c=0; c<metadata_4d_coef_diff.template dims<3>(); ++c)
        for(int d=0; d<metadata_4d_coef_diff.template dims<4>(); ++d)      
          coef_diff.get_value<0>(i,j,1,c,d)=(i+j+c+d)*2.;
      // note: for coef_intp dims<2> is killed (i.e. no K dimension)
      for(int a=0; a<metadata_5d.template dims<3>(); ++a)
         for(int d=0; d<metadata_5d.template dims<4>(); ++d)      
            for(int t=0; t<metadata_5d.template dims<5>(); ++t)      
              coef_intp.get_value<0>(i,j,1,a,d,t)=(i+j+a+t)/2.;
    }  
    for(int i=0; i<metadata_3d_k_kh.template dims<0>(); ++i)
       for(int j=0; j<metadata_3d_k_kh.template dims<1>(); ++j)
        for(int k=0; k<metadata_3d_k_kh.template dims<2>(); ++k)
          kh.get_value<0>(i,j,k)=(i+j+k)*2.;
  }  // l loop

    // flux functor
    typedef arg<0, storage_type_4d>                                   p_vt;
    typedef arg<1, field <storage_type_3d_K,ADM_lall>::type>          p_scl;
    typedef arg<2, field <storage_type_5d,ADM_lall>::type>            p_coef_intp;
    
    // consume fucntor
    typedef arg<3, field <storage_type_3d_K,ADM_lall>::type>          p_dscl;
    typedef arg<4, field <storage_type_4d_coef_diff,ADM_lall>::type>  p_coef_diff;
    typedef arg<5, field <storage_type_3d_K_kh,ADM_lall>::type>       p_kh;
    
    // accessor list of all place holders used
    typedef boost::mpl::vector<p_vt, p_scl, p_coef_intp, p_dscl, p_coef_diff, p_kh> accessor_list;

    gridtools::domain_type<accessor_list> domain ((p_vt() = vt), (p_scl() = scl), (p_coef_intp() = coef_intp), (p_dscl() = dscl ), (p_coef_diff() = coef_diff), (p_kh() = kh));

    // i-1,j-1 ---- i+1, j+1
    // extra lower halo added to account for operation dependecy
    uint_t halo_size = 1; 
    uint_t di[5] = {halo_size*2, halo_size, halo_size*2,(uint_t) ADM_iall-halo_size-1,(uint_t) ADM_iall};
    uint_t dj[5] = {halo_size*2, halo_size, halo_size*2,(uint_t) ADM_jall-halo_size-1,(uint_t) ADM_jall};
    gridtools::grid<axis> grid(di,dj);
    

    // single splitter
    grid.value_list[0] = 0;
    grid.value_list[1] = ADM_kall-1;

    auto diffusion = make_computation<backend_t>
        (
         domain, grid,
         make_mss
         (
          execute<forward>(), // parallel means no vertical dependency, forward from 0 to K-1, otherwise it is backward
          make_esf<flux_function>(p_vt(), p_scl(), p_coef_intp() ),//define in same order as accessors unique identifier          
          make_esf<consume_function_1>(p_dscl(), p_vt(), p_coef_diff(), p_kh()),// in case of several DO functions due to dependinceis, add make_esf for each functor
          make_esf<consume_function_2>(p_dscl(), p_vt(), p_coef_diff(), p_kh()),// in case of several DO functions due to dependinceis, add make_esf for each functor
          make_esf<consume_function_3>(p_dscl(), p_vt(), p_coef_diff(), p_kh()),// in case of several DO functions due to dependinceis, add make_esf for each functor
          make_esf<consume_function_4>(p_dscl(), p_vt(), p_coef_diff(), p_kh()),// in case of several DO functions due to dependinceis, add make_esf for each functor
          make_esf<consume_function_5>(p_dscl(), p_vt(), p_coef_diff(), p_kh()),// in case of several DO functions due to dependinceis, add make_esf for each functor
          make_esf<consume_function_6>(p_dscl(), p_vt(), p_coef_diff(), p_kh())// in case of several DO functions due to dependinceis, add make_esf for each functor
          
          )
         );

    // allocate and init storages
    diffusion->ready();

    // copy to device
    diffusion->steady();

    start = clock();
    for(int iter=0; iter < SET_iteration; ++iter) {  
      diffusion->run(); 
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
    diffusion->finalize();

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
        std::cout << diffusion->print_meter() << std::endl;
#endif

    ASSERT_TRUE(result);
}
*/

    return 0;
}


/**
@}
*/
