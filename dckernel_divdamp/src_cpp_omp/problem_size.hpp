#ifndef PROBLEM_SIZE_H
#define PROBLEM_SIZE_H

//############dump I/O#############//
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstring>

#define H_SHORT  16
#define H_MID    64
#define H_LONG  256

/* action type */
#define H_FREAD   0
#define H_FWRITE  1
#define H_FAPPEND 2

/* return type */
#define ERROR_CODE   -1
#define SUCCESS_CODE  1

using namespace std;

typedef float (real32_t);
typedef double(real64_t);

/* file structure */
typedef struct{
  char fname[H_LONG];
  int32_t opened;
  FILE *fp;
} fileinfo_t;

extern fileinfo_t *dumpio_finfo;
//############dump I/O#############//


//-------------------------------------------------------------------------------
//
// Problem size and global parameters
//
//-------------------------------------------------------------------------------

  int ADM_NSYS     = 32;
  int ADM_MAXFNAME = 1024;
  int ADM_LOG_FID  = 6;

  int IO_FREAD = 0;

  //--- Identifier of triangle element (i-axis-side or j-axis side)
  int TI  = 1;
  int TJ  = 2;

  //--- Identifier of line element (i-axis-side, ij-axis side, or j-axis side)
  int AI  = 1;
  int AIJ = 2;
  int AJ  = 3;

  //--- Identifier of 1 variable
  int K0  = 1;

  int ADM_nxyz = 3; // dimension of the spacial vector

  //--- region
  static const int ADM_lall      = 1;     // number of regular region per process
  int ADM_lall_pl   = 2;     // number of pole    region per process

  //--- horizontal grid
  int ADM_gall      = 16900; // number of horizontal grid per regular region
  int ADM_gall_1d    = 130;   // number of horizontal grid (1D)
  int ADM_gmin       = 2;     // start index of 1D horizontal grid
  int ADM_gmax       = 129;   // end   index of 1D horizontal grid

  int ADM_iall      = 130;
  int ADM_imin      = 2;
  int ADM_imax      = 129;
  int ADM_jall      = 130;
  int ADM_jmin      = 2;
  int ADM_jmax      = 129;

  int ADM_gall_pl   = 6;     // number of horizontal grid for pole region
  int ADM_gslf_pl   = 1;     // index for pole point
  int ADM_gmin_pl   = 2;     // start index of grid around the pole point
  int ADM_gmax_pl   = 6;     // end   index of grid around the pole point
  int ADM_vlink     = 5;     // number of grid around the pole point

  //--- vertical grid
  int ADM_vlayer    = 40;    // number of vertical layer
  int ADM_kall      = 42;    // number of vertical grid
  int ADM_kmin      = 2 ;    // start index of vertical grid
  int ADM_kmax      = 41;    // end   index of vertical grid

  // NOTE: to run with a pole region and pentagon
  // set the following values to TRUE
  bool ADM_have_pl     = false; // this ID manages pole region?
  bool ADM_have_sgp[1] = {false}; // region have singlar point?

  //--- constant parameters
  float PI    = 3.14159265358979; // pi
  float EPS   = 1.E-16;           // small number
  float GRAV  = 9.80616;  // Gravitational accerlaration of the Earth [m/s2]
  float Rdry  =   287.0;   // Gas constant of air
  float CVdry =   717.5;   // Specific heat of air (consant volume)
  float NON_HYDRO_ALPHA = 1.0; // Nonhydrostatic/hydrostatic flag

  //--- mod_grd
  int XDIR = 1;
  int YDIR = 2;
  int ZDIR = 3;

  int GRD_LAT = 1;
  int GRD_LON = 2;

  float GRD_rscale = 6.37122E+6; // radius of the planet [m]

  std::string vgrid_fname ="./vgrid40_600m_24km.dat";
  std::string GRD_grid_type ="ON_SPHERE";

  //--- mod_gmtr
  std::string GMTR_polygon_type ="ON_SPHERE";

  int GMTR_p_nmax = 8;

  int P_AREA  = 1;
  int P_RAREA = 2;
  int P_IX    = 3;
  int P_IY    = 4;
  int P_IZ    = 5;
  int P_JX    = 6;
  int P_JY    = 7;
  int P_JZ    = 8;

  int GMTR_t_nmax = 5;

  int T_AREA  = 1;
  int T_RAREA = 2;
  int W1      = 3;
  int W2      = 4;
  int W3      = 5;

  int GMTR_a_nmax    = 12;
  int GMTR_a_nmax_pl = 18;

  int HNX  = 1;
  int HNY  = 2;
  int HNZ  = 3;
  int HTX  = 4;
  int HTY  = 5;
  int HTZ  = 6;
  int TNX  = 7;
  int TNY  = 8;
  int TNZ  = 9;
  int TTX  = 10;
  int TTY  = 11;
  int TTZ  = 12;

  int TN2X = 13;
  int TN2Y = 14;
  int TN2Z = 15;
  int TT2X = 16;
  int TT2Y = 17;
  int TT2Z = 18;

  int SET_iteration = 1;
  int SET_prc_me    = 1;
  double SET_dt     = 60.0;

//############dump I/O#############//
inline int index2D(int i, int j){ return i * ADM_jall + j;}    
inline int index3D(int i, int j, int k){ return i * ADM_kall * ADM_jall + j * ADM_lall + k;}    
inline int index4D(int i, int j, int k, int l){ return i * ADM_lall * ADM_kall * ADM_jall + j * ADM_lall * ADM_kall + k * ADM_lall + l;}    
inline int index5D(int i, int j, int k, int l, int m){ return i * ADM_jall * 6 * ADM_nxyz * ADM_lall + j * 6 * ADM_nxyz * ADM_lall + k * ADM_nxyz * ADM_lall + l * ADM_lall + m;}
inline int index6D(int i, int j, int k, int l, int m, int n){ return i * ADM_jall * 3 * ADM_nxyz * TJ * ADM_lall + j * 3 * ADM_nxyz * TJ * ADM_lall + k * ADM_nxyz * TJ * ADM_lall + l * TJ * ADM_lall + m * ADM_lall + n;}

/******************************************************************************/
/* C functions                                                                */
/******************************************************************************/

/** endian change *****************************************************/
extern void dumpio_ednchg( void* avp_pointer,
                           const int32_t ai_size,
                           const int32_t ai_num   );

/** filename generator ************************************************/
extern void dumpio_mk_fname( char *fname,
                             const char *base,
                             const char *ext,
                             int32_t i,
                             int32_t y );

/** string preprocess *************************************************/
extern void dumpio_set_str( char *_str,
                            const char *str,
                            const int str_len );

/** string postprocess *************************************************/
extern void dumpio_trim_str( char *_str,
                             const char *str,
                             const int str_len );

/** check system & initialze ******************************************/
extern int32_t dumpio_syscheck( void );

/** open file IO stream ***********************************************/
extern int32_t dumpio_fopen( const char *vname_in, int32_t mode );

/** close file IO stream **********************************************/
extern int32_t dumpio_fclose( int32_t fid );

/** write data array **************************************************/
extern int32_t dumpio_write_data( int32_t fid,
                                  int32_t idxsize,
                                  void *data );

/** read data array (full size) ***************************************/
extern int32_t dumpio_read_data( int32_t fid,
                                 int32_t idxsize,
                                 void *data );

/* file ID counter */
int32_t dumpio_num_of_file = 0;

/* package+data+status container */
fileinfo_t *dumpio_finfo = NULL;

/* system information */
int32_t dumpio_system_ednchg = 0;

/** endian change *****************************************************/
void dumpio_ednchg( void* avp_pointer,
                    const int32_t ai_size,
                    const int32_t ai_num   )
{
  int ai_csize, ai_cnum;
  char ac_buf[16];
  char *acp_tmp;
  char *acp_local;

  memset(ac_buf,'\0',sizeof(ac_buf));

  acp_tmp   = (char*) avp_pointer;
  acp_local = (char*) avp_pointer;
  for( ai_cnum=0; ai_cnum<ai_num; ai_cnum++ ) {
    memcpy(ac_buf, acp_local, ai_size);
    for( ai_csize=0; ai_csize<ai_size; ai_csize++ ) {
      *acp_local = ac_buf[ai_size-ai_csize-1];
      acp_local++;
    }
    acp_tmp += ai_size;
  }
}

/** filename generator ************************************************/
void dumpio_mk_fname( char *fname,
                      const char *base,
                      const char *ext,
                      int32_t i,
                      int32_t y )
{
  char _fname[H_LONG];

  switch (y) {
  case 4 :
    sprintf(_fname,"%s.%s%04d",base,ext,i);
    break;
  case 5 :
    sprintf(_fname,"%s.%s%05d",base,ext,i);
    break;
  case 6 :
    sprintf(_fname,"%s.%s%06d",base,ext,i);
    break;
  default :
    break;
  }

  dumpio_trim_str( fname, _fname, H_LONG-1 );
}

/** string preprocess *************************************************/
void dumpio_set_str( char *_str,
                     const char *str,
                     const int str_len )
{
  int i;

  strncpy(_str, str, str_len);

  _str[str_len] = '\0'; /* [fix] H.Yashiro 20120621 */
  for( i=str_len-1; i>=0; i-- ) {
    if( _str[i] == ' ' ) {
      _str[i] = '\0';
    } else {
      break;
    }
  }
}

/** string postprocess *************************************************/
void dumpio_trim_str( char *_str,
                      const char *str,
                      const int str_len )
{
  int i;

  strncpy(_str, str, str_len);

  _str[str_len] = ' ';
  for( i=str_len-1; i>=0; i-- ) {
    if( _str[i] == '\0' ) {
      _str[i] = ' ';
    } else {
      break;
    }
  }
}

/** check system & initialze ******************************************/
int32_t dumpio_syscheck( void )
{
  int32_t i=1;

  if ( (sizeof(real32_t)!=4) || (sizeof(real64_t)!=8) ) {
    printf("Data type (real) is inconsistent!\n");
    exit(1);
  }

  if ( *(char*)&i ) {
    dumpio_system_ednchg = 1;
  } else {
    dumpio_system_ednchg = 0;
  }

  return(SUCCESS_CODE);
}

/** open file IO stream ***********************************************/
int32_t dumpio_fopen( const char *fname, int32_t mode )
{
  int32_t fid;

  /* get file ID */
  fid = dumpio_num_of_file++;

  /* memory re-allocation (expand by new dumpio_num_of_file) */
  dumpio_finfo=(fileinfo_t *)realloc(dumpio_finfo,sizeof(fileinfo_t)*(dumpio_num_of_file));

  /* intitialize */
  dumpio_set_str(dumpio_finfo[fid].fname,fname,H_LONG-1);
  dumpio_finfo[fid].opened = 0;
  dumpio_finfo[fid].fp     = NULL;

  if ( mode==H_FWRITE ) {
    if ( (dumpio_finfo[fid].fp=fopen(dumpio_finfo[fid].fname,"wb"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",dumpio_finfo[fid].fname);
      exit(1);
    }
  } else if( mode==H_FREAD ) { /* [mod] H.Yashiro 20110907 avoid overwrite action */
    if ( (dumpio_finfo[fid].fp=fopen(dumpio_finfo[fid].fname,"rb"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",dumpio_finfo[fid].fname);
      exit(1);
    }
  } else if( mode==H_FAPPEND ) { /* [add] H.Yashiro 20110907 overwrite mode */
    if ( (dumpio_finfo[fid].fp=fopen(dumpio_finfo[fid].fname,"r+b"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",dumpio_finfo[fid].fname);
      exit(1);
    }
  }
  dumpio_finfo[fid].opened = 1;

  return(fid);
}

/** close file IO stream **********************************************/
int32_t dumpio_fclose( int32_t fid )
{
  fclose(dumpio_finfo[fid].fp);
  dumpio_finfo[fid].opened = 0;

  return(SUCCESS_CODE);
}

/** write data array **************************************************/
int32_t dumpio_write_data( int32_t fid,
                           int32_t idxsize,
                           void *data )
{
  int64_t datasize;
  void *_data;

  datasize = idxsize * 8;

  /* data */
  if(dumpio_system_ednchg){
    _data = malloc(datasize);
    memcpy( _data, data, datasize);
    dumpio_ednchg(_data,8,idxsize);
  }else{
    _data = data;
  }

  fwrite(_data,datasize,1,dumpio_finfo[fid].fp);

  if(dumpio_system_ednchg) { free(_data); }

  return(SUCCESS_CODE);
}

/** read data array (full size) ***************************************/
int32_t dumpio_read_data( int32_t fid,
                          int32_t idxsize,
                          void *data )
{
  int64_t datasize;

  datasize = idxsize * 8;

  fread(data,datasize,1,dumpio_finfo[fid].fp);
  if(dumpio_system_ednchg){
    dumpio_ednchg(data,8,idxsize);
  }

  return(SUCCESS_CODE);
}

/******************************************************************************/
/* Fortran interface                                                          */
/******************************************************************************/

extern void dumpio_mk_fname_( char *fname,
                              const char *base,
                              const char *ext,
                              int32_t *i,
                              int32_t *y,
                              int32_t fname_len,
                              int32_t base_len,
                              int32_t ext_len    );

extern void dumpio_syscheck_( void );

extern void dumpio_fopen_( int32_t *fid, const char *fname, int32_t *mode, int32_t fname_len );

extern void dumpio_fclose_( int32_t *fid );

extern void dumpio_write_data_( int32_t *fid, int32_t *idxsize, void *data );

extern void dumpio_read_data_( int32_t *fid, int32_t *idxsize, void *data );

/** filename generator ************************************************/
void dumpio_mk_fname_( char *fname,
                    const char *base,
                    const char *ext,
                    int32_t *i,
                    int32_t *y,
                    int32_t fname_len,
                    int32_t base_len,
                    int32_t ext_len    )
{
  char _base[H_LONG];
  char _ext[H_SHORT];
  int32_t _base_len = base_len;
  int32_t _ext_len  = ext_len;
  if( _base_len > H_LONG-1  ){ _base_len = H_LONG-1;  }
  if( _ext_len  > H_SHORT-1 ){ _ext_len  = H_SHORT-1; }

  dumpio_set_str( _base, base, _base_len );
  dumpio_set_str( _ext,  ext,  _ext_len  );

  dumpio_mk_fname( fname, _base, _ext, *i, *y );
}

/** check system & initialze ******************************************/
void dumpio_syscheck_( void )
{
  int32_t ierr;

  ierr = dumpio_syscheck();
}

/** register new file *************************************************/
void dumpio_fopen_( int32_t *fid, const char *fname, int32_t *mode, int32_t fname_len )
{
  int32_t ierr;

  char _fname[H_LONG];
  int32_t _fname_len = fname_len;
  if( _fname_len > H_LONG-1 ){ _fname_len = H_LONG-1; }

  dumpio_set_str( _fname, fname, _fname_len );

  *fid = dumpio_fopen( _fname, *mode );
}

/** close file IO stream **********************************************/
void  dumpio_fclose_( int32_t *fid )
{
  int32_t ierr;

  ierr = dumpio_fclose( *fid );
}

/** write data array **************************************************/
void dumpio_write_data_( int32_t *fid, int32_t *idxsize, void *data )
{
  int32_t ierr;

  ierr = dumpio_write_data( *fid, *idxsize, data );
}

/** read data array (full size) ***************************************/
void dumpio_read_data_( int32_t *fid, int32_t *idxsize, void *data )
{
  int32_t ierr;

  ierr = dumpio_read_data( *fid, *idxsize, data );
}
//############dump I/O#############//

#endif


