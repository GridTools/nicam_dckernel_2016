!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (divdamp operator)
!
!-------------------------------------------------------------------------------
program dckernel_divdamp
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_oprt3d, only: &
     OPRT3D_divdamp_putcoef, &
     OPRT3D_divdamp
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_ddivdx           (:,:,:)
  real(DP), allocatable :: ORG_ddivdx_pl        (:,:,:)
  real(DP), allocatable :: ORG_ddivdy           (:,:,:)
  real(DP), allocatable :: ORG_ddivdy_pl        (:,:,:)
  real(DP), allocatable :: ORG_ddivdz           (:,:,:)
  real(DP), allocatable :: ORG_ddivdz_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogvx           (:,:,:)
  real(DP), allocatable :: ORG_rhogvx_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogvy           (:,:,:)
  real(DP), allocatable :: ORG_rhogvy_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogvz           (:,:,:)
  real(DP), allocatable :: ORG_rhogvz_pl        (:,:,:)
  real(DP), allocatable :: ORG_rhogw            (:,:,:)
  real(DP), allocatable :: ORG_rhogw_pl         (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAM        (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAM_pl     (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAMH       (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGAMH_pl    (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGSQRTH     (:,:,:)
  real(DP), allocatable :: ORG_VMTR_RGSQRTH_pl  (:,:,:)
  real(DP), allocatable :: ORG_VMTR_C2WfactGz   (:,:,:,:)
  real(DP), allocatable :: ORG_VMTR_C2WfactGz_pl(:,:,:,:)
  real(DP), allocatable :: ORG_cinterp_TN       (:,:,:,:)
  real(DP), allocatable :: ORG_cinterp_TN_pl    (:,:,:,:)
  real(DP), allocatable :: ORG_cinterp_HN       (:,:,:,:)
  real(DP), allocatable :: ORG_cinterp_HN_pl    (:,:  ,:)
  real(DP), allocatable :: ORG_cinterp_TRA      (:,:,:)
  real(DP), allocatable :: ORG_cinterp_TRA_pl   (:,:)
  real(DP), allocatable :: ORG_cinterp_PRA      (:,:)
  real(DP), allocatable :: ORG_cinterp_PRA_pl   (:,:)

  real(RP), allocatable :: ddivdx           (:,:,:)
  real(RP), allocatable :: ddivdx_pl        (:,:,:)
  real(RP), allocatable :: ddivdy           (:,:,:)
  real(RP), allocatable :: ddivdy_pl        (:,:,:)
  real(RP), allocatable :: ddivdz           (:,:,:)
  real(RP), allocatable :: ddivdz_pl        (:,:,:)
  real(RP), allocatable :: rhogvx           (:,:,:)
  real(RP), allocatable :: rhogvx_pl        (:,:,:)
  real(RP), allocatable :: rhogvy           (:,:,:)
  real(RP), allocatable :: rhogvy_pl        (:,:,:)
  real(RP), allocatable :: rhogvz           (:,:,:)
  real(RP), allocatable :: rhogvz_pl        (:,:,:)
  real(RP), allocatable :: rhogw            (:,:,:)
  real(RP), allocatable :: rhogw_pl         (:,:,:)
  real(RP), allocatable :: VMTR_RGAM        (:,:,:)
  real(RP), allocatable :: VMTR_RGAM_pl     (:,:,:)
  real(RP), allocatable :: VMTR_RGAMH       (:,:,:)
  real(RP), allocatable :: VMTR_RGAMH_pl    (:,:,:)
  real(RP), allocatable :: VMTR_RGSQRTH     (:,:,:)
  real(RP), allocatable :: VMTR_RGSQRTH_pl  (:,:,:)
  real(RP), allocatable :: VMTR_C2WfactGz   (:,:,:,:)
  real(RP), allocatable :: VMTR_C2WfactGz_pl(:,:,:,:)
  real(RP), allocatable :: cinterp_TN       (:,:,:,:)
  real(RP), allocatable :: cinterp_TN_pl    (:,:,:,:)
  real(RP), allocatable :: cinterp_HN       (:,:,:,:)
  real(RP), allocatable :: cinterp_HN_pl    (:,:  ,:)
  real(RP), allocatable :: cinterp_TRA      (:,:,:)
  real(RP), allocatable :: cinterp_TRA_pl   (:,:)
  real(RP), allocatable :: cinterp_PRA      (:,:)
  real(RP), allocatable :: cinterp_PRA_pl   (:,:)

  real(RP), allocatable :: check_ddivdx     (:,:,:)
  real(RP), allocatable :: check_ddivdx_pl  (:,:,:)
  real(RP), allocatable :: check_ddivdy     (:,:,:)
  real(RP), allocatable :: check_ddivdy_pl  (:,:,:)
  real(RP), allocatable :: check_ddivdz     (:,:,:)
  real(RP), allocatable :: check_ddivdz_pl  (:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dckernel_divdamp"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( ddivdx           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ddivdx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ddivdy           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ddivdy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ddivdz           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ddivdz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( rhogvx           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( rhogvx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( rhogvy           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( rhogvy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( rhogvz           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( rhogvz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( rhogw            (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( rhogw_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( VMTR_RGAM        (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( VMTR_RGAMH       (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( VMTR_RGSQRTH     (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   ) )
  allocate( VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl) )
  allocate( cinterp_TN       (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( cinterp_TN_pl    (ADM_gall_pl,ADM_lall_pl,2    ,ADM_nxyz) )
  allocate( cinterp_HN       (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( cinterp_HN_pl    (ADM_gall_pl,ADM_lall_pl,      ADM_nxyz) )
  allocate( cinterp_TRA      (ADM_gall   ,ADM_lall   ,TI:TJ         ) )
  allocate( cinterp_TRA_pl   (ADM_gall_pl,ADM_lall_pl               ) )
  allocate( cinterp_PRA      (ADM_gall   ,ADM_lall                  ) )
  allocate( cinterp_PRA_pl   (ADM_gall_pl,ADM_lall_pl               ) )
  ! Todo : first touch with considering NUMA
  ddivdx           (:,:,:)   = 0.0_RP
  ddivdx_pl        (:,:,:)   = 0.0_RP
  ddivdy           (:,:,:)   = 0.0_RP
  ddivdy_pl        (:,:,:)   = 0.0_RP
  ddivdz           (:,:,:)   = 0.0_RP
  ddivdz_pl        (:,:,:)   = 0.0_RP
  rhogvx           (:,:,:)   = 0.0_RP
  rhogvx_pl        (:,:,:)   = 0.0_RP
  rhogvy           (:,:,:)   = 0.0_RP
  rhogvy_pl        (:,:,:)   = 0.0_RP
  rhogvz           (:,:,:)   = 0.0_RP
  rhogvz_pl        (:,:,:)   = 0.0_RP
  rhogw            (:,:,:)   = 0.0_RP
  rhogw_pl         (:,:,:)   = 0.0_RP
  VMTR_RGAM        (:,:,:)   = 0.0_RP
  VMTR_RGAM_pl     (:,:,:)   = 0.0_RP
  VMTR_RGAMH       (:,:,:)   = 0.0_RP
  VMTR_RGAMH_pl    (:,:,:)   = 0.0_RP
  VMTR_RGSQRTH     (:,:,:)   = 0.0_RP
  VMTR_RGSQRTH_pl  (:,:,:)   = 0.0_RP
  VMTR_C2WfactGz   (:,:,:,:) = 0.0_RP
  VMTR_C2WfactGz_pl(:,:,:,:) = 0.0_RP
  cinterp_TN       (:,:,:,:) = 0.0_RP
  cinterp_TN_pl    (:,:,:,:) = 0.0_RP
  cinterp_HN       (:,:,:,:) = 0.0_RP
  cinterp_HN_pl    (:,:  ,:) = 0.0_RP
  cinterp_TRA      (:,:,:  ) = 0.0_RP
  cinterp_TRA_pl   (:,:    ) = 0.0_RP
  cinterp_PRA      (:,:    ) = 0.0_RP
  cinterp_PRA_pl   (:,:    ) = 0.0_RP

  allocate( check_ddivdx     (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_ddivdx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( check_ddivdy     (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_ddivdy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( check_ddivdz     (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_ddivdz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  check_ddivdx     (:,:,:)   = 0.0_RP
  check_ddivdx_pl  (:,:,:)   = 0.0_RP
  check_ddivdy     (:,:,:)   = 0.0_RP
  check_ddivdy_pl  (:,:,:)   = 0.0_RP
  check_ddivdz     (:,:,:)   = 0.0_RP
  check_ddivdz_pl  (:,:,:)   = 0.0_RP

  !###############################################################################

  !---< read input data >---
  allocate( ORG_ddivdx           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_ddivdx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_ddivdy           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_ddivdy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_ddivdz           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_ddivdz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_rhogvx           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_rhogvx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_rhogvy           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_rhogvy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_rhogvz           (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_rhogvz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_rhogw            (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_rhogw_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_VMTR_RGAM        (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_VMTR_RGAMH       (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_VMTR_RGSQRTH     (ADM_gall   ,ADM_kall,ADM_lall   )   )
  allocate( ORG_VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)   )
  allocate( ORG_VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   ) )
  allocate( ORG_VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl) )
  allocate( ORG_cinterp_TN       (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( ORG_cinterp_TN_pl    (ADM_gall_pl,ADM_lall_pl,2    ,ADM_nxyz) )
  allocate( ORG_cinterp_HN       (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( ORG_cinterp_HN_pl    (ADM_gall_pl,ADM_lall_pl,      ADM_nxyz) )
  allocate( ORG_cinterp_TRA      (ADM_gall   ,ADM_lall   ,TI:TJ         ) )
  allocate( ORG_cinterp_TRA_pl   (ADM_gall_pl,ADM_lall_pl               ) )
  allocate( ORG_cinterp_PRA      (ADM_gall   ,ADM_lall                  ) )
  allocate( ORG_cinterp_PRA_pl   (ADM_gall_pl,ADM_lall_pl               ) )

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.dc_divdamp3d','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_ddivdx         (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_ddivdx_pl      (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_ddivdy         (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_ddivdy_pl      (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_ddivdz         (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_ddivdz_pl      (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_rhogvx         (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_rhogvx_pl      (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_rhogvy         (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_rhogvy_pl      (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_rhogvz         (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_rhogvz_pl      (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_rhogw          (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_rhogw_pl       (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_VMTR_RGAM      (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_VMTR_RGAM_pl   (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_VMTR_RGAMH     (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_VMTR_RGAMH_pl  (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*ADM_lall   , ORG_VMTR_RGSQRTH   (:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*ADM_lall_pl, ORG_VMTR_RGSQRTH_pl(:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_kall*6*ADM_lall   , ORG_VMTR_C2WfactGz   (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_kall*6*ADM_lall_pl, ORG_VMTR_C2WfactGz_pl(:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_lall   *3*3, ORG_cinterp_TN    (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_lall_pl*2*3, ORG_cinterp_TN_pl (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_lall   *3*3, ORG_cinterp_HN    (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_lall_pl*3  , ORG_cinterp_HN_pl (:,:,  :) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_lall   *2  , ORG_cinterp_TRA   (:,:,:  ) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_lall_pl    , ORG_cinterp_TRA_pl(:,:    ) )
  call dumpio_read_data( EX_fid, ADM_gall   *ADM_lall       , ORG_cinterp_PRA   (:,:    ) )
  call dumpio_read_data( EX_fid, ADM_gall_pl*ADM_lall_pl    , ORG_cinterp_PRA_pl(:,:    ) )

  call dumpio_fclose(EX_fid)

  ddivdx           (:,:,:)   = real( ORG_ddivdx           (:,:,:)  , kind=RP )
  ddivdx_pl        (:,:,:)   = real( ORG_ddivdx_pl        (:,:,:)  , kind=RP )
  ddivdy           (:,:,:)   = real( ORG_ddivdy           (:,:,:)  , kind=RP )
  ddivdy_pl        (:,:,:)   = real( ORG_ddivdy_pl        (:,:,:)  , kind=RP )
  ddivdz           (:,:,:)   = real( ORG_ddivdz           (:,:,:)  , kind=RP )
  ddivdz_pl        (:,:,:)   = real( ORG_ddivdz_pl        (:,:,:)  , kind=RP )
  rhogvx           (:,:,:)   = real( ORG_rhogvx           (:,:,:)  , kind=RP )
  rhogvx_pl        (:,:,:)   = real( ORG_rhogvx_pl        (:,:,:)  , kind=RP )
  rhogvy           (:,:,:)   = real( ORG_rhogvy           (:,:,:)  , kind=RP )
  rhogvy_pl        (:,:,:)   = real( ORG_rhogvy_pl        (:,:,:)  , kind=RP )
  rhogvz           (:,:,:)   = real( ORG_rhogvz           (:,:,:)  , kind=RP )
  rhogvz_pl        (:,:,:)   = real( ORG_rhogvz_pl        (:,:,:)  , kind=RP )
  rhogw            (:,:,:)   = real( ORG_rhogw            (:,:,:)  , kind=RP )
  rhogw_pl         (:,:,:)   = real( ORG_rhogw_pl         (:,:,:)  , kind=RP )
  VMTR_RGAM        (:,:,:)   = real( ORG_VMTR_RGAM        (:,:,:)  , kind=RP )
  VMTR_RGAM_pl     (:,:,:)   = real( ORG_VMTR_RGAM_pl     (:,:,:)  , kind=RP )
  VMTR_RGAMH       (:,:,:)   = real( ORG_VMTR_RGAMH       (:,:,:)  , kind=RP )
  VMTR_RGAMH_pl    (:,:,:)   = real( ORG_VMTR_RGAMH_pl    (:,:,:)  , kind=RP )
  VMTR_RGSQRTH     (:,:,:)   = real( ORG_VMTR_RGSQRTH     (:,:,:)  , kind=RP )
  VMTR_RGSQRTH_pl  (:,:,:)   = real( ORG_VMTR_RGSQRTH_pl  (:,:,:)  , kind=RP )
  VMTR_C2WfactGz   (:,:,:,:) = real( ORG_VMTR_C2WfactGz   (:,:,:,:), kind=RP )
  VMTR_C2WfactGz_pl(:,:,:,:) = real( ORG_VMTR_C2WfactGz_pl(:,:,:,:), kind=RP )
  cinterp_TN       (:,:,:,:) = real( ORG_cinterp_TN       (:,:,:,:), kind=RP )
  cinterp_TN_pl    (:,:,:,:) = real( ORG_cinterp_TN_pl    (:,:,:,:), kind=RP )
  cinterp_HN       (:,:,:,:) = real( ORG_cinterp_HN       (:,:,:,:), kind=RP )
  cinterp_HN_pl    (:,:  ,:) = real( ORG_cinterp_HN_pl    (:,:  ,:), kind=RP )
  cinterp_TRA      (:,:,:  ) = real( ORG_cinterp_TRA      (:,:,:  ), kind=RP )
  cinterp_TRA_pl   (:,:    ) = real( ORG_cinterp_TRA_pl   (:,:    ), kind=RP )
  cinterp_PRA      (:,:    ) = real( ORG_cinterp_PRA      (:,:    ), kind=RP )
  cinterp_PRA_pl   (:,:    ) = real( ORG_cinterp_PRA_pl   (:,:    ), kind=RP )

  deallocate( ORG_ddivdx            )
  deallocate( ORG_ddivdx_pl         )
  deallocate( ORG_ddivdy            )
  deallocate( ORG_ddivdy_pl         )
  deallocate( ORG_ddivdz            )
  deallocate( ORG_ddivdz_pl         )
  deallocate( ORG_rhogvx            )
  deallocate( ORG_rhogvx_pl         )
  deallocate( ORG_rhogvy            )
  deallocate( ORG_rhogvy_pl         )
  deallocate( ORG_rhogvz            )
  deallocate( ORG_rhogvz_pl         )
  deallocate( ORG_rhogw             )
  deallocate( ORG_rhogw_pl          )
  deallocate( ORG_VMTR_RGAM         )
  deallocate( ORG_VMTR_RGAM_pl      )
  deallocate( ORG_VMTR_RGAMH        )
  deallocate( ORG_VMTR_RGAMH_pl     )
  deallocate( ORG_VMTR_RGSQRTH      )
  deallocate( ORG_VMTR_RGSQRTH_pl   )
  deallocate( ORG_VMTR_C2WfactGz    )
  deallocate( ORG_VMTR_C2WfactGz_pl )
  deallocate( ORG_cinterp_TN        )
  deallocate( ORG_cinterp_TN_pl     )
  deallocate( ORG_cinterp_HN        )
  deallocate( ORG_cinterp_HN_pl     )
  deallocate( ORG_cinterp_TRA       )
  deallocate( ORG_cinterp_TRA_pl    )
  deallocate( ORG_cinterp_PRA       )
  deallocate( ORG_cinterp_PRA_pl    )

  !###############################################################################

  call GRD_setup ! allocate GRD_rdgz

  !---< operator module setup >---
  call OPRT3D_divdamp_putcoef( VMTR_RGAM     (:,:,:)  , VMTR_RGAM_pl     (:,:,:)  , & ![IN]
                               VMTR_RGAMH    (:,:,:)  , VMTR_RGAMH_pl    (:,:,:)  , & ![IN]
                               VMTR_RGSQRTH  (:,:,:)  , VMTR_RGSQRTH_pl  (:,:,:)  , & ![IN]
                               VMTR_C2WfactGz(:,:,:,:), VMTR_C2WfactGz_pl(:,:,:,:), & ![IN]
                               cinterp_TN    (:,:,:,:), cinterp_TN_pl    (:,:,:,:), & ![IN]
                               cinterp_HN    (:,:,:,:), cinterp_HN_pl    (:,:,:)  , & ![IN]
                               cinterp_TRA   (:,:,:)  , cinterp_TRA_pl   (:,:)    , & ![IN]
                               cinterp_PRA   (:,:)    , cinterp_PRA_pl   (:,:)      ) ![IN]

  check_ddivdx   (:,:,:) = ddivdx   (:,:,:)
  check_ddivdx_pl(:,:,:) = ddivdx_pl(:,:,:)
  check_ddivdy   (:,:,:) = ddivdy   (:,:,:)
  check_ddivdy_pl(:,:,:) = ddivdy_pl(:,:,:)
  check_ddivdz   (:,:,:) = ddivdz   (:,:,:)
  check_ddivdz_pl(:,:,:) = ddivdz_pl(:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"
  call DEBUG_rapstart('DC_divdamp_kernel')

  do iteration = 1, SET_iteration



     call OPRT3D_divdamp( ddivdx(:,:,:), ddivdx_pl(:,:,:), &
                          ddivdy(:,:,:), ddivdy_pl(:,:,:), &
                          ddivdz(:,:,:), ddivdz_pl(:,:,:), &
                          rhogvx(:,:,:), rhogvx_pl(:,:,:), &
                          rhogvy(:,:,:), rhogvy_pl(:,:,:), &
                          rhogvz(:,:,:), rhogvz_pl(:,:,:), &
                          rhogw (:,:,:), rhogw_pl (:,:,:)  )



     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'check_ddivdx   ', check_ddivdx   (:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdx_pl', check_ddivdx_pl(:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdy   ', check_ddivdy   (:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdy_pl', check_ddivdy_pl(:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdz   ', check_ddivdz   (:,:,:) )
     call DEBUG_valuecheck( 'check_ddivdz_pl', check_ddivdz_pl(:,:,:) )
     call DEBUG_valuecheck( 'rhogvx         ', rhogvx         (:,:,:) )
     call DEBUG_valuecheck( 'rhogvx_pl      ', rhogvx_pl      (:,:,:) )
     call DEBUG_valuecheck( 'rhogvy         ', rhogvy         (:,:,:) )
     call DEBUG_valuecheck( 'rhogvy_pl      ', rhogvy_pl      (:,:,:) )
     call DEBUG_valuecheck( 'rhogvz         ', rhogvz         (:,:,:) )
     call DEBUG_valuecheck( 'rhogvz_pl      ', rhogvz_pl      (:,:,:) )
     call DEBUG_valuecheck( 'rhogw          ', rhogw          (:,:,:) )
     call DEBUG_valuecheck( 'rhogw_pl       ', rhogw_pl       (:,:,:) )

     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'ddivdx         ', ddivdx         (:,:,:) )
     call DEBUG_valuecheck( 'ddivdx_pl      ', ddivdx_pl      (:,:,:) )
     call DEBUG_valuecheck( 'ddivdy         ', ddivdy         (:,:,:) )
     call DEBUG_valuecheck( 'ddivdy_pl      ', ddivdy_pl      (:,:,:) )
     call DEBUG_valuecheck( 'ddivdz         ', ddivdz         (:,:,:) )
     call DEBUG_valuecheck( 'ddivdz_pl      ', ddivdz_pl      (:,:,:) )
  enddo

  write(ADM_LOG_FID,*) '### Varidation : grid-by-grid diff ###'
  check_ddivdx   (:,:,:) = check_ddivdx   (:,:,:) - ddivdx   (:,:,:)
  check_ddivdx_pl(:,:,:) = check_ddivdx_pl(:,:,:) - ddivdx_pl(:,:,:)
  check_ddivdy   (:,:,:) = check_ddivdy   (:,:,:) - ddivdy   (:,:,:)
  check_ddivdy_pl(:,:,:) = check_ddivdy_pl(:,:,:) - ddivdy_pl(:,:,:)
  check_ddivdz   (:,:,:) = check_ddivdz   (:,:,:) - ddivdz   (:,:,:)
  check_ddivdz_pl(:,:,:) = check_ddivdz_pl(:,:,:) - ddivdz_pl(:,:,:)
  call DEBUG_valuecheck( 'check_ddivdx   ', check_ddivdx   (:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdx_pl', check_ddivdx_pl(:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdy   ', check_ddivdy   (:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdy_pl', check_ddivdy_pl(:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdz   ', check_ddivdz   (:,:,:) )
  call DEBUG_valuecheck( 'check_ddivdz_pl', check_ddivdz_pl(:,:,:) )

  call DEBUG_rapend('DC_divdamp_kernel')
  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dckernel_divdamp
!-------------------------------------------------------------------------------
