!**********************************************************************************************************************************
!  calchr.f90 for bcc-Fe,V,Cr,W,Ta,Nb,Mo,Ba  20012.8.27
!    2012. 1.30  F90 Version 1.00  written by M.Inukai
!  chrbe.f90
!    2012. 3.10  H. Sato version written by M. Inukai
!    2013. 6.11  H. Sato version written by M. Inukai
!    2013. 7.19  monoclinic crystral
!
! Compile option for gfortran
! gfortran -o makethhrinput -Wall -pedantic -std=f95 -fbounds-check -O
! -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
! -g makethhrinput.f90 -free -m64
!**********************************************************************************************************************************
!**********************************************************************************************************************************
! file data ----------------------------------------------------------
! file number :   1    2    3    4     5     6    7     8     9     10
! file name   :  f01      f03         key  dis  f07   f08   f09    f10
! file type   :  inp      inp         inp  out  out   out   out    out
! file data ----------------------------------------------------------
! file data ----------------------------------------------------------
! file number :  11   12   13   14    15    16    17    18    19    20
! file name   :           f13  f14   f15   f16   f17
! file type   :           out  out   inp   out   inp
! file data ----------------------------------------------------------
! file data ----------------------------------------------------------
! file number :  21   22   23   24    25    26    27    28    29    30
! file name   :                                  f27
! file type   :                            inp   inp
! file data ----------------------------------------------------------
! file number :  41   42   43   44    45    46    47    48    49    50
! file name   :                                  f47   f48   f49
! file type   :                                  out   out   out
! file data ----------------------------------------------------------
! file number :  51   52   53   54    55    56    57    58    59    60
! file name   :
! file type   :
! file data ----------------------------------------------------------
! file number :  61   62   63   64    65    66    67    68    69    70
! file name   :                      f65                           f70
! file type   :                      inp                           out
! file data ----------------------------------------------------------
! file number :  71   72   73   74    75    76    77    78    79    80
! file name   :                                                    f80
! file type   :                                                    out
! file data ----------------------------------------------------------
! file number :  81   82   83   84    85    86    87    88    89    90
! file name   :                                                    f90
! file type   :                                                    out
! file data ----------------------------------------------------------
! file number :  91   92   93   94    95    96    97    98    99
! file name   :
! file type   :
! file data ----------------------------------------------------------
!=======================================================================
program main
!**********************************************************************************************************************************
!module ff_cons
  implicit none
  ! 0: noting, 1: console output, 2: f90-98 output 3: more information
  ! CHOP is check option
  integer,parameter :: check_file_option  = 0
  !
  ! 0:no, 1:show bond and anti-bond
  integer,parameter :: bond_state_switch = 0
  !
  ! cut eigenvalute at threshold, NEED cut_Ry_threshold_cons < 0
  real(8),parameter :: cut_Ry_threshold_cons = -10.0
  !
  !
! H. Sato program line start
  ! CNRYEV = 13.6058
  real(8),parameter :: convert_Ry2eV_cons = 13.60569805
  real(8),parameter :: PI=3.141592653
! H. Sato program line end
  !
!end module ff_cons
!**********************************************************************************************************************************
! HR plot
  ! 0:off, 1: 3D plot
  integer,parameter :: hr_3D_plot_switch = 0
  !
!  real(8),parameter :: PI = 3.14159265
  !
  ! 0:off, 1: on
  integer,parameter :: default_setting_switch = 1
  !
  ! 0:off, 1: on
!  integer,parameter print_switch == 0
  integer print_switch
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!module ff_read_data
!  implicit none
  !
  integer energy_array_number ! ENNM = energy_array_number >= energy_line_end_num
  integer ff_array_number     ! FFNM = ff_array_number
  !
  character(11), allocatable :: energy_num_chara(:) ! energy(ENNM)
  !
  real(8), allocatable :: eigenvalue_array_Ry(:)                       ! eigenv(ENNM)
  !
  integer wave_num(3)                                                  ! wavenum=wavenu(3)
  real(8), allocatable :: wave_K2(:)                                   ! wave_K2(FFNM)
  real(8), allocatable :: wave_coef_real(:,:)                          ! wavecr(FFNM, ENNM)
  real(8), allocatable :: wave_coef_imag(:,:)                          ! waveci(FFNM, ENNM)
  real(8), allocatable :: wave_coef_2(:,:)                             ! wavcr2(FFNM, ENNM)
  real(8), allocatable :: wave_coef_2_bond(:,:)                        ! wavec2b(FFNM, ENNM)
  real(8), allocatable :: wave_coef_2_anti_bond(:,:)                   ! wavec2a(FFNM, ENNM)
  !
!  integer, allocatable :: degenerate_pw(:)                             ! degapw(FFNM)
  !
  real(8), allocatable :: shelter_array(:)                             ! shelt4(ENNM)
  real(8), allocatable :: shelter_array_bond(:)                        ! shelt4b(ENNM)
  real(8), allocatable :: shelter_array_anti_bond(:)                   ! shelt4a(ENNM)
  !
  real(8), allocatable :: eigenvalue_array_eV(:)                       ! eigeve(ENNM)
  real(8), allocatable :: eigenvalue_array_sub_ef_eV(:)                ! neigev(ENNM)
  !
!end module ff_read_data
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!module character_cons
!  implicit none
  character(45) warning_chara
    DATA warning_chara / 'Warning, find out ********* in PW coefficient' /
  !
  character(29) eigenvalue_start_line_chara
    DATA eigenvalue_start_line_chara / '     EIGENVALUES ARE:        ' /
  !
  character(29) eigenvalue_end_line_chara
    DATA eigenvalue_end_line_chara / '       **********************' /
  !
  character(29) ff_start_end_line_chara
    DATA ff_start_end_line_chara / '   RECIPROCAL LATTICE VECTORS' /
  !
  character(29) imaginary_part_chara
    DATA imaginary_part_chara / 'IMAGPART                     ' /
  !
  character(11) empty_chara
    DATA empty_chara / '           ' /
  !
  character(7)  read_k_pre_chara
    DATA read_k_pre_chara / '     K=' /
!end module character_cons
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!module character_set
!  implicit none
  character(130) reador
  character(40) pre_line_nlo_chara
  character(29) dum ! dump (read "arbitrary character array")
  character(10) read_k_point_chara, read_now_k_point_chara
! read EF automatically
  character(43) read_dump ! Sato Auto
! syesno is select yes or no
!  character(1) axis_select
!end module character_set
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!module integer_set
!  implicit none
  integer wavevector_line_set_up_twice
  ! 2:real part, 3:imaginary part
  integer read_data_line
  integer eigenvalue_start_line_num
  integer eigenvalue_end_line_num
  integer ff_coefficient_line_num
  integer ff_start_line_num
  integer ff_end_line_num
  integer energy_line_num
  integer now_num_eigenvalue_line
  integer true_eigenvalue_start_line_num
  integer true_eigenvalue_end_line_num
  integer ff_line_set_num
  integer now_search_ff_line_num
  integer now_search_ff_total_line_num
  integer end_line_num
  integer num_eigenvalue_line
  ! energy_line_end_num <= energy_array_number
  integer energy_line_end_num
! get fflend, when read start
  integer num_of_miller_index
  integer new_pw_data_num
  integer new_ff_line_end_num
  integer ff_count
  integer i, j, k
  integer int_energy_line_num_divided
!  integer NL1, NL5, NL6
  integer temp_eigen_line_calc_num
!  integer degeneracy_num
  integer real_imaginary_switch
  integer eigen_value_num
! c reduce zone wave h, k, l
  integer num_matrix
  integer matrix_num
  integer nlo_num
  integer read_nlo
!  integer num_of_k_point
!  integer hcp_on
  integer ff_start_line
  integer now_matrix_num
  integer cut_energy_line_num
  integer num_eigenvalue_array_Ry
  integer ef_read_type ! Sato Auto
!  integer ff_list_max_num ! Sato Auto
  integer II ! Sato Auto
  integer new_ff_line_end_num_min ! Sato Auto
  integer new_ff_line_end_num_temp ! Sato Auto
  integer low_num_eigenvale ! Sato Auto
!end module integer_set
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!module real_set
!  implicit none
  real(8) energy_line_num_divided_nine
  real(8) true_energy_line_num
  real(8) ef_energy_Ry
  real(8) miller_index_h, miller_index_k, miller_index_l
  real(8) shelter
  real(8) axis_a, axis_b, axis_c, axis_n
  real(8) angle_a, angle_b, angle_c
  real(8) angle_a_radian, angle_b_radian, angle_c_radian
!  real(8) energy_eV_min, energy_eV_max
  real(8) num_eigenvalue_sub_ef_Ry
  !
!end module real_set
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!program main
  !
!  use ff_cons, only : convert_Ry2eV_cons
!  use character_cons
!  use character_set
!  use integer_set
!  use real_set
  !
!  implicit none
  ! -----------------------
  !
!  integer k_point_num
!  integer energy_array_number ! ENNM = energy_array_number >= energy_line_end_num
!  integer ff_array_number     ! FFNM = ff_array_number
  !
!  character(11), allocatable :: energy_num_chara(:)  ! energy(ENNM)
  !
!  real(8), allocatable :: eigenvalue_array_Ry(:)      ! eigenv(ENNM)
  !
!  integer wave_num(3)                                  ! wavenum=wavenu(3)
!  real(8), allocatable :: wave_K2(:)                  ! wave_K2(FFNM)
!  real(8), allocatable :: wave_coef_real(:,:)         ! wavecr(FFNM, ENNM)
!  real(8), allocatable :: wave_coef_imag(:,:)         ! waveci(FFNM, ENNM)
!  real(8), allocatable :: wave_coef_2(:,:)            ! wavcr2(FFNM, ENNM)
!  real(8), allocatable :: wave_coef_2_bond(:,:)       ! wavec2b(FFNM, ENNM)
!  real(8), allocatable :: wave_coef_2_anti_bond(:,:)  ! wavec2a(FFNM, ENNM)
  !
!  integer, allocatable :: degenerate_pw(:)            ! degapw(FFNM)
  !
!  real(8), allocatable :: shelter_array(:)            ! shelt4(ENNM)
!  real(8), allocatable :: shelter_array_bond(:)       ! shelt4b(ENNM)
!  real(8), allocatable :: shelter_array_anti_bond(:)  ! shelt4a(ENNM)
  !
!  real(8), allocatable :: eigenvalue_array_eV(:)      ! eigeve(ENNM)
!  real(8), allocatable :: eigenvalue_array_sub_ef_eV(:)  ! neigev(ENNM)
  ! -----------------------
  !
!**********************************************************************************************************************************
!  call read_line_and_data(ff_array_number)
!
! HR calculation continue
!-------------- -------------------------------------------------------
!  character(16) name_of_column_chara(12) ! COMPNE
!  character(3) end_chara
  !
  integer total_k_point          ! TWEGIH
  integer now_num_of_k_point     ! NNUOKP
  integer num_of_energy_step     ! NOSTEP
!  integer dis_end_k_point_line
!  integer energy_line_end_num    ! enlend
  integer end_of_k_point_line    ! ENDOKL
  integer ix, iy, iz, idv
  integer HR_i, HR_j, HR_k, HR_l, HR_m
  integer num_total_K2_step_for_3D_plot
  !
!  real(8) N
  real(8) correct_volume_factor
  real(8) correct_volume_factor_a
  real(8) fac
  real(8) K2_at_C2_1st_max        ! maxwg2
  real(8) C2_1st_max              ! maxwc2c
  real(8) K2_at_C2_2nd_max        ! maxwg2_2nd
  real(8) C2_2nd_max              ! maxwc2c_2nd
  real(8) K2_at_C2_max_limited_K  ! maxwg2_Kl
  real(8) C2_max_limited_K        ! maxwc2c_Kl
  real(8) First_C2_max_limited_K  ! FPF_K_limit
!  real(8) start_energy            ! STROEN
!  real(8) energy_step             ! ENESTE
!  real(8) K_limit
!  real(8) C2_limit                ! X_limit
!  real(8) E_limit
!  real(8) K2_step_for_3D_plot     ! K2STE_Elimit
!  real(8) sigma                   ! SIGMA
!  real(8) sigma_max               ! SIGMA_MAX
!  real(8) sigma2_max              ! SIGMA2_MAX
!  real(8) temp_sigma2_max         ! RSIGMA2
!  real(8) sigma_percent           ! SIGMAP
!  real(8) sigma2_percent          ! SIGMA2P
!  real(8) corrected_total_K_max   ! KGT
!  real(8) corrected_total_K2_max  ! KGT2
!  real(8) corr_total_K_max_limi_K ! KGT_KL
!  real(8) corr_total_K2_max_limi_K! KGT2_KL
!  real(8) num_K2_step             ! K2STEP
!  real(8) energy_for_hr_plot      ! ENERHR
!  real(8) total_ave_C2_max_dw     ! CKT
  !
  character(10), allocatable :: num_of_k_point_chara(:)
  real(8), allocatable :: C2_max_at_k_point(:,:)    ! C2_Kmax(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: weight_at_C2_max(:,:)     ! CTWEGIA_Ckmax(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: total_K_max(:,:)          ! KGA(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: total_K_max_limited_K(:,:)! KGT_K_limit(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: count_weight(:,:)         ! CTWEGIA(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: total_C2_for_3D_plot(:)   ! K2_Elimit_T(num_total_K2_step_for_3D_plot)
  real(8), allocatable :: total_C2_for_3D_plot_dw(:)! K2_Elimit(num_total_K2_step_for_3D_plot)
  real(8), allocatable :: weight_for_3D_plot(:)     ! CTWEGIA_Elimit(num_total_K2_step_for_3D_plot)
  real(8), allocatable :: weight(:)                 ! WEIGHT(end_of_k_point_line)
!  real(8), allocatable :: K2_distribution_at_E(:)   ! DENERGY(dis_end_k_point_line)
!  real(8), allocatable :: total_averaged_K_max(:)   ! KGTA(num_of_energy_step)
!  real(8), allocatable :: true_weight_at_k_point(:) ! CTWETA(num_of_energy_step)
!  real(8), allocatable :: true_weight_at_C2_max(:)  ! CTWETA_Ckmax(num_of_energy_step)
!  real(8), allocatable :: tot_ave_K_max_limit(:)    ! KGT_K_limit_T(num_of_energy_step)
!  real(8), allocatable :: total_ave_C2_max(:)       ! C2_Kmax_T(num_of_energy_step)
!  real(8), allocatable :: sigma2(:)                 ! SIGMA2(num_of_energy_step)
!
! H. Sato program line start
  real(8), allocatable :: bk2(:)
  real(8), allocatable :: fk2(:)
!  real(8) AA
!  real(8) CC
  real(8) BXX
  real(8) BXY
  real(8) BXZ
  real(8) BYX
  real(8) BYY
  real(8) BYZ
  real(8) BZX
  real(8) BZY
  real(8) BZZ
  !
  real(8) BG1X
  real(8) BG1Y
  real(8) BG1Z
  real(8) BG2X
  real(8) BG2Y
  real(8) BG2Z
  real(8) BG3X
  real(8) BG3Y
  real(8) BG3Z
   !
  real(8) WVX
  real(8) WVY
  real(8) WVZ
  integer HR_l2
  integer k2
  real(8) EBOT
  real(8) AVV1
  real(8) AKK1
  real(8) fftt
  real(8) w2k
  real(8) w2f
  integer(4) kcount ! SO
  real(8) wave_K2_old ! SO
  integer(4) LO_search_flag ! SO
  integer(4) new_ff_array_number ! SO
!  real(8) AKN
! H. Sato program line end
  integer(4) IATO
  integer(4) NATO
  integer(4), allocatable ::  MIF(:)
  integer(4) MULT,ISPLIT
  integer(4) j_k, j_l, j_k_No
  real(8) end_l, start_k_maxC
!  integer(4) HR_init
!
!-------------- -------------------------------------------------------
!
!**********************************************************************************************************************************
!  call read_cons_data(k_point_chara_set)
!
! bag fix
! integer
ff_start_line = 0
now_num_of_k_point = 1
i = 0
j = 0
k = 0
HR_i = 0
HR_j = 0
HR_k = 0
HR_l = 0
HR_m = 0
ff_end_line_num = 0
ff_start_line_num = 0
II = 0 ! Sato Auto
new_ff_line_end_num_min = 0 ! Sato Auto
new_ff_line_end_num_temp = 0 ! Sato Auto
kcount = 0 ! SO
read_nlo = 0 ! SO
LO_search_flag = 0 ! SO
new_ff_array_number = 0 ! SO
! bag fix
! real
EBOT = 0.0
wave_K2_old = 0.0 ! SO
!
! H. Sato program line start
num_of_energy_step = 1
print_switch = 0
! H. Sato program line end
!
! H. Sato program line start
!      open( 3, file = 'f03' ) ! Sfe2.klist
      open( 4, file = 'input_K2.dat' )
      open(10, file = 'input_C2.dat' )
      ! nefene=0.49239
!      ef_energy_Ry = 0.57495 ;Auto read from f15
      !
!      AA = 5.03811976
!      CC = (PI/AA)**2 * convert_Ry2eV_cons
      !
! H. Sato program line end
!
! -----------------------
! f10 : FF EF input data from .dos1(f10)
!  open (10, file = 'f10')
  !
!  read(10, '(A43)') read_dump
!  read(10, '(4x,f10.5)') ef_energy_Ry
!  write(6, * ) 'Auto Read Fermi Energy: EF=', ef_energy_Ry, ' (Ry)'
  !
!  close (10)
!
! -----------------------
!  open (65, file = 'f65')
!    do i=1,14
!      read(65, '(A43)') read_dump
      !write(6,*) read_dump
!    end do
!    read(65, '(7x,i5,16x,i5)') IATO,NATO
    !write(6, '(7x,i5,16x,i5)') IATO,NATO
!  close(65)
  open (60, file = 'f60')
    IATO=0
    NATO=0
    read(60, '(A43)') read_dump
    i = 1
    read(60, '(27x,i3)') IATO
    write(6,'(a5,i3)') "IATO=",IATO
    i=i+1
    do while ( read_dump(11:39)/="NUMBER OF SYMMETRY OPERATIONS" )
      read(60, '(A43)') read_dump
      !write(6,*) read_dump
      i=i+1
    end do
    allocate(MIF(i+1))
    !
    do j=1,i
      MIF(i)=0
    end do
    !
    rewind(60)
    read(60, '(A43)') read_dump
    i = 1
    do while ( read_dump(11:39)/="NUMBER OF SYMMETRY OPERATIONS" )
      read(60, '(A43)') read_dump
      !write(6,*) read_dump
      i=i+1
      if( read_dump(11:15)=="MULT=" )then
        MIF(i)=1
      end if
    end do
    !
    rewind(60)
    read(60, '(A43)') read_dump
    i = 1
    j = 1
    do while ( read_dump(11:39)/="NUMBER OF SYMMETRY OPERATIONS" )
      read(60, '(A43)') read_dump
      !write(6,*) read_dump
      i=i+1
      if( MIF(i+1)==1 )then
        read(60,'(15x,i2,17x,i2)') MULT,ISPLIT
        write(6,'(a4,i3,a6,i2,a8,i2)') 'ATOM',j,',MULT=',MULT,',ISPLIT=',ISPLIT
        !NATO = NATO + abs(MULT*ISPLIT)
        NATO = NATO + abs(MULT)
        i=i+1
        j=j+1
      end if
    end do
    deallocate(MIF)
    !
    rewind(60)
    read(60, '( a130 )') reador
    !write(6, '( a45 )') reador(1:45)
    read(60, '( a130 )') reador
    !write(6, '( a45 )') reador(1:45)
    if (reador(1:1) == 'F') then
      NATO = NATO * 4
    else if (reador(1:1) == 'B') then
      NATO = NATO * 2
    else if (reador(1:1) == 'P') then
      NATO = NATO * 1
    else if (reador(1:1) == 'H') then
      NATO = NATO * 1
    else
      NATO = NATO * 1
    end if
    !
    rewind(60)
    write(6, '(a5,i5,a6,i5)') 'IATO=',IATO,',NATO=',NATO
!    write(6, '(a5,i5)') 'IATO=',IATO
  close(60)
  i=0
  j=0
  !
!  open(65, file = 'f65')
!    read(65, '(A43)') read_dump
!    i = 1
!    do while ( read_dump(1:16)/="number of atoms:" )
!      read(65, '(A43)') read_dump
      !write(6,*) read_dump
!      i=i+1
!    end do
    !
!    rewind(65)
!    read(65, '(A43)') read_dump
!    j = 1
!    do while ( j /= i )
!      read(65, '(A43)') read_dump
      !write(6,*) read_dump
!      j=j+1
!      if( j == (i-1) )then
!        read(65, '(16x,i4)') NATO
!        write(6,'(a5,i5)') 'NATO=', NATO
!      end if
!    end do
!  close(65)
!  i=0
!  j=0
  !
  open (70, file = 'f70')
    write(70,'(a5,i5)') 'IATO=', IATO
    write(70,'(a5,i5)') 'NATO=', NATO
  close(70)
! -----------------------
  open (15, file = 'f15')
  !
  i = 0
  do while( 1==1 )
    read(15, '(A43)') read_dump
    i = i + 1
    ! version 9 if ( read_dump(1:18) == '   FERMI ENERGY AT') exit
    if ( read_dump(1:18) == '   FERMI ENERGY AT') then
      ef_read_type = 0
      exit
    end if
    if ( read_dump(1:38) == ':FER  : F E R M I - ENERGY(FERMI-SM.)=') then
      ef_read_type = 1
      exit
    end if
  end do
  rewind(15)
  do j=1, i-1
    read(15, '(A43)') read_dump
  !  write(6,*) read_dump
  end do
  if ( ef_read_type == 0 ) then
    read(15, '(18x,f20.15)') ef_energy_Ry
  endif
  if ( ef_read_type == 1 ) then
    read(15, '(38x,f10.5)') ef_energy_Ry
  end if
  write(6, * ) 'Auto Read Fermi Energy: EF=', ef_energy_Ry, ' (Ry)'
  !
  close (15)
!
! -----------------------
! read input form .inputd
!  open (2, file = 'f02')
  !
!  read( 2, '(A3)' ) dum
!  read( 2, '(2f10.5, i10, i5)' ) start_energy, energy_step, num_of_energy_step, print_switch
!  if( print_switch == 0 ) then
!    write( 6, *) 'sigma 1 calculation version'
!  end if
!  if( print_switch == 1 ) then
!    write( 6, *) 'sigma 1 calculation version, and print all data'
!  end if
  !
!  close(2)
!
! -----------------------
! f60 : read structure data from .struct(f60)
  open(60, file = 'f60')
  rewind(60)
  !
  read(60, '( a130 )') reador
    write(6, '( a45 )') reador(1:45)
  read(60, '( a130 )') reador

  write(6, '( a45 )') reador(1:45)
!  hcp_on = 0
!  if (reador(1:1) == 'H') then
!    hcp_on = 1
!  endif
  read(60, '( a130 )') reador
  read(60, '( 6f10.6 )') axis_a, axis_b, axis_c, angle_a, angle_b, angle_c
  !
  write(6, * ) ' axis_a  ,  axis_b  ,  axis_c  , angle_a  , angle_b  , angle_c  '
  !
  !axis_a = axis_a * 0.529177249
  !axis_b = axis_b * 0.529177249
  !axis_c = axis_c * 0.529177249
  !
  write(6, '( 6f10.6 )') axis_a, axis_b, axis_c, angle_a, angle_b, angle_c
  !
  close (60)
!
! Normalized by axis
!  write(6, * ) 'Normalized by axis : a, b, c or n(=normalization by a axis)?'
!  if ( default_setting_switch == 1 ) then
!    write(6, * ) 'normalization by a axis'
!    axis_n = axis_a
!  else
!    read(5, '( a1 )') axis_select
!    if (axis_select == 'a') then
!      axis_n = axis_a
!      write(6, * ) 'Normalized by ', axis_select, ' axis'
!    elseif(axis_select == 'b') then
!      axis_n = axis_b
!      write(6, * ) 'Normalized by ', axis_select, ' axis'
!    elseif(axis_select == 'c') then
!      axis_n = axis_c
!      write(6, * ) 'Normalized by ', axis_select, ' axis'
!    else
!      write(6, * ) 'normalization by a axis'
!      axis_n = axis_a
!    hcp_on = -1
!    endif
!  end if
  !
! change radian
  write(6,*) 'This code can calculate FCC, BCC, HCP, SC, orthorombic, only'
  write(6,*) 'Range (0.001 < gamma < 180), alpha=90.0 and beta=90.0, only'
  angle_a_radian = angle_a*PI/180.0
  angle_b_radian = angle_b*PI/180.0
  angle_c_radian = angle_c*PI/180.0
!      write(6,*) sin(angle_a_radian), sin(angle_b_radian), sin(angle_c_radian)
!      write(6,*) cos(angle_a_radian), cos(angle_b_radian), cos(angle_c_radian)
! http://en.wikipedia.org/wiki/Bravais_lattice
! Triclinic (This formula is volume calculation for all structure)
  correct_volume_factor = ( ( axis_a * axis_b * axis_c * &
    ( 1.0 - (cos(angle_a_radian))**2.0 - (cos(angle_b_radian))**2.0 - (cos(angle_c_radian))**2.0 + &
    2.0*cos(angle_a_radian)*cos(angle_b_radian)*cos(angle_c_radian) )**(1.0/2.0) ))
  correct_volume_factor_a = axis_a * axis_a * axis_a
  fac = ( 1.0 - (cos(angle_a_radian))**2.0 - (cos(angle_b_radian))**2.0 - (cos(angle_c_radian))**2.0 + &
    2.0*cos(angle_a_radian)*cos(angle_b_radian)*cos(angle_c_radian) )**(1.0/2.0)
  open (14, file = 'f14')
    write(6,'(a8,f10.5)')  'Vabc/Va=', correct_volume_factor/correct_volume_factor_a
    write(6,'(a8,f10.5)')  'V/Vcub =', fac
    write(14,'(a8,f10.5)') 'Vabc/Va=', correct_volume_factor/correct_volume_factor_a
    write(14,'(a8,f10.5)') 'V/Vcub =', fac
  close(14)
! read EF input data from .struct : end
!
!
! check start --------------------------------------------------------
!  if (check_file_option >= 1) then
!    write(6, * ) 'read line start'
!  endif
! check end ----------------------------------------------------------
!
! f17 : read G1, G2, G3
  open (17, file = 'f17')
  i = 0
  do while(i==0)
    read (17, '(a26)') reador
    if (reador(1:26)  == '    G1        G2        G3') then
      !read(17, '(1x,3f10.6)') BXX, BYX, BZX
      !read(17, '(1x,3f10.6)') BXY, BYY, BZY
      !read(17, '(1x,3f10.6)') BXZ, BYZ, BZZ
      read(17, '(1x,3f10.6)') BG1X, BG2X, BG3X
      read(17, '(1x,3f10.6)') BG1Y, BG2Y, BG3Y
      read(17, '(1x,3f10.6)') BG1Z, BG2Z, BG3Z
      i = 1
    endif
  end do
  close (17)
  write(6,*)  'G1        G2        G3'
  !write(6,*) BXX, BYX, BZX
  !write(6,*) BXY, BYY, BZY
  !write(6,*) BXZ, BYZ, BZZ
  write(6,*) BG1X, BG2X, BG3X
  write(6,*) BG1Y, BG2Y, BG3Y
  write(6,*) BG1Z, BG2Z, BG3Z
  !
! -----------------------
! f03 : read k point data from .klist(f03)
  open (3, file = 'f03')
  !
  i = 0
  do while(1==1)
    read (3, '(A10)') read_k_point_chara
    if (read_k_point_chara == 'END       ') exit
    if (read_k_point_chara /= '          ') then
      i = i + 1
    endif
  end do
  !
  end_of_k_point_line = i
  allocate( num_of_k_point_chara(end_of_k_point_line) )
  allocate( weight(end_of_k_point_line) )
  !
  rewind(3)
  total_k_point = 0
  do i = 1, end_of_k_point_line
! format (A10. 4I10, 3F5.2) name, ix, iy, iz, idv, weight
    read( 3, '( A10, 4i10, f5.2)' ) num_of_k_point_chara(i), ix, iy, iz, idv, weight(i)
    total_k_point = total_k_point + int(weight(i))
  end do
  !
  write(6,*) 'total k point =', total_k_point
  !
  close (3)
  !
!
! -----------------------
! f29 : data from .inputd(f29)
!  open (29, file = 'f29')
!  i = 0
!  do while(1)
!    read( 29, '(A3)' ) end_chara
!    if( end_chara == 'END' ) exit
!    i = i + 1
!  end do
! set number of k point
!  dis_end_k_point_line = i
!  allocate ( K2_distribution_at_E(dis_end_k_point_line) )                   ! DENERGY(dis_end_k_point_line)
!  rewind(29)
!  do i = 1, dis_end_k_point_line
! format (F5.2) ENERGY for distribution
!    read( 29, '( f5.2 )' ) K2_distribution_at_E(i)
!     K2_distribution_at_E(i) = K2_distribution_at_E(i) / CNRYEV
!  end do
  !
!  close(29)
!
!  write( 6,*) '-----------------------'
! f27 : read G vs coefficients data from .gcene(f27)
!  open(27, file = 'f27')
!  rewind(27)
!  read(27, '( a130 )') reador
!  read(27, '( a130 )') reador
!  read(27, '(2f10.5)') energy_eV_min, energy_eV_max
!  write(6, * ) 'K vs coefficients tabel calculation range'
!  write(6, * ) 'energy_eV_min=', energy_eV_min
!  write(6, * ) 'energy_eV_max=', energy_eV_max
!  write(6, * )
!  close (27)
  !
!  write( 6,*) '-----------------------'
! read X input data : tart
!  write(6,*) 'Please, input X.X %. This code wast Ckmax1st,'
!  write(6,*) 'if (Ckmax1st - Ckmax2nd) / Ckmax1st * 100 < X.X %'
!  if ( default_setting_switch == 1 ) then
!    C2_limit = 0.01
!  else
!    read(5, '(f5.2)' ) C2_limit
!    if( C2_limit == 0.0 ) then
!      C2_limit = 5.0
!    end if
!  end if
!  write(6, *) 'SET= ', C2_limit,'%'
! read X input data : end
!
!  write( 6,*) '-----------------------'
! read K input data : start
!  write(6,*) 'Please, input |2K|^2 for |2*(k+G)|^2 < |2K|^2 (e.g. 2.01)'
!  if ( default_setting_switch == 1 ) then
!    K_limit = 0.0
!    write( 6,*) 'First |C(Kmax)|^2 Peak Plot'
!  else
!    read(5, '(f10.5)' ) K_limit
!    if( K_limit == 0.0 ) then
!      write( 6,*) 'First |C(Kmax)|^2 Peak Plot'
!    else
!      write( 6,*) 'SET = ', K_limit
!    end if
!  end if
! read K input data : end
!
!  write( 6,*) '-----------------------'
! read E input data : start
!  write(6,*) 'Please, input E. |2K|^2 vs |C|^2 plot at E eV (e.g. -1.0)'
!  if ( default_setting_switch == 1 ) then
!    E_limit = 0.0
!  else
!    read(5, '(f10.5)' ) E_limit
!  end if
!  write( 6,*) 'SET = ', E_limit, 'eV'
!      E_limit = E_limit/CNRYEV + nefene
! read E input data : end
  write( 6,*) '-----------------------'
!
!
!**********************************************************************************************************************************
!  call read_line_and_data(hr_array_number)
!
  allocate ( C2_max_at_k_point(end_of_k_point_line,num_of_energy_step) )    ! C2_Kmax(end_of_k_point_line, num_of_energy_step)
  allocate ( weight_at_C2_max(end_of_k_point_line,num_of_energy_step) )     ! CTWEGIA_Ckmax(end_of_k_point_line, num_of_energy_step)
  allocate ( total_K_max(end_of_k_point_line,num_of_energy_step) )          ! KGA(end_of_k_point_line, num_of_energy_step)
  allocate ( total_K_max_limited_K(end_of_k_point_line,num_of_energy_step) )! KGT_K_limit(end_of_k_point_line, num_of_energy_step)
  allocate ( count_weight(end_of_k_point_line,num_of_energy_step) )         ! CTWEGIA(end_of_k_point_line, num_of_energy_step)
  if ( (hr_3D_plot_switch == 1) .and. (print_switch == 1) ) then
    allocate ( total_C2_for_3D_plot(num_total_K2_step_for_3D_plot) )          ! K2_Elimit_T(num_total_K2_step_for_3D_plot)
    allocate ( total_C2_for_3D_plot_dw(num_total_K2_step_for_3D_plot) )       ! K2_Elimit(num_total_K2_step_for_3D_plot)
    allocate ( weight_for_3D_plot(num_total_K2_step_for_3D_plot) )            ! CTWEGIA_Elimit(num_total_K2_step_for_3D_plot)
  else
    allocate ( total_C2_for_3D_plot(1) )          ! K2_Elimit_T(num_total_K2_step_for_3D_plot)
    allocate ( total_C2_for_3D_plot_dw(1) )       ! K2_Elimit(num_total_K2_step_for_3D_plot)
    allocate ( weight_for_3D_plot(1) )            ! CTWEGIA_Elimit(num_total_K2_step_for_3D_plot)
  end if
!  allocate ( total_averaged_K_max(num_of_energy_step) )                     ! KGTA(num_of_energy_step)
!  allocate ( true_weight_at_k_point(num_of_energy_step) )                   ! CTWETA(num_of_energy_step)
!  allocate ( true_weight_at_C2_max(num_of_energy_step) )                    ! CTWETA_Ckmax(num_of_energy_step)
!  if ( print_switch == 1 ) then
!    allocate ( tot_ave_K_max_limit(num_of_energy_step) )                      ! KGT_K_limit_T(num_of_energy_step)
!  else
!    allocate ( tot_ave_K_max_limit(1) )                      ! KGT_K_limit_T(num_of_energy_step)
!  end if
!  allocate ( total_ave_C2_max(num_of_energy_step) )                         ! C2_Kmax_T(num_of_energy_step)
!  allocate ( sigma2(num_of_energy_step) )                                   ! SIGMA2(num_of_energy_step)
  ! initialize array
  do HR_i = 1,end_of_k_point_line
    do HR_j = 1,num_of_energy_step
      C2_max_at_k_point(HR_i,HR_j)     = 0.0
      weight_at_C2_max(HR_i,HR_j)      = 0.0
      total_K_max(HR_i,HR_j)           = 0.0
      total_K_max_limited_K(HR_i,HR_j) = 0.0
      count_weight(HR_i,HR_j)          = 0.0
    end do
  end do
  !
!  do HR_i = 1,num_of_energy_step
!    total_averaged_K_max(HR_i)   = 0.0
!    true_weight_at_k_point(HR_i) = 0.0
!    true_weight_at_C2_max(HR_i)  = 0.0
!    if ( print_switch == 1 ) then
!      tot_ave_K_max_limit(HR_i)    = 0.0
!    else
!      tot_ave_K_max_limit(1)    = 0.0
!    end if
!    total_ave_C2_max(HR_i)       = 0.0
!    sigma2(HR_i)                 = 0.0
!  end do
  !
  !
  do HR_i = 1, end_of_k_point_line
    now_num_of_k_point = HR_i
!**********************************************************************************************************************************
!  call read_line_and_data(ff_array_number)
!
! f01 : read G vs coefficients data from cutting file(cut 26)
!**********************************************************************************************************************************
!  call correct_data()
!  integer ff_start_line, now_matrix_num
  !
! f01 : read G vs coefficients data from cutting file(cut 26)
  open (1, file = 'f01')
  rewind (1)
! f26 : read G vs coefficients data from .output1(f26)
!
! H. Sato program line start
  open (26, file = 'f26')
!  open(26, file = 'Sfe.output1' )
! H. Sato program line end
!
  rewind (26)
  !
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write (6, * ) 'read line start'
  endif
! check end ----------------------------------------------------------
  !
!  write( 6, * ) 'Now number of k point =', now_num_of_k_point, '/', end_of_k_point_line
!  num_of_k_point = 1
  do while(1==1)
! version dependence -------------------------------
    if (reador(1:28) == 'Time for iouter   (hamilt) :') then
      read(26, '( 40x, i9 )') read_nlo
      if (now_num_of_k_point == 1) then
        write(6, * ) 'WIEN2k version 9'
      endif
      pre_line_nlo_chara = 'Time for iouter   (hamilt) :'
    endif
! version dependence -------------------------------
    if (reador(1:28) == 'Time for distrib  (hamilt, c') then
      read(26, '( 40x, i9 )') read_nlo
      if (now_num_of_k_point == 1) then
        write (6, * ) 'WIEN2k version 10 and later'
      endif
      pre_line_nlo_chara = 'Time for distrib  (hamilt, c'
    endif
! version dependence -------------------------------
    read(26, '( a130 )') reador
!    if ( (reador(1:7) == read_k_pre_chara) .and. (reador(41:50) == input_k_point) ) exit
!    if ( (reador(1:7) == read_k_pre_chara)) then
!      write(6, '( a130 )') reador
!    end if
    !if ( (reador(1:7) == read_k_pre_chara) .and. (reador(41:50) == num_of_k_point_chara(now_num_of_k_point)) ) then
!      write(6, '( a130 )') reador
    !  exit
    !end if
! lapwso  start1
    if ( (reador(1:7) == read_k_pre_chara) ) then
      kcount = kcount + 1 
!      write(6, *) kcount
    end if
    if ( kcount == now_num_of_k_point) then
!      write(6, *) kcount
      kcount = 0
      exit
    end if
! lapwso  start1
  end do
  !
  read(26, '( 17x, i5 )') matrix_num
  !
  read_data_line = 2
  !
  do while(1==1)
    read(26, '( a130 )') reador
    if (reador(123:130) == 'IMAGPART') then
      read_data_line = 3
    endif
    if (reador(1:28) == pre_line_nlo_chara(1:28) ) then
      read (26, '( 40x, i9 )') read_nlo
      exit
    endif
    if (reador(1:26) == '       NUMBER OF K-POINTS:') exit
  end do
  !
  num_of_miller_index = matrix_num - read_nlo
! ----------
  rewind (26)
  !
  wavevector_line_set_up_twice = 0
  !
  do while(1==1)
    read(26, '( a130 )') reador
!    if ( (reador(1:7) == read_k_pre_chara) .and. (reador(41:50) == input_k_point) ) exit
!    if ( (reador(1:7) == read_k_pre_chara) .and. (reador(41:50) == num_of_k_point_chara(now_num_of_k_point)) ) exit
! lapwso start2
    if ( (reador(1:7) == read_k_pre_chara) ) then
      kcount = kcount + 1 
    end if
    if ( kcount == now_num_of_k_point) then
!      write(6, *) kcount
      kcount = 0
      exit
    end if
! lapwso end2
  end do
  write(1, '( a130 )') reador
  !
  ! check file
  if (check_file_option >= 2) then
    ! file f80 is check file ( case.output1 for ********* )
    open (80, file = 'f80')
  end if
  ! check file
  !
  i = 1 ! i = line number
  do while(1==1)
    read(26, '( a130 )') reador
! add replace ********* by 0.000000----------------------------------
    i = i + 1 ! i = line number
    j = j + 1 ! j = start near APW coefficient
    !
    if (reador(1:29) == eigenvalue_end_line_chara) then
      ff_start_line = i
    endif
    if (reador(1:29) == ff_start_end_line_chara) then
      j = 0
    endif
    now_matrix_num = j / read_data_line
    !
    if (i > ff_start_line) then
! 1
      if (reador(17:27) == '  *********') then
        reador(17:27)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  17, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
! warning_chara = 'Warning, find out ********* in PW coefficient'
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 2
      if (reador(28:38) == '  *********') then
        reador(28:38)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  28, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 3
      if (reador(39:49) == '  *********') then
        reador(39:49)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  39, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 4
      if (reador(50:60) == '  *********') then
        reador(50:60)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  50, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 5
      if (reador(61:71) == '  *********') then
        reador(61:71)  = '   0.000000'
        if (check_file_option >= 2) then
          write (80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write (80,  * ) '********* existed at line ', i, ', row  61, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 6
      if (reador(72:82) == '  *********') then
        reador(72:82)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  72, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 7
      if (reador(83:93) == '  *********') then
        reador(83:93)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  83, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 8
      if (reador (94:104) == '  *********') then
        reador(94:104)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  94, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
! 9
      if (reador(105:115) == '  *********') then
        reador(105:115)  = '   0.000000'
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row 105, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
      endif
    endif
    ! --------------------------
    write (1, '( a130 )') reador
    if ((reador (1:7) == read_k_pre_chara) .or. (reador (1:26) == '       NUMBER OF K-POINTS:')) exit
  end do
  !
  close (26)
  close (1)
  !
  !
  if ( HR_i == 1 ) then
    write( 6, * ) 'READ |G|^2 DATA'
  end if
!  write( 6, *) ' # k point, # Energy, # APW    , Max|2(G+k)|'
  !
  open (1, file = 'f01')
  rewind (1)
  !
  read (1, '( 7x, 3f10.5, 3x, 10a  )') miller_index_h, miller_index_k, miller_index_l, read_now_k_point_chara
  read (1, '( 17x, i5 )') num_matrix
  !
  i = 3
  do
    ! read needless data -------------------------------------------------
    read (1, '( a29 )') dum
    ! eigenvalue start line ----------------------------------------------
    if (dum == eigenvalue_start_line_chara) then
      eigenvalue_start_line_num = i
    endif
    ! eignevalue end line ------------------------------------------------
    if (dum == eigenvalue_end_line_chara) then
      eigenvalue_end_line_num = i
    endif
    ! ff start line ------------------------------------------------------
    if (dum == ff_start_end_line_chara) then
      if (wavevector_line_set_up_twice == 0) then
        ff_start_line_num = i
        wavevector_line_set_up_twice = wavevector_line_set_up_twice + 1
        !
        read(1, '( 122x, a29 )') dum
        read(1, '( 122x, a29 )') dum
        read(1, '( 122x, a29 )') dum
        read(1, '( 122x, a29 )') dum
        read(1, '( 122x, a29 )') dum
        !
        if (dum == imaginary_part_chara) then
          if ( HR_i == 1 ) then
            write(6, * ) 'INCLUDE IMAGINARY PART CALCULATION'
          end if
          read_data_line = 3
        else
          if ( HR_i == 1 ) then
            write(6, * ) 'ONLY REAL PART CALCULATION'
          end if
          read_data_line = 2
      endif
      ! ff end line --------------------------------------------------------
      elseif (wavevector_line_set_up_twice == 1) then
        ff_end_line_num = i + 5
        wavevector_line_set_up_twice = wavevector_line_set_up_twice + 1
      endif
    endif
    !
    if ( (dum(1:20) == '         allocate HS') .and. (wavevector_line_set_up_twice == 1) ) then
      ff_end_line_num = i + 5
      wavevector_line_set_up_twice = wavevector_line_set_up_twice + 1
    endif
    !
    if (dum (1:28) == pre_line_nlo_chara(1:28) ) then
      read (1, '( 40x, i9 )') nlo_num
    endif
    if (dum (1:26) == '       NUMBER OF K-POINTS:') then
      nlo_num = read_nlo
    endif
    if ((dum (1:28) == pre_line_nlo_chara(1:28) ) .or. (dum (1:26) == '       NUMBER OF K-POINTS:')) exit
  i = i + 1
  end do
!---- read line end-----------------------------------------------------
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write (6, * ) 'read line end'
  endif
! check end ----------------------------------------------------------
  !
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'set data calculation start 1'
  endif
! check end ----------------------------------------------------------
  !
!---- set data calculation start 1--------------------------------------
  energy_line_num = eigenvalue_end_line_num - eigenvalue_start_line_num - 2
  eigen_value_num = energy_line_num * 5
  ff_coefficient_line_num = ff_end_line_num - ff_start_line_num
  true_eigenvalue_start_line_num = eigenvalue_start_line_num + 1
  true_eigenvalue_end_line_num = eigenvalue_end_line_num - 2
! initialize
  ff_line_set_num = 0
  now_search_ff_line_num = 0
! initialize
  energy_line_num_divided_nine = energy_line_num / 9
  int_energy_line_num_divided = NINT(energy_line_num_divided_nine)
  temp_eigen_line_calc_num = 5 * (energy_line_num - int_energy_line_num_divided)
  if (MOD(temp_eigen_line_calc_num, 9) /= 0) then
    true_energy_line_num = REAL(temp_eigen_line_calc_num) / 9 + 0.5
  else
    true_energy_line_num = REAL(temp_eigen_line_calc_num) / 9
  endif
  end_line_num = ff_coefficient_line_num * NINT(true_energy_line_num) + ff_start_line_num
  ff_array_number = num_matrix - nlo_num
! check start --------------------------------------------------------
  if (check_file_option >= 2) then
    ! file f90 is check file
    open (90, file = 'f90')
    write(90, * ) 'basic data'
    write(90, * ) '----------------------------------------------'
    write(90, * ) 'eigen value line             : ', energy_line_num
    write(90, * ) 'eigen value number         : <=', eigen_value_num
    write(90, * ) '(h,k,l) and coefficient line : ', ff_coefficient_line_num
    write(90, * ) 'read eige nvalue start line  : ', true_eigenvalue_start_line_num
    write(90, * ) 'read eige nvalue end line    : ', true_eigenvalue_end_line_num
    write(90, * ) 'int_energy_line_num_divided  : ', int_energy_line_num_divided
    write(90, * ) 'energy_line_num_divided_nine : ', energy_line_num_divided_nine
    write(90, * ) 'true_energy_line_num         : ', true_energy_line_num
    write(90, * ) 'end line                     : ', end_line_num
    write(90, * ) 'matrix size                  : ', num_matrix
    write(90, * ) 'nlo (hamilt)                 : ', nlo_num
    write(90, * ) 'matrix size - nlo (hamilt)   : ', ff_array_number
    write(90, * ) 'K point                      : ', read_now_k_point_chara
    write(90, * ) 'h miller index               : ', miller_index_h
    write(90, * ) 'k miller index               : ', miller_index_k
    write(90, * ) 'l miller index               : ', miller_index_l
    write(90, * ) '----------------------------------------------'
    close (90)
  endif
  ! write(6,*) "Matrix:", num_matrix, " - # of LO", nlo_num, " = # of PW:", ff_array_number
  energy_array_number = NINT(true_energy_line_num) * 9
! check end ----------------------------------------------------------
  num_eigenvalue_line = 1
!---- set data calculation end
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'set data calculation end 1'
  endif
  !
  close (1)
  write( 6, * ) 'Now number of k point =', now_num_of_k_point, '/', end_of_k_point_line
! check end ----------------------------------------------------------
!
!
  energy_array_number = energy_array_number + 100
  ff_array_number = ff_array_number + 100
  allocate ( eigenvalue_array_Ry(energy_array_number) )                    ! eigenv(ENNM)
!
!
! f01 : read G vs coefficients data from cutting file(cut 26)
  open (1, file = 'f01')
  !
  rewind (1)
! read eigenvalues, ff coefficients and energies
  now_num_eigenvalue_line = 0
  do i = 1, true_eigenvalue_end_line_num
! read needless data -------------------------------------------------
    if ( i < true_eigenvalue_start_line_num ) then
      read(1, '( a29 )') dum
    endif
! read needless data -------------------------------------------------
! read eigenvalue ----------------------------------------------------
    if ( (i >= true_eigenvalue_start_line_num) .and. (i <= true_eigenvalue_end_line_num) ) then
      num_eigenvalue_line = i - true_eigenvalue_start_line_num + 1
      if (MOD(num_eigenvalue_line, 9) /= 0) then
        now_num_eigenvalue_line = now_num_eigenvalue_line + 1
! read eigenvalue --------------------------------------------------
        read(1, '( 2x, 5f13.7 )') (eigenvalue_array_Ry(5 * (now_num_eigenvalue_line - 1) + j) , j = 1, 5)
      else
! read empty -------------------------------------------------------
        read(1, '( a29 )') dum
      endif
    endif
  end do
  !
  cut_energy_line_num = 0
  num_eigenvalue_array_Ry = NINT(true_energy_line_num) - 1
  do i = 1, num_eigenvalue_array_Ry
    if ( eigenvalue_array_Ry(i * 9) <= cut_Ry_threshold_cons ) then
      cut_energy_line_num = i
    end if
  end do
  energy_array_number = energy_array_number - (cut_energy_line_num * 9)
  ff_array_number = ff_array_number + 100
  !
!**********************************************************************************************************************************
  ! -----------------------
!  allocate ( eigenvalue_array_Ry(energy_array_number) )                    ! eigenv(ENNM)
  allocate ( energy_num_chara(energy_array_number) )                       ! energy(ENNM)
  allocate ( wave_K2(ff_array_number) )                                    ! wave_K2(FFNM)
  allocate ( wave_coef_real(ff_array_number, energy_array_number) )        ! wavecr(FFNM, ENNM)
  allocate ( wave_coef_imag(ff_array_number, energy_array_number) )        ! waveci(FFNM, ENNM)
  allocate ( wave_coef_2(ff_array_number, energy_array_number) )           ! wavcr2(FFNM, ENNM)
  if( bond_state_switch == 1 ) then
    allocate ( wave_coef_2_bond(ff_array_number, energy_array_number) )      ! wavec2a(FFNM, ENNM)
    allocate ( wave_coef_2_anti_bond(ff_array_number, energy_array_number) ) ! wavec2b(FFNM, ENNM)
  allocate ( shelter_array_bond(energy_array_number) )                     ! shelt4b(ENNM)
  allocate ( shelter_array_anti_bond(energy_array_number) )                ! shelt4a(ENNM)
  else
    allocate ( wave_coef_2_bond(1,1) )      ! wavec2a(FFNM, ENNM)
    allocate ( wave_coef_2_anti_bond(1,1) ) ! wavec2b(FFNM, ENNM)
    allocate ( shelter_array_bond(1) )                     ! shelt4b(ENNM)
    allocate ( shelter_array_anti_bond(1) )                ! shelt4a(ENNM)
  endif
!  allocate ( degenerate_pw(ff_array_number) )
  allocate ( shelter_array(energy_array_number) )                          ! shelt4(ENNM)
  allocate ( eigenvalue_array_eV(energy_array_number) )                    ! eigeve(ENNM)
  allocate ( eigenvalue_array_sub_ef_eV(energy_array_number) )             ! neigev(ENNM)
!
! H. Sato program line start
  allocate ( bk2(ff_array_number) )
  allocate ( fk2(ff_array_number) )
! H. Sato program line end
!
  ! -----------------------
  !
  energy_array_number = energy_array_number - 100
  ff_array_number = ff_array_number - 100
!**********************************************************************************************************************************
!  call initialize_data(energy_array_number,energy_num_chara, shelter_array, eigenvalue_array_Ry, &
!  eigenvalue_array_eV, eigenvalue_array_sub_ef_eV, wave_num, ff_array_number, wave_K2, &
!  wave_coef_real, wave_coef_imag, wave_coef_2, wave_coef_2_bond, wave_coef_2_anti_bond)
  !
!
!  write( 6, * ) 'Now number of k point =', now_num_of_k_point, '/', end_of_k_point_line
!
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write (6, * ) 'initialize array start'
  endif
! check end ----------------------------------------------------------
  !
  do i = 1, energy_array_number
    energy_num_chara(i) = 'XXXXXXXXXXX'
    shelter_array(i) = 0.0
    eigenvalue_array_Ry(i) = 0.0
    eigenvalue_array_eV(i) = 0.0
    eigenvalue_array_sub_ef_eV(i) = 0.0
  end do
  !
  wave_num(1) = 0
  wave_num(2) = 0
  wave_num(3) = 0
  !
  do j = 1, ff_array_number
    wave_K2(j) = 0.0
!
! H. Sato program line start
    bk2(j) = 0.0
    fk2(j) = 0.0
! H. Sato program line end
!
    do k = 1, energy_array_number
      wave_coef_real(j, k) = 0.0
      wave_coef_imag(j, k) = 0.0
      wave_coef_2(j, k) = 0.0
      if( bond_state_switch == 1 ) then
        wave_coef_2_bond(j, k) = 0.0
        wave_coef_2_anti_bond(j, k) = 0.0
      else
        wave_coef_2_bond(1, 1) = 0.0
        wave_coef_2_anti_bond(1, 1) = 0.0
      endif
    end do
  end do
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write (6, * ) 'initialize array end'
  endif
! check end ----------------------------------------------------------
!**********************************************************************************************************************************
!  call read_g_coefficients_data(energy_array_number, ff_array_number, &
!  eigenvalue_array_Ry, energy_num_chara, wave_num, wave_coef_real, wave_coef_imag, wave_coef_2 &
!  )
  ! read_g_coefficients_data() call K_calculation()
!  call K_calculation(ff_array_number, wave_K2, wave_num)
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write (6, * ) 'read start'
  endif
! check end ----------------------------------------------------------
!
! f01 : read G vs coefficients data from cutting file(cut 26)
  open (1, file = 'f01')
  !
  rewind (1)
! read eigenvalues, ff coefficients and energies
  now_num_eigenvalue_line = 0
  ff_line_set_num = -cut_energy_line_num
  do i = 1, end_line_num
! read needless data -------------------------------------------------
    if ( i < true_eigenvalue_start_line_num ) then
      read(1, '( a29 )') dum
    endif
! read needless data -------------------------------------------------
! read eigenvalue ----------------------------------------------------
    if ( (i >= true_eigenvalue_start_line_num) .and. (i <= true_eigenvalue_end_line_num) ) then
      num_eigenvalue_line = i - true_eigenvalue_start_line_num + 1
      if (MOD(num_eigenvalue_line, 9) /= 0) then
        now_num_eigenvalue_line = now_num_eigenvalue_line + 1
! read eigenvalue --------------------------------------------------
        read(1, '( 2x, 5f13.7 )') (eigenvalue_array_Ry(5 * (now_num_eigenvalue_line - 1) + j) , j = 1, 5)
      else
! read empty -------------------------------------------------------
        read(1, '( a29 )') dum
      endif
    endif
! read eigenvalue ----------------------------------------------------
! read needless data -------------------------------------------------
    if ( (i > true_eigenvalue_end_line_num) .and. (i < ff_start_line_num) ) then
      read(1, '( a29 )') dum
!      write(6, '( a29 )') dum
    endif
! read needless data -------------------------------------------------
! read ff coefficients and energies ----------------------------------
    if (i >= ff_start_line_num) then
! reset ff coefficient and energy line -------------------------------
      now_search_ff_total_line_num = i - ff_start_line_num
      now_search_ff_line_num = MOD(now_search_ff_total_line_num, ff_coefficient_line_num)
      ! through needless line
      if ( (now_search_ff_line_num == 0) .or. &
           (now_search_ff_line_num == 1) .or. &
           (now_search_ff_line_num == (ff_coefficient_line_num - 3) ) .or. &
           (now_search_ff_line_num == (ff_coefficient_line_num - 2) ) .or. &
           (now_search_ff_line_num == (ff_coefficient_line_num - 1) ) &
      ) then
! read empty -------------------------------------------------------
        read(1, '( a29 )') dum
        if ( (now_search_ff_line_num == 0) .and. (dum /= ff_start_end_line_chara) ) exit
        if (dum(1:26) == '       NUMBER OF K-POINTS:') exit
        !
      else ! get g and coefficient data
        if (now_search_ff_line_num == 2) then
          ff_line_set_num = ff_line_set_num + 1
          if (ff_line_set_num > 0) then
! read energies ------------------------------------------------------
            read(1, '( 17x, 9a11 )') (energy_num_chara(9 * (ff_line_set_num - 1) + j) , j = 1, 9)
!          WRITE(6, '( 17x, 9a11 )') (energy_num_chara(9 * (ff_line_set_num - 1) + j) , j = 1, 9)
          else
            read(1, '( a29 )') dum
          end if
        endif
        !
! wave number line start from "now_search_ff_line_num = 3"
        if (read_data_line == 2) then
          real_imaginary_switch = 1
        else
          real_imaginary_switch = 0
        endif
        !
        if ( (MOD(now_search_ff_line_num, read_data_line) == real_imaginary_switch) .and. (now_search_ff_line_num >= 3) ) then
          ff_count = now_search_ff_line_num / read_data_line
          if ( (ff_line_set_num > 0) .and. (ff_count <= ff_array_number) ) then
! read wave value ----------------------------------------------------
            read(1, '( 3i4 )') (wave_num(j) , j = 1, 3)
! K calculate --------------------------------------------------------
!**********************************************************************************************************************************
!            call K_calculation()
! check start --------------------------------------------------------
  if (check_file_option == 1) then
    write(6, * ) 'G**2 calculation start'
  endif
! check end ----------------------------------------------------------
            !
!  i = ff_count
! calculate squared K ( K = G + k ) once
!
  if (ff_line_set_num == 1) then
    !
! H. Sato program line start
! Bravais lattice (monoclinic crystal {0.0 < angle_c < 180})
! C.J.Bradley and A.P.Cracknell., "The mathematical theroy of symmetry in solids", Oxford Univ. Press, (1972) 85.
      !
      if ( (angle_a == 90 ) .and. (angle_b == 90) .and. (angle_c ==90 ) ) then
      BXX=  1.0
      BXY=  0.0
      BXZ=  0.0
      !
      BYX=  0.0
      BYY=  1.0
      BYZ=  0.0
      !
      BZX=  0.0
      BZY=  0.0
      BZZ=  1.0
          !
      else if ( (angle_a == 90) .and. (angle_b ==90) .and. (angle_c == 120) ) then
      BXX = 2.0/sqrt(3.0)
      BXY = 1.0/sqrt(3.0)*(axis_b/axis_a)
      BXZ = 0.0
      !
      BYX=  0.0
      BYY=  1.0
      BYZ=  0.0
      !
      BZX=  0.0
      BZY=  0.0
      BZZ=  1.0
          !
      else
        BXX = BG1X * axis_a
        BXY = BG2X * axis_b
        BXZ = BG3X * axis_c
        BYX = BG1Y * axis_a
        BYY = BG2Y * axis_b
        BYZ = BG3Y * axis_c
        BZX = BG1Z * axis_a
        BZY = BG2Z * axis_b
        BZZ = BG3Z * axis_c
      end if
      !
      axis_n = axis_a
! calculate squared G ( G = G + K )
        WVX=BXX*( float( wave_num(1) ) + miller_index_h ) * (2 * axis_n / axis_a) &
           +BXY*( float( wave_num(2) ) + miller_index_k ) * (2 * axis_n / axis_b) &
           +BXZ*( float( wave_num(3) ) + miller_index_l ) * (2 * axis_n / axis_c)
        !
        WVY=BYX*( float( wave_num(1) ) + miller_index_h ) * (2 * axis_n / axis_a) &
           +BYY*( float( wave_num(2) ) + miller_index_k ) * (2 * axis_n / axis_b) &
           +BYZ*( float( wave_num(3) ) + miller_index_l ) * (2 * axis_n / axis_c)
        !
        WVZ=BZX*( float( wave_num(1) ) + miller_index_h ) * (2 * axis_n / axis_a) &
           +BZY*( float( wave_num(2) ) + miller_index_k ) * (2 * axis_n / axis_b) &
           +BZZ*( float( wave_num(3) ) + miller_index_l ) * (2 * axis_n / axis_c)
        !
        wave_K2(ff_count) = WVX**2 + WVY**2 + WVZ**2
        ! LO search
        !if( (wave_K2_old > wave_K2(ff_count) ) .and. (LO_search_flag == 0) ) then
        !  LO_search_flag = 1
        !  new_ff_array_number = (ff_count - 1)/2*1.90
        !  write(6,*) "find # of LO:",int ((ff_count - 1)/2*1.90)
        !end if
        !wave_K2_old = wave_K2(ff_count)
!
! test all symmetry group program
!        WVX=BG1X*( float( wave_num(1)*2 ) + miller_index_h*2 ) &
!           +BG1Y*( float( wave_num(2)*2 ) + miller_index_k*2 ) &
!           +BG1Z*( float( wave_num(3)*2 ) + miller_index_l*2 )
        !
!        WVY=BG2X*( float( wave_num(1)*2 ) + miller_index_h*2 ) &
!           +BG2Y*( float( wave_num(2)*2 ) + miller_index_k*2 ) &
!           +BG2Z*( float( wave_num(3)*2 ) + miller_index_l*2 )
        !
!        WVZ=BG3X*( float( wave_num(1)*2 ) + miller_index_h*2 ) &
!           +BG3Y*( float( wave_num(2)*2 ) + miller_index_k*2 ) &
!           +BG3Z*( float( wave_num(3)*2 ) + miller_index_l*2 )
        !
!        wave_K2(ff_count) = ( WVX**2 + WVY**2 + WVZ**2 ) * (correct_volume_factor**(1.0/3.0))
! H. Sato program line end
!
!    if (hcp_on == 0) then
!      wave_K2(ff_count) = &
!      ( (float(wave_num(1) ) + miller_index_h) * (2 * axis_n / axis_a) ) **2 + &
!      ( (float(wave_num(2) ) + miller_index_k) * (2 * axis_n / axis_b) ) **2 + &
!      ( (float(wave_num(3) ) + miller_index_l) * (2 * axis_n / axis_c) ) **2
!    endif
!    if (hcp_on == -1) then
!      wave_K2(ff_count) = &
!      ( (float(wave_num(1) ) + miller_index_h) * 2) **2 + &
!      ( (float(wave_num(2) ) + miller_index_k) * 2) **2 + &
!      ( (float(wave_num(3) ) + miller_index_l) * 2) **2
!    endif
! ----normalized for hcp
!    if (abs(angle_c - 120.0*PI/180.0) <= 0.001 ) then
!      wave_K2(ff_count) = &
!      ( ( (float(wave_num(1) ) + miller_index_h) * 2 + (float(wave_num(2) ) + miller_index_k) ) * &
!      (1 / 1.732050808) * (2 * axis_n / axis_a) ) **2 + &
!      ( ( (float(wave_num(2) ) + miller_index_k) ) * (2 * axis_n / axis_b) ) **2 + &
!      ( ( (float(wave_num(3) ) + miller_index_l) ) * (2 * axis_n / axis_c) ) **2
!    endif
!    write (6, *) ff_count, wave_K2(ff_count)
! check end ----------------------------------------------------------
  if (check_file_option == 1) then
    write(6, * ) 'G**2 calculation end'
  endif
! check end ----------------------------------------------------------
  !
  endif
!
!**********************************************************************************************************************************
! read ff coefficients at real part ----------------------------------
            read(1, '( 16x, 9f11.6 )') (wave_coef_real(ff_count, ( (ff_line_set_num - 1) * 9 + j) ) , j = 1, 9)
!            WRITE(6, '( 16x, 9f11.6 )') (wave_coef_real(ff_count, ( (ff_line_set_num - 1) * 9 + j) ) , j = 1, 9)
! read ff coefficients at imaginary part -----------------------------
            if (read_data_line == 3) then
              read(1, '( 16x, 9f11.6 )') (wave_coef_imag(ff_count, ( (ff_line_set_num - 1) * 9 + j) ) , j = 1, 9)
!              WRITE(6, '( 16x, 9f11.6 )') (wave_coef_imag(ff_count, ( (ff_line_set_num - 1) * 9 + j) ) , j = 1, 9)
            endif
          else
            read(1, '( a29 )') dum
            read(1, '( a29 )') dum
            if (read_data_line == 3) then
              read(1, '( a29 )') dum
            endif
          endif
        endif
        !
      endif
    endif
!    WRITE(6, *) i, end_line_num, j, ff_count,ff_array_number, &
!wave_K2(ff_count), wave_coef_real(ff_count,10), now_search_ff_line_num
  end do
!---- read end----------------------------------------------------------
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'read end'
  endif
! check end ----------------------------------------------------------
  !
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'set data calculation start 2'
  endif
! check end ----------------------------------------------------------
! when include imaginary part, repeat by 3 line (read_data_line = 3)
! When not include imagirnay part, repeat by 2 line (read_data_line = 2)
!      ff_end_line_num = ( ff_coefficient_line_num - 5 ) / read_data_line
! or ff_end_line_num = ff_count is OK, also.
!      ff_end_line_num = ff_count
  do i = 1, energy_array_number
    eigenvalue_array_Ry(i) = eigenvalue_array_Ry(i + (cut_energy_line_num * 9))
    if ((energy_num_chara(i) == empty_chara) .or. (energy_num_chara(i) == 'XXXXXXXXXXX')) then
      energy_line_end_num = i - 1
      exit
    endif
!    write(6, * ) energy_num_chara(i), eigenvalue_array_Ry(i)
  enddo
!---- set data calculation end2-----------------------------------------
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'set data calculation end 2'
  endif
! check end ----------------------------------------------------------
  !
  close (1)
!**********************************************************************************************************************************
!  call squared_coefficient_calculate(energy_array_number, ff_array_number, wave_coef_2, &
!    wave_coef_real, wave_coef_imag, wave_coef_2_bond, wave_coef_2_anti_bond)
!
  !if( new_ff_array_number <= 0 )then
  !  ff_array_number = new_ff_array_number
  !end if
  !ff_array_number = you select number of PW
  do i = 1, ff_array_number
    do j = 1, energy_array_number
      wave_coef_2(i, j) = wave_coef_real(i, j) **2 + wave_coef_imag(i, j) **2
      !
! ----bond, anti-bond calculation
      if( bond_state_switch == 1 ) then
        if (read_data_line == 2) then
          if (wave_coef_real(i, j) >= 0) then
            wave_coef_2_bond(i, j) = wave_coef_2(i, j)
          endif
          if (wave_coef_real(i, j) < 0) then
            wave_coef_2_anti_bond(i, j) = -1 * wave_coef_2(i, j)
          endif
        endif
      endif
! ----bond, anti-bond calculation
      !
    end do
  end do
!**********************************************************************************************************************************
!  call bubble_sort(wave_K2, wave_coef_2, wave_coef_real, wave_coef_imag, &
!  wave_coef_2_bond, wave_coef_2_anti_bond, shelter_array)
!
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'bubble sort start'
    write(6, * ) 'about', ff_end_line_num, '/2 second'
  endif
! check end ----------------------------------------------------------
  !
!---- bubble sort start-------------------------------------------------
! very easy sort program
! for example, from bottom array to top array,  arrange like
! ----------------------------------------------------------------------
! G**2 small
! ..........
! G**2 large
! ----------------------------------------------------------------------
  do i = 2, ff_array_number
    do j = ff_array_number, i, - 1
      if (wave_K2(j - 1) > wave_K2(j) ) then
! main; sort G**2
        shelter = wave_K2(j - 1)
        wave_K2(j - 1) = wave_K2(j)
        wave_K2(j) = shelter
! sub: sort wave coefficient according to G**2

        do k = 1, energy_line_end_num
          shelter_array(k) = wave_coef_2(j - 1, k)
          wave_coef_2(j - 1, k) = wave_coef_2(j, k)
          wave_coef_2(j, k) = shelter_array(k)
          !
! ----bond, anti-bond calculation
          if( bond_state_switch == 1 ) then
            shelter_array(k) = wave_coef_2_bond(j - 1, k)
            wave_coef_2_bond(j - 1, k) = wave_coef_2_bond(j, k)
            wave_coef_2_bond(j, k) = shelter_array(k)
            !
            shelter_array(k) = wave_coef_2_anti_bond(j - 1, k)
            wave_coef_2_anti_bond(j - 1, k) = wave_coef_2_anti_bond(j, k)
            wave_coef_2_anti_bond(j, k) = shelter_array(k)
          endif
! ----bond, anti-bond calculation
          !
        end do
      endif
    end do
!    write(6,*) i, ff_array_number, wave_K2(i)
  end do
!  do i = 1, ff_array_number
!    write(6,*) wave_K2(i)
!  end do
  !
! check start --------------------------------------------------------
!  if (check_file_option >= 1) then
!    write(6, * ) 'bubble sort end'
!  endif
! check end ----------------------------------------------------------
!**********************************************************************************************************************************
!  call assort(wave_K2, degenerate_pw, &
!  wave_coef_2, shelter_array, shelter_array_bond, wave_coef_real, shelter_array_anti_bond, &
!  wave_coef_imag, wave_coef_2_bond, wave_coef_2_anti_bond)
  !
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'assort start'
  endif
! check end ----------------------------------------------------------
  !
!---- assort start------------------------------------------------------
!---- initialize start------------------------------------------------
! for G^2, shelt3 change from i type to f type
  shelter = wave_K2(1)
!  do k = 1, ff_array_number
!    degenerate_pw(k) = 1
!  end do

  do k = 1, energy_line_end_num
    shelter_array(k) = wave_coef_2(1, k)
    if( bond_state_switch == 1 ) then
      shelter_array_bond(k) = wave_coef_2_bond(1, k)
      shelter_array_anti_bond(k) = wave_coef_2_anti_bond(1, k)
    endif
  end do
  new_pw_data_num = 1
  new_ff_line_end_num = 0
!  degeneracy_num = 1
!---- initialize end--------------------------------------------------
! assort  ------------------------------------------------------------
  do i = 1, ff_array_number - 1
    if (wave_K2(i) == wave_K2(i + 1) ) then
      shelter = wave_K2(i + 1)
!      degeneracy_num = degeneracy_num + 1
      do k = 1, energy_line_end_num
        shelter_array(k) = shelter_array(k) + wave_coef_2(i + 1, k)
        if( bond_state_switch == 1 ) then
          shelter_array_bond(k) = shelter_array_bond(k) + wave_coef_2_bond(i + 1, k)
          shelter_array_anti_bond(k) = shelter_array_anti_bond(k) + wave_coef_2_anti_bond(i + 1, k)
        endif
      end do
    else
      new_ff_line_end_num = new_ff_line_end_num + 1
      wave_K2(new_pw_data_num) = shelter
!      degenerate_pw(new_pw_data_num) = degeneracy_num
!      degeneracy_num = 0
      do k = 1, energy_line_end_num
        wave_coef_2(new_pw_data_num, k) = shelter_array(k)
        if( bond_state_switch == 1 ) then
          wave_coef_2_bond(new_pw_data_num, k) = shelter_array_bond(k)
          wave_coef_2_anti_bond(new_pw_data_num, k) = shelter_array_anti_bond(k)
        endif
      end do
      shelter = wave_K2(i + 1)
      do k = 1, energy_line_end_num
        shelter_array(k) = wave_coef_2(i + 1, k)
        if( bond_state_switch == 1 ) then
          shelter_array_bond(k) = wave_coef_2_bond(i + 1, k)
          shelter_array_anti_bond(k) = wave_coef_2_anti_bond(i + 1, k)
        endif
      end do
!write(6, *) wave_K2(new_pw_data_num)
      new_pw_data_num = new_pw_data_num + 1
    endif
  end do
  !
!---- next initial data-----------------------------------------------
  do i = new_pw_data_num, ff_array_number
    wave_K2(i) = 0.0
    do k = 1, energy_line_end_num
      wave_coef_2(i, k) = 0.0
      if( bond_state_switch == 1 ) then
        wave_coef_2_bond(i, k) = 0.0
        wave_coef_2_anti_bond(i, k) = 0.0
      endif
    end do
  end do
  !
!---- next initial data-----------------------------------------------
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'assort end'
  endif
! check end ----------------------------------------------------------
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!subroutine eV_calculate()
!  implicit none
!  use ff_cons, only : check_file_option, convert_Ry2eV_cons
!  use ff_read_data, only : eigenvalue_array_eV, eigenvalue_array_Ry, eigenvalue_array_sub_ef_eV, &
!    ff_array_number, energy_array_number
!  use character_cons, only :
!  use character_set, only :
!  use integer_set, only : new_pw_data_num, new_ff_line_end_num, degeneracy_num, ff_end_line_num, energy_line_end_num
!  use real_set, only : num_eigenvalue_sub_ef_Ry
  !
!  integer i, j, k ! use do loop & read line
  !
! check start --------------------------------------------------------
  if (check_file_option == 1) then
    write(6, * ) 'energy calculate start'
  endif
! check end ----------------------------------------------------------
  !
  do i = 1, energy_line_end_num
    eigenvalue_array_eV(i) = convert_Ry2eV_cons * eigenvalue_array_Ry(i)
    num_eigenvalue_sub_ef_Ry = eigenvalue_array_Ry(i) - ef_energy_Ry
    eigenvalue_array_sub_ef_eV(i) = convert_Ry2eV_cons * num_eigenvalue_sub_ef_Ry
  end do
!---- energy output start-----------------------------------------------
!---- energy output end-------------------------------------------------
!  return
!end eV_calculate!***************************************************************************************************
!
!**************************************************************************************************
!  call write_K_vs_coefficients(eigenvalue_array_sub_ef_eV, degenerate_pw, &
!  wave_K2, wave_coef_2, ff_array_number)
!
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'write_K_vs_coefficients start'
  endif
! check end ----------------------------------------------------------
  !
!---- write_K_vs_coefficients ------------------------------------------
! for example, print like ( NOTE!: G^2 = G**2 ) (igor type)
! ----------------------------------------------------------------------
!      l       energy 1      l     energy 2   (eV) l.........   (eV) l
! G^2  lcoefficient for G^2  l .........           l                 l
! ----------------------------------------------------------------------
! use; eigenvalue_array_Ry(energy_array_number), i, energy_line_end_num, &
!ff_end_line_num, wave_coef_2(ff_array_number,energy_array_number),j
! write energy data
!---- write_K_vs_coefficients ------------------------------------------
  !
  ! file f08 is output file

!  open (8, file = 'f08')
  !
!  write(8, 810) (eigenvalue_array_sub_ef_eV(i), i = 1, energy_line_end_num)
!  810 FORMAT( 'degeneracy fact', ',', '   G^2 / eV    ', ',', 999( 'E=', f13.8, ',') )
!  810 FORMAT( '   G^2 / eV    ', ',', 999( 'E=', f13.8, ',') )
  ! write G**2 and coefficient data
!  do i = 1, new_ff_line_end_num
!    write(8, 820) degenerate_pw(i), wave_K2(i), (wave_coef_2(i, j), j = 1, energy_line_end_num)
!  820 FORMAT( i15, ',', f15.8, ',', 999( f15.8, ',' ) )
!    write(8, 820) wave_K2(i), (wave_coef_2(i, j), j = 1, energy_line_end_num)
!  820 FORMAT( f15.8, ',', 999( f15.8, ',' ) )
!  end do
  !
!  close (8)
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'write_K_vs_coefficients end'
  endif
! check end ----------------------------------------------------------
!**********************************************************************************************************************************
!  call write_coefficients_vs_K(wave_K2, wave_coef_2, ff_array_number, degenerate_pw, &
!  wave_coef_2_anti_bond, wave_coef_2_bond)
!
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'write_coefficients_vs_K start'
  endif
! check end ----------------------------------------------------------
  !
!---- write_coefficients_vs_K ------------------------------------------
! for example, print like ( NOTE!: G^2 = G**2 )
! ----------------------------------------------------------------------
!           l         G^2 =        l        G^2 =        l
! energy(eV)l coefficient for G^2  l .........           l
! ----------------------------------------------------------------------
!---- write_coefficients_vs_K ------------------------------------------
  ! file f09 is output file ( change line and row )
!  open (9, file = 'f09')
  !
!  write(9, 900) (wave_K2(i), i = 1, new_ff_line_end_num)
!  900 FORMAT('  ENERGY / eV  ', ',', 999('G^2=', f9.4, '  ,') )
! write G**2 and coefficient data
!  do i = 1, energy_line_end_num
!    write(9, 910) eigenvalue_array_sub_ef_eV(i), (wave_coef_2(j, i), j = 1, new_ff_line_end_num)
!  910 FORMAT( f15.10, ',', 999( f15.8, ',') )
!  end do
!    write(9, 920) (degenerate_pw(i), i = 1, new_ff_line_end_num)
!  920 FORMAT( 15x, ',', 999( i15, ',') )
  !
!  close (9)
  !
!---- write_coefficients_vs_K ------------------------------------------
if( bond_state_switch == 1 ) then
!---- write_coefficients_vs_K ------------------------------------------
! for example, print like ( NOTE!: G^2 = G**2 ) bond calculation results
! ----------------------------------------------------------------------
!           l         G^2 =        l        G^2 =        l

! energy(eV)l coefficient for G^2  l .........           l
! ----------------------------------------------------------------------
!---- write_coefficients_vs_K ------------------------------------------
  ! file f47 is bond and bond calculation results
!  open (47, file = 'f47')
  !
!  write(47, 4700) (wave_K2(i), i = 1, new_ff_line_end_num)
! 4700 FORMAT('  ENERGY / eV  ', ',', 999('G^2=', f9.4, '  ,') )
! write G**2 and coefficient data
!  do i = 1, energy_line_end_num
!    write(47, 4710) eigenvalue_array_sub_ef_eV(i), (wave_coef_2_bond(j, i), j = 1, new_ff_line_end_num)
! 4710 FORMAT( f15.10, ',', 999( f15.8, ',') )
!  end do
!  write(47, 4720) (degenerate_pw(i), i = 1, new_ff_line_end_num)
! 4720 FORMAT( 15x, ',', 999( i15, ',') )
 !
! close (47)
!---- write_coefficients_vs_K ------------------------------------------
! for example, print like ( NOTE!: G^2 = G**2 ) anti-bond calculation re
! ----------------------------------------------------------------------
!           l         G^2 =        l        G^2 =        l
! energy(eV)l coefficient for G^2  l .........           l
! ----------------------------------------------------------------------
!---- write_coefficients_vs_K ------------------------------------------
  ! file f48 is bond and anti-bond calculation results
!  open (48, file = 'f48')
  !
!  write(48, 4800) (wave_K2(i), i = 1, new_ff_line_end_num)
! 4800 FORMAT('  ENERGY / eV  ', ',', 999('G^2=', f9.4, '  ,') )
! write G**2 and coefficient data
!  do i = 1, energy_line_end_num
!    write(48, 4810) eigenvalue_array_sub_ef_eV(i), (wave_coef_2_anti_bond(j, i), j = 1, new_ff_line_end_num)
! 4810 FORMAT( f15.10, ',', 999( f15.8, ',') )
!  end do
!  write(48, 4820) (degenerate_pw(i), i = 1, new_ff_line_end_num)
! 4820 FORMAT( 15x, ',', 999( i15, ',') )
 !
! close (48)
endif
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write(6, * ) 'write_coefficients_vs_K end'
  endif
! check end ----------------------------------------------------------
!**********************************************************************************************************************************
  ! -----------------------
  deallocate ( eigenvalue_array_Ry )                                   ! eigenv(ENNM)
  deallocate ( energy_num_chara )                                      ! energy(ENNM)
!  deallocate ( wave_K2 )                                               ! wave_K2(FFNM)
  deallocate ( wave_coef_real )                                        ! wavecr(FFNM, ENNM)
  deallocate ( wave_coef_imag )                                        ! waveci(FFNM, ENNM)
!  deallocate ( wave_coef_2 )                                           ! wavcr2(FFNM, ENNM)
  deallocate ( wave_coef_2_bond )                                      ! wavec2a(FFNM, ENNM)
  deallocate ( wave_coef_2_anti_bond )                                 ! wavec2b(FFNM, ENNM)
!  deallocate ( degenerate_pw )
  deallocate ( shelter_array )                                      ! shelt4(ENNM)
  deallocate ( shelter_array_bond )                                 ! shelt4b(ENNM)
  deallocate ( shelter_array_anti_bond )                            ! shelt4a(ENNM)
  deallocate ( eigenvalue_array_eV )                                ! eigeve(ENNM)
!  deallocate ( eigenvalue_array_sub_ef_eV )                         ! neigev(ENNM)
  ! -----------------------
  !
!=======================================================================
!**********************************************************************************************************************************
!  call read_line_and_data(ff_array_number)
!
! HR calculation continue
!-------------- -------------------------------------------------------
!  character(16) name_of_column_chara(12) ! COMPNE
!  character(3) end_chara
  !
!  integer num_of_energy_step     ! NOSTEP
!  integer energy_line_end_num    ! enlend
!  integer end_of_k_point_line    ! ENDOKL
!  integer HR_i, HR_j, HR_k, HR_l, HR_m
!  integer N
  !
!  real(8) K2_at_C2_1st_max        ! maxwg2
!  real(8) C2_1st_max              ! maxwc2c
!  real(8) K2_at_C2_2nd_max        ! maxwg2_2nd
!  real(8) C2_2nd_max              ! maxwc2c_2nd
!  real(8) K2_at_C2_max_limited_K  ! maxwg2_Kl
!  real(8) C2_max_limited_K        ! maxwc2c_Kl
!  real(8) First_C2_max_limited_K  ! FPF_K_limit
!  real(8) start_energy            ! STROEN
!  real(8) energy_step             ! ENESTE
!  real(8) K_limit
!  real(8) C2_limit                ! X_limit
!  real(8) E_limit
!  real(8) K2_step_for_3D_plot     ! K2STE_Elimit
!  real(8) num_total_K2_step_for_3D_plot ! TK2STE_Elimit
!  real(8) sigma                   ! SIGMA
!  real(8) sigma_max               ! SIGMA_MAX
!  real(8) sigma2_max              ! SIGMA2_MAX
!  real(8) temp_sigma2_max         ! RSIGMA2
!  real(8) sigma_percent           ! SIGMAP
!  real(8) sigma2_percent          ! SIGMA2P
!  real(8) corrected_total_K_max   ! KGT
!  real(8) corrected_total_K2_max  ! KGT2
!  real(8) corr_total_K_max_limi_K ! KGT_KL
!  real(8) corr_total_K2_max_limi_K! KGT2_KL
!  real(8) num_K2_step             ! K2STEP
!  real(8) energy_for_hr_plot      ! ENERHR
!  real(8) total_ave_C2_max_dw     ! CKT
  !
!  real(8), allocatable :: C2_max_at_k_point(:,:)    ! C2_Kmax(end_of_k_point_line, num_of_energy_step)
!  real(8), allocatable :: weight_at_C2_max(:,:)     ! CTWEGIA_Ckmax(end_of_k_point_line, num_of_energy_step)
!  real(8), allocatable :: total_K_max(:,:)          ! KGA(end_of_k_point_line, num_of_energy_step)
!  real(8), allocatable :: total_K_max_limited_K(:,:)! KGT_K_limit(end_of_k_point_line, num_of_energy_step)
!  real(8), allocatable :: count_weight(:,:)         ! CTWEGIA(end_of_k_point_line, num_of_energy_step)
!  real(8), allocatable :: total_C2_for_3D_plot(:)   ! K2_Elimit_T(num_total_K2_step_for_3D_plot)
!  real(8), allocatable :: total_C2_for_3D_plot_dw(:)! K2_Elimit(num_total_K2_step_for_3D_plot)
!  real(8), allocatable :: weight_for_3D_plot(:)     ! CTWEGIA_Elimit(num_total_K2_step_for_3D_plot)
!  real(8), allocatable :: weight(:)                 ! WEIGHT(end_of_k_point_line)
!  real(8), allocatable :: K2_distribution_at_E(:)   ! DENERGY(dis_end_k_point_line)
!  real(8), allocatable :: total_averaged_K_max(:)   ! KGTA(num_of_energy_step)
!  real(8), allocatable :: true_weight_at_k_point(:) ! CTWETA(num_of_energy_step)
!  real(8), allocatable :: true_weight_at_C2_max(:)  ! CTWETA_Ckmax(num_of_energy_step)
!  real(8), allocatable :: tot_ave_K_max_limit(:)    ! KGT_K_limit_T(num_of_energy_step)
!  real(8), allocatable :: total_ave_C2_max(:)       ! C2_Kmax_T(num_of_energy_step)
!  real(8), allocatable :: sigma2(:)                 ! SIGMA2(num_of_energy_step)
!-------------- -------------------------------------------------------
!    do HR_j = 1, num_of_energy_step
      do HR_k = 1, energy_line_end_num
        K2_at_C2_1st_max       = 0.0
        C2_1st_max             = 0.0
        K2_at_C2_2nd_max       = 0.0
        C2_2nd_max             = 0.0
        K2_at_C2_max_limited_K = 0.0
        C2_max_limited_K       = 0.0
        First_C2_max_limited_K = 1
!
! H. Sato program line start
! Inukai start
!        if( (HR_i == 1) .and. (HR_k == 5) )then
        if( (HR_i == 1) .and. (HR_k == 1) )then
          write(6,*) 'eignvalue list (eV) at gamma point'
          !write(6,*) (j,':', eigenvalue_array_sub_ef_eV(j), j=1, energy_line_end_num)
          if (energy_line_end_num >= 800) then
            write(6,*) (j,':', eigenvalue_array_sub_ef_eV(j), j=1, 800)
          else
           write(6,*) (j,':', eigenvalue_array_sub_ef_eV(j), j=1, energy_line_end_num)
          end if
          write(6,*) 'input eigenvlue number at bottom in Valence band'
                 end_l = 1.0
                 start_k_maxC = 0.0
                              j_k_No = 0
                 
                   do j_l = 1, new_ff_line_end_num 
                     if ( (wave_K2( j_l ) < end_l) .and. (wave_K2( j_l ) >= 0.0) ) then
                       end_l = wave_K2( j_l )
                                        do j_k = 1, energy_line_end_num
                         if ( wave_coef_2( j_l, j_k ) >= start_k_maxC ) then
                           start_k_maxC = wave_coef_2( j_l, j_k )
                           EBOT = eigenvalue_array_sub_ef_eV( j_k )
                                              j_k_No = j_k
                         end if
                                        end do
                     end if
                   end do
              write(6,*) 'program recommend the eignvalue'
              write(6,*) ' No. and eV at gamma point, No.', j_k_No, ':', EBOT
                        k2 = j_k_No
              write(6,*) 'Please input start eigenvalue No.'
              read(5,*) k2
              EBOT = eigenvalue_array_sub_ef_eV(k2)
        end if
        low_num_eigenvale = energy_line_end_num - k2 + 1
!        AKN = ( eigenvalue_array_sub_ef_eV(HR_k) - EBOT ) / CC
! H. Sato program line end
!
!        if( (eigenvalue_array_sub_ef_eV(HR_k) > (start_energy + energy_step*(HR_j - 1) - energy_step/2)) .and. &
!        (eigenvalue_array_sub_ef_eV(HR_k) <= (start_energy + energy_step*(HR_j - 1) + energy_step/2)) ) then
!
! Inukai end
! H. Sato program line start
        do HR_l = 1, new_ff_line_end_num
!      write(6,*) ' e(',k,')=',neigev(k),' AKK=',AKN
          bk2(HR_l) = wave_K2(HR_l)
          fk2(HR_l) = wave_coef_2(HR_l,HR_k)
          if( wave_coef_2(HR_l,HR_k) > C2_1st_max ) then
!            K2_at_C2_2nd_max = K2_at_C2_1st_max
!            C2_2nd_max       = C2_1st_max
            K2_at_C2_1st_max = wave_K2(HR_l)
            C2_1st_max       = wave_coef_2(HR_l,HR_k)
          end if
! H. Sato program line end
!
!-------------- -------------------------------------------------------
!            if( ( sqrt(wave_K2(HR_l)) < K_limit) .and. (K_limit /= 0.0 ) )then
!              if( wave_coef_2(HR_l,HR_k) > C2_max_limited_K ) then
!                K2_at_C2_max_limited_K = wave_K2(HR_l)
!                C2_max_limited_K       = wave_coef_2(HR_l,HR_k)
!              end if
!            end if
!-------------- -------------------------------------------------------
!            if( (wave_coef_2(HR_l,HR_k) > C2_max_limited_K) .and. (K_limit == 0.0 ) .and. (First_C2_max_limited_K == 1) ) then
!              if( wave_coef_2(HR_l,HR_k) > C2_max_limited_K ) then
!                K2_at_C2_max_limited_K = wave_K2(HR_l)
!                C2_max_limited_K       = wave_coef_2(HR_l,HR_k)
!              end if
!            else
!              First_C2_max_limited_K = 0
!            end if
        end do
!
! H. Sato program line start
        do HR_l = 1, new_ff_line_end_num
          do HR_l2 = (HR_l + 1), new_ff_line_end_num
            AVV1 = fk2(HR_l)
            AKK1 = bk2(HR_l)
            if(fk2(HR_l2) > AVV1) then
              fk2(HR_l)  = fk2(HR_l2)
              bk2(HR_l)  = bk2(HR_l2)
              fk2(HR_l2) = AVV1
              bk2(HR_l2) = AKK1
            end if
          end do
        end do
! test remove H.Sato pro by inukai: start
        fftt = 0.0
        do HR_l = 1, new_ff_line_end_num
          w2k = bk2(HR_l)
          w2f = fk2(HR_l)
          fftt = fftt + w2f
          if(HR_l.gt.1.and.abs(bk2(1)-w2k).lt.0.0001) then
          fk2(1)=fk2(1) +w2f
          fk2(HR_l)=0.0
          endif
        end do
! test remove H.Sato pro by inukai: end
!
        if(eigenvalue_array_sub_ef_eV(HR_k) >= EBOT) then
          if(eigenvalue_array_sub_ef_eV(HR_k) == EBOT) then
            k2=HR_k
          end if
! Inukai start
!          if(HR_k <= (k2+21)) then
!            write(4,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
!              bk2(1),bk2(2),bk2(3),bk2(4),bk2(5),  &
!              bk2(6),bk2(7),bk2(8),bk2(9),bk2(10), &
!              bk2(11),bk2(12),bk2(13),bk2(14),bk2(15), &
!              bk2(16),bk2(17),bk2(18),bk2(19),bk2(20), &
!              bk2(21),bk2(22),bk2(23),bk2(24),bk2(25), &
!              bk2(26),bk2(27),bk2(28),bk2(29),bk2(30), &
!              bk2(31),bk2(32),bk2(33),bk2(34),bk2(35), &
!              bk2(36),bk2(37),bk2(38)
!,bk2(39),bk2(40)
!              bk2(41),bk2(42),bk2(43),bk2(44),bk2(45), &
!              bk2(46),bk2(47),bk2(48),bk2(49),bk2(50)
!          end if
          if(HR_k <= energy_line_end_num) then

!            write(4,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
!            (bk2(ff_list_max_num),ff_list_max_num=1, new_ff_line_end_num)
!            write(4,'(i4,f13.8,i8)') HR_k,eigenvalue_array_sub_ef_eV(HR_k),new_ff_line_end_num
!            do ff_list_max_num=1, new_ff_line_end_num
!              write(4, '(f13.8)') bk2(ff_list_max_num)
!            end do
            write(4,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
              bk2(1),bk2(2),bk2(3),bk2(4),bk2(5),  &
              bk2(6),bk2(7),bk2(8),bk2(9),bk2(10), &
              bk2(11),bk2(12),bk2(13),bk2(14),bk2(15), &
              bk2(16),bk2(17),bk2(18),bk2(19),bk2(20), &
              bk2(21),bk2(22),bk2(23),bk2(24),bk2(25), &
              bk2(26),bk2(27),bk2(28),bk2(29),bk2(30)
          II = II + 1
          end if
!          if(HR_k <= (k2+21)) then
!            write(10,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
!              fk2(1),fk2(2),fk2(3),fk2(4),fk2(5),  &
!              fk2(6),fk2(7),fk2(8),fk2(9),fk2(10), &
!              fk2(11),fk2(12),fk2(13),fk2(14),fk2(15), &
!              fk2(16),fk2(17),fk2(18),fk2(19),fk2(20), &
!              fk2(21),fk2(22),fk2(23),fk2(24),fk2(25), &
!              fk2(26),fk2(27),fk2(28),fk2(29),fk2(30), &
!              fk2(31),fk2(32),fk2(33),fk2(34),fk2(35), &
!              fk2(36),fk2(37),fk2(38)
!,fk2(39),fk2(40)
!              fk2(41),fk2(42),fk2(43),fk2(44),fk2(45), &
!              fk2(46),fk2(47),fk2(48),fk2(49),fk2(50)
!          endif
          if(HR_k <= energy_line_end_num) then
!            write(10,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
!            (fk2(ff_list_max_num),ff_list_max_num=1, new_ff_line_end_num)
!            write(10,'(i4,f13.8,i8)') HR_k,eigenvalue_array_sub_ef_eV(HR_k),new_ff_line_end_num
!            do ff_list_max_num=1, new_ff_line_end_num
!              write(10,'(f13.8)') fk2(ff_list_max_num)
!            end do
            write(10,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
              fk2(1),fk2(2),fk2(3),fk2(4),fk2(5),  &
              fk2(6),fk2(7),fk2(8),fk2(9),fk2(10), &
              fk2(11),fk2(12),fk2(13),fk2(14),fk2(15), &
              fk2(16),fk2(17),fk2(18),fk2(19),fk2(20), &
              fk2(21),fk2(22),fk2(23),fk2(24),fk2(25), &
              fk2(26),fk2(27),fk2(28),fk2(29),fk2(30)
          endif
!          if( (HR_k == 5) .or. (HR_k == 6) ) then
!            write(6,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k),bk2(1),fk2(1)
!            write(6,'(i4,3f13.8)') HR_k,eigenvalue_array_sub_ef_eV(HR_k),bk2(1),fk2(1)
!          end if
 889      format(i4,31f13.8)
! 889      format(i4,f13.8)
        end if
      end do
      new_ff_line_end_num_temp = new_ff_line_end_num
      if( new_ff_line_end_num_temp >= new_ff_line_end_num_min) then
        new_ff_line_end_num_min = new_ff_line_end_num_temp
      end if
! Inukai end
!      write(6,*) 'N(eig)=',energy_line_end_num,' N(ff)=',new_ff_line_end_num,' P=',fftt
!    end do
! HR calculation end
    !
! H. Sato program line end
!
!------------ ---------------------------------------------------------
!          if( ( (C2_1st_max - C2_2nd_max) * 100 / C2_1st_max ) >= C2_limit ) then
!            C2_max_at_k_point(HR_i,HR_j) = C2_max_at_k_point(HR_i,HR_j) + C2_1st_max
!            weight_at_C2_max(HR_i,HR_j)  = weight_at_C2_max(HR_i,HR_j) + 1.0
!          end if
!------------ ---------------------------------------------------------
!          total_K_max(HR_i,HR_j) = total_K_max(HR_i,HR_j) + SQRT(K2_at_C2_1st_max)
!          total_K_max_limited_K(HR_i,HR_j) = total_K_max_limited_K(HR_i,HR_j) + SQRT(K2_at_C2_max_limited_K)
!          count_weight(HR_i,HR_j) = count_weight(HR_i,HR_j) + 1.0
!          write( 6, * ) total_K_max(HR_i,HR_j),total_K_max_limited_K(HR_i,HR_j),count_weight(HR_i,HR_j)
!------------ ---------------------------------------------------------
          ! 3D plot
!          if ( (hr_3D_plot_switch == 1) .and. (print_switch == 1)) then
            !
!            if( (E_limit > (start_energy + energy_step*(HR_j - 1) - energy_step/2)) .and. &
!            (E_limit <= (start_energy + energy_step*(HR_j - 1) + energy_step/2)) ) then
!              K2_step_for_3D_plot = 0.05
!              num_total_K2_step_for_3D_plot = new_ff_line_end_num * int(1 / K2_step_for_3D_plot)
!              do HR_l = 1, new_ff_line_end_num
!                do HR_m = 1, num_total_K2_step_for_3D_plot
!                  if( wave_K2(HR_l) > (K2_step_for_3D_plot*(HR_m - 1) - K2_step_for_3D_plot/2) .and. &
!                    (wave_K2(HR_l) <= (K2_step_for_3D_plot*(HR_m - 1) + K2_step_for_3D_plot/2)) ) then
                    !
!                    total_C2_for_3D_plot(HR_m) = total_C2_for_3D_plot(HR_m) + wave_coef_2(HR_l,HR_k)
!                    weight_for_3D_plot(HR_m)   = weight_for_3D_plot(HR_m) + 1
!                  end if
!                end do
!              end do
!            end if
            !
!          end if
!------------ ---------------------------------------------------------
!          if ( print_switch == 1 ) then
!            if( (K2_distribution_at_E(1) > (start_energy + energy_step*(HR_j - 1) - energy_step/2)) .and. &
!            (K2_distribution_at_E(1) <= (start_energy + energy_step*(HR_j - 1) + energy_step/2)) ) then
!              write(30, '(3i10, f20.5)') HR_i, HR_k, HR_l, SQRT(K2_at_C2_1st_max)
!            end if
!          end if
!------------ ---------------------------------------------------------
!        end if
!      end do
!    end do
! HR calculation end
!------------ ---------------------------------------------------------
    deallocate ( wave_K2 )
    deallocate ( wave_coef_2 )
    deallocate ( eigenvalue_array_sub_ef_eV )
!
! H. Sato program line start
    deallocate ( bk2 )
    deallocate ( fk2 )
! H. Sato program line end
!
!------------ ---------------------------------------------------------
  end do
!--------Sato start
  open(16, file="f16")
    write(16,'(a4,i6)') "MMIN=",(energy_line_end_num + 10)
    write(16,'(a4,i6)') "II= ",II
    write(16,'(a4,i6)') "NTK ",end_of_k_point_line
  close(16)
!--------Sato end
  close(4)
  close(10)
  !
! read FF data and HR calculation end
!  write( 6, * ) 'AVERAGE |G|^2 CALCULATION'
!  do HR_i = 1, num_of_energy_step
!    do HR_j = 1, end_of_k_point_line
!      total_averaged_K_max(HR_i) = total_averaged_K_max(HR_i) + total_K_max(HR_j,HR_i) * weight(HR_j)
!      true_weight_at_k_point(HR_i) = true_weight_at_k_point(HR_i) + count_weight(HR_j,HR_i) * weight(HR_j)
!      true_weight_at_C2_max(HR_i) = true_weight_at_C2_max(HR_i) + weight_at_C2_max(HR_j,HR_i) * weight(HR_j)
!      total_ave_C2_max(HR_i) = total_ave_C2_max(HR_i) + C2_max_at_k_point(HR_j,HR_i) * weight(HR_j)
!      if ( print_switch == 1 ) then
!        tot_ave_K_max_limit(HR_i) = tot_ave_K_max_limit(HR_i) + total_K_max_limited_K(HR_j,HR_i) * weight(HR_j)
!      end if
!    end do
!  end do
  !
!------------ ---------------------------------------------------------
!  deallocate ( total_K_max )
!  deallocate ( weight )
!  deallocate ( count_weight )
!  deallocate ( weight_at_C2_max )
!  deallocate ( C2_max_at_k_point )
!  deallocate ( total_K_max_limited_K )
!------------ ---------------------------------------------------------
! old program
!  write( 6, * ) 'SIGMA2 CALCULATION (old definition)'
!  do HR_i = 1, num_of_energy_step
!    do HR_j = 1, end_of_k_point_line
!      if( (true_weight_at_k_point(HR_i) /= 0.0) .and. (count_weight(HR_j,HR_i) /= 0.0)) then

!        sigma2(HR_i) = sigma2(HR_i) + count_weight(HR_j,HR_i) * &
!weight(HR_j) * ( ( total_K_max(HR_j,HR_i) / count_weight(HR_j,HR_i) - &
!        ( total_averaged_K_max(HR_i) / true_weight_at_k_point(HR_i) ) ) ** 2 ) / true_weight_at_k_point(HR_i)
!      end if
!    end do
!  end do
!-----------------------------------------------------------------------------
  !
! --------write : start
!  open(4, file = 'f04')
!  write( 6, * ) 'SIGMA CALCULATION  AND  DATA WRITING'
!  if ( print_switch == 1 ) then
    ! --------SIGMA_MAX, SIGMA2_MAX, SIGMAP, SIGMA2P : start
!    sigma_max = 0.0
!    sigma2_max = 0.0
!    do HR_i = 1, num_of_energy_step
!      sigma = 0.0
!      temp_sigma2_max = sigma2(HR_i)
!      sigma = SQRT(temp_sigma2_max)
!      if( sigma >= sigma_max ) then
!        sigma_max = sigma
!      end if
!      if( sigma2(HR_i) >= sigma2_max ) then
!        sigma2_max = sigma2(HR_i)
!      end if
!    end do
! --------SIGMA_MAX, SIGMA2_MAX, SIGMAP, SIGMA2P : end
    !
!    name_of_column_chara(1)(1:16)= ' ENERGY_eV     ,'
!    name_of_column_chara(2)(1:16)= ' K             ,'
!    name_of_column_chara(3)(1:16)= ' K2            ,'
!    name_of_column_chara(4)(1:16)= ' SIGMA         ,'
!    name_of_column_chara(5)(1:16)= ' SIGMA2        ,'
!    name_of_column_chara(6)(1:16)= ' SIGMA_per     ,'
!    name_of_column_chara(7)(1:16)= ' SIGMA2_per    ,'
!    name_of_column_chara(8)(1:16)= ' N             ,'
!    name_of_column_chara(9)(1:16)= ' NW            ,'
!    name_of_column_chara(10)(1:16)=' CKmax*CKmax/NW,'
!    name_of_column_chara(11)(1:16)=' CKmax*CKmax   ,'
!    name_of_column_chara(12)(1:16)=' limited_K*K   ,'
!    write( 4, '(12(a16))' ) (name_of_column_chara(HR_i)(1:16), HR_i=1, 12)
    !
!    do HR_i = 1, num_of_energy_step
!      energy_for_hr_plot = start_energy + energy_step*(HR_i - 1)
      !
!      corrected_total_K_max = 0.0
!      if( true_weight_at_k_point(HR_i) /= 0.0 ) then
!        corrected_total_K_max = total_averaged_K_max(HR_i) / true_weight_at_k_point(HR_i) * sqrt(correct_volume_factor)
!      end if
      !
!      corrected_total_K2_max = 0.0
!      if( true_weight_at_k_point(HR_i) /= 0.0 ) then
!        corrected_total_K2_max = (total_averaged_K_max(HR_i) * total_averaged_K_max(HR_i)) / &
!        (true_weight_at_k_point(HR_i) * true_weight_at_k_point(HR_i)) * correct_volume_factor
!      end if
      !
!      temp_sigma2_max = sigma2(HR_i)
!      sigma = SQRT(temp_sigma2_max)
      !
!      sigma_percent = 0.0
!      if( sigma_max /= 0.0 ) then
!        sigma_percent = sigma / sigma_max * 100
!      end if
      !
!      sigma2_percent = 0.0
!      if( sigma2_max /= 0.0 ) then
!        sigma2_percent = sigma2(HR_i) / sigma2_max * 100
!      end if
      !
!      total_ave_C2_max_dw = 0.0
!      if( true_weight_at_C2_max(HR_i) /= 0.0 ) then
!        total_ave_C2_max_dw = total_ave_C2_max(HR_i) / true_weight_at_C2_max(HR_i)
!      end if
      !
!      corr_total_K_max_limi_K = 0.0
!      if( true_weight_at_k_point(HR_i) /= 0.0 ) then
!        corr_total_K_max_limi_K = tot_ave_K_max_limit(HR_i) / true_weight_at_k_point(HR_i) * sqrt(correct_volume_factor)
!      end if
!      corr_total_K2_max_limi_K = corr_total_K_max_limi_K * corr_total_K_max_limi_K
      !
!      N =  0.0
!      do HR_j = 1, end_of_k_point_line
!        N = N  + count_weight(HR_j,HR_i)
!      end do
      !
!      write( 4, 6200) energy_for_hr_plot, corrected_total_K_max, corrected_total_K2_max, &
!        sigma, temp_sigma2_max, sigma_percent, sigma2_percent, N, true_weight_at_k_point(HR_i), &
!        total_ave_C2_max_dw, total_ave_C2_max(HR_i), corr_total_K2_max_limi_K
! 6200 format( f12.7, ',', 11( f20.5, ',' ) )
      !
!    end do
    !
!  else
!    name_of_column_chara(1)(1:16)= ' ENERGY_eV     ,'
!    name_of_column_chara(2)(1:16)= ' k+G_2         ,'
!    name_of_column_chara(3)(1:16)= ' SIGMA2        ,'
!    name_of_column_chara(4)(1:16)= ' CKmax*CKmax/NW,'
!    write( 4, '(4(a16))' ) (name_of_column_chara(HR_i)(1:16), HR_i=1, 4)
    !
    ! ENERGY_eV, K_G_2, SIGMA2, CKmax_CKmax_NW
!    do HR_i = 1, num_of_energy_step
!      energy_for_hr_plot = start_energy + energy_step*(HR_i - 1)
      !
!      corrected_total_K2_max = 0.0
!      if( true_weight_at_k_point(HR_i) /= 0.0 ) then
!        corrected_total_K2_max = (total_averaged_K_max(HR_i) * total_averaged_K_max(HR_i)) / &
!        (true_weight_at_k_point(HR_i) * true_weight_at_k_point(HR_i)) * correct_volume_factor
!      end if
      !
!      temp_sigma2_max = sigma2(HR_i)
      !
!      total_ave_C2_max_dw = 0.0
!      if( true_weight_at_C2_max(HR_i) /= 0.0 ) then
!        total_ave_C2_max_dw = total_ave_C2_max(HR_i) / true_weight_at_C2_max(HR_i)
!      end if
      !
!      write( 4, 6250) energy_for_hr_plot, corrected_total_K2_max, temp_sigma2_max, total_ave_C2_max_dw
! 6250 format( f12.7, ',', 4( f20.5, ',' ) )
!    end do
!  end if
!  close(4)
  !
! f70 : read structure data from .struct(f70)
!  open(70, file = 'f70')
!  if ( hr_3D_plot_switch == 1 ) then
!    rewind(70)
!    write(70,*) '  |K|^2  ,     |C|^2     '
!    do HR_i=1, num_total_K2_step_for_3D_plot
!      if( weight_for_3D_plot(HR_i) /= 0.0 ) then
!        total_C2_for_3D_plot_dw(HR_i) = total_C2_for_3D_plot(HR_i) / weight_for_3D_plot(HR_i)
!      end if
!      num_K2_step = K2_step_for_3D_plot*(HR_i - 1)
!      write(70, 6300) num_K2_step, total_C2_for_3D_plot_dw(HR_i)
! 6300 format( f10.5, ',', f15.8 )
!    end do
!  end if
!  close(70)
! --------write : end
 !
  deallocate ( num_of_k_point_chara )
  deallocate ( C2_max_at_k_point )       ! C2_Kmax(end_of_k_point_line, num_of_energy_step)
  deallocate ( weight_at_C2_max )        ! CTWEGIA_Ckmax(end_of_k_point_line, num_of_energy_step)
  deallocate ( total_K_max )             ! KGA(end_of_k_point_line, num_of_energy_step)
  deallocate ( total_K_max_limited_K )   ! KGT_K_limit(end_of_k_point_line, num_of_energy_step)
  deallocate ( count_weight )            ! CTWEGIA(end_of_k_point_line, num_of_energy_step)
    deallocate ( total_C2_for_3D_plot )    ! K2_Elimit_T(num_total_K2_step_for_3D_plot)
    deallocate ( total_C2_for_3D_plot_dw ) ! K2_Elimit(num_total_K2_step_for_3D_plot)
    deallocate ( weight_for_3D_plot )      ! CTWEGIA_Elimit(num_total_K2_step_for_3D_plot)
  deallocate ( weight )                  ! WEIGHT(end_of_k_point_line)
!  deallocate ( K2_distribution_at_E )    ! DENERGY(dis_end_k_point_line)
!  deallocate ( total_averaged_K_max )    ! KGTA(num_of_energy_step)
!  deallocate ( true_weight_at_k_point )  ! CTWETA(num_of_energy_step)
!  deallocate ( true_weight_at_C2_max )   ! CTWETA_Ckmax(num_of_energy_step)
!  deallocate ( tot_ave_K_max_limit )     ! KGT_K_limit_T(num_of_energy_step)
!  deallocate ( total_ave_C2_max )        ! C2_Kmax_T(num_of_energy_step)
!  deallocate ( sigma2 )                  ! SIGMA2(num_of_energy_step)
  !
end program main
!**********************************************************************************************************************************
