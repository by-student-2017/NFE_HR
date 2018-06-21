!**********************************************************************************************************************************
!  calchr.f90 for bcc-Fe,V,Cr,W,Ta,Nb,Mo,Ba  20012.8.27
!    2012. 1.30  F90 Version 1.00  written by M.Inukai
!    2012. 3.10  H. Sato version written by M. Inukai
!    2013. 6.11  H. Sato version written by M. Inukai
!    2013. 7.19  tried to make monoclinic crystral version, but failed
!    2014. 2.10  Refined by M. Inukai
!
! Compile option for gfortran
! gfortran -o makethhrinput -Wall -pedantic -std=f95 -fbounds-check -O
! -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace
! -g makethhrinput.f90 -free -m64
!
! Compile option for ifort
! ifort -check all -warn all -std -fpe0 -traceback -g -o
! makethhrinput makethhrinput.f90
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
! file name   :                                               f99
! file type   :                                               out
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
  ! H. Sato program line start
  ! CNRYEV = 13.6058
  real(8),parameter :: convert_Ry2eV_cons = 13.60569805
  real(8),parameter :: PI=3.141592653
  ! H. Sato program line end
  !
!end module ff_cons
!**********************************************************************************************************************************
! HR plot
  ! 0:off, 1: on
  integer,parameter :: default_setting_switch = 1
  !
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
  integer temp_eigen_line_calc_num
  integer real_imaginary_switch
  integer eigen_value_num
  integer num_matrix
  integer matrix_num
  integer nlo_num
  integer read_nlo
  integer ff_start_line
  integer now_matrix_num
  integer cut_energy_line_num
  integer num_eigenvalue_array_Ry
  integer ef_read_type ! Sato Auto
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
  real(8) num_eigenvalue_sub_ef_Ry
  !
!end module real_set
!**********************************************************************************************************************************
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
  integer end_of_k_point_line    ! ENDOKL
  integer ix, iy, iz, idv
  integer HR_i, HR_j, HR_k, HR_l, HR_m
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
  !
  character(10), allocatable :: num_of_k_point_chara(:)
  real(8), allocatable :: C2_max_at_k_point(:,:)    ! C2_Kmax(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: weight_at_C2_max(:,:)     ! CTWEGIA_Ckmax(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: total_K_max(:,:)          ! KGA(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: total_K_max_limited_K(:,:)! KGT_K_limit(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: count_weight(:,:)         ! CTWEGIA(end_of_k_point_line, num_of_energy_step)
  real(8), allocatable :: weight(:)                 ! WEIGHT(end_of_k_point_line)
  !
  ! H. Sato program line start
  real(8), allocatable :: bk2(:)
  real(8), allocatable :: fk2(:)
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
  integer(4) kcount ! For SO calculation
  ! H. Sato program line end
  integer(4) IATO
  integer(4) NATO
  integer(4), allocatable ::  MIF(:)
  integer(4) MULT,ISPLIT
  integer(4) j_k, j_l, j_k_No, j_k_start_No, j_k_end_No
  real(8) end_l, start_k_maxC
  real(8) upper_range_eV
  integer(4) upper_range_num
  integer(4) HR_ke
!
!-------------- -------------------------------------------------------
!
!**********************************************************************************************************************************
! bag fix, integer
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
  energy_line_end_num = 0
  ff_start_line_num = 0
  ff_end_line_num = 0
  eigenvalue_start_line_num = 0
  eigenvalue_end_line_num = 0
  II = 0 ! Sato Auto
  new_ff_line_end_num_min = 0 ! Sato Auto
  new_ff_line_end_num_temp = 0 ! Sato Auto
  kcount = 0 ! SO
! bag fix, real
  EBOT = 0.0
!**********************************************************************************************************************************
  ! H. Sato program line start
  num_of_energy_step = 1
  print_switch = 0
  ! H. Sato program line end
  !
  ! H. Sato program line start
  open( 4, file = 'input_K2.dat' )
  open(10, file = 'input_C2.dat' )
  ! SG No.2
  ! AA = 5.03811976
  ! CC = (PI/AA)**2 * convert_Ry2eV_cons
  ! H. Sato program line end
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
      ! write(6,*) read_dump
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
      ! write(6,*) read_dump
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
      ! write(6,*) read_dump
      i=i+1
      if( MIF(i+1)==1 )then
        read(60,'(15x,i2,17x,i2)') MULT,ISPLIT
        write(6,'(a4,i3,a6,i2,a8,i2)') 'ATOM',j,',MULT=',MULT,',ISPLIT=',ISPLIT
        ! NATO = NATO + abs(MULT*ISPLIT)
        NATO = NATO + abs(MULT)
        i=i+1
        j=j+1
      end if
    end do
    deallocate(MIF)
    !
    rewind(60)
    read(60, '( a130 )') reador
    ! write(6, '( a45 )') reador(1:45)
    read(60, '( a130 )') reador
    ! write(6, '( a45 )') reador(1:45)
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
  close(60)
  i=0
  j=0
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
    ! write(6,*) read_dump
  end do
  if ( ef_read_type == 0 ) then
    read(15, '(18x,E27.16)') ef_energy_Ry
  endif
  if ( ef_read_type == 1 ) then
    read(15, '(38x,E25.15)') ef_energy_Ry
  end if
  write(6, * ) 'Auto Read Fermi Energy: EF=', ef_energy_Ry, ' (Ry)'
  !
  close (15)
! -----------------------
! f60 : read structure data from .struct(f60)
  open(60, file = 'f60')
  rewind(60)
  !
  read(60, '( a130 )') reador
    write(6, '( a45 )') reador(1:45)
  read(60, '( a130 )') reador

  write(6, '( a45 )') reador(1:45)
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
  ! change theta to radian
  write(6,*) 'This code can calculate FCC, BCC, HCP and SC only'
  angle_a_radian = angle_a*PI/180.0
  angle_b_radian = angle_b*PI/180.0
  angle_c_radian = angle_c*PI/180.0
  ! write(6,*) sin(angle_a_radian), sin(angle_b_radian), sin(angle_c_radian)
  ! write(6,*) cos(angle_a_radian), cos(angle_b_radian), cos(angle_c_radian)
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
! -----------------------
! f17 : read G1, G2, G3
  open (17, file = 'f17')
  i = 0
  do while(i==0)
    read (17, '(a26)') reador
    if (reador(1:26)  == '    G1        G2        G3') then
      read(17, '(1x,3f10.6)') BG1X, BG2X, BG3X
      read(17, '(1x,3f10.6)') BG1Y, BG2Y, BG3Y
      read(17, '(1x,3f10.6)') BG1Z, BG2Z, BG3Z
      i = 1
    endif
  end do
  close (17)
  write(6,*)  'G1        G2        G3'
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
  write( 6,*) '-----------------------'
!
! -----------------------
! f66 : read k point data from upper_range(f66)
  open (66, file = 'f66')
  read (66, '(f10.5)') upper_range_eV
  write( 6,*) '-----------------------'
  write( 6,*) 'upper range: ', upper_range_eV
  write( 6,*) '-----------------------'
  close (66)
!**********************************************************************************************************************************
!  call read_line_and_data(hr_array_number)
!
  allocate ( C2_max_at_k_point(end_of_k_point_line,num_of_energy_step) )    ! C2_Kmax(end_of_k_point_line, num_of_energy_step)
  allocate ( weight_at_C2_max(end_of_k_point_line,num_of_energy_step) )     ! CTWEGIA_Ckmax(end_of_k_point_line, num_of_energy_step)
  allocate ( total_K_max(end_of_k_point_line,num_of_energy_step) )          ! KGA(end_of_k_point_line, num_of_energy_step)
  allocate ( total_K_max_limited_K(end_of_k_point_line,num_of_energy_step) )! KGT_K_limit(end_of_k_point_line, num_of_energy_step)
  allocate ( count_weight(end_of_k_point_line,num_of_energy_step) )         ! CTWEGIA(end_of_k_point_line, num_of_energy_step)
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
  do HR_i = 1, end_of_k_point_line
    now_num_of_k_point = HR_i
!**********************************************************************************************************************************
! call read_line_and_data(ff_array_number)
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
  open(99, file = 'nlo' )
  do while(1==1)
! version dependence -------------------------------
    if (reador(1:28) == 'Time for iouter   (hamilt) :') then
      read(26, '( 40x, i9 )') read_nlo
      write(99, *) read_nlo
      if (now_num_of_k_point == 1) then
        write(6, * ) 'WIEN2k version 9'
        pre_line_nlo_chara = 'Time for iouter   (hamilt) :'
      endif
    endif
! version dependence -------------------------------
    if (reador(1:28) == 'Time for distrib  (hamilt, c') then
      read(26, '( 40x, i9 )') read_nlo
      if (now_num_of_k_point == 1) then
        write (6, * ) 'WIEN2k version 10, 11, 12'
        write (6, * ) 'Number of local orbitals =', read_nlo
        write(99, *) read_nlo
      endif
      pre_line_nlo_chara = 'Time for distrib  (hamilt, c'
    endif
! version dependence -------------------------------
    read(26, '( a130 )') reador
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
  read(26, '( 17x, i6 )') matrix_num
  write (6, * ) 'MATRIX SIZE = ', matrix_num
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
    if (reador(1:26) == '       NUMBER OF K-POINTS:') then
       rewind(99)
       read (99, *) read_nlo
      exit
    endif
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
      ! write(6, *) kcount
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
        !
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
        !
      endif
! 2
      if (reador(28:38) == '  *********') then
        reador(28:38)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  28, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
      endif
! 3
      if (reador(39:49) == '  *********') then
        reador(39:49)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  39, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
      endif
! 4
      if (reador(50:60) == '  *********') then
        reador(50:60)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  50, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
      endif
! 5
      if (reador(61:71) == '  *********') then
        reador(61:71)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write (80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write (80,  * ) '********* existed at line ', i, ', row  61, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
      endif
! 6
      if (reador(72:82) == '  *********') then
        reador(72:82)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  72, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
      endif
! 7
      if (reador(83:93) == '  *********') then
        reador(83:93)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  83, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
      endif
! 8
      if (reador (94:104) == '  *********') then
        reador(94:104)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row  94, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
      endif
! 9
      if (reador(105:115) == '  *********') then
        reador(105:115)  = '   0.000000'
        !
        if (check_file_option >= 2) then
          write(80, * ) 'NMATRIX NUMBER ', now_matrix_num
          write(80,  * ) '********* existed at line ', i, ', row 105, ', reador(118:130)
          if (now_matrix_num <= num_of_miller_index) then
            write(80,  * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
            write(6, * ) warning_chara, '( Line ', now_matrix_num, ' in ', num_of_miller_index, ')'
          endif
          write(80,  * ) '--------------------------------------------------------------'
        endif
        !
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
  ! write( 6, *) ' # k point, # Energy, # APW    , Max|2(G+k)|'
  !
  open (1, file = 'f01')
  rewind (1)
  !
  read (1, '( 7x, 3f10.5, 3x, 10a  )') miller_index_h, miller_index_k, miller_index_l, read_now_k_point_chara
  !read (1, '( 17x, i5 )') num_matrix
  read (1, '( 17x, i6 )') num_matrix
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
    if ((dum (1:28) == pre_line_nlo_chara(1:28) ) .or. (dum (1:26) == '       NUMBER OF K-POINTS:')) then
       rewind(99)
       read (99, *) read_nlo
       nlo_num = read_nlo
      exit
    endif
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
  energy_array_number = NINT(true_energy_line_num) * 9
  ff_array_number     = ff_array_number
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
energy_array_number = energy_array_number + 30
  allocate ( eigenvalue_array_Ry(energy_array_number) )                    ! eigenv(ENNM)
energy_array_number = energy_array_number - 30
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
  !
energy_array_number = energy_array_number + 30
ff_array_number = ff_array_number + 30
!**********************************************************************************************************************************
  ! -----------------------
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
  allocate ( shelter_array(energy_array_number) )                          ! shelt4(ENNM)
  allocate ( eigenvalue_array_eV(energy_array_number) )                    ! eigeve(ENNM)
  allocate ( eigenvalue_array_sub_ef_eV(energy_array_number) )             ! neigev(ENNM)
  !
  ! H. Sato program line start
  allocate ( bk2(ff_array_number) )
  allocate ( fk2(ff_array_number) )
  ! H. Sato program line end
  !
!**********************************************************************************************************************************
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
  !
energy_array_number = energy_array_number - 30
ff_array_number = ff_array_number - 30
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write (6, * ) 'initialize array end'
  endif
! check end ----------------------------------------------------------
!**********************************************************************************************************************************
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
            !write(6, '( 3i4 )') (wave_num(j) , j = 1, 3)
! K calculate --------------------------------------------------------
!**********************************************************************************************************************************
! check start --------------------------------------------------------
  if (check_file_option == 1) then
    write(6, * ) 'G**2 calculation start'
  endif
! check end ----------------------------------------------------------
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
    WVX=BXX*( float( wave_num(1) ) + miller_index_h ) * (2.0 * axis_n / axis_a) &
       +BXY*( float( wave_num(2) ) + miller_index_k ) * (2.0 * axis_n / axis_b) &
       +BXZ*( float( wave_num(3) ) + miller_index_l ) * (2.0 * axis_n / axis_c)
    !
    WVY=BYX*( float( wave_num(1) ) + miller_index_h ) * (2.0 * axis_n / axis_a) &
       +BYY*( float( wave_num(2) ) + miller_index_k ) * (2.0 * axis_n / axis_b) &
       +BYZ*( float( wave_num(3) ) + miller_index_l ) * (2.0 * axis_n / axis_c)
    !
    WVZ=BZX*( float( wave_num(1) ) + miller_index_h ) * (2.0 * axis_n / axis_a) &
       +BZY*( float( wave_num(2) ) + miller_index_k ) * (2.0 * axis_n / axis_b) &
       +BZZ*( float( wave_num(3) ) + miller_index_l ) * (2.0 * axis_n / axis_c)
    !
    wave_K2(ff_count) = WVX**2 + WVY**2 + WVZ**2
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
            !WRITE(6, '( 16x, 9f11.6 )') (wave_coef_real(ff_count, ( (ff_line_set_num - 1) * 9 + j) ) , j = 1, 9)
! read ff coefficients at imaginary part -----------------------------
            if (read_data_line == 3) then
              read(1, '( 16x, 9f11.6 )') (wave_coef_imag(ff_count, ( (ff_line_set_num - 1) * 9 + j) ) , j = 1, 9)
              !WRITE(6, '( 16x, 9f11.6 )') (wave_coef_imag(ff_count, ( (ff_line_set_num - 1) * 9 + j) ) , j = 1, 9)
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
    ! WRITE(6, *) i, end_line_num, j, ff_count,ff_array_number, &
    ! wave_K2(ff_count), wave_coef_real(ff_count,10), now_search_ff_line_num
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
  energy_line_end_num = energy_array_number
  do i = 1, energy_array_number
    eigenvalue_array_Ry(i) = eigenvalue_array_Ry(i + (cut_energy_line_num * 9))
    if ((energy_num_chara(i) == empty_chara) .or. (energy_num_chara(i) == 'XXXXXXXXXXX')) then
      energy_line_end_num = i - 1
      exit
    endif
    ! write(6, * ) energy_num_chara(i), eigenvalue_array_Ry(i)
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
!
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
    ! write(6,*) i, ff_array_number, wave_K2(i)
  end do
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
  !
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
      ! write(6, *) wave_K2(new_pw_data_num)
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
!end eV_calculate
!**************************************************************************************************
!
!**************************************************************************************************
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
    write(6, * ) 'write_coefficients_vs_K end', HR_i
  endif
! check end ----------------------------------------------------------
!**********************************************************************************************************************************
  ! -----------------------
  deallocate ( eigenvalue_array_Ry )                                   ! eigenv(ENNM)
  deallocate ( energy_num_chara )                                      ! energy(ENNM)
  deallocate ( wave_coef_real )                                        ! wavecr(FFNM, ENNM)
  deallocate ( wave_coef_imag )                                        ! waveci(FFNM, ENNM)
  deallocate ( wave_coef_2_bond )                                      ! wavec2a(FFNM, ENNM)
  deallocate ( wave_coef_2_anti_bond )                                 ! wavec2b(FFNM, ENNM)
  deallocate ( shelter_array )                                      ! shelt4(ENNM)
  deallocate ( shelter_array_bond )                                 ! shelt4b(ENNM)
  deallocate ( shelter_array_anti_bond )                            ! shelt4a(ENNM)
  deallocate ( eigenvalue_array_eV )                                ! eigeve(ENNM)
!**********************************************************************************************************************************
!-------------- -------------------------------------------------------
      !write(6, * ) HR_i, HR_k, energy_line_end_num
      do HR_k = 1, energy_line_end_num
        if (check_file_option >= 1) then
          write(6, * ) 'write_NFE_DATA', HR_i
        endif
        !write(6, * ) HR_i, HR_k
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
        if( (HR_i == 1) .and. (HR_k == 1) )then
          write(6,*) 'eigenvalue list (eV) at gamma point'
          !write(6,*) (j,':', eigenvalue_array_sub_ef_eV(j), j=1, energy_line_end_num)
          ! find data above -25 eV 
          j_k_start_No = 0
          do j_k = 1, energy_line_end_num
            if ( eigenvalue_array_sub_ef_eV( j_k ) <= -25 ) then
              j_k_start_No = j_k
            end if
          end do
          j_k_start_No = j_k_start_No + 1
          ! find data below 1 eV
          j_k_end_No = 0
          do j_k = 1, energy_line_end_num
            if ( eigenvalue_array_sub_ef_eV( j_k ) <= 1 ) then
              j_k_end_No = j_k
            end if
          end do
          j_k_end_No = j_k_end_No + 1
          !
          !if (energy_line_end_num <= 800) then
            if(j_k_start_No >= 3)then
              write(6,*) (j,':', eigenvalue_array_sub_ef_eV(j), j=(j_k_start_No-2), j_k_end_No)
            else
              write(6,*) (j,':', eigenvalue_array_sub_ef_eV(j), j=j_k_start_No, j_k_end_No)
            end if
          !else
          !  write(6,*) (j,':', eigenvalue_array_sub_ef_eV(j), j=j_k_start_No, energy_line_end_num)
          !end if
          write(6,*) 'input eigenvlue number at bottom in Valence band'
          ! auto
          end_l = 1.0
          start_k_maxC = 0.0
          j_k_No = 0
          do j_l = 1, new_ff_line_end_num 
            if ( (wave_K2( j_l ) < end_l) .and. (wave_K2( j_l ) >= 0.0) ) then
              end_l = wave_K2( j_l )
              do j_k = j_k_start_No, energy_line_end_num
                if ( wave_coef_2( j_l, j_k ) >= start_k_maxC ) then
                  start_k_maxC = wave_coef_2( j_l, j_k )
                  EBOT = eigenvalue_array_sub_ef_eV( j_k )
                  j_k_No = j_k
                end if
              end do
            end if
          end do
          ! auto
          write(6,*) 'program recommend the eignvalue'
          write(6,*) ' No. and eV at gamma point, No.', j_k_No, ':', EBOT, "eV"
          k2 = j_k_No
          write(6,*) 'Please, input start eigenvalue No., (0:auto)'
          read(5,*) k2
          if( k2 == 0 ) then
            k2 = j_k_No
            write(6,*) "use auto start eigenvalue No.", k2
          end if
          EBOT = eigenvalue_array_sub_ef_eV(k2)
          !
          do HR_ke = 1, energy_line_end_num
            if( eigenvalue_array_sub_ef_eV(energy_line_end_num) < upper_range_eV ) then
              upper_range_num = HR_ke
            end if
            upper_range_num = upper_range_num
          end do
          write(6,*) "use end eigenvalue No.", upper_range_num
          write(6,*) 'calc. range:', EBOT, '-', upper_range_eV, 'eV'
        end if
        low_num_eigenvale = energy_line_end_num - k2 + 1
        ! H. Sato program line end
        ! Inukai end
        !
        ! H. Sato program line start
        do HR_l = 1, new_ff_line_end_num
          bk2(HR_l) = wave_K2(HR_l)
          fk2(HR_l) = wave_coef_2(HR_l,HR_k)
          if( wave_coef_2(HR_l,HR_k) > C2_1st_max ) then
            K2_at_C2_1st_max = wave_K2(HR_l)
            C2_1st_max       = wave_coef_2(HR_l,HR_k)
          end if
        ! H. Sato program line end
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
          !if(HR_k <= energy_line_end_num) then
          if(HR_k <= upper_range_num) then
          ! write(4,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
          ! (bk2(ff_list_max_num),ff_list_max_num=1, new_ff_line_end_num)
          ! write(4,'(i4,f13.8,i8)') HR_k,eigenvalue_array_sub_ef_eV(HR_k),new_ff_line_end_num
          ! do ff_list_max_num=1, new_ff_line_end_num
          !   write(4, '(f13.8)') bk2(ff_list_max_num)
          ! end do
            write(4,'(i4,31f13.8)') HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
              bk2(1),bk2(2),bk2(3),bk2(4),bk2(5),  &
              bk2(6),bk2(7),bk2(8),bk2(9),bk2(10), &
              bk2(11),bk2(12),bk2(13),bk2(14),bk2(15), &
              bk2(16),bk2(17),bk2(18),bk2(19),bk2(20), &
              bk2(21),bk2(22),bk2(23),bk2(24),bk2(25), &
              bk2(26),bk2(27),bk2(28),bk2(29),bk2(30)
          II = II + 1
          end if
          !
          !if(HR_k <= energy_line_end_num) then
          if(HR_k <= upper_range_num) then
          ! write(10,889) HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
          ! (fk2(ff_list_max_num),ff_list_max_num=1, new_ff_line_end_num)
          ! write(10,'(i4,f13.8,i8)') HR_k,eigenvalue_array_sub_ef_eV(HR_k),new_ff_line_end_num
          ! do ff_list_max_num=1, new_ff_line_end_num
          !   write(10,'(f13.8)') fk2(ff_list_max_num)
          ! end do
            write(10,'(i4,31f13.8)') HR_k,eigenvalue_array_sub_ef_eV(HR_k), &
              fk2(1),fk2(2),fk2(3),fk2(4),fk2(5),  &
              fk2(6),fk2(7),fk2(8),fk2(9),fk2(10), &
              fk2(11),fk2(12),fk2(13),fk2(14),fk2(15), &
              fk2(16),fk2(17),fk2(18),fk2(19),fk2(20), &
              fk2(21),fk2(22),fk2(23),fk2(24),fk2(25), &
              fk2(26),fk2(27),fk2(28),fk2(29),fk2(30)
          endif
          !
        end if
      end do
      new_ff_line_end_num_temp = new_ff_line_end_num
      if( new_ff_line_end_num_temp >= new_ff_line_end_num_min) then
        new_ff_line_end_num_min = new_ff_line_end_num_temp
      end if
      ! Inukai end
    ! HR calculation end
    !
    ! H. Sato program line end
!
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
    !write(16,'(a4,i6)') "MMIN=",(energy_line_end_num + 10)
    write(16,'(a4,i6)') "MMIN=",(upper_range_num - k2 + 1)
    write(16,'(a4,i6)') "II= ",II
    write(16,'(a4,i6)') "NTK ",end_of_k_point_line
  close(16)
  !--------Sato end
  close(4)
  close(10)
  close(99)
  !
  deallocate ( num_of_k_point_chara )
  deallocate ( C2_max_at_k_point )       ! C2_Kmax(end_of_k_point_line, num_of_energy_step)
  deallocate ( weight_at_C2_max )        ! CTWEGIA_Ckmax(end_of_k_point_line, num_of_energy_step)
  deallocate ( total_K_max )             ! KGA(end_of_k_point_line, num_of_energy_step)
  deallocate ( total_K_max_limited_K )   ! KGT_K_limit(end_of_k_point_line, num_of_energy_step)
  deallocate ( count_weight )            ! CTWEGIA(end_of_k_point_line, num_of_energy_step)
  deallocate ( weight )                  ! WEIGHT(end_of_k_point_line)
  !
end program main
!**********************************************************************************************************************************
