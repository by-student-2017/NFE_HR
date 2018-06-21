!**********************************************************************************************************************************
!  calcff.f90
!    2012. 1.26  F90 Version 1.00  written by M.Inukai
!**********************************************************************************************************************************
!**********************************************************************************************************************************
! file data ----------------------------------------------------------  
! file number :   1    2    3    4     5     6    7     8     9     10  
! file name   :  f01      f03         key  dis  f07   f08   f09    f10  
! file type   :  inp      inp         inp  out  out   out   out     in  
! file data ----------------------------------------------------------  
! file data ----------------------------------------------------------  
! file number :  11   12   13   14    15    16    17    18    19    20  
! file name   :                                                         
! file type   :                                                         
! file data ----------------------------------------------------------  
! file data ----------------------------------------------------------  
! file number :  21   22   23   24    25    26    27    28    29    30  
! file name   :                                                           
! file type   :                                                         
! file data ----------------------------------------------------------  
! file number :  41   42   43   44    45    46    47    48    49    50  
! file name   :                                  f47   f48   f49        
! file type   :                                  out   out   out        
! file data ----------------------------------------------------------  
! file number :  51   52   53   54    55    56    57    58    59    60  
! file name   :                                                         
! file type   :                                                         
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
  ! CNRYEV
  real(8),parameter :: convert_Ry2eV_cons = 13.6058 
  real(8),parameter :: PI=3.141592653
  !
  ! cut eigenvalute at threshold, NEED cut_Ry_threshold_cons < 0
  real(8),parameter :: cut_Ry_threshold_cons = -1.5
  !
  real(8) wave_K2_diff
!end module ff_cons
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
  character(10), allocatable :: k_point_chara_set(:) ! KPOINT(k_point_num)
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
  character(10) input_k_point, read_k_point_chara, read_now_k_point_chara
! read EF automatically
  character(43) read_dump
! syesno is select yes or no
  character(1) axis_select
!end module character_set
  character(30) filename
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
  integer num_of_k_point
  integer hcp_on
  integer ff_start_line
  integer now_matrix_num
  integer cut_energy_line_num
  integer num_eigenvalue_array_Ry
  integer while_infinite_loop
  !integer start_K2_E
  !integer end_K2_E
  integer ef_read_type
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
  real(8) threshold_energy
  real(8) sigma_C2
  real(8) sigma_CG
  real(8) CG
  real(8) correct_volume_factor, correct_volume_factor_n, COK
  real(8) calculation_start_eV, calculation_end_eV, calculation_C2_threshold, times
!end module real_set
!**********************************************************************************************************************************
!**********************************************************************************************************************************
!  call read_cons_data(k_point_chara_set)
!
! bag fix
ff_start_line = 0
while_infinite_loop = 1
wave_K2_diff = 0.0
  !
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
  open (10, file = 'f10')
  !
  i = 0 
  do while( while_infinite_loop == 1 )
    read(10, '(A43)') read_dump
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
  rewind(10)
  do j=1, i-1
    read(10, '(A43)') read_dump
!    write(6,*) read_dump
  end do
  if ( ef_read_type == 0 ) then
    read(10, '(18x,E25.15)') ef_energy_Ry
  endif
  if ( ef_read_type == 1 ) then
    read(10, '(38x,E25.5)') ef_energy_Ry
  end if
  write(6, * ) 'Auto Read Fermi Energy: EF=', ef_energy_Ry, ' (Ry)' 
  !
  close (10)
!
! -----------------------
! K2 plot
  open (20, file = 'f20')
  !
  read(20, '(A10)') read_dump
  read(20, '(f10.5)') threshold_energy
  write(6, * ) 'ABS(K2 plot) <= ', threshold_energy, ' (eV)' 
  !
  close (20)
!
! -----------------------
! f60 : read structure data from .struct(f60)
  open(60, file = 'f60')
  !
  read(60, '( a130 )') reador 
  write(6, '( a45 )') reador(1:45) 
  read(60, '( a130 )') reador 
  write(6, '( a45 )') reador(1:45) 
  hcp_on = 0 
  if (reador(1:1) == 'H') then 
    hcp_on = 1 
  endif
  read(60, '( a130 )') reador 
  read(60, '( 6f10.6 )') axis_a, axis_b, axis_c, angle_a, angle_b, angle_c
  !
  write(6, * ) ' axis_a  ,  axis_b  ,  axis_c  , angle_a  , angle_b  , angle_c  '
  !
  axis_a = axis_a * 0.529177249
  axis_b = axis_b * 0.529177249
  axis_c = axis_c * 0.529177249
  !
  write(6, '( 6f10.6 )') axis_a, axis_b, axis_c, angle_a, angle_b, angle_c
  !
  close (60) 
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
!
! Normalized by axis
  write(6, * ) 'Normalized by axis : a, b, c or n(=No normalization)? : default = a axis'
  read(5, '( a1 )') axis_select
  if (axis_select == 'a') then
    axis_n = axis_a
    write(6, * ) 'Normalized by ', axis_select, ' axis'
  elseif(axis_select == 'b') then 
    axis_n = axis_b
    write(6, * ) 'Normalized by ', axis_select, ' axis'
  elseif(axis_select == 'c') then 
    axis_n = axis_c
    write(6, * ) 'Normalized by ', axis_select, ' axis'
  elseif(axis_select == 'n') then 
    write(6, * ) 'No normalization'
    hcp_on = -1
    COK = 1
  else
    axis_n = axis_a
    write(6, * ) 'Normalized by a axis'
  endif
  if(hcp_on /= -1) then
    correct_volume_factor_n = axis_n * axis_n * axis_n
    COK = correct_volume_factor / correct_volume_factor_n
    COK = COK**(2.0/3.0)
  end if
  write(6,*) "volume correct = ", COK, "times"
!
!
! check start --------------------------------------------------------  
  if (check_file_option >= 1) then 
    write(6, * ) 'read line start' 
  endif
! check end ----------------------------------------------------------  
!
! f03 : read k point data from .klist(f03)
  open (3, file = 'f03')
  !
  i = 0 
  do while( while_infinite_loop == 1 )
    read (3, '(A10)') read_k_point_chara 
    if (read_k_point_chara == 'END       ') exit
    if (read_k_point_chara /= '          ') then
      i = i + 1
    endif
  end do
  !
  allocate( k_point_chara_set(i) )
  !
  rewind (3)
  i = 0 
  do while( while_infinite_loop == 1 )
    read (3, '(A10)') read_k_point_chara 
    if (read_k_point_chara == 'END       ') exit
    if (read_k_point_chara /= '          ') then
      i = i + 1
      k_point_chara_set(i) = read_k_point_chara
    endif
  end do
  !
  write(6, * )
  write(6, * ) 'Please, input k-point, e.g., '
  write(6, '(5a10)') (k_point_chara_set(j) , j = 1, i)
  write(6, * ) ', then push enter.'
  read(5, '( a10 )') input_k_point
  write(6, * )
  !
  deallocate( k_point_chara_set )
  !
  close (3)
!
!
! f27 : read G vs coefficients data from .gcene(f27)
!  open(27, file = 'f27')
!  read(27, '( a130 )') reador
!  read(27, '( a130 )') reador
!  read(27, '(2f10.5)') energy_eV_min, energy_eV_max
!  write(6, * ) 'K vs coefficients tabel calculation range' 
!  write(6, * ) 'energy_eV_min=', energy_eV_min
!  write(6, * ) 'energy_eV_max=', energy_eV_max
!  write(6, * )
!  close (27)
!**********************************************************************************************************************************
!  call correct_data()
!  integer ff_start_line, now_matrix_num
  !
! f01 : read G vs coefficients data from cutting file(cut 26)
  open (1, file = 'f01')
  rewind (1)
! f26 : read G vs coefficients data from .output1(f26)
  open (26, file = 'f26')
  rewind (26)
  !
! check start --------------------------------------------------------
  if (check_file_option >= 1) then
    write (6, * ) 'read line start'
  endif
! check end ----------------------------------------------------------
  !
  num_of_k_point = 1
  do while( while_infinite_loop == 1 )
! version dependence -------------------------------
    if (reador(1:28) == 'Time for iouter   (hamilt) :') then
      read(26, '( 40x, i9 )') read_nlo
      if (num_of_k_point == 1) then
        write(6, * ) 'WIEN2k version 9'
        num_of_k_point = 0
      endif
      pre_line_nlo_chara = 'Time for iouter   (hamilt) :'
    endif

! version dependence -------------------------------
    if (reador(1:28) == 'Time for distrib  (hamilt, c') then
      read(26, '( 40x, i9 )') read_nlo
      if (num_of_k_point == 1) then
        write (6, * ) 'WIEN2k version 10 or later'
        num_of_k_point = 0
      endif
      pre_line_nlo_chara = 'Time for distrib  (hamilt, c'
    endif
! version dependence -------------------------------
    read(26, '( a130 )') reador
    if ( (reador(1:7) == read_k_pre_chara) .and. (reador(41:50) == input_k_point) ) exit
  end do
  !
  read(26, '( 17x, i5 )') matrix_num
  !
  read_data_line = 2
  !
  do while( while_infinite_loop == 1 )
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
  do while( while_infinite_loop == 1 )
    read(26, '( a130 )') reador                                   
    if ( (reador(1:7) == read_k_pre_chara) .and. (reador(41:50) == input_k_point) ) exit
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
  do while( while_infinite_loop == 1 )
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
!**********************************************************************************************************************************
!  call read_line_and_data(ff_array_number)
! 
! f01 : read G vs coefficients data from cutting file(cut 26)
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
          write(6, * ) 'INCLUDE IMAGINARY PART CALCULATION' 
          read_data_line = 3 
        else
          write(6, * ) 'ONLY REAL PART CALCULATION' 
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
! check end ----------------------------------------------------------
!
!
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
  ! -----------------------
  !
!**********************************************************************************************************************************
!  call initialize_data(energy_array_number,energy_num_chara, shelter_array, eigenvalue_array_Ry, &
!  eigenvalue_array_eV, eigenvalue_array_sub_ef_eV, wave_num, ff_array_number, wave_K2, &
!  wave_coef_real, wave_coef_imag, wave_coef_2, wave_coef_2_bond, wave_coef_2_anti_bond)
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
  if (ff_line_set_num == 1) then
    !
    if (hcp_on == 0) then 
      wave_K2(ff_count) = &
      ( (float(wave_num(1) ) + miller_index_h) * (2.0 * axis_n / axis_a) ) **2 + &
      ( (float(wave_num(2) ) + miller_index_k) * (2.0 * axis_n / axis_b) ) **2 + &
      ( (float(wave_num(3) ) + miller_index_l) * (2.0 * axis_n / axis_c) ) **2
    endif
    if (hcp_on == -1) then 
      wave_K2(ff_count) = &
      ( (float(wave_num(1) ) + miller_index_h) * 2.0) **2 + &
      ( (float(wave_num(2) ) + miller_index_k) * 2.0) **2 + &
      ( (float(wave_num(3) ) + miller_index_l) * 2.0) **2
    endif
! ----normalized for hcp
    if (hcp_on == 1) then 
      wave_K2(ff_count) = &
      ( ( (float(wave_num(1) ) + miller_index_h) * 2.0 + (float(wave_num(2) ) + miller_index_k) ) * &
      (1.0 / sqrt(3.0)) * (2.0 * axis_n / axis_a) ) **2 + &
      ( ( (float(wave_num(2) ) + miller_index_k) ) * (2.0 * axis_n / axis_b) ) **2 + &
      ( ( (float(wave_num(3) ) + miller_index_l) ) * (2.0 * axis_n / axis_c) ) **2
    endif
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
!    WRITE(6, *) i, end_line_num, j, ff_count,ff_array_number, wave_K2(ff_count), &
!  wave_coef_real(ff_count,10), now_search_ff_line_num
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
!    if (wave_K2(i) == wave_K2(i + 1) ) then 
    wave_K2_diff = wave_K2(i) - wave_K2(i + 1)
    if ( abs(wave_K2_diff) < 0.001) then 
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
!end eV_calculate
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
!  call write_K_vs_coefficients(eigenvalue_array_sub_ef_eV, degenerate_pw, &
!  wave_K2, wave_coef_2, ff_array_number)
!
! check start --------------------------------------------------------  
  if (check_file_option >= 1) then 
    write(6, * ) 'write_K_vs_coefficients start' 
  endif
! check end ----------------------------------------------------------  
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
  ! real calculation range
  open (39, file = "calcff_start_end_eV" )
    read(39, '(a)') read_dump
    read(39, '(9x,E25.15)') calculation_start_eV
    read(39, '(9x,E25.15)') calculation_end_eV
    read(39, '(9x,E25.15)') calculation_C2_threshold
    read(39, '(9x,E25.15)') times
    write(6,*)
    write(6,*) read_dump
    write(6,*) "Start eV:", calculation_start_eV
    write(6,*) "End   eV:", calculation_end_eV
    write(6,*) "|C|^2 >=:", calculation_C2_threshold
    write(6,*) "times   :", times
  close (39)
  !
  filename(1:20)  = 'calcff_energy_vs_C2_'
  filename(21:26) = input_k_point
  filename(27:30) = '.txt'
  open (9, file = filename )
  !
  write(9,*) "eV_",input_k_point,",","C2_",input_k_point
  do j = 1, new_ff_line_end_num
    do i = 1, energy_line_end_num
      if( (eigenvalue_array_sub_ef_eV(i) >= calculation_start_eV) .and. (eigenvalue_array_sub_ef_eV(i) <= calculation_end_eV) )then
        if( wave_coef_2(j, i) > calculation_C2_threshold ) then
          write(9,'(f14.8,",",a4)') eigenvalue_array_sub_ef_eV(i),"NaN"
          write(9,'(f14.8,",",f14.8)') eigenvalue_array_sub_ef_eV(i), wave_K2(j)*COK
          write(9,'(f14.8,",",f14.8)') eigenvalue_array_sub_ef_eV(i), wave_K2(j)*COK + wave_coef_2(j, i)*times
        end if
      end if
    end do
  end do
  close (9)
  !
  filename(1:20)  = 'calcff_energy_vs_CG_'
  filename(21:26) = input_k_point
  filename(27:30) = '.txt'
  open (29, file = filename )
  !
  write(29,*) "energy_",input_k_point,",","CG_",input_k_point                                      
  do j = 1, new_ff_line_end_num
    sigma_C2 = 0.0
    sigma_CG = 0.0
    CG = 0.0
    do i = 1, energy_line_end_num
      if( (eigenvalue_array_sub_ef_eV(i) >= calculation_start_eV) .and. (eigenvalue_array_sub_ef_eV(i) <= calculation_end_eV) )then
        if( wave_coef_2(j, i) > calculation_C2_threshold ) then
          sigma_C2 = sigma_C2 + wave_coef_2(j, i)
          sigma_CG = sigma_CG + eigenvalue_array_sub_ef_eV(i)*wave_coef_2(j, i)
        end if
      end if
    end do
    CG = sigma_CG / sigma_C2
    write(29,'(f14.8,",",f14.8)') CG, wave_K2(j)*COK
  end do
  close (29)
  !
! check start --------------------------------------------------------  
  if (check_file_option >= 1) then 
    write(6, * ) 'write_coefficients_vs_K end' 
  endif
! check end ----------------------------------------------------------  
!**********************************************************************************************************************************
  ! -----------------------
  deallocate ( eigenvalue_array_Ry )                                   ! eigenv(ENNM)
  deallocate ( energy_num_chara )                                      ! energy(ENNM)
  deallocate ( wave_K2 )                                               ! wave_K2(FFNM)
  deallocate ( wave_coef_real )                                        ! wavecr(FFNM, ENNM)
  deallocate ( wave_coef_imag )                                        ! waveci(FFNM, ENNM)
  deallocate ( wave_coef_2 )                                           ! wavcr2(FFNM, ENNM)
  deallocate ( wave_coef_2_bond )                                      ! wavec2a(FFNM, ENNM)
  deallocate ( wave_coef_2_anti_bond )                                 ! wavec2b(FFNM, ENNM)
!  deallocate ( degenerate_pw )
  deallocate ( shelter_array )                                      ! shelt4(ENNM)
  deallocate ( shelter_array_bond )                                 ! shelt4b(ENNM)
  deallocate ( shelter_array_anti_bond )                            ! shelt4a(ENNM)
  deallocate ( eigenvalue_array_eV )                                ! eigeve(ENNM)
  deallocate ( eigenvalue_array_sub_ef_eV )                         ! neigev(ENNM)
  ! -----------------------
  !
!=======================================================================
  !
end program main
!**********************************************************************************************************************************
