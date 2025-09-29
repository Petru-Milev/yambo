!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
function TDDFT_ALDA_eh_space_G_kernel(BS_n_g_fxc, BS_T_grp_ip, i_Tp, BS_T_grp_ik, i_Tk) result(H_x)
 !
 ! Calculates the F_xc scattering 
 !
 !  (ic(1),ik(1),is(1)) --<--:...:--<-- (ic(2),ik(2),is(2))
 !                              :Fxc:
 !  (iv(1),ik(3),is(3)) -->--:...:-->-- (iv(2),ik(4),is(4))
 !
 use pars,           ONLY:SP,pi,cZERO, cI_DP
 use FFT_m,          ONLY:fft_size
 use BS,             ONLY:BS_T_group
 use wrapper_omp,    ONLY:Vstar_dot_V_omp, V_dot_V_omp
 use electrons,      ONLY:n_spinor,spin_occ
 use R_lattice,      ONLY:g_vec
 use xc_functionals, ONLY:n_xc, f_xc_mat_xyz_dim
 use electrons,      ONLY:n_spin,n_spinor
 !
 implicit none
 !
 integer,                  intent(in) :: BS_n_g_fxc,i_Tp,i_Tk
 type(BS_T_group), target, intent(in) :: BS_T_grp_ip,BS_T_grp_ik
 !
 complex(SP) :: H_x
 complex(SP) :: tmp_coll(BS_n_g_fxc,i_Tk,n_spinor,n_spinor,n_xc,f_xc_mat_xyz_dim)
 complex(SP) :: tmp_coll_arr(BS_n_g_fxc)  
 !%O_tddft_L(:,i_T_el,i_spinor,j_spinor, n_xc, f_xc_mat_xyz_dim)
 !TO CHECK THE USAGE OF i_Tk here
 !
 integer     :: i_sp_xc,o_sp_xc,i_sp_c,i_sp_v, i_dir 
 !
 i_sp_c=BS_T_grp_ik%table(i_Tk,4)
 i_sp_v=BS_T_grp_ik%table(i_Tk,5)
 !
 !!%O_tddft_L(:,i_T_el,i_spinor,j_spinor, n_xc, f_xc_mat_xyz_dim) 6 ind
 !
 tmp_coll_arr = cmplx(0.0_SP, kind=SP)
 if (n_spinor==1) then
   !LDA PART
   H_x=Vstar_dot_V_omp(BS_n_g_fxc,BS_T_grp_ik%O_tddft_L(:,i_Tk,1,1,1,1),&
   &                              BS_T_grp_ip%O_tddft_R(:,i_Tp,i_sp_c,i_sp_v,1,1))
   !GGA part  
   !Multiplic left collision w iG 
   do i_dir=1,3
     tmp_coll_arr(:) = tmp_coll_arr(:) + cI_DP * g_vec(1:BS_n_g_fxc, i_dir)*conjg(BS_T_grp_ik%O_tddft_L(:,i_Tk,1,1,2,i_dir)) 
   enddo
   H_x=H_x+(2/fft_size)*V_dot_V_omp(BS_n_g_fxc, tmp_coll_arr(:), BS_T_grp_ip%O_tddft_R(:,i_Tp,i_sp_c,i_sp_v,2,1))


 else if (n_spinor==2) then
   H_x=cZERO
   do i_sp_xc=1,n_spinor
     do o_sp_xc=1,n_spinor
       H_x=H_x+Vstar_dot_V_omp(BS_n_g_fxc,BS_T_grp_ik%O_tddft_L(:,i_Tk,i_sp_xc,o_sp_xc,1,1),&
       &                                  BS_T_grp_ip%O_tddft_R(:,i_Tp,i_sp_xc,o_sp_xc,1,1))
     enddo
   enddo
 endif
 !
 ! tddft_alda_r_space should be mutiplied by X, it is mutiplied by Co in K
 ! X  =      spin_occ/DL_vol/Nq 
 ! Co = 4*pi*spin_occ/DL_vol/Nq
 ! --> X/Co = 1/4/pi
 !
 H_x=H_x/4._SP/pi !TO CHECK
 !
end function TDDFT_ALDA_eh_space_G_kernel
