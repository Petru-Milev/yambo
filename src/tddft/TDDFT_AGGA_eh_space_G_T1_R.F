
! License-Identifier: GPL
!
! Copyright (C) 2025 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
#include<dev_defs.h>
!
subroutine TDDFT_AGGA_eh_space_G_T1_R(iq,Xk,i_T_grp,NG,l_bs_tddft_wf_in_loop,tddft_wf,mode)

! Calculates the first gga term of Fxc for GGA.
! Term in LaTeX 
! T^{GGA1}_{(n1,n2),(n3,n4)} &=
! `\int dr_1 dr_2\,  \phi_{n1}(r_1) \phi^*_{n2}(r_1)\,
! f_2(r_1) \frac{d \sigma(\mathbf{r_1})}{d n(\mathbf{r_2})}\,
! \phi^*_{n3}(r_2) \phi_{n4}(r_2)
! \notag \\


 use pars,           ONLY:SP,DP,pi,cZERO
 use FFT_m,          ONLY:fft_size
 use wave_func,      ONLY:WF
 use interfaces,     ONLY:WF_load,WF_free
 use R_lattice,      ONLY:qindx_X,bz_samp,minus_G
 use BS,             ONLY:l_BS_magnons,BS_T_grp,BS_bands
 use electrons,      ONLY:n_spinor,n_sp_pol,n_spin
 use timing_m,       ONLY:timing
 use openmp,         ONLY:OPENMP_update,master_thread
 use interfaces,     ONLY:WF_apply_symm
 use TDDFT,          ONLY:tddft_wf_t
 use xc_functionals, ONLY:xc_fnctl_t, xc_setup_fnctl, XC_FAMILY_LDA, &
         & XC_FAMILY_GGA, xc_spin 
 use xc_f03_lib_m,   ONLY:xc_f03_gga_fxc

 implicit none 

 type(bz_samp),intent(in) :: Xk
 integer,      intent(in) :: iq,i_T_grp,NG(2)
 logical,      intent(in) :: l_bs_tddft_wf_in_loop
 character(1), intent(in) :: mode
 type(tddft_wf_t), target, intent(inout) :: tddft_wf

 logical     :: l_load_WFs
 integer     :: is1, is1, i1        !is - iter over spins 
 integer     :: i_spinor,j_spinor,ip_spinor,jp_spinor,ifft,irhotw,is_yambo,os_yambo,i_spinor_y,j_spinor_y
 integer     :: i_T_el,i_T_el_p,N_T_el_p,i_c,i_v,i_sp_c,i_sp_v,i_k_bz,i_k,i_s,i_g0,i_g1,i_g2,i_p_bz,i_g_p,i_p,&
&               i_sp,i_T_grp_p,is(4),os(4),qs(3),NK(2), FUNCTIONAL, ixc 
 integer(IPL) :: size_IPL
 !
 ! XC routines arrays 
 !
 !to check the usage os xc_spin(1) here 
 real(DP) :: spin_rho_DP(n_spin, fft_size), v2rho2(xc_spin(1),fft_size),v2rhosigma(xc_spin(1)*n_spin,fft_size),&
 & v2sigma2(xc_spin(1)*n_spin,fft_size), f1(xc_spin(1),fft_size), f2(xc_spin(1)*n_spin,fft_size), f3((xc_spin(1)*n_spin,fft_size)
 real(DP) :: sigma(xc_spin(1), fft_size)
 real(SP) :: drho(fft_size,n_spin,3) 
 type(xc_fnctl_t) :: fnctl(2)
 !
 complex(DP), pointer DEV_ATTR :: rhotwr_p(:)
 complex(SP), pointer :: WF_symm1_p(:,:)
 complex(SP), pointer :: WF_symm2_p(:,:)
 !
 call timing('T_space ALDA Osc.',OPR='start')

 !put them under if gga func condition when alloc in the main stuff
 v2rho2=0._DP
 v2rhosigma=0._DP
 v2sigma2=0._DP
 f1=0._DP !f1 = v2rho2 ... etc 
 f2=0._DP !see notation in the notes   
 f3=0._DP 

 !Getting information about the functional to be used by libxc 

 call xc_setup_functional(fnctl,FUNCTIONAL,n_spin)

  
 WF_symm1_p => tddft_wf%WF_symm1
 WF_symm2_p => tddft_wf%WF_symm2
 rhotwr_p   => tddft_wf%rhotwr_DP
 !
 WF_symm1_p=cZERO
 WF_symm2_p=cZERO
 rhotwr_p=  cZERO

 size_IPL=int(fft_size,kind=IPL) !What is IPL?

 !drho - gradient of the density, sigma - density gradient in libxc formalism. 
 do i1 = 1,fft_size
   do is1 = 1,n_spin
     do is2 = is1,n_spin
       is_xc = xc_spin_map1(is1,is2)
       sigma(is_xc,i1) = real(drho(i1,is1,1)*drho(i1,is2,1)+drho(i1,is1,2)*drho(i1,is2,2)+drho(i1,is1,3)*drho(i1,is2,3),SP)
     end do
   end do
 end do

 !sum over exch and corr components 
 do ixc = 1,2 
   if (fnctl(ixc)%id == 0 .or. fnctl(ixc)%family == XC_FAMILY_LDA) cycle
   if (fnctl(ixc)%family==XC_FAMILY_GGA)
     v2rho2=0._DP
     v2rhosigma=0._DP
     v2sigma2=0._DP
     !check that it is actually returning the proper value, and not the what was before + value 
     call xc_f03_gga_fxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma(1,1), v2rho2(1,1), v2rhosigma(1,1), v2sigma2(1,1))
     f1(:,:) = f1(:,:) + v2rho2(:,:)     !LDA term 
     f2(:,:) = f2(:,:) + v2rhosigma(:,:) !GGA1 and GGA2 term 
     f3(:,:) = f3(:,:) + v2sigma2(:,:)   !GGA3 term 
   end if 
 end do 

 !What does this do stand for? 
 do i_T_el=1,BS_T_grp(i_T_grp)%size
   !
   i_k_bz=BS_T_grp(i_T_grp)%table(i_T_el,1)
   i_k   =Xk%sstar(i_k_bz,1)
   i_s   =Xk%sstar(i_k_bz,2)
   !
   i_p_bz=qindx_X(iq,i_k_bz,1)
   i_g_p =qindx_X(iq,i_k_bz,2)
   !
   i_p   =Xk%sstar(i_p_bz,1)
   i_sp  =Xk%sstar(i_p_bz,2)
   !
   i_v   =BS_T_grp(i_T_grp)%table(i_T_el,2)
   i_c   =BS_T_grp(i_T_grp)%table(i_T_el,3)
   i_sp_c=BS_T_grp(i_T_grp)%table(i_T_el,4)
   i_sp_v=BS_T_grp(i_T_grp)%table(i_T_el,5)
   !
   is=(/0,i_k,i_s ,i_sp_c/)
   os=(/0,i_p,i_sp,i_sp_v/)
   qs=(/minus_G(i_g_p),iq,1/)
   !
   if (mode=="R".or.mode=="C".or.mode=="F") then
     is(1)=i_c
     os(1)=i_v
   else if (mode=="A".or.mode=="Q") then
     ! Same momenta and spin indexes, inverted bands
     is(1)=i_v
     os(1)=i_c
   endif
   !
   l_load_WFs= l_bs_tddft_wf_in_loop .and. (NK(1)/=min(i_k,i_p).or.NK(2)/=max(i_k,i_p))
   if (l_load_WFs) then
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp critical
#endif
     if (NK(2)/=-1) call WF_free(WF,keep_fft=.true.,keep_states_to_load=.true.)
     NK=(/min(i_k,i_p),max(i_k,i_p)/)
     call WF_load(WF,NG(1),NG(2),BS_bands,NK,k_extrema_only=.true.,quiet=.true.,&
     &            space='R',title="Kernel exch",keep_states_to_load=.true.)
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp end critical
#endif
   endif
   !
   !
   call WF_apply_symm(is,WF_symm1_p)
   call WF_apply_symm(os,WF_symm2_p)
   !
   if( n_spinor==1 ) then
     rhotwr_p(:)=conjg(WF_symm1_p(:,1)*cmplx(WF_symm2_p(:,1)),kind=DP)
     call perform_fft_3d(qs,rhotwr_p,BS_T_grp(i_T_grp)%O_tddft_gga_1_R(:,i_T_el,1,1))
   endif
   !
   if(n_spinor==2) then
     !
     do i_spinor=1,n_spinor
       do j_spinor=1,n_spinor
         !
     do ifft=1,fft_size
       rhotwr_p(ifft)=cmplx(WF_symm1_p(ifft,i_spinor)*conjg(WF_symm2_p(ifft,j_spinor)),kind=DP)
         enddo
         call perform_fft_3d(qs,rhotwr_p,BS_T_grp(i_T_grp)%O_tddft_L(:,i_T_el,i_spinor,j_spinor))
         !
       enddo
     enddo
   endif
   !
 enddo
 !

 !
 if(l_bs_tddft_wf_in_loop) call WF_free(WF,keep_fft=.true.,keep_states_to_load=.true.)
 !
 nullify(WF_symm1_p)
 nullify(WF_symm2_p)
 nullify(rhotwr_p)
 !
 call timing('T_space ALDA Osc.',OPR='stop')

end subroutine TDDFT_AGGA_eh_space_G_T1_R 

