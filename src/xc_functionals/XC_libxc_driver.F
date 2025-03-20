!
! License-Identifier: GPL
!
! Copyright (C) 2010 The Yambo Team
!
! Authors (see AUTHORS file for details): MG DS DV 
!
subroutine XC_libxc_driver(en,Xk,FUNCTIONAL,ORDER,EXX_FRACTION,EXX_SCREENING)
 !
 ! This routine evaluates the xc functional in yambo 
 !   - 03/05/12: LDA OK, GGA OK => EXTENSIVE TESTS NEEDED) 
 !   - 21/11/13: for HYB_GGA the GGA part calculated => TO TEST
 !   - 16/06/21: GGA part improved with more accurate calculation of drho for order=2 and n_spin=1
 !        To do: compute a more precise drho also in the case n_spin>1 (see comment around line 200)
 !
 use pars,           ONLY:SP,DP,IPL,pi,pi_DP,RZERO
 use com,            ONLY:msg
 use parser_m,       ONLY:parser
 use R_lattice,      ONLY:bz_samp
 use D_lattice,      ONLY:DL_vol,alat
 use drivers,        ONLY:l_real_time,l_nl_optics
#if defined _SC || defined _RT
 use hamiltonian,    ONLY:rho_n,magn_n
#endif
 use electrons,      ONLY:levels,n_spin,n_spinor
 use FFT_m,          ONLY:fft_size
 use wrapper,        ONLY:V_dot_V
 use BS,             ONLY:l_BS_magnons
 use TDDFT,          ONLY:l_Fxc_from_Vxc
 use xc_functionals, ONLY:E_xc,E_xc_val,magn,F_xc,V_xc,xc_setup_hyb,&
&                         XC_LDA_C_KP,xc_fnctl_t,xc_setup_fnctl,xc_f03_func_end, &
&                         xc_spin,xc_spin_map1, &
&                         XC_HYB_GGA_XC_GAUPBE,XC_FACTOR,&
&                         XC_FAMILY_LDA,XC_FAMILY_GGA,XC_FAMILY_HYB_GGA,&
&                         XC_GGA_X_PBE, XC_GGA_C_PBE
 use xc_f03_lib_m,   ONLY:xc_f03_lda_exc,xc_f03_lda_vxc,xc_f03_lda_fxc, &
&                         xc_f03_gga_exc,xc_f03_gga_vxc,xc_f03_gga_fxc
 use interfaces,     ONLY:el_density_and_current,el_magnetization
 use pseudo,         ONLY:pp_rho_nlcc,pp_has_nlcc,PP_nlcc_free
 !
#include<y_memory.h>
 !
 type(levels) ::en
 type(bz_samp)::Xk
 integer      ::FUNCTIONAL,FUNCTIONALTMP,ORDER
 real(SP), optional :: EXX_FRACTION, EXX_SCREENING
 !
 ! Work Space
 !
 integer  :: i1,i2,is_xc,is1,is2,ixc,ifft,ORDER_NOW, spin_iter, i3 
 integer(IPL) :: size_IPL
 real(SP) :: mod_magn
 real(SP) :: local_exx_fraction, PERT_FAC_SIGMA
 real(DP) :: exsr,PERT_FAC
 !
 ! XC routines arrays
 !
 real(SP) :: spin_rho_SP(fft_size,n_spin),rho(fft_size)
 real(DP) :: spin_rho_DP(n_spin,fft_size), spin_rho_DP_tmp(n_spin,fft_size), rspts(fft_size)
 real(SP),allocatable::   drho_tmp(:,:),V_xc_pert(:,:,:), tmp(:), tmp_2dim (:,:), drho(:,:,:), drho_pert(:,:,:),&
  &                                   vxc_tmp(:,:)! ,vxc_tmp(:,:), V_xc_pert_tmp(:,:) 
 real(DP),allocatable::   sigma(:,:), sigma_pert(:,:,:,:), V_xc_pert_cd(:,:,:,:), sigma_tmp(:,:)
 complex(SP), allocatable:: drho_cmplx(:,:) !complex(SP), allocatable
 complex(SP), allocatable:: spin_rho_SP_cmplx(:), spin_rho_DP_cmplx(:)
 !
 ! note that multidimensional fxc & vxc is allowed only in xcspol
 !
 real(DP),allocatable :: tmp_exc(:)
 real(DP),allocatable :: v1rho(:,:),v1rho_pert(:,:),v1sigma(:,:),v1sigma_pert(:,:), vxc_rho_pert(:,:,:,:), &
 & vxc_sigma_pert(:,:,:,:), f1rho(:,:), f1sigma(:,:)
 real(DP),allocatable :: vx1(:,:),vx2(:,:)
 real(DP),allocatable :: v2rho2(:,:),v2rhosigma(:,:),v2sigma2(:,:)
 logical              :: EvaluateDGradient,l_Fxc_Libxc
 !
 ! Functional and Functional infos
 !
 type(xc_fnctl_t) :: fnctl(2) 
 !
 !=============================== 
 ! Convert and Initialize
 !===============================
 !
 call parser('FxcLibxc',l_Fxc_Libxc)
 l_Fxc_from_Vxc=.not.l_Fxc_Libxc
 !
 l_Fxc_from_Vxc=.TRUE.
 l_Fxc_Libxc=.FALSE.
 exsr=0._DP
 FUNCTIONALTMP=0
!Set the functional as a PBE to get Exc, next calculate the GAU-PBE Exc as
!Exc=Exc(PBE)-alfa*Ex(PBE,SR) next add the HF part.
 !
 if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) then 
   FUNCTIONALTMP=XC_HYB_GGA_XC_GAUPBE*XC_FACTOR
   FUNCTIONAL=XC_GGA_X_PBE*XC_FACTOR+XC_GGA_C_PBE  !101130   This is semi-local PBE  
 endif
 !
 call xc_setup_fnctl(fnctl,FUNCTIONAL,n_spin)
 if (FUNCTIONALTMP==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) FUNCTIONAL=XC_HYB_GGA_XC_GAUPBE*XC_FACTOR
 !
 local_exx_fraction = RZERO
 if (present(EXX_FRACTION).and.fnctl(1)%family==XC_FAMILY_HYB_GGA) then 
   call xc_setup_hyb(fnctl,EXX_FRACTION,EXX_SCREENING)
   local_exx_fraction = fnctl(1)%exx_fraction
 end if
 if (FUNCTIONALTMP==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) local_exx_fraction=EXX_FRACTION
 !
 EvaluateDGradient = any((fnctl(:)%family==XC_FAMILY_GGA).or.any(fnctl(:)%family==XC_FAMILY_HYB_GGA))
 !
 !==========================
 ! Electronic/spin density 
 !==========================
 !
#if defined _SC || defined _RT
 if (l_real_time.or.l_nl_optics) then !.or.l_sc_run) then
   rho=rho_n
   if(n_spin> 1) magn=magn_n
 endif
#endif
 !
 if (.not.(l_real_time.or.l_nl_optics)) then !.or.l_sc_run)) then
   call el_density_and_current(en,Xk,rho=rho)
   if (n_spin> 1) call el_magnetization(en,Xk,magn)
 endif
 !
 ! non-linear core corrections (NLCC)
 ! to be improved for yambo_rt, otherwise this is done at each time step
 !
 call XC_nlcc_setup()
 if (pp_has_nlcc) rho(:)=rho(:)+real(pp_rho_nlcc(:),SP)
 call PP_nlcc_free()
 !
 if(n_spin==1) spin_rho_SP(:,1)=rho(:)
 if(n_spin> 1) then
   if (n_spinor==1) then
     !$omp parallel do default(shared), private(ifft)
     do ifft=1,fft_size
       spin_rho_SP(ifft,1)=0.5*(rho(ifft)+magn(ifft,3))
       spin_rho_SP(ifft,2)=0.5*(rho(ifft)-magn(ifft,3))
     enddo
     !$omp end parallel do
   else
     !$omp parallel do default(shared), private(ifft,mod_magn)
     do ifft=1,fft_size
       mod_magn=sqrt(magn(ifft,1)**2+magn(ifft,2)**2+magn(ifft,3)**2)
       spin_rho_SP(ifft,1)=0.5*(rho(ifft)+mod_magn)
       spin_rho_SP(ifft,2)=0.5*(rho(ifft)-mod_magn)
     enddo
     !$omp end parallel do
   endif
 endif
 !
 do i1=1,fft_size
   do i2=1,n_spin
     spin_rho_DP(i2,i1)=real(spin_rho_SP(i1,i2),DP)*real(fft_size,DP)/real(DL_vol,DP)
     ! DEBUG <
     ! This was needed for absorption using old libxc
     ! Keeping it with magnons gives numberical noise
     !if (spin_rho_DP(i2,i1)<1.E-10_DP) spin_rho_DP(i2,i1)=0._DP
     ! DEBUG >
   enddo
 enddo

 open(10,file='spin_rho_DP.dat',status='unknown')
  write(10,*) "spin_rho_DP"
  do ifft=1,fft_size
    write(10,*) ifft, (spin_rho_DP(i1,ifft), i1=1,n_spin)
  enddo
  close(10)
 !
 !=============================
 ! Electronic density gradient
 !=============================
 if (ORDER==0) then
   YAMBO_ALLOC(tmp_exc,(fft_size))
   tmp_exc=0._SP
 endif
 if (ORDER==1.or.ORDER==2) then
   YAMBO_ALLOC(v1rho,(n_spin,fft_size))
   v1rho=0._DP
 endif
 if (ORDER==2.and..not.l_Fxc_from_Vxc) then
   YAMBO_ALLOC(v2rho2,(xc_spin(1),fft_size))
   v2rho2=0._DP
 endif
 !
 if (EvaluateDGradient) then 
   !
   YAMBO_ALLOC(drho,(fft_size,n_spin,3))
   YAMBO_ALLOC(sigma,(xc_spin(1),fft_size))
   if (ORDER==1.or.ORDER==2) then
     YAMBO_ALLOC(v1sigma,(xc_spin(1),fft_size))
     v1sigma=0._DP
   endif
   if (ORDER==2.and..not.l_Fxc_from_Vxc) then
     YAMBO_ALLOC(v2rhosigma,(xc_spin(1)*n_spin,fft_size))
     YAMBO_ALLOC(v2sigma2,(xc_spin(1)*n_spin,fft_size))
     v2rhosigma=0._DP
     v2sigma2=0._DP
   endif
   if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) then
     YAMBO_ALLOC(vx1,(n_spin,fft_size))
     YAMBO_ALLOC(vx2,(xc_spin(1),fft_size))
     vx1=0._DP
     vx2=0._DP
   endif 
   !
   YAMBO_ALLOC(drho_cmplx,(fft_size,3))
   drho_cmplx=0._SP
   if (n_spin>1 .or. order<2 .or. l_Fxc_from_Vxc) then
     YAMBO_ALLOC(spin_rho_SP_cmplx,(fft_size))  
     do is1 = 1,n_spin
       spin_rho_SP_cmplx(:)=cmplx(spin_rho_SP(:,is1),0.0_SP,kind=SP)

       call msg('s',' Here')
       call eval_Gradient(spin_rho_SP_cmplx,drho_cmplx,1,"density")
       drho(:,is1,:)=real(drho_cmplx)
     end do
     !Y(MBO_FREE(drho_cmplx)
     !Y_AMBO_FREE(spin_rho_SP_cmplx)
   else
     ! DS: This procedure computes a more precise drho (much less noise)
     !     This is particularly important for the case order=2
     ! To do: Extend the procedure to the case n_spin>1
     YAMBO_ALLOC(drho_tmp,(fft_size,3))
     call el_density_and_current(en,Xk,drho=drho_tmp)
     drho(:,1,:) = drho_tmp(:,:)
     YAMBO_FREE(drho_tmp)
   endif
   !
   forall (is1=1:n_spin,i1=1:fft_size) drho(i1,is1,:) = drho(i1,is1,:)*fft_size/DL_vol*2._SP*pi/alat(:)
   !
   do i1 = 1,fft_size
     do is1 = 1,n_spin
       do is2 = is1,n_spin
         is_xc = xc_spin_map1(is1,is2)
         sigma(is_xc,i1) = real(drho(i1,is1,1)*drho(i1,is2,1)+drho(i1,is1,2)*drho(i1,is2,2)+drho(i1,is1,3)*drho(i1,is2,3),DP)
       end do
     end do
   end do 

 endif
 !
 !===========================
 ! Evaluate the xc functional
 !===========================
 !
 if (order==0)                                  E_xc=0._SP
 if (order==1.or.(order==2.and.l_Fxc_from_Vxc)) V_xc=0._SP
 if (order==2)                                  F_xc=0._SP
 !
 
 !TODO put some if stataments here
  YAMBO_ALLOC(f1rho,(n_spin,fft_size))
  f1rho=0._DP
  YAMBO_ALLOC(f1sigma,(xc_spin(1),fft_size))
  f1sigma=0._DP
  YAMBO_ALLOC(vxc_sigma_pert,(fft_size, n_spin, xc_spin(1) ,2))  
  !For vxc_sigma_pert. Second index - values of Vxc up or down, 
  !third (xc_spin) - index of which component of sigma is perturbed
  !fourth - 1 plus perturbation, 2 minus perturbation
  vxc_sigma_pert = 0._DP
  YAMBO_ALLOC(vxc_rho_pert,(fft_size,n_spin,n_spin,2))
  vxc_rho_pert = 0._DP
  YAMBO_ALLOC(vxc_tmp,(fft_size, n_spin))
  vxc_tmp = 0._DP
  YAMBO_ALLOC(sigma_tmp,(xc_spin(1),fft_size))
  sigma_tmp = 0._DP
 
if (EvaluateDGradient) then
do is1 = 1,n_spin
  spin_rho_SP_cmplx(:)=cmplx(spin_rho_SP(:,is1),0.0_SP,kind=SP)

  call msg('s',' Here')
  call eval_Gradient(spin_rho_SP_cmplx,drho_cmplx,1,"density")
  drho(:,is1,:)=real(drho_cmplx)
end do
endif

 if (FUNCTIONAL == XC_LDA_C_KP) then ! Special case, XC_LDA_C_KP not available from libxc
   select case(ORDER)
     case(0)
       forall (i1=1:fft_size) rspts(i1)=(3._DP/4._DP/pi_DP/real(spin_rho_SP(i1,1)/DL_vol,DP))**(1._DP/3._DP)
       call xc_rpa_kp(rspts,v1rho)
       E_xc = real(tmp_exc,SP) 
       E_xc_val = V_dot_V(fft_size,rho,E_xc(:))
     case(1)
       call error(' Vxc not coded for LDA C KP')
     case(2)
       call error(' Fxc not coded for LDA C KP')
   end select
   return
 end if
 !
 size_IPL=int(fft_size,kind=IPL)
 !
 ORDER_NOW=ORDER
 if (l_Fxc_from_Vxc.and.ORDER==2) then
   call msg('s',' Evaluating Fxc as numerical functional derivative of Vxc')
   !TODO test this better, and put if condition for case of LDA and GGA separately
   !TODO put the perturbation factor in the input file?
   !TODO CHANGE IT NOW FOR THE TESTS
   !PERT_FAC=0.0000001_SP
   PERT_FAC=0.1D-2
   PERT_FAC_SIGMA=0.0000001_SP
   ORDER_NOW=1
   YAMBO_ALLOC(V_xc_pert,(fft_size,n_spin,n_spin)) 
   ! DS: first n_spin is for the diagonal elements of v_xc, the second is for the spin_up and spin_down perturbation
   YAMBO_ALLOC(V_xc_pert_cd,(fft_size,n_spin,n_spin,2)) 
   V_xc_pert=0._SP
   V_xc_pert_cd=0._DP
   !last n_spin(3) is for spin up/down perturbations, and 2 corresponts to plus and minus perturbations 
   !Y_AM__B_O_ALL_O_C(V_xc_pert_tmp, (fft_size,n_spin)) to delete
   YAMBO_ALLOC(spin_rho_DP_cmplx, (fft_size))
   YAMBO_ALLOC(v1rho_pert,(n_spin,fft_size))
   YAMBO_ALLOC(v1sigma_pert,(xc_spin(1),fft_size))
   YAMBO_ALLOC(drho_pert,(fft_size,n_spin,3))
 endif
  !
 do ixc = 1,2
   if (fnctl(ixc)%id == 0) cycle
   if (fnctl(ixc)%family==XC_FAMILY_LDA) then 
     select case(ORDER_NOW)
       case(0)
         call xc_f03_lda_exc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), tmp_exc(1))
         E_xc = E_xc + real(tmp_exc,SP)    
       case(1)
         call xc_f03_lda_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), v1rho(1,1))
         forall(i1=1:fft_size,i2=1:n_spin) V_xc(i1,i2) = V_xc(i1,i2) + real(v1rho(i2,i1),SP)
         if (l_Fxc_from_Vxc.and.ORDER==2) then
           !spin_rho_DP=spin_rho_DP*PERT_FAC
           YAMBO_ALLOC(v1rho_pert,(n_spin,fft_size))
           do is1=1,n_spin
              spin_rho_DP(is1,:)=spin_rho_DP(is1,:)+PERT_FAC
              call xc_f03_lda_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), v1rho_pert(1,1))
              !forall(i1=1:fft_size, i2=1:n_spin) V_xc_pert(i1,i2,is1) = V_xc_pert(i1,i2,is1) + real(v1rho_pert(i2,i1),SP) !WHY SP?
              forall(i1=1:fft_size, i2=1:n_spin) V_xc_pert_cd(i1,i2,is1,1) = V_xc_pert_cd(i1,i2,is1,1) + real(v1rho_pert(i2,i1),DP) 
              !spin_rho_DP(is1,:)=spin_rho_DP(is1,:)-PERT_FAC
              spin_rho_DP(is1,:)=spin_rho_DP(is1,:)-(2*PERT_FAC)
              call xc_f03_lda_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), v1rho_pert(1,1))
              forall(i1=1:fft_size, i2=1:n_spin) V_xc_pert_cd(i1,i2,is1,2) = V_xc_pert_cd(i1,i2,is1,2) + real(v1rho_pert(i2,i1),DP) 
              !WHY SP?
              spin_rho_DP(is1,:)=spin_rho_DP(is1,:)+PERT_FAC
           enddo  
           !call xc_f03_lda_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), v1rho(1,1))
           !forall(i1=1:fft_size,i2=1:n_spin) V_xc_pert(i1,i2) = V_xc_pert(i1,i2) + real(v1rho(i2,i1),SP)
           !YAM_BO_FR_EE(v1rho_pert)
         endif
       case(2)
         if( l_BS_magnons) then
           call xc_f03_lda_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), v1rho(1,1))
           forall(i1=1:fft_size,i2=1:n_spin) V_xc(i1,i2) = V_xc(i1,i2) + real(v1rho(i2,i1),SP)
         endif
         ! I could avoid this, however I need it due to numerical noise
         ! See also comment in Build_F_xc_mat
         !if( (.not.l_BS_magnons) .or. n_spinor==2) then
           call xc_f03_lda_fxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), v2rho2(1,1))
           ! Yambo internal for Kxc
           call XC_eval_lda_kernel(v1rho,v2rho2)
         !endif
     end select
          !GGA PART     
   elseif (fnctl(ixc)%family==XC_FAMILY_GGA.or.fnctl(ixc)%family==XC_FAMILY_HYB_GGA) then
        !GAUPBE short range x-energy  and potentials part
        call msg('s','Here 1 ')
        if(FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) call XC_eval_pbegaux(spin_rho_DP(1,1), sigma(1,1), exsr,vx1(1,1),vx2(1,1))
     select case(ORDER_NOW)
       case(0)
         call xc_f03_gga_exc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma(1,1), tmp_exc(1))
         E_xc = E_xc + real(tmp_exc,SP)    
       case(1)
         call xc_f03_gga_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma(1,1), v1rho(1,1), v1sigma(1,1))
         if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR.and.ixc==1) then
           !GAUPBE short energy part of the x potential,factor 2 in vx2 for
           !consistency qith QE
           v1rho=v1rho-EXX_FRACTION*vx1
           v1sigma=v1sigma-(EXX_FRACTION*vx2)/2._DP
         endif
         ! Yambo internal for V_xc
         call XC_eval_gga_potential(v1rho,v1sigma,drho,V_xc)
         open (13, file='vxc_default.dat', status='unknown')
         write(13,*) 'vxc_sigma_pert'
         do ifft=1,fft_size
           write(13,*) ifft, V_xc(ifft,1), V_xc(ifft,2)
         enddo
         close(13)


         if (l_Fxc_from_Vxc.and.ORDER==2) then
           !TO DELETE .or.n_spin==1 when the n_spin==1 case is implemented
             ! Here we peturb only density, without perturbing the sigma
            do is1=1,n_spin ! \partial V_xc^{alpha} / \partial rho^{is1} - notation used from the notes
              call msg('s','Here 4')
              f1rho = 0._DP
              f1sigma = 0._DP
              spin_rho_DP_tmp(:,:)=spin_rho_DP(:,:)
              spin_rho_DP_tmp(is1,:)=spin_rho_DP_tmp(is1,:)+PERT_FAC
              !spin_rho_DP(is1,:)=spin_rho_DP(is1,:)+PERT_FAC
              call msg('s','Here 4.25')
              call xc_f03_gga_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP_tmp(1,1), sigma(1,1),& 
                                  v1rho_pert(1,1),v1sigma_pert(1,1))
              call msg('s','Here 4.05')     
              f1rho(:,:)=f1rho(:,:)+real(v1rho_pert(:,:), DP)
              f1sigma(:,:)=f1sigma(:,:)+real(v1sigma_pert(:,:), DP)
              call msg('s','Here 4.1')
              !call XC_eval_gga_potential(v1rho_pert,v1sigma_pert,drho_pert,vxc_tmp(:,:)) 
              !returned \partial V_xc^{up} / \partial rho^{is1} and partial V_xc^{down} / \partial rho^{is1}
              !vxc_rho_pert(:,:,is1,1) = vxc_rho_pert(:,:,is1,1) + real(vxc_tmp(:,:))
              !forall (i1=1:fft_size, i2=1:n_spin) vxc_rho_pert(i1,i2,is1,1) = &
              !&                                   vxc_rho_pert(i1,i2,is1,1) + real(v1rho_pert(i2,i1),DP)
              spin_rho_DP_tmp(:,:)=spin_rho_DP(:,:)
              spin_rho_DP_tmp(is1,:)=spin_rho_DP_tmp(is1,:)-PERT_FAC
              !spin_rho_DP(is1,:)=spin_rho_DP(is1,:)-2*PERT_FAC
              
              
              call xc_f03_gga_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP_tmp(1,1), sigma(1,1),& 
                                  v1rho_pert(1,1),v1sigma_pert(1,1))
              !forall (i2=1:n_spin) f1rho(i2,:) = real(real(vxc_rho_pert(:,i1,i2,1) - vxc_rho_pert(:,i1,i2,2),DP)/(2*PERT_FAC),DP)
              f1rho(:,:)=real(f1rho(:,:)-real(v1rho_pert(:,:),DP),DP)/(2*PERT_FAC)
              f1sigma(:,:)=real(f1sigma(:,:)-real(v1sigma_pert(:,:),DP),DP)/(2*PERT_FAC)
              call msg('s','Here 4.15')
              call XC_eval_gga_potential(f1rho,f1sigma,drho,vxc_tmp(:,:))
              vxc_rho_pert(:,:,is1,1)=vxc_rho_pert(:,:,is1,1) + real(vxc_tmp(:,:))


              !call XC_eval_gga_potential(v1rho_pert,v1sigma_pert,drho_pert,vxc_tmp(:,:)) 
              call msg('s','Here 4.2')
              !forall (i1=1:fft_size, i2=1:n_spin) vxc_rho_pert(i1,i2,is1,2) = &
              !&                                   vxc_rho_pert(i1,i2,is1,2) + real(v1rho_pert(i2,i1),DP)
              !vxc_rho_pert(:,:,is1,2)= vxc_rho_pert(:,:,is1,2) + real(vxc_tmp(:,:)) !CHECK AND MAYBE TO CHANGE INDECES
              !spin_rho_DP(is1,:)=spin_rho_DP(is1,:)+PERT_FAC
            end do 
            
            
            ! Here we perturn only the sigma, without perturbing the density
            !do spin_iter = 1,n_spin ! \partial V_xc^{spin_iter} / \partial sigma^{is1,is2} - notation used from the notes
              call msg('s','Here 5')
              do is1=1,n_spin
                do is2=is1,n_spin
                    is_xc = xc_spin_map1(is1,is2)
                    f1rho = 0._DP
                    f1sigma = 0._DP
                    !if (is_xc == 1) then 
                    !  spin_rho_SP(:,1)=spin_rho_SP(:,1)+PERT_FAC_SIGMA
                    !else if (is_xc == 2) then
                    !  spin_rho_SP(:,1)=spin_rho_SP(:,1)+PERT_FAC_SIGMA
                    !  spin_rho_SP(:,2)=spin_rho_SP(:,2)+PERT_FAC_SIGMA
                    !else if (is_xc == 3) then
                    !  spin_rho_SP(:,2)=spin_rho_SP(:,2)+PERT_FAC_SIGMA
                    !end if
                    !call msg('s','Here 5.1')
                    !do i3 = 1,n_spin
                    ! spin_rho_SP_cmplx(:)=cmplx(spin_rho_SP(:,i3),0.0_SP,kind=SP)
                    !  call eval_Gradient(spin_rho_SP_cmplx,drho_cmplx,1,"density")
                    !  drho_pert(:,i3,:)=real(drho_cmplx)
                    !end do
                    !call msg('s','Here 5.2')
                    !forall (i3=1:n_spin,i1=1:fft_size) drho_pert(i1,i3,:) = drho_pert(i1,i3,:)*fft_size/DL_vol*2._SP*pi/alat(:)
                    !sigma_tmp(:,:) = sigma_pert(:,:,is_xc,1)
                    !call msg('s','Here 5.3')
                    sigma_tmp(:,:) = sigma(:,:)
                    sigma_tmp(is_xc,:) = sigma_tmp(is_xc,:) + PERT_FAC_SIGMA
                    call msg('s','Here 5.31')
                    call xc_f03_gga_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma_tmp(1,1),& 
                                    v1rho_pert(1,1),v1sigma_pert(1,1))
                    f1rho(:,:)=f1rho(:,:)+real(v1rho_pert(:,:))
                    f1sigma(:,:)=f1sigma(:,:)+real(v1sigma_pert(:,:))
                    !call xc_f03_gga_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma(1,1),& 
                    !                v1rho_pert(1,1),v1sigma_pert(1,1))
                    !call XC_eval_gga_potential(v1rho_pert,v1sigma_pert,drho,vxc_tmp(:,:)) !wil return Vxc_tmp(fft, n_spin)
                    !call msg('s','Here 5.4')
                    !vxc_sigma_pert(:,:,is_xc,1)=vxc_sigma_pert(:,:,is_xc,1) + vxc_tmp(:,:)

                    !if (is_xc == 1) then 
                    !  spin_rho_SP(:,1)=spin_rho_SP(:,1)-2*PERT_FAC_SIGMA
                    !else if (is_xc == 2) then
                    !  spin_rho_SP(:,1)=spin_rho_SP(:,1)-2*PERT_FAC_SIGMA
                    !  spin_rho_SP(:,2)=spin_rho_SP(:,2)-2*PERT_FAC_SIGMA
                    !else if (is_xc == 3) then
                    !  spin_rho_SP(:,2)=spin_rho_SP(:,2)-2*PERT_FAC_SIGMA
                    !end if

                    !do i3 = 1,n_spin
                    !  spin_rho_SP_cmplx(:)=cmplx(spin_rho_SP(:,i3),0.0_SP,kind=SP)
                    !  call eval_Gradient(spin_rho_SP_cmplx,drho_cmplx,1,"density")
                    !  drho_pert(:,i3,:)=real(drho_cmplx)
                    !end do
                    !forall (i3=1:n_spin,i1=1:fft_size) drho_pert(i1,i3,:) = drho_pert(i1,i3,:)*fft_size/DL_vol*2._SP*pi/alat(:)
                    sigma_tmp(:,:) = sigma(:,:)
                    sigma_tmp(is_xc,:) = sigma_tmp(is_xc,:) - PERT_FAC_SIGMA

                    call msg('s','Here 5.5')
                    call xc_f03_gga_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma_tmp(1,1),& 
                                    v1rho_pert(1,1),v1sigma_pert(1,1))
                    f1rho(:,:)=real(f1rho(:,:)-real(v1rho_pert(:,:)), DP)/(2*PERT_FAC_SIGMA)
                    f1sigma(:,:)=real(f1sigma(:,:)-real(v1sigma_pert(:,:), DP))/(2*PERT_FAC_SIGMA)
                    !call XC_eval_gga_potential(v1rho_pert,v1sigma_pert,drho,vxc_tmp(:,:)) !wil return Vxc_tmp(fft, n_spin)\
                    call XC_eval_gga_potential(f1rho,f1sigma,drho,vxc_tmp(:,:))
                    vxc_sigma_pert(:,:,is_xc,1)=vxc_sigma_pert(:,:,is_xc,1) + vxc_tmp(:,:)
                    !vxc_sigma_pert(:,:,is_xc,2)=vxc_sigma_pert(:,:,is_xc,2) + vxc_tmp(:,:)
                    !call msg('s','Here 5.6')
                    !if (is_xc == 1) then 
                    !  spin_rho_SP(:,1)=spin_rho_SP(:,1)+PERT_FAC_SIGMA
                    !else if (is_xc == 2) then
                    !  spin_rho_SP(:,1)=spin_rho_SP(:,1)+PERT_FAC_SIGMA
                    !  spin_rho_SP(:,2)=spin_rho_SP(:,2)+PERT_FAC_SIGMA
                    !else if (is_xc == 3) then
                    !  spin_rho_SP(:,2)=spin_rho_SP(:,2)+PERT_FAC_SIGMA
                    !end if

                    call msg('s','Here 6')
                end do 
              end do 
            !end do
           
           !call XC_eval_gga_potential(v1rho,v1sigma,drho,V_xc_pert)
           !spin_rho_DP=spin_rho_DPPERT_FAC
         end if
       case(2)
         if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR) call error(' Fxc not coded for GAU-PBEP')
         if( l_BS_magnons) then
           call xc_f03_gga_vxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma(1,1), v1rho(1,1), v1sigma(1,1))
           ! Yambo internal for V_xc
           call XC_eval_gga_potential(v1rho,v1sigma,drho,V_xc)
         endif
         ! I could avoid this, however I need it due to numerical noise
         ! See also comment in Build_F_xc_mat
         !if( (.not.l_BS_magnons) .or. n_spinor==2) then
           call xc_f03_gga_fxc(fnctl(ixc)%conf, size_IPL, spin_rho_DP(1,1), sigma(1,1), v2rho2(1,1), v2rhosigma(1,1), v2sigma2(1,1))
           ! Yambo internal for Kxc
           call XC_eval_gga_kernel(v1rho,v1sigma,v2rho2,v2rhosigma,v2sigma2)
         !endif
     end select
   end if
   call xc_f03_func_end(fnctl(ixc)%conf)
 enddo
 !


 if (order==0) then
   E_xc_val = V_dot_V(fft_size,rho,E_xc(:))
   if (FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR)  E_xc_val= E_xc_val-EXX_FRACTION*real(exsr,SP)
 endif
 !


 if (l_Fxc_from_Vxc.and.ORDER==2) then
   if (any(fnctl(:)%family==XC_FAMILY_LDA)) then
     call msg('s','Inside the numder loop. LDA')
     do i1=1,n_spin
       do i2=1,n_spin
         do ifft=1, fft_size
           F_xc(ifft,i1,i2)=real(real(V_xc_pert_cd(ifft,i1,i2,1)-V_xc_pert_cd(ifft,i1,i2,2),DP)/(2*PERT_FAC),SP)
         enddo 
       enddo
     enddo
     if (n_spin==2) then
       YAMBO_ALLOC(tmp, (fft_size))
       !12 and 21 components should be equal, here we try to reduce numerical error 
       !Can there exist systems or cases when they are not the same?
       !TO take comment out 
       !tmp=(F_xc(:,1,2)+F_xc(:,2,1))/2._DP !to check  
       !F_xc(:,1,2)=tmp 
       !F_xc(:,2,1)=tmp 
       YAMBO_FREE(tmp)
     endif 
   else if (any((fnctl(:)%family==XC_FAMILY_GGA).or.any(fnctl(:)%family==XC_FAMILY_HYB_GGA))) then
     call msg('s','Inside the numder loop. GGA')
     !vxc_rho_pert is (fft_size, n_spin, n_spin, 2) and vxc_sigma_pert is (
     !f1rho is (n_spin, fft_size) and f1sigma is (xc_spin, fft_size)
     !XC_eval_gga_potential return vxc_tmp(fft_size, n_spin)
     !Here we used the XC_eval_gga because to calc Fxc we arrived to the same structure as for the Vxc
     !If it works it will be explained in the notes
      do i1=1,n_spin
        !do i2=1,n_spin
          forall (ifft=1:fft_size, i3=1:n_spin) f1rho(i3,ifft) = real(vxc_rho_pert(ifft,i1,i3,1))
          forall (i3=1:xc_spin(1), ifft=1:fft_size) f1sigma(i3,ifft) = real(vxc_sigma_pert(ifft,i1,i3,1))
          call XC_eval_gga_potential(f1rho,f1sigma,drho,vxc_tmp)
          F_xc(:,i1,:) = vxc_tmp(:,:)          
        !enddo
      enddo

      open (13, file='vxc_sigma_pert.dat', status='unknown')
      write(13,*) 'vxc_sigma_pert'
      do ifft=1,fft_size
        write(13,*) ifft, vxc_sigma_pert(ifft,1,1,1), vxc_sigma_pert(ifft,1,2,1), vxc_sigma_pert(ifft,1,3,1)
        write(13,*) ifft, vxc_sigma_pert(ifft,2,1,1), vxc_sigma_pert(ifft,2,2,1), vxc_sigma_pert(ifft,2,3,1)
      enddo
      close(13)

      do i1=1,n_spin
        !forall (i2=1:n_spin) f1rho(i2,:) = real(real(vxc_rho_pert(:,i1,i2,1) - vxc_rho_pert(:,i1,i2,2),DP)/(2*PERT_FAC),DP)
        !forall (i2=1:n_spin) f1rho(i2,:) = real(real(vxc_rho_pert(:,i1,i2,1) - vxc_rho_pert(:,i1,i2,2),DP)/(2*PERT_FAC),DP)
        !TO CHECK THIS
        forall (i2=1:n_spin, ifft=1:fft_size) f1rho(i2,ifft) = real(vxc_rho_pert(ifft,i2,i1,1))
        !forall (i2=1:xc_spin(1)) f1sigma(i2,:) = real(real((vxc_sigma_pert(:,i1,i2,1) - &
        !&                                          vxc_sigma_pert(:,i1,i2,2)),DP)/(2*PERT_FAC_SIGMA),DP)
        !forall (i2=1:xc_spin(1)) f1sigma(i2,:) = real(real((vxc_sigma_pert(:,i1,i2,1) - v1sigma(i2,:)),DP)/(2*PERT_FAC_SIGMA),DP)
        forall (i2=1:xc_spin(1), ifft=1:fft_size) f1sigma(i2,ifft) = real(vxc_sigma_pert(ifft,i1,i2,1))
        call msg('s','Here post f1sigma')
        !f1rho(:,:) = vxc_rho_pert(:,:,is1,1)



      if (i1==1) then 
        open (13, file='vxc_rho_pert_i1.dat', status='unknown')
        write(13,*) 'vxc_rho_pert'
        do i2=1,fft_size
          write(13,*) i2, vxc_rho_pert(i2,1,1,1), vxc_rho_pert(i2,1,2,1)
          write(13,*) i2, vxc_rho_pert(i2,2,1,1), vxc_rho_pert(i2,2,2,1)
          !write(13,*) i2, vxc_rho_pert(i2,i1,1,2), vxc_rho_pert(i2,i1,2,2)
        enddo
        close(13)
      elseif (i1==2) then
        open(14, file='vxc_rho_pert_i2.dat', status='unknown')
        write(14,*) 'vxc_rho_pert'
        do i2=1,fft_size
          write(14,*) i2, vxc_rho_pert(i2,i1,1,1), vxc_rho_pert(i2,i1,2,1)
          !write(14,*) i2, vxc_rho_pert(i2,i1,1,2), vxc_rho_pert(i2,i1,2,2)
        enddo 
        close(14)
      endif

          
        call XC_eval_gga_potential(f1rho,f1sigma,drho,vxc_tmp) 
        !call XC_eval_gga_potential(f1rho,f1sigma,drho,vxc_tmp) 
        !
        !
        !
        !F_xc(:,:,i1) = vxc_tmp(:,:) !tmp_2dim(:,:) vxc_tmp(:,:)
        !
        !
        !
        if (i1==1) then 
          open (13, file='file_f1rho.dat', status='unknown')
          write(13,*) 'F1rho'
          do i2=1,fft_size
            write(13,*) i2, f1rho(1,i2), f1rho(2,i2)
          enddo
          close(13)
  
          open (14, file='file_f1sigma.dat', status='unknown')
          write(14,*) 'F1sigma'
          do i2=1,fft_size
            write(14,*) i2, f1sigma(1,i2), f1sigma(2,i2), f1sigma(3,i2)
          enddo
          close(14)
          open(11, file='vxc_tmp_1.dat', status='unknown')
          write(11,*) 'Vxc_tmp_1'
          do i2=1,fft_size  
            write(11,*) i2, vxc_tmp(i2,1), vxc_tmp(i2,2)
          enddo
          close(11)
        elseif (i1==2) then
          open(12, file='vxc_tmp_2.dat', status='unknown')
          write(12,*) 'Vxc_tmp_2'
          do i2=1,fft_size
            write(12,*) i2, vxc_tmp(i2,1), vxc_tmp(i2,2)
          enddo
          close(12)
        endif 
      enddo
   else 
     call error('Functional not implemented for numerical derivative')
   endif      
end if !if (l_Fxc_from_Vxc.and.ORDER==2) 

! to Delete Start
if (ORDER==2) then
  ! This is a testing statament to print the elements of Fxc
  !     It is not needed for the code
  open(10,file='Fxc.dat',status='unknown')
  write(10,*) 'Fxc'
  do i1=1,fft_size
    write(10,*) i1, F_xc(i1,1,1), F_xc(i1,1,2), F_xc(i1,2,1), F_xc(i1,2,2)
  enddo
  close(10)
endif

open(10,file='spin_rho_DP_final.dat',status='unknown')
write(10,*) "spin_rho_DP"
do ifft=1,fft_size
  write(10,*) ifft, (spin_rho_DP(i1,ifft), i1=1,n_spin)
enddo
close(10)
!to Delete End
call msg('s','Here 7.04')

 if (EvaluateDGradient) then
   YAMBO_FREE(sigma)
   YAMBO_FREE(drho)
 endif
 if (ORDER==1) then
   YAMBO_FREE(tmp_exc)
 endif
 if (EvaluateDGradient.and.(ORDER==1)) then
   YAMBO_FREE(v1rho)
   YAMBO_FREE(v1sigma)
 endif
 if (EvaluateDGradient.and.(ORDER==2)) then
   YAMBO_FREE(v2rho2)
   YAMBO_FREE(v2rhosigma)
   YAMBO_FREE(v2sigma2)
 endif
 if (EvaluateDGradient.and.(FUNCTIONAL==XC_HYB_GGA_XC_GAUPBE*XC_FACTOR)) then
   YAMBO_FREE(vx1)
   YAMBO_FREE(vx2)
 endif
 !
 if (present(EXX_FRACTION)) EXX_FRACTION = local_exx_fraction 
 !
 call msg('s','Here 7.10')
 if (l_Fxc_from_Vxc) then 
   YAMBO_FREE(V_xc_pert)
   YAMBO_FREE(V_xc_pert_cd)
   YAMBO_FREE(spin_rho_SP_cmplx)
   YAMBO_FREE(sigma_pert)
   YAMBO_FREE(v1rho_pert)
   YAMBO_FREE(v1sigma_pert)
   YAMBO_FREE(vxc_rho_pert)
   YAMBO_FREE(vxc_sigma_pert)
   YAMBO_FREE(vxc_tmp)
   YAMBO_FREE(drho_pert)
   YAMBO_FREE(sigma_tmp) !To delete
   YAMBO_FREE(f1rho)
   YAMBO_FREE(f1sigma)  
 endif 
 call msg('s','Here 7.11')
end subroutine XC_libxc_driver

