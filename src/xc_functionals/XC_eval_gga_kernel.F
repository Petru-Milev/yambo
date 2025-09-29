!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine XC_eval_gga_kernel(drho,v2rho2,v2rhosigma,v2sigma2)
 !
 ! Expression to be derived ...
 ! 
 use pars,           ONLY:SP,DP
 use electrons,      ONLY:n_spin
 use FFT_m,          ONLY:fft_size
 use xc_functionals, ONLY:xc_spin, F_xc !F_xc(fft, n_spin, n_spin, n_xc, xyz) 
 !
 implicit none
 !
 real(DP),intent(in) :: drho(fft_size,n_spin,3)
 real(DP),intent(in) :: v2rho2(xc_spin(1),fft_size),v2rhosigma(xc_spin(1)*n_spin,fft_size),v2sigma2(xc_spin(1)*n_spin,fft_size)
 !
 integer :: is1,is2,is_xc,ifft, i 
 !
 if (n_spin==1) then 
   do ifft=1,fft_size
     F_xc(ifft, 1, 1, 1, 1) =  F_xc(ifft, 1, 1, 1, 1) + real(v2rho2(1,ifft),SP) ! LDA Terms  
     do i=1,3
       F_xc(ifft, 1, 1, 2, i) = drho(ifft,1,i) * v2rhosigma(1,ifft) !First gga term  
     end do 
   end do 

 endif 

end subroutine XC_eval_gga_kernel


