c-----------------------------------------------------------------------
c     Inside cem_drift_newton.F we use Jacobi-free Newton Krylov method
c     (JFNK) to solve the steady state of PNP equation
c     By Ping Hsuan 06/27/2017
c-----------------------------------------------------------------------
c     This subroutine aim to initial variables for JFNK
c     Those variables are defined in NEWTON
      subroutine cem_drift_newton_init
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'      
      include 'POISSON'
      include 'NEWTON'

      call rzero(sp_k  ,npts)
      call rzero(sn_k  ,npts)
      call rzero(gp_k  ,npts)
      call rzero(gn_k  ,npts)
      call rzero(cpi   ,npts)
      call rzero(cni   ,npts)
      call rzero(cpo   ,npts)
      call rzero(cno   ,npts)
      call rzero(cp_k  ,npts)
      call rzero(cn_k  ,npts)
      call rzero(fp    ,npts)
      call rzero(fn    ,npts)

      return
      end
c---------------------------------------------------------------------
c     This subroutine is the main core of JFNK method which involves
c     using variables cN and cP in DRIFT
      subroutine cem_drift_newton
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'      
      include 'POISSON'
      include 'NEWTON'
      include 'CTIMER'
      include 'RTIMER'

      real alpha, tt
      real new_tol, tol
      real glmax,glsc2,rnorm,rnorm1,rnorm2
      real rinorm1, rinorm2, rinorm
      real fnorm,fnorm1,fnorm2
      integer i, j, k, n
      integer maxit
      real glsum
      real f_pre, f_now

      alpha = 1      ! relaxation parameter \alpha \in (0,1]             
      nsteps = 200   ! steps for pseudo-transient continuation
      maxit = 20     ! iteration for Newton method
      jaceps = 1e-5  ! perturbation parameter for Jacobi-free
      epsinv = 1./jaceps   ! eps for approximation for jacobian
      tol = 1e-2     ! tolerance for GMRES
      new_tol = 1e-5 ! tolerance for Newton method

      dt_newton = param(1)  ! pseudo-transient time step
                            ! it can grow as istep increase 
                            
      dtinv = 1./dt
      dt_newton_inv = 1./dt_newton

      write(*,*)'Jaceps=',jaceps

      cpu_t = 0.0
      cpu_chk = 0.0

      call cem_drift_newton_init

      call copy (cp_k,cP,npts) ! get the initial condition 
      call copy (cn_k,cN,npts) ! get the initial condition

      call copy(cpi,cp_k,npts)
      call copy(cni,cn_k,npts)

      call compute_gk      ! compute rhs g_k for GMRES

      fnorm1 = glsc2(fp,fp,npts)  
      fnorm2 = glsc2(fn,fn,npts)
      fnorm1 = sqrt(fnorm1)
      fnorm2 = sqrt(fnorm2)
      fnorm  = max(fnorm1,fnorm2)

      f_now = fnorm

      call chsign (gp_k,npts)
      call chsign (gn_k,npts)

      call compute_energy  ! check the energy of the system in 
                           ! each pseudo-time step    

      do istep = 1, nsteps 

         cpu_dtime = dclock()

         call rzero(sp_k,npts)   ! set sp_k = 0 
         call rzero(sn_k,npts)   ! set sn_k = 0

         if ((istep .eq. 1)) then
            dt_newton = param(1)
c            dt = dt_newton
         else 
            dt_newton = dt_newton * f_pre/f_now
c            dt = dt_newton
         end if

         dtinv = 1./dt
         dt_newton_inv = 1./dt_newton

c         rinorm1 = glsc2(gp_k,gp_k,npts)
c         rinorm2 = glsc2(gn_k,gn_k,npts)
c         rinorm1 = sqrt(rinorm1)
c         rinorm2 = sqrt(rinorm2)
c         rinorm  = max(rinorm1,rinorm2)
c         write(*,*)'rinorm=',rinorm
         
         do newton_iter = 1, maxit
            
            call drift_hmh_gmres_newton   ! solve sn_k
     $      (sn_k,gn_k,cn_k,cp_k,fn,npts,tol,2)
            call chk_amax('sn_k',sn_k,npts)
            
            call drift_hmh_gmres_newton   ! solve sp_k
     $      (sp_k,gp_k,cp_k,cn_k,fp,npts,tol,1)
            call chk_amax('sp_k',sp_k,npts)
            
            call add2s2(cp_k,sp_k,alpha,npts)   ! cp_k = cp_k + alpha * sp_k
            call add2s2(cn_k,sn_k,alpha,npts)   ! cn_k = cn_k + alpha * sn_k

            write(*,*)'newton_iter=',newton_iter            
            write(*,*)'cpi and cni inside newton'
            call chk_amax  ('cpn',cpi,npts)
            call chk_amax  ('cnn',cni,npts)

            call compute_gk   ! compute rhs g_k for checking tolerance
                              ! and also store for next Newton iteration
            call chsign (gp_k,npts)
            call chsign (gn_k,npts)

            if (newton_iter .eq. 1 ) then

               rinorm1 = glsc2(gp_k,gp_k,npts)
               rinorm2 = glsc2(gn_k,gn_k,npts)
               rinorm1 = sqrt(rinorm1)
               rinorm2 = sqrt(rinorm2)
               rinorm  = max(rinorm1,rinorm2)
               
               if (nid .eq. 0) then
                  write(*,*)'newton',istep,newton_iter,rinorm,rinorm1
     $                              ,rinorm2,dt_newton,dt
               end if 
            else 
               rnorm1 = glsc2(gp_k,gp_k,npts)
               rnorm2 = glsc2(gn_k,gn_k,npts)
               rnorm1 = sqrt(rnorm1)
               rnorm2 = sqrt(rnorm2)
               rnorm  = max(rnorm1,rnorm2)

            if (nid .eq. 0) then
               write(*,*)'newton',istep,newton_iter,rnorm/rinorm,rnorm1
     $                           ,rnorm2,dt_newton,dt
            end if 

c`            if (rnorm .gt. 1e5) call exitt

            if (rnorm/rinorm .lt. new_tol) goto 900

            endif
         enddo
900   continue
      

      call copy(cpi,cp_k,npts)   ! update cp_i for compute gp
      call copy(cni,cn_k,npts)   ! update cn_i for compute gn

      write(*,*)'cpi and cni outside newton'
      call chk_amax('cpo',cpi,npts)
      call chk_amax('cno',cni,npts)

      call copy(cP,cp_k,npts) ! just for cem_out
      call copy(cN,cn_k,npts) ! just for cem_out
      
      fnorm1 = glsc2(fp,fp,npts)  
      fnorm2 = glsc2(fn,fn,npts)
      fnorm1 = sqrt(fnorm1)
      fnorm2 = sqrt(fnorm2)
      fnorm  = max(fnorm1,fnorm2)
      f_pre = f_now
      f_now = fnorm

      cpu_dtime = dclock()-cpu_dtime
      cpu_t = cpu_t+cpu_dtime
      cpu_t_step = cpu_t/istep
      cpu_p_t = glsum(cpu_t_step /npts,1)/np
      
      time = time + dt
  
c      dt_newton = dt_newton * 1.5

c      dt = dt_newton 

      call compute_energy
      call userchk
      call cem_out
      
      enddo
      
      call cem_end

      return

      end
c-----------------------------------------------------------------------
c     This routine computes gk for each Newton iteration
c
c          g^{n}_k = u^{n}_k - \delta t f(u^{n}_k) - u^{n}_0
c
c     Input    cpi,cni denote the initial of each Newton iteration 
c     Input    cp_k,cn_k
c     Output   gp_k,gn_k
      subroutine compute_gk
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'

      call rzero(fp,npts)
      call rzero(fn,npts)

      call compute_f ! compute f(u^{n}_k)
     
      call rzero(gp_k,npts)
      call rzero(gn_k,npts)

      call add2s2(gp_k,cpi,-1.0,npts)
      call add2s2(gn_k,cni,-1.0,npts)

      call add2s2(gp_k,fp,-1.0*dt_newton,npts)
      call add2s2(gn_k,fn,-1.0*dt_newton,npts)

      call add2s2(gp_k,cp_k,1.0,npts)
      call add2s2(gn_k,cn_k,1.0,npts)
      
      write(*,*)'dtnewtoninv_gk',dt_newton_inv
      call cmult(gp_k,dt_newton_inv,npts)   ! case with divided by dt_newton
      call cmult(gn_k,dt_newton_inv,npts)   ! case with divided by dt_newton
      return
      
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      
      call rzero(fp,npts)
      call rzero(fn,npts)

      call copy(cP,cp_k,npts) 
      call copy(cN,cn_k,npts)
      
      if (nid .eq. 0)   write(*,*)'compute_fb_dt',dt 
      call cem_drift_source
      call cem_drift_poisson
      call cem_drift_sem_bdf1
      if (nid .eq. 0)   write(*,*)'compute_fa_dt',dt

      call copy(cpo,cP,npts)
      call copy(cno,cN,npts)

      call sub3(fp,cpo,cp_k,npts)
      call sub3(fn,cno,cn_k,npts)

      call cmult(fp,dtinv,npts)
      call cmult(fn,dtinv,npts)
      if (nid .eq. 0)    write(*,*)'compute_f_dtinv',dtinv
      return
      end
c-----------------------------------------------------------------------
c     This routine computes Jp where J is the Jacobian matrix and p 
c     is a vector. Note that we do not construct the Jacobian matrix 
c     exactly but use the following formula to compute Jp
c   
c        J_k s_k = s_k - (dt /eps) *  ( f(u_k + eps*s_k) - f(u_k) )
c   
c     where f(u_k) has been store in each Newton iteration
c           f(u_k + eps*s_k) has to be computed in each GMRES iteration
      subroutine Jacobi_free(Jp,u0,u1,f,p,iflag)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'

      real     Jp(npts), p(npts)
      real     u0(npts), u1(npts)
      real     tmp(npts), tmp1(npts)
      real     fo(npts), f(npts)
      integer  iflag
      real     eps
      real     pnorm,unorm,glsc2
     
      
      pnorm = glsc2(p,p,npts)
      unorm = glsc2(u0,u0,npts)
      pnorm = sqrt(pnorm)
      unorm = sqrt(unorm)
      eps = (1+unorm)*1e-14
      eps = sqrt(eps)
      eps = eps/pnorm
      epsinv = 1./eps
      write(*,*)'1_epsinv',epsinv,1./eps
      if (nid .eq. 0) write(*,*)istep,newton_iter,'eps=',eps

      call rzero  (Jp,npts)
      call copy   (tmp,u0,npts)
      call add2s2 (tmp,p,eps,npts) !jaceps -> eps
      
      if (iflag .eq. 1) then
      
      call copy   (cP,tmp,npts)
      call copy   (cN,u1 ,npts)
      
      else if (iflag .eq. 2) then
      
      call copy   (cP,u1,npts)
      call copy   (cN,tmp,npts)

      endif
      
      call cem_drift_source
      call cem_drift_poisson
      call cem_drift_sem_bdf1

      if (iflag .eq. 1) then
      call sub3   (fo,cP,tmp,npts)     ! fo = cP - tmp
      call cmult  (fo,dtinv,npts)      ! fo = fo / dt
      call sub3   (tmp1,fo,f,npts)     ! tmp1 = fo - f
      call cmult  (tmp1,epsinv,npts)   ! tmp1 = (fo-f)/eps
      call cmult  (tmp1,dt_newton,npts)       ! tmp1 = (fo-f)*dt_newton/eps
      call sub3   (Jp,p,tmp1,npts)     ! Jp = p - tmp1
      call cmult  (Jp,dt_newton_inv,npts)        ! Jp = p/dt_newton - (fo-f)/eps     
      write(*,*)'2_epsinv',epsinv                  ! case with divided by dt_newton
      else if (iflag .eq. 2) then
      call sub3   (fo,cN,tmp,npts)     ! fo = cN - tmp
      call cmult  (fo,dtinv,npts)      ! fo = fo / dt
      call sub3   (tmp1,fo,f,npts)     ! tmp1 = fo - f
      call cmult  (tmp1,epsinv,npts)   ! tmp1 = (fo-f)/eps
      call cmult  (tmp1,dt_newton,npts)       ! tmp1 = (fo-f)*dt_newton/eps
      call sub3   (Jp,p,tmp1,npts)     ! Jp = p - tmp1
      call cmult  (Jp,dt_newton_inv,npts)        ! Jp = p/dt_newton - (fo-f)/eps
                                       ! case with divided by dt_newton
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_newton
     $           (phi,res,u0,u1,f,n,tol,iflag)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'

      integer  n,outer, iflag
      real     phi(n),res(n),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp
      real     eps, u0(n), u1(n), f(n)
      real*8   etime1,dnekclock

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      
      call rzero(x,n)
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r,res,n)                  ! r = res
         else
            write(*,*)'restart'
            !update residual
            call copy   (r,res,n)                  ! r = res
            call Jacobi_free (w,u0,u1,f,x,iflag)! w = A x
            call add2s2 (r,w,-1.,n) ! r = r - w
         endif

         gamma(1) = glsc2(r,r,n)                   ! gamma  = (r,r)
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         !write(6,*) 'initial residule',gamma(1)
                                                   !      1

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1
            write(*,*)'iter=',iter
            call chk_amax('vj1',v(1,j),npts)   
            call Jacobi_free (w,u0,u1,f,v(1,j),iflag) ! w = A v
            call chk_amax('wjf',w,npts)

c           !modified Gram-Schmidt
            do i=1,j
               call chk_amax('wgs',w,npts)
               call chk_amax('v1i',v(1,i),npts)
               h(i,j)=glsc2(w,v(1,i),n)           ! h    = (w,v )
                                                  ! i,j       i
               write(*,*)'h_ij',h(i,j)
               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            call chk_amax('w2',w,npts)

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            call chk_amax('w3',w,npts)
            alpha = sqrt(glsc2(w,w,n))        ! alpha =  \/ (w,w)
            write(*,*)'alpha=',alpha
            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))

            if ((nid.eq.0).and.(istep.le.5))
     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
            write(*,*)'temp',temp            !  j+1
            call chk_amax('vjf',v(1,j+1),npts)
         enddo
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,v(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

       call copy(phi,x,n)

c      call ortho   (res) ! Orthogonalize wrt null space, if present

      if (nid.eq.0) then
          write(*,*)'istep=',istep
          if ((mod(istep,iocomm).eq.0).or.(istep.le.3))
     $    write(6,9999) istep,newton_iter,iter,tolpss
      endif
 9999 format(' ',' ',i9,i5,' gmres_newton_iteration#',i5,1p3e12.4)

      return
      end
c=======================================================================         
