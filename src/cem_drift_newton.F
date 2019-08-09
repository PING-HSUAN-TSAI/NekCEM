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
      include 'DRIFT'  ! lcdim
      include 'NEWTON'

      integer ic

      do ic = 1,lcdim
        call rzero(fn  (1,ic),npts)
        call rzero(gn_n(1,ic),npts)
        call rzero(c_k0(1,ic),npts)
        call rzero(cn_k(1,ic),npts)
        call rzero(sn_k(1,ic),npts)
      enddo
      call izero(BDFcnt,lcdim+1)  ! count the BDF1 usage


c     Set default ParamNT
      if (nid.eq.0) write(*,*)'ifparamNT',ifparamNT

      if (ifparamNT) then
        ! turn this on in user file and use parameters set in uservp
      else

        maxit     = 20       ! iteration for Newton method
        dtNT = param(1)      ! pseudo-transient time step
                             ! it can grow as istep increase 

        ifsep     = 0        ! -1=fully couple, 0=couple, 1=separate
        inexactNT_flag = 0   ! 0: fixed 0.1  1: r/r_pre  2: abs(r-res)/r_pre
        dt_flag   = 1        ! 0: fixed dt  1: dt = const * dtNT
        dtNT_flag = 2        ! 1: dtNT=c*dtNT  2:SER  -2:SER + constraint

        alphaNT   = 1.0      ! relaxation parameter of NT solver \alpha \in (0,1]
        jaceps    = 1e-5     ! perturbation parameter for Jacobi-free
        dt_const  = 1e-4     ! constant c for the formula dt = c * dtNT
                             ! c = 1, 0.1, 0.01, 0.001, 0.0001
                             ! c = 0.0001 is the best case up to now

        tolNT     = 1e-5     ! tolerance for Newton method
        tolNTstep = 1e-9     ! to avoid NT loop stalling
        max_dtNT  = 1E20     ! Max. of dtNT, avoid NaN
        min_ftol  = 1E-12    ! Min. of f norm, avoid NaN
        max_dt    = 1E-1     ! Max. of dt, should be satisfied CFL

      endif

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

      real glmax,glsc3,glsum
      real rnorm,rnorms(lcdim),rnorms_pre(lcdim),res_gmres(lcdim)
      real rinorm,rinorms(lcdim),ratio
      real fnorm,fnorms(lcdim)
      real norms_sk(lcdim),norms_ck0(lcdim),sk_ratio(lcdim),NTstep_norm
      real norm_ck0,norm_sk
      real temp(lpts1)

      integer i,j,k,n,ic

      integer irestart ! print restart file for NT
      integer ifirst,ilast

      irestart = param(83)

      call cem_drift_newton_init ! paramNT

c     adjust restart time step
      if (.not.usesteps) then
         if (ifrestart) then
            nsteps = int((fintim-restart_time)/dt)
            if (nid.eq.0) then
               write(6,*) 'fintim',fintim
               write(6,*) 'restart_time',restart_time
               write(6,*) 'fintim-restart_time',fintim-restart_time
            endif
         else
            nsteps = int((fintim)/dt)
            if (nid.eq.0) write(6,*) 'fintim',fintim
         endif
         if (nid.eq.0) write(6,*) 'adjusted nsteps',nsteps
      endif

      ifirst = irststep+1
      ilast = nsteps+irststep

      if (nid.eq.0) then
         write(6,*) '============================'
         write(6,*) '========  Begin Run ========'
         write(6,*) '============================'
         write(6,3) ifirst-1,ilast
      endif
    3 format(' istep =', i10,'  to',i10)

c     the initial step
      istep = ifirst-1
      call userchk
      call cem_out

      if (nsteps.eq.0) call cem_end

c     compute initials
      epsinv = 1./jaceps   ! eps for approximation for jacobian

      call get_dt(dt,max_dt,dtNT,dt_const,dt_flag)
      dtinv = 1./dt
      
      do ic = 1,lcdim
         call copy (cn_k(1,ic),cN(1,ic),npts) ! get the initial condition
         call copy (c_k0(1,ic),cn_k(1,ic),npts)
      enddo

      call compute_f_couple(fn,cn_k,npts)

      fnorm = 0.
      do ic = 1,lcdim
         fnorms(ic) = glsc3(fn(1,ic),fn(1,ic),mult,npts)
c        fnorms(ic) = sqrt(fnorms(ic))
         fnorm = fnorm + fnorms(ic)
      enddo
      fnorm = sqrt(fnorm)
c     fnorm = glmax(fnorms,lcdim)
      f_now = fnorm

c     Start timer
      cpu_t = 0.0
      cpu_chk = 0.0

c     Start pseudo-time step
      istep = ifirst
      do while(.true.)

         cpu_dtime = dclock()

c        determine dtNT, dt
         if ((istep.eq.1)) then
            dtNT = param(1)
         else
            call get_dtNT(dtNT,max_dtNT,f_pre,f_now,min_ftol,dtNT_flag)
            call get_dt(dt,max_dt,dtNT,dt_const,dt_flag)
         endif
         dtNTinv = 1./dtNT
         dtinv = 1./dt

         call compute_cfl(npts)

         if (dt_flag.ne.0) call compute_f_couple(fn,c_k0,npts) ! update f_k0 if dt changed
         call compute_gk(gn_n,fn,cn_k,c_k0,npts) ! compute rhs g_k for GMRES
         do ic = 1,lcdim
            call chsign (gn_n(1,ic),npts)
         enddo
        
c        if (mod(istep,10).eq.0) then
c        if (istep.le.15)then
c           call get_J_new(cn_k,jaceps,dtNTinv,istep,npts,2)
c        endif


c        Start Newton iteration
         do iterNT=1,maxit

            norm_ck0 = 0.
            do ic = 1,lcdim
               norms_ck0(ic) = glsc3(c_k0(1,ic),c_k0(1,ic),mult,npts)
c              norms_ck0(ic) = sqrt(norms_ck0(ic))
               norm_ck0 = norm_ck0 + norms_ck0(ic)
            enddo
            norm_ck0 = sqrt(norm_ck0)

c           determine inexactNT
            do ic = 1,lcdim
               call get_inexactNT(inexactNT(ic),rnorms(ic),
     $         rnorms_pre(ic),res_gmres(ic),iterNT,inexactNT_flag)
            enddo

            if (ifsteric) then
c              call drift_hmh_gmres_NTsteric(sn_k,gn_n,cn_k,fn
c     $           ,mult,lcdim,npts,tolGMRES)
            else
               do ic = 1,lcdim
                  call drift_hmh_gmres_newton(sn_k(1,ic),gn_n(1,ic)
     $               ,cn_k,fn(1,ic),mult,lcdim,npts
     $               ,res_gmres(ic),inexactNT(ic),ic,ifsep)
               enddo
            endif
            do ic = 1,lcdim
               call add2s2(cn_k(1,ic),sn_k(1,ic),alphaNT,npts) ! cp_k=cp_k+alpha*sp_k
            enddo

            call compute_f_couple(fn,cn_k,npts)
            call compute_gk(gn_n,fn,cn_k,c_k0,npts)  ! checking tol. + next NT iter
            do ic = 1,lcdim
               call chsign(gn_n(1,ic),npts)
            enddo

            rnorm = 0.
            do ic = 1,lcdim
               rnorms(ic) = glsc3(gn_n(1,ic),gn_n(1,ic),mult,npts)
               rnorm = rnorm + rnorms(ic)
c              rnorms(ic) = sqrt(rnorms(ic))
            enddo
            rnorm = sqrt(rnorm)
c           rnorm = glmax(rnorms,lcdim)

            if (iterNT.eq.1) then
               rinorm=rnorm
               if (nid.eq.0) write(6,*)istep,iterNT,rinorm,'rinorm'
            endif
            ratio=rnorm/rinorm 

c           check rel. update step size of NT to avoid stalling
            norm_sk = 0.
            do ic = 1,lcdim
               norms_sk(ic) = glsc3(sn_k(1,ic),sn_k(1,ic),mult,npts)
               norm_sk = norm_sk + norms_sk(ic)
c              norms_sk(ic) = sqrt(norms_sk(ic))
c              sk_ratio(ic) = norms_sk(ic) / norms_ck0(ic)
            enddo
            norm_sk = sqrt(norm_sk)
c           NTstep_norm = glmax(sk_ratio,lcdim)
            NTstep_norm = norm_sk/norm_ck0

            if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) write(6,90)
     $         istep,iterNT,ratio,rnorm,rinorm,dtNT,dt,NTstep_norm

            if ((ratio.lt.tolNT).OR.(NTstep_norm.lt.tolNTstep)) goto 900

         enddo ! newton iteration
 90      format('newton iter',2i6,1p6e12.4)   
900      continue

         do ic = 1,lcdim
            call copy(c_k0(1,ic),cn_k(1,ic),npts)   ! update c_k0 for compute gn
         enddo                           
c        for next pseudo-time step
         do ic = 1,lcdim
            call copy(cN(1,ic),cn_k(1,ic),npts) ! just for cem_out
         enddo

c        compute the norm of f(c_p)^{istep-1} and f(c_n)^{istep-1}
         fnorm = 0.
         do ic = 1,lcdim
            fnorms(ic) = glsc3(fn(1,ic),fn(1,ic),mult,npts)  
            fnorm = fnorm + fnorms(ic)
c           fnorms(ic) = sqrt(fnorms(ic))
         enddo
         fnorm = sqrt(fnorm)
c        fnorm  = glmax(fnorms,lcdim)
         f_pre = f_now  ! set up old fnorm to f_pre
         f_now = fnorm  ! set up new fnorm to f_now

         if (nid.eq.0) write(6,*)istep,iterNT,f_pre,f_now,'fnorm'
         if (nid.eq.0) write(6,*)istep,iterNT,ratio,rinorm,rnorm,'gnorm'
         if (nid.eq.0) write(6,*)istep,iterNT,NTstep_norm,'stepnorm'
         if (nid.eq.0) write(6,*)istep,iterNT,dt,dtNT,'dtNT'

c        compute the CPU_time
         cpu_dtime = dclock()-cpu_dtime
         cpu_t = cpu_t+cpu_dtime
         cpu_t_step = cpu_t/istep
         cpu_p_t = glsum(cpu_t_step /npts,1)/np
        
         time = time + dtNT
    
         cpu_chk = dclock()
         call userchk
         cpu_chk = dclock()-cpu_chk

         call cem_out

         if (nid.eq.0.AND.mod(istep,iocomm).eq.0) then ! conut #BDF1
            write(6,*)'newton BDFcnt: istep = ',istep
            write(6,*)'newton BDFcnt: #Poisson solver = ',BDFcnt(1)
            do ic = 1,lcdim
               write(6,991)ic,BDFcnt(ic+1)
            enddo
         endif

         if (irestart.gt.0) then ! write newton.restart
            if (mod(istep,irestart).eq.0) then
               if (nid.eq.0) then
                  write(6,*)'Writing newton.restart! at istep=',istep
                  call write_paramNT
               endif
            endif
         endif 

         istep = istep+1
         if (lastep.eq.1) goto 901
         if (usesteps) then
            if (istep.gt.ilast) goto 901
         else
            if (time.ge.fintim) goto 901
         endif

      enddo ! pseudo-time step

  901 continue
  991 format(' newton BDFcnt: #cN',I2,' solver = ',i10)
      
      call cem_end

      return
      end
c-----------------------------------------------------------------------
c     This routine computes gk for each Newton iteration
c
c          g^{n}_k = u^{n}_k - \delta t f(u^{n}_k) - u^{n}_0
c
c     Input    c_k0 denote the initial of each Newton iteration 
c     Input    cn_k
c     Output   gn_n
      subroutine compute_gk(g_out,f,ck,c0,n)
c-----------------------------------------------------------------------      
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' !lcdim
      include 'NEWTON'

      integer n,ic
      real g_out(lpts1,lcdim),
     $     ck   (lpts1,lcdim),
     $     c0   (lpts1,lcdim),
     $     f    (lpts1,lcdim)

      do ic = 1,lcdim
         call add3s3(g_out(1,ic),ck(1,ic),f(1,ic),c0(1,ic)
     $              ,1.0,-1.0*dtNT,-1.0,n)
         call cmult(g_out(1,ic),dtNTinv,n) ! case with divided by dt_newton
      enddo

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f_sep(f_out,c_in,iflag,n)
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'

      integer n,ic,iflag
      real f_out(lpts1) ! f output
      real c_in (lpts1,lcdim) ! c input
      real c_out(lpts1,lcdim) ! c output
      
      call cem_drift_sem_bdf1_newton(c_in,c_out,n,iflag)

      call sub3(f_out,c_out(1,abs(iflag)),c_in(1,abs(iflag)),n)
      call cmult(f_out,dtinv,n)

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f_couple(f_out,c_in,n)
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'

      integer n,ic
      real f_out(lpts1,lcdim) ! f output
      real c_in (lpts1,lcdim) ! c input
      real c_out(lpts1,lcdim) ! c output

      call cem_drift_sem_bdf1_newton(c_in,c_out,n,0)

      do ic = 1,lcdim
         call sub3(f_out(1,ic),c_out(1,ic),c_in(1,ic),n)
         call cmult(f_out(1,ic),dtinv,n)
      enddo

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
      subroutine JacobiMatVec(Jp,p,uc,fi,n,iflag,isep)
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON'

      real     Jp(lpts1),p(lpts1)
      real     uc(lpts1,lcdim)      ! variable 
      real     fi(lpts1)            ! parameters
      real     uep(lpts1,lcdim)     ! uc + eps p
      real     foi(lpts1),fon(lpts1,lcdim)
      integer  n,iflag,ic,isep
      real     eps,pnorm,unorm,glsc3 

      pnorm = glsc3(p,p,mult,n)
      unorm = glsc3(uc(1,iflag),uc(1,iflag),mult,n)
      pnorm = sqrt(pnorm)
      unorm = sqrt(unorm)

      eps = (1+unorm)*1e-14   ! formula for varing eps
      eps = sqrt(eps)         ! for using the case using  
      eps = eps/pnorm         ! "ub" filled with previous
      epsinv = 1./eps         ! step info first and Dirichelt 
                              ! info on the boundary this
                              ! eps will blow up becasue of 
                              ! pnorm = 0
      do ic = 1,lcdim
         call copy(uep(1,ic),uc(1,ic),n)
      enddo
      call add3s2 (uep(1,iflag),uc(1,iflag),p,1.0,eps,n) ! uep = u + eps*p

      if (isep.eq.0) then
         call compute_f_couple(fon,uep,n)
      elseif (isep.eq.1) then
         call compute_f_sep(fon(1,ic),uep,ic,n)
      endif
      call sub3   (foi,fon(1,iflag),fi,n)

      call cmult  (foi,epsinv*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
      call sub3   (Jp,p,foi,n)
      call cmult  (Jp,dtNTinv,n)      ! Jp = p/dt_newton - (fo-f)/eps
                                      ! case with divided by dt_newton
      return
      end
c-----------------------------------------------------------------------
      subroutine JacobiMatVec_fc(Jp,p,uc,fi,n)
c     Fully-coupled Jacobian
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON'

      real     Jp(lpts1,lcdim),p(lpts1,lcdim)
      real     uc(lpts1,lcdim)      ! variable
      real     fi(lpts1,lcdim)      ! parameters
      real     uep(lpts1,lcdim)     ! uc + eps p
      real     foi(lpts1,lcdim),fon(lpts1,lcdim)
      integer  n,ic
      real     eps,pnorm,unorm,glsc3

      pnorm = 0.0
      unorm = 0.0
      do ic=1,lcdim
         pnorm = pnorm + glsc3(p(1,ic),p(1,ic),mult,n)
         unorm = unorm + glsc3(uc(1,ic),uc(1,ic),mult,n)
      enddo
      pnorm = sqrt(pnorm)
      unorm = sqrt(unorm)

      eps = (1+unorm)*1e-14   ! formula for varing eps
      eps = sqrt(eps)         ! for using the case using
                              ! "ub" filled with previous
                              !  step info first and
                              !  Dirichlet info on the boundary
      eps = eps/pnorm         ! this eps will blow up because of
                              ! pnorm = 0
      epsinv = 1./eps

      do ic = 1,lcdim
         call add3s2(uep(1,ic),uc(1,ic),p(1,ic),1.0,eps,n) ! uep = u + eps*p
      enddo

      call compute_f_couple(fon,uep,n)

      do ic=1,lcdim
         call sub3(foi(1,ic),fon(1,ic),fi(1,ic),n)
      enddo

      do ic=1,lcdim
         call cmult(foi(1,ic),epsinv*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
         call sub3 (Jp(1,ic),p(1,ic),foi(1,ic),n)
         call cmult(Jp(1,ic),dtNTinv,n)      ! Jp = p/dt_newton - (fo-f)/eps
      enddo                                  ! case with divided by dt_newton

      return
      end
c-----------------------------------------------------------------------
      subroutine JacobiMatVec_steric(Jp,p,uc,fi,n)! w = A v 
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON'

      real     Jp(lpts1,lcdim),p(lpts1,lcdim)
      real     uc(lpts1,lcdim)      ! variable 
      real     fi(lpts1,lcdim)      ! parameters
      real     uep(lpts1,lcdim)     ! uc + eps p
      real     foi(lpts1,lcdim),fon(lpts1,lcdim)
      integer  n,nion,iflag,ic,isep
      real     pnorms(lcdim),unorms(lcdim)
      real     epss(lcdim),epssinv(lcdim)
      real     glsc3

      n = lx1*ly1*lz1*nelt
      do ic = 1,lcdim
         pnorms(ic) = glsc3(p(1,ic),p(1,ic),mult,n)
         unorms(ic) = glsc3(uc(1,ic),uc(1,ic),mult,n)
         pnorms(ic) = sqrt(pnorms(ic))
         unorms(ic) = sqrt(unorms(ic))
         epss(ic) = (1+unorms(ic))*1e-14 ! formula for varing eps
         epss(ic) = sqrt(epss(ic))       ! for using the case using 
                                         ! "ub" filled with previous
                                         ! step info first and 
                                         ! Dirichlet info on the boundary
         epss(ic) = epss(ic)/pnorms(ic)! this eps will blow up because of 
                                    ! pnorm = 0
         epssinv(ic) = 1./epss(ic)
      enddo

      do ic = 1,lcdim
         call add3s2 (uep(1,ic),uc(1,ic),p(1,ic),1.0,epss(ic),n) ! uep = u + eps*p
      enddo

      call compute_f_couple(fon,uep,n)

      do ic = 1,lcdim
         call sub3   (foi(1,ic),fon(1,ic),fi(1,ic),n)
         call cmult  (foi(1,ic),epssinv(ic)*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
         call sub3   (Jp(1,ic),p(1,ic),foi(1,ic),n)
         call cmult  (Jp(1,ic),dtNTinv,n)     ! Jp = p/dt_newton - (fo-f)/eps
      enddo                                             ! case with divided by dt_newton
      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_newton
     $           (phi,res,uc,f,wt,nion,n,res_out,forcing,iflag,isep)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'DRIFT'
      include 'NEWTON'
      include 'NGMRES'

      integer  n,noin,outer,iflag,isep
      real     phi(lpts1),res(lpts1),wt(lpts1)
      real     res_out,forcing,tol,alpha,l,temp
      real     eps,uc(lpts1,nion),f(lpts1)
      real*8   etime1,dnekclock

      integer  nion,i,j,k,iconv,iter,m
      real     rnorm,tolpss
      real     glsc3

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

c      tolps = tol
c      tolpss= tolps
      iconv = 0
      
      call rzero(x_NT,n)
      call rzero(h_NT,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r_NT,res,n) ! r = res
         else
            !update residual
            call copy   (r_NT,res,n) ! r = res
            call JacobiMatVec(w_NT,x_NT,uc,f,n,iflag,isep)! w = A x
            call add2s2 (r_NT,w_NT,-1.,n) ! r = r - w
         endif

         gamma_NT(1) = glsc3(r_NT,r_NT,wt,n) ! gamma  = (r,r)
         gamma_NT(1) = sqrt(gamma_NT(1))     ! gamma  = sqrt{ (r,r) }

         tolps = forcing*gamma_NT(1)  ! tolerance for gmres                   
         tolpss = tolps        ! by using inexact Newton method
                               ! 0.1 for forcing term and is changable

         ! check for lucky convergence
         rnorm = 0.
         if(gamma_NT(1).eq.0.) goto 9000
         temp = 1./gamma_NT(1)
         call cmult2(v_NT(1,1),r_NT,temp,n) ! v  = r / gamma
                                            !  1            1
         do j=1,m
            iter = iter+1

            call JacobiMatVec(w_NT,v_NT(1,j),uc,f,n,iflag,isep) ! w = A v

            !modified Gram-Schmidt
            do i=1,j
               h_NT(i,j)=glsc3(w_NT,v_NT(1,i),wt,n)        ! h    = (w,v )
                                                           ! i,j       i
               call add2s2(w_NT,v_NT(1,i),-h_NT(i,j),n)    ! w = w - h    v
            enddo                                          !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_NT(i,j)
               h_NT(i  ,j)=  c_NT(i)*temp + s_NT(i)*h_NT(i+1,j)
               h_NT(i+1,j)= -s_NT(i)*temp + c_NT(i)*h_NT(i+1,j)
            enddo
                                                    !            ______
            alpha = sqrt(glsc3(w_NT,w_NT,wt,n))     ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h_NT(j,j)*h_NT(j,j)+alpha*alpha)
            temp = 1./l
            c_NT(j) = h_NT(j,j) * temp
            s_NT(j) = alpha  * temp
            h_NT(j,j) = l
            gamma_NT(j+1) = -s_NT(j) * gamma_NT(j)
            gamma_NT(j)   =  c_NT(j) * gamma_NT(j)

            rnorm = abs(gamma_NT(j+1))

c            if ((nid.eq.0).and.(istep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_NT(1,j+1),w_NT,temp,n)   ! v    = w / alpha
                                                   !  j+1
         enddo
  900    iconv = 1
         res_out = rnorm   ! output |Ax-b|
         if (nid.eq.0) write(6,*) istep,iterNT,
     $                 'forcing_term',rnorm/tolpss
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma_NT(k)
            do i=j,k+1,-1
               temp = temp - h_NT(k,i)*c_NT(i)
            enddo
            c_NT(k) = temp/h_NT(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_NT,v_NT(1,i),c_NT(i),n)     ! x = x + c  z
         enddo                               !          i  i
      enddo
 9000 continue

      call copy(phi,x_NT,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
          write(6,9999) istep,iterNT,iter,tolpss
      endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_newton_fc
     $           (phi_fc,res_fc,uc,f_fc,wt,nion,n,res_out,forcing)
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
c     fully coupled
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'DRIFT'
      include 'NEWTON'
      include 'NMGMRES'

      real     phi_fc(lpts1,nion),res_fc(lpts1,nion),wt(lpts1)
      real     eps,uc(lpts1,nion),f_fc(lpts1,nion)
      real     res_out,forcing,tol,alpha,l,temp
      real     rnorm,tolpss
      real     glsc3
      real*8   etime1,dnekclock

      integer  i,j,k,ic,iconv,iter,m
      integer  n,nion,outer

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

c     tolps = tol
c     tolpss= tolps
      iconv = 0
      
      call rzero(x_NT,lpts1*(LDIMT-1))
      call rzero(h_NT,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            do ic = 1,nion
               call copy(r_NT(1,ic),res_fc(1,ic),n) ! r = res
            enddo
         else
            ! update residual
            do ic=1,nion
               call copy(r_NT(1,ic),res_fc(1,ic),n) ! r = res
            enddo
               call JacobiMatVec_fc(w_NT,x_NT,uc,f_fc,n) ! w = A x
            do ic=1,nion
               call add2s2 (r_NT(1,ic),w_NT(1,ic),-1.,n) ! r = r - w
            enddo
         endif

         gamma_NT(1) = 0.
         do ic=1,nion
            gamma_NT(1) = gamma_NT(1)
     $                  + glsc3(r_NT(1,ic),r_NT(1,ic),wt,n) ! gamma  = (r,r)
         enddo
         gamma_NT(1) = sqrt(gamma_NT(1)) ! gamma  = sqrt{ (r,r) }

         tolps = forcing*gamma_NT(1) ! tolerance for gmres                   
         tolpss = tolps              ! by using inexact Newton method
                                     ! 0.1 for forcing term and is changable
c        tolps = 1E-10 ! exact Newton,Lan

         ! check for lucky convergence
         rnorm = 0.
         if(gamma_NT(1) .eq. 0.) goto 9000

         temp = 1./gamma_NT(1)
         do ic=1,nion
            call cmult2(v_NT(1,ic,1),r_NT(1,ic),temp,n) ! v  = r / gamma
         enddo                                          !  1            1
                                                  
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call JacobiMatVec_fc(w_NT,v_NT(1,1,j),uc,f_fc,n) ! w = A v

            ! modified Gram-Schmidt
            do i=1,j

               h_NT(i,j)=0.0
               do ic=1,nion
                  h_NT(i,j)=h_NT(i,j)
     $                     + glsc3(w_NT(1,ic),v_NT(1,ic,i),wt,n) ! h    = (w,v )
               enddo                                             !  i,j         i

               do ic=1,nion
                  call add2s2(w_NT(1,ic),v_NT(1,ic,i),-h_NT(i,j),n) ! w = w - h    v
               enddo                                                !          i,j  i
            enddo                                 

            ! apply Givens rotations to new column
            do i=1,j-1
               temp = h_NT(i,j)
               h_NT(i  ,j)=  c_NT(i)*temp + s_NT(i)*h_NT(i+1,j)
               h_NT(i+1,j)= -s_NT(i)*temp + c_NT(i)*h_NT(i+1,j)
            enddo
            alpha = 0.0

            do ic=1,nion                                             !            ______
               alpha = alpha + glsc3(w_NT(1,ic),w_NT(1,ic),wt,n)     ! alpha =  \/ (w,w)
            enddo
            alpha = sqrt(alpha)

            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h_NT(j,j)*h_NT(j,j)+alpha*alpha)
            temp = 1./l
            c_NT(j) = h_NT(j,j) * temp
            s_NT(j) = alpha  * temp
            h_NT(j,j) = l
            gamma_NT(j+1) = -s_NT(j) * gamma_NT(j)
            gamma_NT(j)   =  c_NT(j) * gamma_NT(j)

            rnorm = abs(gamma_NT(j+1))

c            if ((nid.eq.0).and.(istep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            do ic=1,nion
               call cmult2(v_NT(1,ic,j+1),w_NT(1,ic),temp,n)   ! v    = w / alpha
            enddo                                              !  j+1
         enddo

c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
         res_out = rnorm   ! output |Ax-b|
         if (nid.eq.0) write(6,*) istep,iterNT,
     $                 rnorm/tolpss,'forcing_term'
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma_NT(k)
            do i=j,k+1,-1
               temp = temp - h_NT(k,i)*c_NT(i)
            enddo
            c_NT(k) = temp/h_NT(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            do ic=1,nion
               call add2s2(x_NT(1,ic),v_NT(1,ic,i),c_NT(i),n) ! x = x + c  z
            enddo                                             !          i  i
         enddo                               
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      do ic=1,nion
         call copy(phi_fc(1,ic),x_NT(1,ic),n)
      enddo

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
          write(6,9999) istep,iterNT,iter,tolpss
      endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end
c------------------------------------------------------------------------
c      subroutine drift_hmh_gmres_NTsteric(uC,fC,uref,fi,wt,nion,n,tol)
cc-----------------------------------------------------------------------
cc     For Steric term, involving cross term linear system 
cc      Hii = B + dt(1 + ste_aii)*A_Di
cc      Hij = dt*ste_aij*A_Di
cc    
cc      H = (Hij), f = [f1;f2;...fp], u = [u1;u2;...;up], p = # of species
cc      
cc     solving H*u = f
c      implicit none
c      include 'SIZE'
c      include 'TOTAL'
c      include 'FDMH1'
c      include 'DRIFT'
c      include 'MGMRES'
cC
c      integer  nion,n,outer
c      integer  iter,iconv
c      integer  i,j,k,m
c      integer  ic ! safe index fo nion
c      real     divex,etime_p,tolpss
c
c      real     uC(lpts1,lcdim)
c     $       , fC(lpts1,lcdim)
c     $       , mask(lpts1,lcdim)
c     $       , wt(lpts1)
c     $       , fi(lpts1,lcdim)
c     $       , uref(lpts1,lcdim)
c      real     tol,alpha,l,temp,rnorm
c      real     glsc3
c      real*8   etime1,dnekclock
cC
c
c      if (nid.eq.0) write(6,*) 'start: hmh_gmres_newton_steric'
c
c      n = nx1*ny1*nz1*nelv
c
c      etime1 = dnekclock()
c      etime_p = 0.
c      divex = 0.
c      iter  = 0
c      m     = lgmres
c
c      tolps = tol
c      tolpss= tolps
c      iconv = 0
c      call rzero(xxC,lpts*(LDIMT-1))
c      call rzero(h,m*m)
c
c      outer = 0
c      do while (iconv.eq.0.and.iter.lt.5000)
c         outer = outer+1
c         if(iter.eq.0) then
c            do ic=1,nion
c            call copy  (rC(1,ic),fC(1,ic),n)      ! r = res
c            enddo
c         else
c            !update residual
c            do ic=1,nion
c            call copy  (rC(1,ic),fC(1,ic),n)      ! r = res
c            enddo
c            call JacobiMatVec_steric(wC,xxC,uref,fi,n)! w = A x
c
c            do ic=1,nion
c            call dssum  (wC(1,ic),nx1,ny1,nz1)
c            call col2   (wC(1,ic),mask(1,ic),n)
c            enddo
c            do ic=1,nion
c            call add2s2 (rC(1,ic),wC(1,ic),-1.,n)   ! r = r - w
c            enddo
c         endif
c         gamma(1) = 0.0
c         do ic=1,nion
c           gamma(1) = gamma(1)+glsc3(rC(1,ic),rC(1,ic),wt,n)! gamma  = (r,r)
c         enddo
c         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
c         write(6,*) 'initial residule',gamma(1)
c                                                   !      1     include 'SIZE'
c
c         !check for lucky convergence
c         rnorm = 0.
c         if(gamma(1) .eq. 0.) goto 9000
c         temp = 1./gamma(1)
c         do ic=1,nion
c         call cmult2(vC(1,ic,1),rC(1,ic),temp,n)            !  v  = r / gamma
c         enddo
c         !write(6,*) 'start form m-th krylov subspace'
c         do j=1,m
c            iter = iter+1
c            call JacobiMatVec_steric(wC,vC(1,1,j),uref,fi,n)! w = A v 
c            do ic=1,nion
c            call dssum  (wC(1,ic),nx1,ny1,nz1)
c            call col2   (wC(1,ic),mask(1,ic),n)
c            enddo
c
cc           !modified Gram-Schmidt
c            do i=1,j
c               h(i,j)=0.0
c               do ic=1,nion
c               h(i,j) = h(i,j) + glsc3(wC(1,ic),vC(1,ic,i),wt,n) ! h    = (w,v )
c               enddo                                             ! i,j       i
c               do ic=1,nion
c               call add2s2(wC(1,ic),vC(1,ic,i),-h(i,j),n)   ! w = w - h    v
c               enddo                                        !         i,j  i
c            enddo
c
c            !apply Givens rotations to new column
c            do i=1,j-1
c               temp = h(i,j)
c               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
c               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
c            enddo
c            alpha=0.0
c            do ic=1,nion
c            alpha = alpha + glsc3(wC(1,ic),wC(1,ic),wt,n)  !            ______
c            enddo                                          ! alpha =  \/ (w,w)
c            alpha=sqrt(alpha)
c
c            if(alpha.eq.0.) goto 900  !converged
c            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
c            temp = 1./l
c            c(j) = h(j,j) * temp
c            s(j) = alpha  * temp
c            h(j,j) = l
c            gamma(j+1) = -s(j) * gamma(j)
c            gamma(j)   =  c(j) * gamma(j)
c
c            rnorm = abs(gamma(j+1))!*norm_fac     
cc           ratio = rnorm/div0
c            if ((nid.eq.0).and.(istep.le.5))
c     $           write (6,66) iter,tolpss,rnorm,istep
c   66       format(i5,1p2e12.5,i8,' Divergence')
c
c            if (rnorm .lt. tolps) goto 900 !converged
c            if (j.eq.m) goto 1000 !not converged, restart
c
c            temp = 1./alpha
c            do ic=1,nion
c            call cmult2(vC(1,ic,j+1),wC(1,ic),temp,n)   ! v    = w / alpha
c            enddo                                       !  j+1
c         enddo
cc        write(6,*) 'end of forming m-th krylov subspace'
c  900    iconv = 1
c 1000    continue
c
cc        back substitution
cc             -1
cc        c = H   gamma
cc        write(6,*) 'start solving least squre problem'
c         do k=j,1,-1
c            temp = gamma(k)
c            do i=j,k+1,-1
c               temp = temp - h(k,i)*c(i)
c            enddo
c            c(k) = temp/h(k,k)
c         enddo
c         !sum up Arnoldi vectors
c         do i=1,j
c            do ic=1,nion
c            call add2s2(xxC(1,ic),vC(1,ic,i),c(i),n)  ! x = x + c  z
c            enddo                                      !  i  i
c         enddo
cc        write(6,*) 'end of solving least squre problem'
c      enddo
c 9000 continue
c
c       do ic=1,nion
c       call copy(uC(1,ic),xxC(1,ic),n)
c       enddo
c
cc      call ortho   (res) ! Orthogonalize wrt null space, if present
c
c      etime1 = dnekclock()-etime1
c      if (nid.eq.0) then
c          if ((mod(istep,iocomm).eq.0).or.(istep.le.10))
c     $    write(6,9999) istep,iter,tolpss,etime_p,etime1
c      endif
c 9999 format(' ',' ',i9,' gmres   : iteration#',i5,1p3e12.4)
c
c      return
c      end
c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_newton(cinput,coutput,n,iflag)
c     Three mode: (switch iflag)
c      0 couple:        update all the fields: cN(s) and potent
c      i separate:      update potent and one of cN(i)
c     -i fully-separate:update potent(-1) or cN(-2) or cP(-3)
c---------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'
      include 'NEWTON' !BDFcnt

      integer iflag,n,ic
      real cinput(lpts1,lcdim)
      real coutput(lpts1,lcdim)

      do ic = 1,lcdim
         call copy(cN(1,ic),cinput(1,ic),n)
      enddo

      if (iflag.eq.0) then ! coupling every thing
         call cem_drift_op_bdf
         do ic = 1,lcdim+1
            BDFcnt(ic) = BDFcnt(ic) + 1
         enddo
      else ! seperate
         call cem_drift_source
         call cem_drift_poisson
         BDFcnt(1) = BDFcnt(1) + 1

         call cem_drift_sem_bdf1_init

         call cem_drift_rhs_bdf1(iflag)
         call cem_drift_lhs_bdf1(iflag)
         BDFcnt(iflag+1) = BDFcnt(iflag+1) + 1
      endif

      do ic = 1,lcdim
         call copy(coutput(1,ic),cN(1,ic),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_energy !FIXME: restr later, Lan
c-----------------------------------------------------------------------
c      implicit none
c      include 'SIZE'
c      include 'TOTAL'
c      include 'DRIFT'
c      include 'POISSON'
c
c      real energy_density(npts),energy
c      integer i,j,k
c      real dummy1(npts),dummy2(npts),dummy3(npts)
c      real dummy4(npts),dummy5(npts)
c      logical :: lexist
c      real glsum
c
c      call rzero(energy_density,npts)
c      call rzero(dummy1,npts)
c      call rzero(dummy2,npts)
c      call rzero(dummy3,npts)
c
c      call copy(dummy1,cN,npts)
c      call copy(dummy2,cP,npts)
c      call copy(dummy3,cP,npts)
c
c      call add2s2(dummy3,cN,-1.0,npts)
c      call add2s2(dummy3,ni, 1.0,npts)
c
c      call col2(dummy3,potent,npts)
c      call col2(dummy1,LOG(cN),npts)
c      call col2(dummy2,LOG(cP),npts)
c
c      call add2s2(energy_density,dummy1,1.0,npts)
c      call add2s2(energy_density,dummy2,1.0,npts)
c      call add2s2(energy_density,dummy3,0.5,npts)
c
c      call col2(energy_density,bmn,npts)
c      energy =  glsum(energy_density,npts)
c
c      if ((nid.eq.0).and.(mod(istep,iocomm).eq.0))
c     $    write(6,10) energy,time
c 10   format('     energy: ',1p1e18.7, ' at time=',1p1e12.4)
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dtNT(dtNT,max_dtNT,f_pre,f_now,min_ftol,iflag)
c     This routine deals with the selection for dtNT
c     - SER
c        dtNT = dtNT * f_pre/f_now
c     - fixed increasing rate 
c        dtNt = dtNT * c
c     input: all, output: dtNT
c-----------------------------------------------------------------------
      implicit none

      real dtNT,max_dtNT,f_pre,f_now,min_ftol
      real const, ratio
      integer iflag

      if (iflag.eq.1) then ! fixed ratio

         dtNT = dtNT * 10

      elseif (iflag.eq.2) then ! SER

         if (f_now.lt.min_ftol) then
            dtNT = max_dtNT
         else
            ratio= f_pre/ f_now
            dtNT = dtNT * ratio
         endif

      elseif (iflag.eq.-2) then ! SER with constraint

         if (f_now.lt.min_ftol) then ! near steady state
            dtNT = max_dtNT
         else
            ratio= f_pre/ f_now
            if (ratio.gt.1E2) ratio = 1E2  ! max rate
  
            dtNT = dtNT * ratio
            if (dtNT.lt.1E-4) dtNT = 1E-4  ! min control
         endif

      endif

      if (dtNT.gt.max_dtNT) dtNT = dtNT

      return
      end
c-----------------------------------------------------------------------
c     This routine deals with the relation between dtNT and dt
c     input: all, output: dt
      subroutine get_dt(dt,max_dt,dtNT,const,iflag)
c-----------------------------------------------------------------------
      implicit none

      real dt,max_dt,dtNT,const
      integer iflag

      if (iflag.eq.0) then
c     do nothing and use original dt
      elseif (iflag.eq.1) then
         dt = dtNT * const
      endif

      if (dt.gt.max_dt) dt = max_dt

      return
      end
c-----------------------------------------------------------------------
c     This routine computes CFL condition for Advection Diffusion equation
c     For now, we only monitor this value in case JFNK diverges
c     Based on the following formulas:
c       - diffusion term
c         C * dt / dx   < 1,         C: max. abs. of drift 
c       - drift term
c         D * dt / dx^2 < 1/2,       D: diffusion coefficient
      subroutine compute_cfl(n) ! ToDO: should be integrated into drift, Lan
c------------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      real dxmin,cond1,cond2
      real tmp1,tmp2,tmp3,tmp4,tmp5,tmpf1,tmpf2
      real glmax,glamax
      integer n

      call get_dxmin(dxmin)
      tmp1 = glamax(gradp(1,1),n)
      tmp2 = glamax(gradp(1,2),n)
      tmp3 = glamax(gradp(1,3),n)
      tmpf1 = max(tmp1,tmp2)
      tmpf1 = max(tmp3,tmpf1)
      write(6,*) 'maximum of electric field',tmpf1

      tmp4 = glmax(mu_n(1,1),n)
      tmp5 = glmax(mu_n(1,2),n)
      tmpf2 = max(tmp4,tmp5)
      write(6,*) 'maximum of diffusion coefficient',tmpf2

      cond1 = dxmin / tmpf1
      cond2 = dxmin**2/(2*tmpf2)
      write(6,*) 'delta t should not greater than',cond1,cond2

      return
      end
c------------------------------------------------------------------------
c     forcing term is usually fixed in inexact Newton method.
c     It can be the trade-off between tol in GMRES and #iter in Newton.
c     input: all, output: forcing, rnorm_pre
      subroutine get_inexactNT(forcing,rnorm,rnorm_pre,res_gmres
     $                        ,iterNT,flag)
c------------------------------------------------------------------------
      implicit none

      real forcing,rnorm,rnorm_pre,res_gmres
      integer iterNT,flag

      if (flag.eq.0) then     ! fixed forcing term

         forcing = 0.1         ! can be 0.1 0.01 0.001 0.0001
         rnorm_pre = rnorm

      elseif (flag.eq.1) then ! update by the ratio of rnorm

         if (iterNT.eq.1) then
            forcing = 0.1
         else 
            if (rnorm/rnorm_pre.ge.1) then
               forcing = 0.1
            else
               forcing = rnorm/rnorm_pre
            endif
         endif
         rnorm_pre = rnorm

      elseif (flag.eq.2) then ! update by relative change of rnorm and res_gmres

         if (iterNT.eq.1) then
            forcing = 0.1
         else
            forcing = abs(rnorm-res_gmres)/rnorm_pre
         endif
         rnorm_pre = rnorm

      endif

      return
      end
c------------------------------------------------------------------------
      subroutine read_paramNT
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'

      real paramNT(20)
      integer nn,i

C     Read paramNT

      open(unit=11,file='newton.restart')
      nn = 0
      do i=1,20
         read(11,*,end=99)paramNT(i)
         nn=nn+1
      enddo
   99 continue

C     Set paramNT

      ! 1~7 Pseudo-timestep
      ifsep    = int(paramNT(3))   ! (int) 0=couple, 1=separate, Lan
      dtNT     = paramNT(4)        ! pseudo-transient time step
      max_dtNT = paramNT(5)        ! Max. of dtNT, avoid NaN
      min_ftol = paramNT(6)        ! Min. of f_pre, avoid NaN
      f_pre    = paramNT(7)        ! zero to set dtNT (or SER)

      ! 8~12 Newton iter
      alphaNT  = paramNT(8)        ! relaxation parameter \alpha \in (0,1]
      maxit    = int(paramNT(9))   ! (int) iteration for Newton method
      jaceps   = paramNT(10)       ! perturbation parameter for Jacobi-free
      tolNT    = paramNT(11)       ! tolerance for Newton method

      ! 13~20 else

      close(11)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_paramNT
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'

      integer i
      real paramNT(20)

      open(unit=11,file='newton.restart')

      paramNT(1)  = 0
      paramNT(2)  = 0
      paramNT(3)  = ifsep
      paramNT(4)  = dtNT
      paramNT(5)  = max_dtNT
      paramNT(6)  = min_ftol
      paramNT(7)  = f_pre
      paramNT(8)  = alphaNT
      paramNT(9)  = maxit
      paramNT(10) = jaceps
      paramNT(11) = tolNT
      paramNT(12) = 0
      paramNT(14) = 0
      paramNT(15) = 0
      paramNT(16) = 0
      paramNT(17) = 0
      paramNT(18) = 0
      paramNT(19) = 0
      paramNT(20) = 0

      write(11,3)int(paramNT(1)),'1: '
      write(11,3)int(paramNT(2)),'2: '
      write(11,3)int(paramNT(3)),
     $           '3: ifsep,(int) 0 = couple, 1 = seperate'
      write(11,7)paramNT(4),'4: dtNT, pseudo-timestep tau'
      write(11,7)paramNT(5),'5: max_dtNT, Max. of dtNT, avoid NaN'
      write(11,7)paramNT(6),'6: min_ftol, Min. of f norm,  avoid NaN'
      write(11,7)paramNT(7),'7: f_pre, f norm from previous step'
      write(11,7)paramNT(8),
     $           '8: alpha, relaxation parameter \alpha \in (0,1]'
      write(11,3)int(paramNT(9)),
     $           '9: maxit, (int) iteration for Newton method'
      write(11,7)paramNT(10),
     $           '10: jaceps, perturbation parameter for Jacobi-free'
      write(11,7)paramNT(11),'11: tolNT, tolerance for Newton method'
      write(11,7)paramNT(12),'12: '
      write(11,7)paramNT(13),'13: tolGMRES, tolerance for GMRES'
      write(11,7)paramNT(14),'14: '
      write(11,7)paramNT(15),'15: '
      write(11,7)paramNT(16),'16: '
      write(11,7)paramNT(17),'17: '
      write(11,7)paramNT(18),'18: '
      write(11,7)paramNT(19),'19: '
      write(11,7)paramNT(20),'20: '

    3 format(I4,'              ',A)
    7 format(E13.4'     ',A)

      close(11)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_J(c_in,pert,invtau,pstep,nn,ncdim)

      implicit none                         

      include 'SIZE'                       
      include 'TOTAL'                     
      include 'DRIFT'                    
                                       
      real J_out(nn*nn,ncdim)          
      real fi(lpts1,ncdim),ff(lpts1,ncdim)  
      real fdiff(lpts1,ncdim)              
      real c_in(lpts1,ncdim) ! c input    
      real c_out(lpts1,ncdim) ! c output  
      real pert, invtau                    
      real p(lpts1,ncdim)                  
      real w1(nn,ncdim)                 
      real pertinv                         
      integer ncdim,nn
      integer pstep                        
      integer ic,i,j,k,k1                  
      integer e,ee
      integer re,ii,jj,global_n,nglob
      real flag(lx1*ly1*lz1,lelt)
      character(20) filename1, filename2   
                                           
                                           
      pertinv = 1./pert                    
      write(6,*)'====================='    
      write(6,*)'Begin computing J1'        
      write(6,*)'Parameter check'          
      write(6,*)'pert',pert,'pertinv',pertinv
      write(6,*)'dtNTinv',invtau            
      write(6,*)'pstep',pstep               
      write(6,*)'nn',nn,lpts1           
      write(6,*)'====================='     

      call compute_f_couple(fi,c_in,nn)   

      k1 = 0                                
      do j=1,nn                           
         k1 =  k1+1                         
         k  = (k1-1)*nn+1                 
         do ic = 1,ncdim                    
            call rzero(p(1,ic),nn)        
         enddo                              
         p(j,1) = 1                         
         call dssum(p(1,1),nx1,ny1,nz1) 
         do ic = 1,ncdim                    
         call add3s2(c_out(1,ic),c_in(1,ic),p(1,ic),1.0,pert,nn)
         enddo                                                   
         call compute_f_couple(ff,c_out,nn)                    
         do ic = 1,ncdim                                         
c           call cmult(ff(1,ic),invtau,nn)                     
            call sub3(fdiff(1,ic),ff(1,ic),fi(1,ic),nn)        
            call cmult(fdiff(1,ic),pertinv,nn)
            call cmult(p(1,ic),invtau,nn)
         enddo                                                   
         call sub3(J_out(k,1),p(1,1),fdiff(1,1),nn)
c        call add3s2(J_out(k,1),p(1,1),fdiff(1,1),invtau,-1.0,nn)
      enddo                                                            

      write(6,*)'first Jacobian is done'                               

      write(6,*)'====================='    
      write(6,*)'Begin computing J2'        
      write(6,*)'Parameter check'          
      write(6,*)'pert',pert,'pertinv',pertinv
      write(6,*)'dtNTinv',invtau            
      write(6,*)'pstep',pstep               
      write(6,*)'nn',nn,lpts1           
      write(6,*)'====================='     
                                                                       
      k1 = 0                                                           
      do j=1,nn                                                      
         k1 =  k1+1                                                    
         k  = (k1-1)*nn+1                                            
         do ic = 1,ncdim                                               
            call rzero(p(1,ic),nn)                                   
         enddo                                                         
         p(j,2) = 1                                                    
         call dssum(p(1,2),nx1,ny1,nz1) 
        
         do ic = 1,ncdim                                               
         call add3s2(c_out(1,ic),c_in(1,ic),p(1,ic),1.0,pert,nn)     
         enddo                                                     
         call compute_f_couple(ff,c_out,nn)                      
         do ic = 1,ncdim                                           
c           call cmult(ff(1,ic),invtau,nn)                       
            call sub3(fdiff(1,ic),ff(1,ic),fi(1,ic),nn)          
         enddo                                                     
         call add3s2(J_out(k,2),p(1,2),fdiff(1,2),invtau,          
     $         -pertinv,nn)                                      
      enddo                                                        
      write(6,*)'second Jacobian is done'                          
                                                                   
      ! print out J                                                
      write(filename1,'(a,i0,a)') 'Jacobian1_',pstep,'.dat'        
      write(filename2,'(a,i0,a)') 'Jacobian2_',pstep,'.dat'        
                                                                   
      open(50,file=filename1)                                      
         do i = 1,nn**2                                          
            write(50,*)J_out(i,1)                                  
         enddo                                                     
      close(50)                                                    
      open(50,file=filename2)                                      
         do i = 1,nn**2                                          
            write(50,*)J_out(i,2)                                  
         enddo                                                     
      close(50)                                                    
                                                                   
      end           
c-----------------------------------------------------------------------
      subroutine get_J_new(c_in,pert,invtau,pstep,nn,ncdim)

      implicit none                         

      include 'SIZE'                       
      include 'TOTAL'                     
      include 'DRIFT'                    
      include 'POISSON'
                                       
      real J_out(lpts1,ncdim)          
      real fi(lpts1,ncdim),ff(lpts1,ncdim),fdiff(lpts1,ncdim)  
      real c_in(lpts1,ncdim),c_out(lpts1,ncdim) 
      real p(lpts1,ncdim)                  
      real flag(lx1*ly1*lz1,lelt)
      real uin(lx1*ly1*lz1,lelt),uout(lx1*ly1*lz1,lelt)
      real pert,invtau,pertinv                    
      real glsc2
      integer ncdim,nn,pstep
      integer ic,i,j,k,k1                  
      integer e,ee
      integer re,ii,jj,global_n,nglob
      character(20) filename1, filename2   
                                           
      pertinv = 1./pert                    
      write(6,*)'====================='    
      write(6,*)'Begin computing J1'        
      write(6,*)'Parameter check'          
      write(6,*)'pert',pert,'pertinv',pertinv
      write(6,*)'dtNTinv',invtau            
      write(6,*)'pstep',pstep               
      write(6,*)'nn',nn,lpts1           
      write(6,*)'====================='     

      call compute_f_couple(fi,c_in,nn)   

      nxyz = lx1*ly1*lz1
      write(6,*)'nxyz',nxyz

      global_n = glsc2(dmask(1,1),mult,nn)
      write(6,*)'global_n',global_n
      nglob    = global_n + 0.4          ! Rounding

      do e=1,nelv
         re = e 
         call cfill(flag(1,e),float(re),nxyz)
      enddo
      call dsop(flag,'min',nx1,ny1,nz1)
                                            
      j=0
      do e=1,nelv
      do jj=1,nxyz
         re = e
         ! Found a real dof
         if (int(flag(jj,e)).eq.re.and.dmask((e-1)*nxyz+jj,1).ne.0) 
     $   then 
            j = j+1
            flag(jj,e) = j
         else
            flag(jj,e) = 0
         endif
      enddo
      enddo

      write(filename1,'(a,i0,a)') 'Jacobian1_',pstep,'.dat'        
      open(50,file=filename1)                                      
      j = 0
      write(6,*)'check degree',(nelv-1)*nxyz+nxyz,nn
      do e=1,nelv
      do jj=1,nxyz

         do ic = 1,ncdim                                               
            call rzero(p(1,ic),nn)                                   
         enddo                                                         

         if (int(flag(jj,e)).ne.0) then  ! Next column of identity
            j=j+1
            p((e-1)*nxyz+jj,1) = 1                                    
            call dssum  (p(1,1) ,nx1,ny1,nz1) ! Copy "1" to all elements
            call col2   (p(1,1),mult,nn) ! do the average

            do ic = 1,ncdim                             
               do ii=1,nn
                  c_out(ii,ic) = c_in(ii,ic)+pert*p(ii,ic)
               enddo
c              call add3s2(c_out(1,ic),c_in(1,ic),p(1,ic),1.0,pert,nn) 
            enddo                                                     

            call compute_f_couple(ff,c_out,nn)                      
            do ic = 1,ncdim                                           
               call sub3(fdiff(1,ic),ff(1,ic),fi(1,ic),nn)          
            enddo                                                     
            do ii=1,nn
               J_out(ii,1) = invtau*p(ii,1)-pertinv*fdiff(ii,1)
            enddo
c           call add3s2(J_out(1,1),p(1,1),fdiff(1,1),invtau,          
c    $           -pertinv,nn)                                      
            call dssum  (J_out(1,1),nx1,ny1,nz1) ! Copy "1" to all elements

c           Now, output the "j"th column of A
            

            i = 0
            do ee=1,nelv
            do ii=1,nxyz
               if (int(flag(ii,ee)).ne.0) then 
                  write(50,77)J_out((ee-1)*nxyz+ii,1)
               endif
            enddo
            enddo
         endif

      enddo
      enddo

      write(filename2,'(a,i0,a)') 'Jacobian2_',pstep,'.dat'        
      open(50,file=filename2)                                      
      j = 0
      do e=1,nelv
      do jj=1,nxyz

         do ic = 1,ncdim                                               
            call rzero(p(1,ic),nn)                                   
         enddo                                                         

         if (int(flag(jj,e)).ne.0) then  ! Next column of identity
            j=j+1
            p((e-1)*nxyz+jj,2) = 1                                                    
            call dssum  (p(1,2) ,nx1,ny1,nz1) ! Copy "1" to all elements
            call col2   (p(1,2),mult,nn) ! do the average
            do ic = 1,ncdim                             
               do ii=1,nn
                  c_out(ii,ic) = c_in(ii,ic)+pert*p(ii,ic)
               enddo
c              call add3s2(c_out(1,ic),c_in(1,ic),p(1,ic),1.0,pert,nn) 
            enddo                                                     

            call compute_f_couple(ff,c_out,nn)                      
            do ic = 1,ncdim                                           
               call sub3(fdiff(1,ic),ff(1,ic),fi(1,ic),nn)          
            enddo                                                     
            do ii=1,nn
               J_out(ii,2) = invtau*p(ii,2)-pertinv*fdiff(ii,2)
            enddo
c           call add3s2(J_out(1,2),p(1,2),fdiff(1,2),invtau,          
c    $           -pertinv,nn)                                      
            call dssum  (J_out(1,2),nx1,ny1,nz1) ! Copy "1" to all elements

c           Now, output the "j"th column of A
            

            i = 0
            do ee=1,nelv
            do ii=1,nxyz
               if (int(flag(ii,ee)).ne.0) then 
                  write(50,77)J_out((ee-1)*nxyz+ii,2)
               endif
   77          format(e24.15)
            enddo
            enddo
         endif

      enddo
      enddo
      end
c-----------------------------------------------------------------------
      subroutine get_full_J(c_in,pert,invtau,pstep,nn,ncdim)

      implicit none                         

      include 'SIZE'                       
      include 'TOTAL'                     
      include 'DRIFT'                    
      include 'POISSON'
                                       
      real J_out(lpts1,ncdim)          
      real fi(lpts1,ncdim),ff(lpts1,ncdim),fdiff(lpts1,ncdim)  
      real c_in(lpts1,ncdim),c_out(lpts1,ncdim) 
      real p(lpts1,ncdim)                  
      real flag(lx1*ly1*lz1,lelt)
      real uin(lx1*ly1*lz1,lelt),uout(lx1*ly1*lz1,lelt)
      real pert,invtau,pertinv                    
      real glsc2
      integer ncdim,nn,pstep
      integer ic,i,j,k,k1                  
      integer e,ee,icc
      integer re,ii,jj,global_n,nglob
      character(20) filename1
                                           
      pertinv = 1./pert                    
      write(6,*)'====================='    
      write(6,*)'Begin computing J1'        
      write(6,*)'Parameter check'          
      write(6,*)'pert',pert,'pertinv',pertinv
      write(6,*)'dtNTinv',invtau            
      write(6,*)'pstep',pstep               
      write(6,*)'nn',nn,lpts1           
      write(6,*)'====================='     

      call compute_f_couple(fi,c_in,nn)   

      nxyz = lx1*ly1*lz1
      write(6,*)'nxyz',nxyz

      global_n = glsc2(dmask(1,1),mult,nn)
      write(6,*)'global_n',global_n
      nglob    = global_n + 0.4          ! Rounding

      do e=1,nelv
         re = e 
         call cfill(flag(1,e),float(re),nxyz)
      enddo
      call dsop(flag,'min',nx1,ny1,nz1)
                                            
      j=0
      do e=1,nelv
      do jj=1,nxyz
         re = e
         ! Found a real dof
         if (int(flag(jj,e)).eq.re.and.dmask((e-1)*nxyz+jj,1).ne.0) 
     $   then 
            j = j+1
            flag(jj,e) = j
         else
            flag(jj,e) = 0
         endif
      enddo
      enddo

      write(filename1,'(a,i0,a)') 'Jacobian_',pstep,'.dat'        
      open(50,file=filename1)                                      
      j = 0
      write(6,*)'check degree',(nelv-1)*nxyz+nxyz,nn

      do icc=1,ncdim
      do e=1,nelv
      do jj=1,nxyz

         do ic = 1,ncdim                                               
            call rzero(p(1,ic),nn)                                   
         enddo                                                         

         if (int(flag(jj,e)).ne.0) then  ! Next column of identity
            j=j+1
            p((e-1)*nxyz+jj,icc) = 1
            call dssum  (p(1,icc) ,nx1,ny1,nz1) ! Copy "1" to all elements
            call col2   (p(1,icc),mult,nn) ! do the average

            do ic = 1,ncdim                             
               do ii=1,nn
                  c_out(ii,ic) = c_in(ii,ic)+pert*p(ii,ic)
               enddo
c              call add3s2(c_out(1,ic),c_in(1,ic),p(1,ic),1.0,pert,nn) 
            enddo                                                     

            call compute_f_couple(ff,c_out,nn)                      
            do ic = 1,ncdim                                           
               call sub3(fdiff(1,ic),ff(1,ic),fi(1,ic),nn)          
            enddo                                                     
            do ic = 1,ncdim
               do ii=1,nn
                  J_out(ii,ic) = invtau*p(ii,ic)-pertinv*fdiff(ii,ic)
               enddo
            enddo
c           call add3s2(J_out(1,1),p(1,1),fdiff(1,1),invtau,          
c    $           -pertinv,nn)                                      
            do ic = 1,ncfim
               call dssum  (J_out(1,icc),nx1,ny1,nz1) ! Copy "1" to all elements
            enddo

c           Now, output the "j"th column of A
            i = 0
            do icc=1,ncdim
            do ee=1,nelv
            do ii=1,nxyz
               if (int(flag(ii,ee)).ne.0) then 
                  write(50,77)J_out((ee-1)*nxyz+ii,icc)
               endif
            enddo
            enddo
            enddo
         endif

      enddo
      enddo
      enddo

      end
