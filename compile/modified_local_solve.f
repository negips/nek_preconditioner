!     Major modifications for local_solve_fdm

!======================================================================
      subroutine local_solves_fdm3(u,v)
c
c     Given an input vector v, this returns the additive Schwarz solution
c
c                       -1
c                  u = M   v
c
c                      T     -1  
c     where M = sum ( R_i A_i    R_i )
c                i
c
c     The R_i's are simply index_set restriction operators.
c
c     The local solves are performed with the fast diagonalization method.
c
      include 'SIZE'
      include 'INPUT'
      include 'DOMAIN'
      include 'PARALLEL'
      include 'SOLN'          ! vmult
c
      include 'TSTEP'
      include 'CTIMER'
c
      real u(lx2,ly2,lz2,lelv),v(lx2,ly2,lz2,lelv)
      common /scrpre/ v1(lx1,ly1,lz1,lelv)
     $               ,w1(lx1,ly1,lz1),w2(lx1,ly1,lz1)
      common /scrover/ ar(lelv)

      parameter(lxx=lx1*lx1, levb=lelv+lbelv)
      common /fastd/  df(lx1*ly1*lz1,levb)
     $             ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)
      integer e,eb,eoff

      real wk1(lx1,ly1,lz1,lelv)

c
      if (icalld.eq.0) tsolv=0.0
      icalld=icalld+1
      nsolv=icalld
c
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv
c
!     Fill interiors
      call fill_interior(v1,v)
      call dface_ext    (v1)
      call dssum        (v1,lx1,ly1,lz1)
      call dface_add1si (v1,-1.)

      call copy(wk1,v1,n)

c     Fill interiors
      call fill_interior(v1,v)
      call dface_ext    (v1)
      call ext_crnr     (v1)

      call dssum        (v1,lx1,ly1,lz1)
      call dface_add1si (v1,-1.)
      call crnr_sub     (v1,wk1)
c
c     Now solve each subdomain problem:
c
      etime1=dnekclock()

      eoff  = 0
      if (ifield.gt.1) eoff  = nelv

      do e = 1,nelv
         eb = e + eoff
         call fastdm1(v1(1,1,1,e),df(1,eb)
     $                           ,sr(1,eb),ss(1,eb),st(1,eb),w1,w2)
      enddo
      tsolv=tsolv+dnekclock()-etime1
c
c     Exchange/add elemental solutions
c
!      call s_face_to_int (v1,-1.)
!      call dssum         (v1,lx1,ly1,lz1)
!      call s_face_to_int (v1, 1.)
!      if(param(42).eq.0) call do_weight_op(v1)

c     Copy back to pressure grid
c     lx1 -> lx2
      call extract_interior(u,v1)

c
      return
      end
c-----------------------------------------------------------------------

      subroutine dd_solver2(u,v,h1,h2,h2inv,nit,wk1,wk2)

      implicit none

      include 'SIZE'
      include 'DOMAIN'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'CTIMER'

      real h1(1),h2(1),h2inv(1)
      real u(1),v(1),wk1(1),wk2(1)
      real uc
      common /scrprc/ uc(lx1*ly1*lz1*lelt)

      real alpha
      integer nit,i
      integer ntot
      integer intype
      real sigma
      real gl2norm2
      real nor
c
      if (icalld.eq.0) then
         tddsl=0.0
         tcrsl=0.0
         nddsl=0
         ncrsl=0
      endif
      icalld = icalld + 1
      nddsl  = nddsl  + 1
      ncrsl  = ncrsl  + 1

      ntot  = lx2*ly2*lz2*nelv
      call rzero(u,ntot)
      call copy(wk1,v,ntot)
      intype = 1         ! explicit
!      nit    = 50       ! Iterations
      etime1=dnekclock()
      sigma  = 1.00
      do i=1,nit
!        call local_solves_fdm(wk2,wk1)          ! wk2 = inv(M)*wk1
        call local_solves_fdm3(wk2,wk1)
        call add2s2(u,wk2,sigma,ntot)           ! u   = u + \sigma*wk2
        call cdabdtp(wk1,u,h1,h2,h2inv,intype)  ! wk1 = A u
        call sub2(wk1,v,ntot)                   ! wk1 = (Au - g)
        call chsign(wk1,ntot)                   ! wk1 = (g - Au)
      enddo  

      nor = gl2norm2(wk1,ntot)
      if (nio.eq.0) write(6,*) 'Residual Norm',nit,nor


      tddsl=tddsl+dnekclock()-etime1

      etime1=dnekclock()
!      call crs_solve_l2 (uc,v)
      tcrsl=tcrsl+dnekclock()-etime1

      alpha = 10.
c     if (param(89).ne.0.) alpha = abs(param(89))
!      call add2s2(u,uc,alpha,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine ext_edge(x)

      implicit none        

      include 'SIZE'
      include 'INPUT'

      real x(lx1,ly1,lz1,1)
      integer ie,ix,iy,iz

      do ie=1,nelv

         if (if3d) then

          do ix=2,lx1-1
            x(ix,1,1,ie)      = x(ix,2,1,ie)
            x(ix,ly1,1,ie)    = x(ix,ly1-1,1,ie)
            x(ix,1,lz1,ie)    = x(ix,2,lz1,ie)
            x(ix,ly1,lz1,ie)  = x(ix,ly1-1,lz1,ie)
          enddo  

          do iy=2,ly1-1
            x(1,iy,1,ie)      = x(1,iy,2,ie)
            x(1,iy,lz1,ie)    = x(1,iy,lz1-1,ie)
            x(lx1,iy,1,ie)    = x(lx1,iy,2,ie)
            x(lx1,iy,lz1,ie)  = x(lx1,iy,lz1-1,ie)
          enddo  

          do iz=2,lz1-1
            x(1,1,iz,ie)      = x(2,1,iz,ie)
            x(lx1,1,iz,ie)    = x(lx1-1,1,iz,ie)
            x(1,ly1,iz,ie)    = x(2,ly1,iz,ie)
            x(lx1,ly1,iz,ie)  = x(lx1-1,ly1,iz,ie)
          enddo  
         
         else
!          No Edge extensions in 2D.
!          Only Faces And Corners           

        endif
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine ext_crnr(x)
c
      include 'SIZE'
      include 'INPUT'
      real x(lx1,ly1,lz1,1)
c
      do ie=1,nelv
c
         if (if3d) then

          x(1,1,1,ie)       = x(2,2,2,ie)
          x(lx1,1,1,ie)     = x(lx1-1,2,2,ie)

          x(1,ly1,1,ie)     = x(2,ly1-1,2,ie)
          x(lx1,ly1,1,ie)   = x(lx1-1,ly1-1,2,ie)

          x(1,1,lz1,ie)     = x(2,2,lz1-1,ie)
          x(1,ly1,lz1,ie)   = x(2,ly1-1,lz1-1,ie)

          x(lx1,1,lz1,ie)   = x(lx1-1,2,lz1-1,ie)
          x(lx1,ly1,lz1,ie) = x(lx1-1,ly1-1,lz1-1,ie)

         else
            
          x(1,1,1,ie) = x(2,2,1,ie)
          x(lx1,1,1,ie) = x(lx1-1,2,1,ie)

          x(1,ly1,1,ie) = x(2,ly1-1,1,ie)
          x(lx1,ly1,1,ie) = x(lx1-1,ly1-1,1,ie)

        endif
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine crnr_sub(x,y)

      implicit none

      include 'SIZE'
      include 'INPUT'

      integer ie

      real x(lx1,ly1,lz1,1)
      real y(lx1,ly1,lz1,1)
      real s,y1,y2,y3

      do ie=1,nelv
c
         if (if3d) then

!          x(1,1,1,ie)       = x(2,2,2,ie)
!          x(lx1,1,1,ie)     = x(lx1-1,2,2,ie)
!
!          x(1,ly1,1,ie)     = x(2,ly1-1,2,ie)
!          x(lx1,ly1,1,ie)   = x(lx1-1,ly1-1,2,ie)
!
!          x(1,1,lz1,ie)     = x(2,2,lz1-1,ie)
!          x(1,ly1,lz1,ie)   = x(2,ly1-1,lz1-1,ie)
!
!          x(lx1,1,lz1,ie)   = x(lx1-1,2,lz1-1,ie)
!          x(lx1,ly1,lz1,ie) = x(lx1-1,ly1-1,lz1-1,ie)

         else
         
          s                   = x(1,1,1,ie)
          y1                  = y(2,1,1,ie)
          y2                  = y(1,2,1,ie)
          y3                  = y(2,2,1,ie)
          x(1,1,1,ie)         = s - (y1+y2+y3)

          s                   = x(lx1,1,1,ie)
          y1                  = y(lx1-1,1,1,ie)
          y2                  = y(lx1,2,1,ie)
          y3                  = y(lx1-1,2,1,ie)
          x(lx1,1,1,ie)       = s - (y1+y2+y3)

          s                   = x(1,ly1,1,ie)
          y1                  = y(2,ly1,1,ie)
          y2                  = y(1,ly1-1,1,ie)
          y3                  = y(2,ly1-1,1,ie)
          x(1,ly1,1,ie)       = s - (y1+y2+y3)

          s                   = x(lx1,ly1,1,ie)
          y1                  = y(lx1-1,ly1,1,ie)
          y2                  = y(lx1,ly1-1,1,ie)
          y3                  = y(lx1-1,ly1-1,1,ie)
          x(lx1,ly1,1,ie)     = s - (y1+y2+y3)

        endif
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine fill_interior(x,y)

!     Put Mesh2 node values on to the interior
!     Mesh1 nodes.

      implicit none

      include 'SIZE'

      real x(lx1,ly1,lz1,lelv)
      real y(lx2,ly2,lz2,lelv)

      integer iz1,e,ix,iy,iz
      integer n

      n = lx1*ly1*lz1*nelv

      call rzero(x,n)
!     Fill interiors
      iz1 = 0
      if (ndim.eq.3) iz1=1
      do e=1,nelv
         do iz=1,lz2
         do iy=1,ly2
         do ix=1,lx2
           x(ix+1,iy+1,iz+iz1,e) = y(ix,iy,iz,e)
         enddo
         enddo
         enddo
      enddo

      return
      end subroutine 
!---------------------------------------------------------------------- 

       subroutine extract_interior(x,y)

      implicit none

      include 'SIZE'

      real y(lx1,ly1,lz1,lelv)
      real x(lx2,ly2,lz2,lelv)

      integer iz1,e,ix,iy,iz
      integer n
c
c     Map back to pressure grid (extract interior values)
c
      do e=1,nelv
         do iz=1,lz2
         do iy=1,ly2
         do ix=1,lx2
            x(ix,iy,iz,e) = y(ix+1,iy+1,iz+iz1,e)
         enddo
         enddo
         enddo
      enddo


      return
      end subroutine 
!---------------------------------------------------------------------- 
     
      subroutine mysemhat(a,b,c,d,z,dgll,dgllt,jgll,bgl,zgl,dgl,jgl,n,w)

      implicit none

c     Generate matrices for single element, 1D operators:
c
c     a        = Laplacian
c     b        = diagonal mass matrix
c     c        = convection operator b*d
c     d        = derivative matrix
c     dgll     = derivative matrix
c     dgllt    = derivative matrix transpose
c     jgll     = interpolation matrix 
c     z        = GLL points
c
c     zgl      = GL points
c     bgl      = diagonal mass matrix on GL
c     dgl      = derivative matrix,    mapping from velocity nodes to pressure
c     jgl      = interpolation matrix, mapping from velocity nodes to pressure
c
c     n        = polynomial degree (velocity space)
c     w        = work array of size 2*n+2


      integer i,j,k,n,np,nm,n2

      real a(0:n,0:n),b(0:n),c(0:n,0:n),d(0:n,0:n),z(0:n)
      real dgll(0:n,0:n),dgllt(0:n,0:n),jgll(0:n,0:n)
c
      real bgl(1:n-1),zgl(1:n-1)
      real dgl(1:n-1,0:n),jgl(1:n-1,0:n)
c
      real w(0:2*n+1)
c
      np = n+1    ! No of Points on M1 (lx1)
      nm = n-1    ! No of Points on M2 (lx2)
      n2 = n-2    ! Polynomial order of M2 (lx2-1)
c
      call zwgll (z,b,np)
c
      do i=0,n
         call fd_weights_full(z(i),z(0),n,1,w)
         do j=0,n
            d(i,j) = w(j+np)                   !  Derivative matrix
         enddo
      enddo

      if (n.eq.1) return                       !  No interpolation for n=1

!     Weights and Interpolation on GLL points      
      do i=0,n
!        Using n instead of n2
!        as polynomial order      
         call fd_weights_full(z(i),z(0),n,1,w(0))
         do j=0,n
            jgll(i,j) = w(j   )                  !  Interpolation matrix
            dgll(i,j) = w(j+np)                  !  Derivative    matrix
         enddo
      enddo

!     prabal
      call transpose(dgllt,np,dgll,np)
c
      call rzero(a,np*np)
      do j=0,n
      do i=0,n
         do k=0,n
            a(i,j) = a(i,j) + d(k,i)*b(k)*d(k,j)
         enddo
         c(i,j) = b(i)*d(i,j)
      enddo
      enddo
c
      call zwgl (zgl,bgl,nm)
c
      do i=1,n-1
         call fd_weights_full(zgl(i),z,n,1,w)
         do j=0,n
            jgl(i,j) = w(j   )                  !  Interpolation matrix
            dgl(i,j) = w(j+np)                  !  Derivative    matrix
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------

      subroutine exchange_m2(x1,x2,wk1,wk2,wk3)

!     Exchange interior values of neighboring elements
!     On the Pressure mesh        

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      integer n,n2

      real x1(lx1,ly1,lz1,lelv)     ! output
      real x2(lx2,ly2,lz2,lelv)     ! input
      real wk1(lx1,ly1,lz1,lelv)
      real wk2(lx1,ly1,lz1,lelv)
      real wk3(lx1,ly1,lz1,lelv)

      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      call rzero(wk1,n)
      call fill_interior(wk1,x2)    ! Put M2 on M1 interior
      call dface_ext    (wk1)       ! Extend to face
      call dssum        (wk1,lx1,ly1,lz1)
      call dface_add1si (wk1,-1.)

      call rzero(wk2,n) 
      call fill_interior(wk2,x2)
      call dface_ext    (wk2)       ! Extend to face
      call ext_crnr     (wk2)       ! Extend to corners

      call dssum        (wk2,lx1,ly1,lz1)
      call dface_add1si (wk2,-1.)
      call crnr_sub     (wk2,wk1)

      call copy(x1,wk2,n)

      return
      end subroutine 
!---------------------------------------------------------------------- 
      subroutine swap_lengths_mesh2

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WZ'

      real l
      common /swaplengths/ l(lx1,ly1,lz1,lelv)

      real lr, ls, lt   ! not used by swap lengths
      real llr(lelt)    ! length of left element along "r" 
      real lls(lelt)    ! length of left element along "s"
      real llt(lelt)    ! length of left element along "t"
      real lmr(lelt)    ! length of this element along "r"
      real lms(lelt)    ! length of this element along "s"
      real lmt(lelt)    ! length of this element along "t"
      real lrr(lelt)    ! length of right element along "r"
      real lrs(lelt)    ! length of right element along "s"
      real lrt(lelt)    ! length of right element along "t"
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr,lls,llt
     $              , lmr,lms,lmt
     $              , lrr,lrs,lrt

      integer e,i,j,k,f,nfaces
      integer n,n2,nz0,nzn,nx

!     Coordinates
!     M2 -> M1 (Interior Points just copied) 
!     Face points from M1 mesh      
      real x21(lx1,ly1,lz1,lelt)
      real y21(lx1,ly1,lz1,lelt)
      real z21(lx1,ly1,lz1,lelt)

      common /scrsf/ x21,y21,z21


!     Fill interior values      
      call fill_interior(x21,xm2)
      call fill_interior(y21,ym2)
      call rzero(z21,lx1*ly1*lz1*nelv)
      if (if3d) call fill_interior(z21,zm2)

!     Copy all face values      
      nfaces = 2*ndim
      do e=1,nelv
      do f=1,nfaces
        call facec(x21,xm1,e,f,nx1,ny1,nz1,nelv)
        call facec(y21,ym1,e,f,nx1,ny1,nz1,nelv)
        if (if3d) call facec(z21,zm1,e,f,nx1,ny1,nz1,nelv)
      enddo
      enddo

      n2 = lx1-1
      nz0 = 1
      nzn = 1
      nx  = lx1-2
      if (if3d) then
         nz0 = 0
         nzn = n2
      endif
      call plane_space(lmr,lms,lmt,0,n2,wxm1,x21,y21,z21,nx,n2,nz0,nzn)

      n=n2+1
      if (if3d) then
         do e=1,nelv
         do j=2,n2
         do k=2,n2
            l(1,k,j,e) = lmr(e)
            l(n,k,j,e) = lmr(e)
            l(k,1,j,e) = lms(e)
            l(k,n,j,e) = lms(e)
            l(k,j,1,e) = lmt(e)
            l(k,j,n,e) = lmt(e)
         enddo
         enddo
         enddo
         call dssum(l,n,n,n)
         do e=1,nelv
            llr(e) = l(1,2,2,e)-lmr(e)
            lrr(e) = l(n,2,2,e)-lmr(e)
            lls(e) = l(2,1,2,e)-lms(e)
            lrs(e) = l(2,n,2,e)-lms(e)
            llt(e) = l(2,2,1,e)-lmt(e)
            lrt(e) = l(2,2,n,e)-lmt(e)
         enddo
      else
         do e=1,nelv
         do j=2,n2
            l(1,j,1,e) = lmr(e)
            l(n,j,1,e) = lmr(e)
            l(j,1,1,e) = lms(e)
            l(j,n,1,e) = lms(e)
         enddo
         enddo
         call dssum(l,n,n,1)
         do e=1,nelv
            llr(e) = l(1,2,1,e)-lmr(e)
            lrr(e) = l(n,2,1,e)-lmr(e)
            lls(e) = l(2,1,1,e)-lms(e)
            lrs(e) = l(2,n,1,e)-lms(e)
         enddo
      endif
      return
      end
c----------------------------------------------------------------------



