!---------------------------------------------------------------------- 
      subroutine set_overlap_again
c
c     Set up arrays for overlapping Schwartz algorithm *for pressure solver*
c
      include 'SIZE'
      include 'DOMAIN'
      include 'ESOLV'
      include 'INPUT'
      include 'TSTEP'
c
      REAL*8 dnekclock,t0
c
      parameter (          n_tri = 7*ltotd )
      common /scrns/  tri (n_tri)
      integer         tri,elem
c
      common /screv/ x(2*ltotd)
      common /scrvh/ y(2*ltotd)
      common /scrch/ z(2*ltotd)
c
      common /ctmp0/ nv_to_t(2*ltotd)
c
      parameter (lia = ltotd - 2 - 2*lelt)
      common /scrcg/ ntri(lelt+1),nmask(lelt+1)
     $             , ia(lia)
c
      common /scruz/ color   (4*ltotd)
      common /scrmg/ ddmask  (4*ltotd)
      common /ctmp1/ mask    (4*ltotd)

      parameter(lxx=lx1*lx1, levb=lelv+lbelv)
      common /fastd/  df(lx1*ly1*lz1,levb)
     $             ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)


      integer e

      t0 = dnekclock()

      if (lx1.eq.2) param(43)=1.
      if (lx1.eq.2.and.nid.eq.0) write(6,*) 'No mgrid for lx1=2!'

      if (ifaxis) ifmgrid = .false.
      if (param(43).ne.0) ifmgrid = .false.

      npass = 1
      if (ifmhd) npass = 2
      do ipass=1,npass
         ifield = 1

         if (ifsplit.and.ifmgrid) then

            if (ipass.gt.1) ifield = ifldmhd

            call swap_lengths
            call gen_fast_spacing(x,y,z)
 
            call hsmg_setup
            call h1mg_setup

         elseif (.not.ifsplit) then ! Pn-Pn-2

            if (ipass.gt.1) ifield = ifldmhd

            if (param(44).eq.1) then !  Set up local overlapping solves 
               call set_fem_data_l2(nel_proc,ndom,n_o,x,y,z,tri)
            else
               call swap_lengths
            endif
 
            e = 1
            if (ifield.gt.1) e = nelv+1

            call gen_fast_spacing(x,y,z)

!           This generates the local element eigenvalues/eigenvectors
!           for the pressure poisson problem (FDM).            
            call gen_fast_again(df(1,e),sr(1,e),ss(1,e),st(1,e),x,y,z)

            call init_weight_op
            if (param(43).eq.0) call hsmg_setup
         endif

!         call set_up_h1_crs

      enddo

      t0 = dnekclock()-t0
      if (nio.eq.0) then
         write(6,*) 'Local FDM (reset) ',
     $              t0, ' sec'
      endif

 
      return
      end
c-----------------------------------------------------------------------

      subroutine gen_fast_again(df,sr,ss,st,x,y,z)
c
c     Generate fast diagonalization matrices for each element
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'WZ'
c
      parameter(lxx=lx1*lx1)
      real df(lx1*ly1*lz1,1),sr(lxx*2,1),ss(lxx*2,1),st(lxx*2,1)
c
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt 
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
c
      integer lbr,rbr,lbs,rbs,lbt,rbt,e
c
      real x(lx1,ly1,lz1,nelv)
      real y(lx1,ly1,lz1,nelv)
      real z(lx1,ly1,lz1,nelv)
      real axwt(lx2)

      real fldr(lx1,lelv),flds(lx1,lelv),fldt(lx1,lelv)
      real fld(lx1,ly1,lz1,lelv)
      common /ctmpf_rho/ fldr,flds,fldt,fld

      integer ifld
      integer n


      ierr = 0

      if (param(44).eq.1) then
        if (nio.eq.0) then
          write(6,*) 'Density in FEM local solves not implemented'
          write(6,*) 'Exitting in gen_fast_again()'
        endif  

        call exitt
      endif

      ifld = 1
      n = lx1*ly1*lz1*nelv
      call copy(fld,vtrans(1,1,1,1,ifld),n)
      call get_1D_fld(fldr,flds,fldt,fld)   
!      call rone(fldr,lx1*lelv)
!      call rone(flds,lx1*lelv)
!      call rone(fldt,lx1*lelv)

      do e=1,nelv
c
         if (param(44).eq.1) then
           call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,2,ierr)
         else
           call get_fast_bc_again(lbr,rbr,lbs,rbs,lbt,rbt,e,3,ierr)
         endif
c
c        Set up matrices for each element.
c
         if (param(44).eq.1) then
           call set_up_fast_1D_fem( sr(1,e),lr,nr ,lbr,rbr
     $                      ,llr(e),lmr(e),lrr(e),zgm2(1,1),lx2,e)
         else
           call set_up_fast_1D_sem_again( sr(1,e),lr,nr ,lbr,rbr
     $                      ,llr(e),lmr(e),lrr(e),fldr(1,e),e)
         endif
         if (ifaxis) then
            xsum = vlsum(wxm2,lx2)
            do i=1,ly2
               yavg = vlsc2(y(1,i,1,e),wxm2,lx2)/xsum
               axwt(i) = yavg
            enddo
            call set_up_fast_1D_fem_ax( ss(1,e),ls,ns ,lbs,rbs
     $                 ,lls(e),lms(e),lrs(e),zgm2(1,2),axwt,ly2,e)
         else
            if (param(44).eq.1) then
               call set_up_fast_1D_fem( ss(1,e),ls,ns ,lbs,rbs
     $                      ,lls(e),lms(e),lrs(e),zgm2(1,2),ly2,e)
            else
               call set_up_fast_1D_sem_again( ss(1,e),ls,ns ,lbs,rbs
     $                      ,lls(e),lms(e),lrs(e),flds(1,e),e)
            endif
         endif
         if (if3d) then
            if (param(44).eq.1) then
               call set_up_fast_1D_fem( st(1,e),lt,nt ,lbt,rbt
     $                      ,llt(e),lmt(e),lrt(e),zgm2(1,3),lz2,e)
            else
               call set_up_fast_1D_sem_again( st(1,e),lt,nt ,lbt,rbt
     $                      ,llt(e),lmt(e),lrt(e),fldt(1,e),e)
            endif
         endif
c
c        Set up diagonal inverse
c
         if (if3d) then
            eps = 1.e-5 * (vlmax(lr(2),nr-2)
     $                  +  vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
            l   = 1
            do k=1,nt
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j) + lt(k)
               if (diag.gt.eps) then
                  df(l,e) = 1.0/diag
               else
c                 write(6,3) e,'Reset Eig in gen fast:',i,j,k,l
c    $                         ,eps,diag,lr(i),ls(j),lt(k)
c   3             format(i6,1x,a21,4i5,1p5e12.4)
                  df(l,e) = 0.0
               endif
               l = l+1
            enddo
            enddo
            enddo
         else
            eps = 1.e-5*(vlmax(lr(2),nr-2) + vlmax(ls(2),ns-2))
            l   = 1
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j)
               if (diag.gt.eps) then
                  df(l,e) = 1.0/diag
               else
c                 write(6,2) e,'Reset Eig in gen fast:',i,j,l
c    $                         ,eps,diag,lr(i),ls(j)
c   2             format(i6,1x,a21,3i5,1p4e12.4)
                  df(l,e) = 0.0
               endif
               l = l+1
            enddo
            enddo
         endif
c
c        Next element ....
c
      enddo

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('E INVALID BC FOUND in genfast$',ierrmx)
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem_again(s,lam,n,lbc,rbc,ll,lm,lr,
     $                                    rho,ie)

      implicit none

      include 'SIZE'
      include 'SEMHAT'
c
      common /fast1dsem/ g(lr2),w(lr2)
c
      real g,w
      real s(1),lam(1),ll,lm,lr
      integer lbc,rbc
      
      integer bb0,bb1,eb0,eb1,n,n1
      logical l,r

      integer ie
      real rho(lx1),dummy(lx1)
      real rhobh(lx1)     ! density*mass
      integer i,j,j0      ! prabal

      n=lx1
      call col3(rhobh,bh,rho,lx1)
     
      n=lx1-1
      !bcs on E are from normal vel component
      if(lbc.eq.2 .or. lbc.eq.3) then !wall,sym - dirichlet velocity
         eb0=1
      else !outflow,element - neumann velocity
         eb0=0
      endif
      if(rbc.eq.2 .or. rbc.eq.3) then !wall,sym - dirichlet velocity
         eb1=n-1
      else !outflow,element - neumann velocity
         eb1=n
      endif
      !bcs on B are from tangent vel component
      if(lbc.eq.2) then !wall - dirichlet velocity
         bb0=1
      else !outflow,element,sym - neumann velocity
         bb0=0
      endif
      if(rbc.eq.2) then !wall - dirichlet velocity
         bb1=n-1
      else !outflow,element,sym - neumann velocity
         bb1=n
      endif
c
      l = (lbc.eq.0)
      r = (rbc.eq.0)
c
c     calculate E tilde operator
      call set_up_fast_1D_sem_op_again(s,eb0,eb1,l,r,ll,lm,lr,bh,dgl,
     $                                 rho,0)

!     Already included rho scaling in rhobh
!      call set_up_fast_1D_sem_op(s,eb0,eb1,l,r,ll,lm,lr,rhobh,dgl,0)

c     call outmat(s,n+1,n+1,'  Et  ',ie)
c     calculate B tilde operator
!      call rone(dummy,lx1)
!      call set_up_fast_1D_sem_op_again(g,bb0,bb1,l,r,ll,lm,lr,bh,jgl,
!     $                                 dummy,1)
      call set_up_fast_1D_sem_op(g,bb0,bb1,l,r,ll,lm,lr,bh,jgl,1)

      n=n+1
      call generalev(s,g,lam,n,w)
      if(.not.l) call row_zero(s,n,n,1)
      if(.not.r) call row_zero(s,n,n,n)
      call transpose(s(n*n+1),n,s,n) ! compute the transpose of s

      return
      end
c-----------------------------------------------------------------------

      subroutine set_up_fast_1D_sem_op_again
     $                  (g,b0,b1,l,r,ll,lm,lr,bh,jgl,rho,jscl)
c                  -1 T
c     G = J (rho*B)  J
c
c     gives the inexact restriction of this matrix to
c     an element plus one node on either side
c
c     g - the output matrix
c     b0, b1 - the range for Bhat indices for the element
c              (enforces boundary conditions)
c     l, r - whether there is a left or right neighbor
c     ll,lm,lr - lengths of left, middle, and right elements
c     bh - hat matrix for B
c     jgl - hat matrix for J (should map vel to pressure)
c     jscl - how J scales
c            0: J = Jh
c            1: J = (L/2) Jh
c
c     result is inexact because:
c        neighbor's boundary condition at far end unknown
c        length of neighbor's neighbor unknown
c        (these contribs should be small for large N and
c         elements of nearly equal size)
c

!     rho - density

      include 'SIZE'
      real g(0:lx1-1,0:lx1-1)
      real bh(0:lx1-1),jgl(1:lx2,0:lx1-1)
      real ll,lm,lr
      integer b0,b1
      logical l,r
      integer jscl
c
      real bl(0:lx1-1),bm(0:lx1-1),br(0:lx1-1)
      real gl,gm,gr,gll,glm,gmm,gmr,grr
      real fac
      integer n

      real rho(0:lx1-1)       ! density
      real rhobh(0:lx1-1)     ! density*mass

      n=lx1-1

      call col3(rhobh(0),bh(0),rho(0),lx1)

c     compute the scale factors for J      
      if (jscl.eq.0) then
         gl=1.
         gm=1.
         gr=1.
      elseif (jscl.eq.1) then
         gl=0.5*ll
         gm=0.5*lm
         gr=0.5*lr
      endif
      gll = gl*gl
      glm = gl*gm
      gmm = gm*gm
      gmr = gm*gr
      grr = gr*gr
c
c     compute the summed inverse mass matrices for
c     the middle, left, and right elements
      do i=1,n-1
         bm(i)=2. /(lm*rhobh(i))
      enddo
      if (b0.eq.0) then
         bm(0)=0.5*lm*rhobh(0)
         if(l) bm(0)=bm(0)+0.5*ll*rhobh(n)
         bm(0)=1. /bm(0)
      endif
      if (b1.eq.n) then
         bm(n)=0.5*lm*rhobh(n)
         if(r) bm(n)=bm(n)+0.5*lr*rhobh(0)
         bm(n)=1. /bm(n)
      endif
c     note that in computing bl for the left element,
c     bl(0) is missing the contribution from its left neighbor
      if (l) then
         do i=0,n-1
            bl(i)=2. /(ll*rhobh(i))
         enddo
         bl(n)=bm(0)
      endif
c     note that in computing br for the right element,
c     br(n) is missing the contribution from its right neighbor
      if (r) then
         do i=1,n
            br(i)=2. /(lr*rhobh(i))
         enddo
         br(0)=bm(n)
      endif
c      
      call rzero(g,(n+1)*(n+1))
      do j=1,n-1
         do i=1,n-1
            do k=b0,b1
               g(i,j) = g(i,j) + gmm*jgl(i,k)*bm(k)*jgl(j,k)
            enddo
         enddo
      enddo
c      
      if (l) then
         do i=1,n-1
            g(i,0) = glm*jgl(i,0)*bm(0)*jgl(n-1,n)
            g(0,i) = g(i,0)
         enddo
c        the following is inexact
c        the neighbors bc's are ignored, and the contribution
c        from the neighbor's neighbor is left out
c        that is, bl(0) could be off as noted above
c        or maybe i should go from 1 to n
         do i=0,n
            g(0,0) = g(0,0) + gll*jgl(n-1,i)*bl(i)*jgl(n-1,i)
         enddo
      else
         g(0,0)=1.
      endif
c      
      if (r) then
         do i=1,n-1
            g(i,n) = gmr*jgl(i,n)*bm(n)*jgl(1,0)
            g(n,i) = g(i,n)
         enddo
c        the following is inexact
c        the neighbors bc's are ignored, and the contribution
c        from the neighbor's neighbor is left out
c        that is, br(n) could be off as noted above
c        or maybe i should go from 0 to n-1
         do i=0,n
            g(n,n) = g(n,n) + grr*jgl(1,i)*br(i)*jgl(1,i)
         enddo
      else
         g(n,n)=1.
      endif
     
      return
      end
c-----------------------------------------------------------------------

      subroutine get_1D_fld(fldr,flds,fldt,fld)

!     Get 1D version of Density to use in local fdm solver

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'WZ'

      integer nxn,nz0
      
      real w(lx1)       ! weight (wxm1)
      real fld(lx1,ly1,lz1,lelv)
      real lr2,ls2,lt2,wsum,weight

      integer nx,ny,nz
      integer ie,i,j,k
      real fldr(lx1,lelv),flds(lx1,lelv),fldt(lx1,lelv)


      nx = lx1
      ny = ly1
      nz = lz1

      call dssum(fld,nx,ny,nz)
      call rzero(fldr,lx1*nelv)
      call rzero(flds,lx1*nelv)
      call rzero(fldt,lx1*nelv)

      call copy(w,wxm1,lx1)

      do ie=1,nelv
        if (if3d) then
           do i=1,nx
             lr2  = 0.
             wsum = 0.
             do k=1,nz
             do j=1,ny
                weight = w(j)*w(k)
                lr2  = lr2  +   weight*fld(i,j,k,ie)
                wsum = wsum + weight
             enddo
             enddo
             lr2       = lr2/wsum
             fldr(i,ie)= lr2
           enddo
c
           do j=1,ny
             ls2 = 0.
             wsum = 0.
             do k=1,nz
             do i=1,nx
                weight = w(i)*w(k)
                ls2  = ls2  +   weight*fld(i,j,k,ie)
                wsum = wsum + weight
             enddo
             enddo
             ls2       = ls2/wsum
             flds(j,ie)= ls2
           enddo  
c
           do k=1,nz
             lt2 = 0.
             wsum = 0.
             do j=1,ny
             do i=1,nx
                weight = w(i)*w(j)
                lt2  = lt2  +   weight*fld(i,j,k,ie)
                wsum = wsum + weight
             enddo
             enddo
             lt2       = lt2/wsum
             fldt(k,ie)= lt2
           enddo

        else              ! 2D

           do i=1,nx
             lr2 = 0.
             wsum = 0.
             do j=1,ny
                weight = w(j)
                lr2  = lr2  + weight*fld(i,j,1,ie)
                wsum = wsum + weight
             enddo
             lr2       = lr2/wsum
             fldr(i,ie)= lr2
           enddo  
c
           do j=1,ny
             ls2 = 0.
             wsum = 0.
             do i=1,nx
                weight = w(i)
                ls2  = ls2  + weight*fld(i,j,1,ie)
                wsum = wsum + weight
             enddo
             ls2       = ls2/wsum
             flds(j,ie)= ls2
           enddo
        endif           ! if3d
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine get_fast_bc_again(lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ierr)

  
      integer                lbr,rbr,lbs,rbs,lbt,rbt,e,bsym
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'TSTEP'
c
      integer fbc(6)
c
c     ibc = 0  <==>  Dirichlet
c     ibc = 1  <==>  Dirichlet, outflow (no extension)
c     ibc = 2  <==>  Neumann,   


      do iface=1,2*ldim
         ied = eface(iface)
         ibc = -1

         if (ifmhd) call mhd_bc_dn(ibc,iface,e) ! can be overwritten by 'mvn'

         if (cbc(ied,e,ifield).eq.'   ') ibc = 0      ! ibc = 0
         if (cbc(ied,e,ifield).eq.'E  ') ibc = 0      ! ibc = 0   ! prabal
         if (cbc(ied,e,ifield).eq.'msi') ibc = 0      ! ibc = 0
         if (cbc(ied,e,ifield).eq.'MSI') ibc = 0      ! ibc = 0
         if (cbc(ied,e,ifield).eq.'P  ') ibc = 0      ! ibc = 0   ! prabal
         if (cbc(ied,e,ifield).eq.'p  ') ibc = 0      ! ibc = 0
         if (cbc(ied,e,ifield).eq.'O  ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'ON ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'o  ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'on ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'MS ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'ms ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'MM ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'mm ') ibc = 1      ! ibc = 1
         if (cbc(ied,e,ifield).eq.'mv ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'mvn') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'v  ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'V  ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'W  ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'SYM') ibc = bsym   ! ibc = bsym
         if (cbc(ied,e,ifield).eq.'SL ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'sl ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'SHL') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'shl') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'A  ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'S  ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'s  ') ibc = 2      ! ibc = 2
         if (cbc(ied,e,ifield).eq.'J  ') ibc = 0      ! ibc = 0
         if (cbc(ied,e,ifield).eq.'SP ') ibc = 0      ! ibc = 0

         fbc(iface) = ibc

         if (ierr.eq.-1) write(6,1) ibc,ied,e,ifield,cbc(ied,e,ifield)
  1      format(2i3,i8,i3,2x,a3,'  get_fast_bc_error')

      enddo

      if (ierr.eq.-1) call exitti('Error A get_fast_bc$',e)

      lbr = fbc(1)
      rbr = fbc(2)
      lbs = fbc(3)
      rbs = fbc(4)
      lbt = fbc(5)
      rbt = fbc(6)

      ierr = 0 
      if (ibc.lt.0) ierr = lglel(e)

c     write(6,6) e,lbr,rbr,lbs,rbs,(cbc(k,e,ifield),k=1,4)
c   6 format(i5,2x,4i3,3x,4(1x,a3),'  get_fast_bc')

      return
      end
c-----------------------------------------------------------------------


