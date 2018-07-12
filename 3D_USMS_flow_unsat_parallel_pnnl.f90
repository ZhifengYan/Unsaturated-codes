!!!**************** discription ***************
!!! A parallel FVM code with fully-stagger arrangment to solve USMS model
!!! double-precision
!!! consider porosity heterogeneity, ensure mass conservation
!!! walls on both side-boundary
!!! developed by Zhifeng Yan
!!! 2015-07-01
!!!
!!! for imbibition, set por(por==1)=0 is okey, becasue only capillary force works.
!!! for drainage, need to set por(por==1) = a (where a < 1) to get stable simu with a acceptable dt; also the bottom two layers need to set por(por==1) = a.
!!!********************************************



Module GLOBAL
!!! for grid 
  Integer::i,j,k,Nx,Ny,Nz,Nt,LX,LY,LZ,ilenx,ileny,ilenz

!!! for vel and pres
  Double precision,Dimension(:,:,:),Allocatable::u_o,v_o,w_o,h_o,h_star,u,v,w,h,u_star,v_star,w_star
  Double precision,Dimension(:,:,:),Allocatable::tmp_u_o,tmp_v_o,tmp_w_o,tmp_w_star,tmp_w,tmp_h_star,tmp_h
  Double precision,Dimension(:,:,:),Allocatable::ne,tmp_ne,tmp_m,halfx_ne,halfy_ne,halfz_ne,halfx_KKs,halfy_KKs,halfz_KKs
  Double precision,Dimension(:,:,:),Allocatable::halfx_KK,halfy_KK,halfz_KK,halfx_KK_inv,halfy_KK_inv,halfz_KK_inv,rx,ry,rz
  Double precision,Dimension(:,:,:),Allocatable::Se,tmp_Se,halfx_Se,halfy_Se,halfz_Se,S,FF,dSdh,Sm
  Double precision,Dimension(:,:,:),Allocatable::theta,tmp_theta,halfx_theta,halfy_theta,halfz_theta,halfx_m,halfy_m,halfz_m
  Double precision,Dimension(:,:,:),Allocatable::q0,q,a,b,c,d,e,f,g
  Double precision,Dimension(:,:,:),Allocatable::m,tau,lambda
  Double precision,Dimension(:),Allocatable::ma,mb,mc,md,me,mf,mg,HH_star,HH,QQ
  logical, allocatable, dimension(:,:,:):: wall,tmp_wall
  logical, allocatable, dimension(:):: WW
  Double precision,Dimension(:,:,:),Allocatable::alpha_prime
  Double precision::beta_prime,cc_eff,sx,h_cc,DD,a0,a1,a2,a3

!!! for input and output
  Double precision,Dimension(3)::gvt
  Integer::vbc_t,vbc_d,hbc_t,hbc_d,Max_step,Print_step,File_step,t0
  Double precision::h_d,h_t,u_d,v_d,w_d,u_t,v_t,w_t
  Double precision::dx,dy,dz,dt,rho,nu,g0,Cm,alpha,beta,Ms,omega,h_eps,HH_eps,alpha_h,m0,d10,hcc_s
  Character(64)::PMname,OUT_DIR,IN_DIR,Rootname


!!! for computation
  Integer::istep,p,hijk,nxx,nyy,nzz,np,nn,pp
  Double precision::SMALL,mean_por,h_err,h_err_local,HH_err,HH_err_local
  Double precision::frlocalsum,azlocalsum,frsum,azsum
  Double precision::Porlocal,Porsum,ulocalmax,umax,vlocalmax,vmax,wlocalmax,wmax,ulocalmin,umin,vlocalmin,vmin,wlocalmin,wmin
  Double precision::ulocalsum,usum,vlocalsum,vsum,wlocalsum,wsum
  Double precision::total_ax,total_ay,total_az
  Double Precision::starttime,endtime,totaltime

!!! for parallel
  Integer,PARAMETER::MEMO=8
  INTEGER::Nstart,error,mype,npes,idex,iprev,inext,nallgrp

  !!  Parameter (SMALL=1.e-15)
  !!  Parameter (LARGE=1.e+10)

End Module GLOBAL




PROGRAM USMSSolver

  Use GLOBAL
  Implicit None
  include 'mpif.h'
  integer::stat(MPI_STATUS_SIZE)

  Call input

  SMALL=1.e-20

  nxx = int(0.5*Nx)+1
  nyy = int(0.5*Ny)+1
  nzz = int(0.5*Nz)+1
  Nt = Nx*Ny*Nz

  ilenx = (Nx+1)*Ny
  ileny = Nx*(Ny+1)
  ilenz = Nx*Ny



!!! memory issue

  Allocate (u_o(Nx+1,Ny,Nz),v_o(Nx,Ny+1,Nz),w_o(Nx,Ny,Nz+1),h_o(Nx,Ny,Nz),h_star(Nx,Ny,Nz),h(Nx,Ny,Nz))
  Allocate (u_star(Nx+1,Ny,Nz),v_star(Nx,Ny+1,Nz),w_star(Nx,Ny,Nz+1))
  Allocate (tmp_u_o(Nx+1,Ny,2),tmp_v_o(Nx,Ny+1,2),tmp_w_o(Nx,Ny,2),tmp_h(Nx,Ny,2))
  Allocate (tmp_h_star(Nx,Ny,2))

  Allocate (u(Nx+1,Ny,Nz),v(Nx,Ny+1,Nz),w(Nx,Ny,Nz+1))
  Allocate (ne(Nx,Ny,Nz),wall(Nx,Ny,Nz),halfx_ne(Nx+1,Ny,Nz),halfy_ne(Nx,Ny+1,Nz),halfz_ne(Nx,Ny,Nz+1))
  Allocate (theta(Nx,Ny,Nz),halfx_theta(Nx+1,Ny,Nz),halfy_theta(Nx,Ny+1,Nz),halfz_theta(Nx,Ny,Nz+1))
  Allocate (Se(Nx,Ny,Nz),halfx_Se(Nx+1,Ny,Nz),halfy_Se(Nx,Ny+1,Nz),halfz_Se(Nx,Ny,Nz+1),S(Nx,Ny,Nz),Sm(Nx,Ny,Nz))
  Allocate (m(Nx,Ny,Nz),halfx_m(Nx+1,Ny,Nz),halfy_m(Nx,Ny+1,Nz),halfz_m(Nx,Ny,Nz+1),tau(Nx,Ny,Nz),lambda(Nx,Ny,Nz))

  Allocate (tmp_wall(Nx,Ny,2),tmp_ne(Nx,Ny,2),tmp_m(Nx,Ny,2),tmp_theta(Nx,Ny,2),tmp_Se(Nx,Ny,2))
  Allocate (halfx_KKs(Nx+1,Ny,Nz),halfy_KKs(Nx,Ny+1,Nz),halfz_KKs(Nx,Ny,Nz+1))
  Allocate (halfx_KK(Nx+1,Ny,Nz),halfy_KK(Nx,Ny+1,Nz),halfz_KK(Nx,Ny,Nz+1))
  Allocate (halfx_KK_inv(Nx+1,Ny,Nz),halfy_KK_inv(Nx,Ny+1,Nz),halfz_KK_inv(Nx,Ny,Nz+1),rx(Nx+1,Ny,Nz),ry(Nx,Ny+1,Nz),rz(Nx,Ny,Nz+1))
  Allocate (FF(Nx,Ny,Nz),dSdh(Nx,Ny,Nz),alpha_prime(Nx,Ny,Nz))

  Allocate (q0(Nx,Ny,Nz),q(Nx,Ny,Nz),a(Nx,Ny,Nz),b(Nx,Ny,Nz),c(Nx,Ny,Nz),d(Nx,Ny,Nz),e(Nx,Ny,Nz),f(Nx,Ny,Nz),g(Nx,Ny,Nz))
  Allocate (ma(Nt),mb(Nt),mc(Nt),md(Nt),me(Nt),mf(Nt),mg(Nt),HH_star(Nt),HH(Nt),QQ(Nt),WW(Nt))

!!! import grids


!!! initilization


  u_star = 0.
  v_star = 0.
  w_star = 0.
  h_star = 0.

  u = 0.
  v = 0.
  w = 0.
  h = 0.


  !  r = 0.

  q = 0.
  a = 0.
  b = 0.
  c = 0.
  d = 0.
  e = 0.
  f = 0.
  g = 0.

  !  MM = 0.
  ma = 0.
  mb = 0.
  mc = 0.
  md = 0.
  me = 0.
  mf = 0.
  mg = 0.

  HH_star = 0.
  HH = 0.
  QQ = 0.

  !  Nstart = 1
  !  MEMO = 8


  starttime = MPI_WTIME()



!!! read porosity
  Call ReadDensity(ne(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
       Trim(IN_DIR)//'/'//Trim(PMname))

  ne(1,:,:) = 0.
  ne(Nx,:,:) = 0.
  ne(:,1,:) = 0.
  ne(:,Ny,:) = 0.

  wall = .false. 
  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if (ne(i,j,k) .eq. 0.) wall(i,j,k)=.true. 
        enddo
     enddo
  enddo


!!! exchange porosity and wall at interface

  call mpi_sendrecv(ne(:,:,Nz),ilenz,MPI_Double_precision,inext,810, &
       tmp_ne(:,:,1),ilenz,MPI_Double_precision,iprev,810,nallgrp,stat,error)
  call mpi_sendrecv(ne(:,:,1),ilenz,MPI_Double_precision,iprev,820, &
       tmp_ne(:,:,2),ilenz,MPI_Double_precision,inext,820,nallgrp,stat,error)

  call mpi_sendrecv(wall(:,:,Nz),ilenz,MPI_LOGICAL,inext,310, &
       tmp_wall(:,:,1),ilenz,MPI_LOGICAL,iprev,310,nallgrp,stat,error)
  call mpi_sendrecv(wall(:,:,1),ilenz,MPI_LOGICAL,iprev,320, &
       tmp_wall(:,:,2),ilenz,MPI_LOGICAL,inext,320,nallgrp,stat,error)

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)





!!$  Porlocal = sum(ne)
!!$  CALL MPI_REDUCE(Porlocal,Porsum,1,MPI_Double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
!!$
!!$  mean_por = Porsum/(LX-2)/(LY-2)/LZ


  Print*,'mype=',mype,'porosity---------','nxx=',nxx,'nyy=',nyy
  Print*,ne(nxx,nyy,:)


!!$  IF (mype .eq. 0 ) THEN
!!$     print*,'##'
!!$     print*,'average porosity  =', mean_por
!!$     Print*,'porosity---------','nxx=',nxx,'nyy=',nyy
!!$     Print*,ne(nxx,nyy,:)
!!$  Endif


!!! read initial velocity and pressure
  if(t0.eq.0) then
     CALL ReadDensity(u_o(:,:,:),Nz,Nx+1,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_u')
     CALL ReadDensity(v_o(:,:,:),Nz,Nx,Ny+1,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_v')
     CALL ReadDensity(w_o(:,:,:),Nz+1,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_w')
     CALL ReadDensity(h_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_h')
     CALL ReadDensity(Se(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_Se')
  else

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_u_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(u_o(:,:,:),Nz,Nx+1,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_v_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(v_o(:,:,:),Nz,Nx,Ny+1,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_w_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(w_o(:,:,:),Nz+1,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_h_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(h_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_Se_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(Se(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))



  endif

  IF (mype .eq. 0 ) THEN
     print*,'##'
     print*,'initilization done'
  Endif


!!! calculate irreducible water saturation for each voxel

  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if (ne(i,j,k) == 0) then
              Sm(i,j,k) = 0.
           else if ((ne(i,j,k)>0).and.(ne(i,j,k)<0.1)) then
              Sm(i,j,k) = 1+(Cm/0.0063-10)*ne(i,j,k) 
           else if ((ne(i,j,k)>=0.1).and.(ne(i,j,k)<0.9)) then
              Sm(i,j,k) = Cm/ne(i,j,k)**1.2
           else 
              Sm(i,j,k) = 0.02+(0.2-Cm/0.088)*(ne(i,j,k)-1) !2% water vapor hold in pore space (Sm=0.02) 
!              Sm(i,j,k) = 0.01+(0.1-Cm/0.088)*(ne(i,j,k)-1) !1% water vapor hold in pore space (Sm=0.01) 
!              Sm(i,j,k) = 0.001+(0.01-Cm/0.088)*(ne(i,j,k)-1)
           endif
        enddo
     enddo
  enddo





!!! calculate parameters for storage coefficients

  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if (wall(i,j,k)) then
              Se(i,j,k) = 0.
              S(i,j,k) = 0.
              theta(i,j,k) = 0.
              alpha_prime(i,j,k) = 0.
           else
              S(i,j,k) = Se(i,j,k)*(1.-Sm(i,j,k)) + Sm(i,j,k)
              theta(i,j,k) = ne(i,j,k)*S(i,j,k)
              alpha_prime(i,j,k) = alpha*(1.-ne(i,j,k))*rho*g0     
           endif
        enddo
     enddo
  enddo


  beta_prime = beta*rho*g0




!!! calculate parameters for van-genuchten model
  
  a0 = 2.35
  a1 = 0.0822
  a2 = -0.497
  a3 = 1.238

  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if (wall(i,j,k)) then
              m(i,j,k) = 0.
              tau(i,j,k) = 0.
              lambda(i,j,k) = 0.
           else
!!$              cc_eff = cc * (1-ne(i,j,k))
!!$              DD = a0 + (1-exp(a1*cc_eff))/(a2*(1+exp(a1*cc_eff))+a3*(1-exp(a1*cc_eff)))
!!$              m(i,j,k) = (3-DD)/(4-DD)
              m(i,j,k) = m0
              tau(i,j,k) = 1./(1-m(i,j,k))
              sx = 0.72 - 0.35*exp(-tau(i,j,k)**4)
              h_cc = 0.45/d10*(1-ne(i,j,k))/ne(i,j,k)*0.1    !bubling pressure (Kpa)
              lambda(i,j,k) = sx**(1./(tau(i,j,k)-1))/(h_cc+hcc_s)*(sx**(-1./m(i,j,k))-1)**(1-m(i,j,k))*10    !unit (1/m), that is why *10 
           endif
        enddo
     enddo
  enddo


!notice 
print*,'m',maxval(m),minval(m)
print*,'tau',maxval(tau),minval(tau)
print*,'lambda',maxval(lambda),minval(lambda)



!!! exchange water content and m at interface

  call mpi_sendrecv(theta(:,:,Nz),ilenz,MPI_Double_precision,inext,910, &
       tmp_theta(:,:,1),ilenz,MPI_Double_precision,iprev,910,nallgrp,stat,error)
  call mpi_sendrecv(theta(:,:,1),ilenz,MPI_Double_precision,iprev,920, &
       tmp_theta(:,:,2),ilenz,MPI_Double_precision,inext,920,nallgrp,stat,error)

  call mpi_sendrecv(m(:,:,Nz),ilenz,MPI_Double_precision,inext,930, &
       tmp_m(:,:,1),ilenz,MPI_Double_precision,iprev,930,nallgrp,stat,error)
  call mpi_sendrecv(m(:,:,1),ilenz,MPI_Double_precision,iprev,940, &
       tmp_m(:,:,2),ilenz,MPI_Double_precision,inext,940,nallgrp,stat,error)


  CALL MPI_BARRIER(MPI_COMM_WORLD,error)





!!! velocity and pressure boundary condition setting

  if(mype.eq.0) then
     if (vbc_d.eq.1) then
        u_o(:,:,1) = u_d
        v_o(:,:,1) = v_d
     endif
     if (vbc_d.eq.2) then
        u_o(:,:,1) = 0.
        v_o(:,:,1) = 0.
     endif
     Do i=1,Nx
        Do j=1,Ny
           if(wall(i,j,1)) then
              w_o(i,j,1) = 0.
           else
              if (vbc_d.eq.1) then
                 w_o(i,j,1) = w_d
                 w_o(i,j,2) = w_d
              endif
              if (vbc_d.eq.2) then
                 w_o(i,j,1) = w_o(i,j,2)*(theta(i,j,1)+theta(i,j,2))/2/theta(i,j,1)
              endif
           endif
        enddo
     enddo

     if (hbc_d.eq.1) then
        h_o(:,:,1) = h_d
     endif
     if (hbc_d.eq.2) then
        h_o(:,:,1) = h_o(:,:,2)
     endif

  endif




  if(mype.eq.npes-1) then
     if (vbc_t.eq.1) then
        u_o(:,:,Nz) = u_t
        v_o(:,:,Nz) = v_t
     endif
     if (vbc_t.eq.2) then
        u_o(:,:,Nz) = 0.
        v_o(:,:,Nz) = 0.
     endif

     Do i=1,Nx
        Do j=1,Ny
           if(wall(i,j,Nz)) then
              w_o(i,j,Nz+1) = 0.
           else
              if (vbc_t.eq.1) then
!notice   this setting make sure no water enter the soil when w_t=0.
                 w_o(i,j,Nz+1) = w_t
                 w_o(i,j,Nz) = w_t  
              endif
              if (vbc_t.eq.2) then
                 w_o(i,j,Nz+1) = w_o(i,j,Nz)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./theta(i,j,Nz)
              endif
           endif
        enddo
     enddo

     if (hbc_t.eq.1) then
        h_o(:,:,Nz) = h_t
     endif
     if (hbc_t.eq.2) then
        h_o(:,:,Nz) =  h_o(:,:,Nz-1)
     endif

  endif






  !! wall location and pore area
  Do i=1,Nx
     Do j=1,Ny
        Do k=1,Nz
           if(wall(i,j,k)) then
              u_o(i,j,k) = 0.
              u_o(i+1,j,k) = 0.
              v_o(i,j,k) = 0.
              v_o(i,j+1,k) = 0.
              w_o(i,j,k) = 0.
              w_o(i,j,k+1) = 0.
              h_o(i,j,k) = 0.
           endif
        Enddo
     Enddo
  Enddo



  IF (mype.eq.0) THEN
     CALL WriteMatlab()
     print*,'begin computation'
  ENDIF





!!!************* begin to compute ************

  Do istep=t0+1,Max_step


     IF(mod(istep,Print_step).eq.0) THEN
        IF (mype .eq. 0 ) THEN
           print*,'-------------------------'
           print*,'istep=', istep
        Endif
     endif



!!! boundary condition setting and
!!! communicate info between bottom and top planes



     call mpi_sendrecv(u_o(:,:,Nz),ilenx,MPI_double_precision,inext,110, &
          tmp_u_o(:,:,1),ilenx,MPI_double_precision,iprev,110,nallgrp,stat,error)

     call mpi_sendrecv(u_o(:,:,1),ilenx,MPI_double_precision,iprev,120, &
          tmp_u_o(:,:,2),ilenx,MPI_double_precision,inext,120,nallgrp,stat,error)

     call mpi_sendrecv(v_o(:,:,Nz),ileny,MPI_double_precision,inext,130, &
          tmp_v_o(:,:,1),ileny,MPI_double_precision,iprev,130,nallgrp,stat,error)

     call mpi_sendrecv(v_o(:,:,1),ileny,MPI_double_precision,iprev,140, &
          tmp_v_o(:,:,2),ileny,MPI_double_precision,inext,140,nallgrp,stat,error)

     call mpi_sendrecv(w_o(:,:,Nz),ilenz,MPI_double_precision,inext,150, &
          tmp_w_o(:,:,1),ilenz,MPI_double_precision,iprev,150,nallgrp,stat,error)

     call mpi_sendrecv(w_o(:,:,2),ilenz,MPI_double_precision,iprev,160, &
          tmp_w_o(:,:,2),ilenz,MPI_double_precision,inext,160,nallgrp,stat,error)


!!$     CALL MPI_BARRIER(MPI_COMM_WORLD,error)






!!!***** solve u_star, v_star, w_star *******


!!! for u_star
     Do i=2,Nx  
        Do j=2,Ny-1
           Do k=1,Nz
              if(wall(i-1,j,k).or.wall(i,j,k)) then
                 u_star(i,j,k) = 0.
              else
                 if(wall(i,j-1,k)) then
                    u_o(i,j-1,k)=-u_o(i,j,k)
                 endif
                 if(wall(i,j+1,k)) then
                    u_o(i,j+1,k)=-u_o(i,j,k)
                 endif

                 if(k.eq.1) then
                    if(tmp_wall(i,j,1)) then
                       tmp_u_o(i,j,1)=-u_o(i,j,1)
                    endif
                    if(wall(i,j,1+1)) then
                       u_o(i,j,1+1)=-u_o(i,j,1)
                    endif
                    u_star(i,j,1) = u_o(i,j,1) - dt/(0.5*(theta(i-1,j,1)+theta(i,j,1)))* &
                         (-0.5*(u_o(i-1,j,1)+u_o(i,j,1))*0.5*(u_o(i-1,j,1)+u_o(i,j,1))*theta(i-1,j,1)/dx + &
                         0.5*(u_o(i,j,1)+u_o(i+1,j,1))*0.5*(u_o(i,j,1)+u_o(i+1,j,1))*theta(i,j,1)/dx - &
                         0.5*(v_o(i-1,j,1)+v_o(i,j,1))*0.5*(u_o(i,j-1,1)+u_o(i,j,1))* &
                         0.25*(theta(i-1,j-1,1)+theta(i-1,j,1)+theta(i,j-1,1)+theta(i,j,1))/dy + &
                         0.5*(v_o(i-1,j+1,1)+v_o(i,j+1,1))*0.5*(u_o(i,j,1)+u_o(i,j+1,1))* &
                         0.25*(theta(i-1,j,1)+theta(i-1,j+1,1)+theta(i,j,1)+theta(i,j+1,1))/dy - &
                         0.5*(w_o(i-1,j,1)+w_o(i,j,1))*0.5*(tmp_u_o(i,j,1)+u_o(i,j,1))* &
                         0.25*(tmp_theta(i-1,j,1)+theta(i-1,j,1)+tmp_theta(i,j,1)+theta(i,j,1))/dz + &
                         0.5*(w_o(i-1,j,1+1)+w_o(i,j,1+1))*0.5*(u_o(i,j,1)+u_o(i,j,1+1))* &
                         0.25*(theta(i-1,j,1)+theta(i-1,j,1+1)+theta(i,j,1)+theta(i,j,1+1))/dz) + &
                         dt/(0.5*(theta(i-1,j,1)+theta(i,j,1)))*nu*(-(u_o(i,j,1)-u_o(i-1,j,1))/dx**2*theta(i-1,j,1) + &
                         (u_o(i+1,j,1)-u_o(i,j,1))/dx**2*theta(i,j,1) - &
                         (u_o(i,j,1)-u_o(i,j-1,1))/dy**2*0.25* &
                         (theta(i-1,j-1,1)+theta(i-1,j,1)+theta(i,j-1,1)+theta(i,j,1)) + &
                         (u_o(i,j+1,1)-u_o(i,j,1))/dy**2*0.25* &
                         (theta(i-1,j,1)+theta(i-1,j+1,1)+theta(i,j,1)+theta(i,j+1,1)) - &
                         (u_o(i,j,1)-tmp_u_o(i,j,1))/dz**2*0.25* &
                         (tmp_theta(i-1,j,1)+theta(i-1,j,1)+tmp_theta(i,j,1)+theta(i,j,1)) + &
                         (u_o(i,j,1+1)-u_o(i,j,1))/dz**2*0.25* &
                         (theta(i-1,j,1)+theta(i-1,j,1+1)+theta(i,j,1)+theta(i,j,1+1))) + &
                         dt*gvt(1)
                 else if(k.eq.Nz) then
                    if(wall(i,j,Nz-1)) then
                       u_o(i,j,Nz-1)=-u_o(i,j,Nz)
                    endif
                    if(tmp_wall(i,j,2)) then
                       tmp_u_o(i,j,2)=-u_o(i,j,Nz)
                    endif
                    u_star(i,j,Nz) = u_o(i,j,Nz) - dt/(0.5*(theta(i-1,j,Nz)+theta(i,j,Nz)))* &
                         (-0.5*(u_o(i-1,j,Nz)+u_o(i,j,Nz))*0.5*(u_o(i-1,j,Nz)+u_o(i,j,Nz))*theta(i-1,j,Nz)/dx + &
                         0.5*(u_o(i,j,Nz)+u_o(i+1,j,Nz))*0.5*(u_o(i,j,Nz)+u_o(i+1,j,Nz))*theta(i,j,Nz)/dx - &
                         0.5*(v_o(i-1,j,Nz)+v_o(i,j,Nz))*0.5*(u_o(i,j-1,Nz)+u_o(i,j,Nz))* &
                         0.25*(theta(i-1,j-1,Nz)+theta(i-1,j,Nz)+theta(i,j-1,Nz)+theta(i,j,Nz))/dy + &
                         0.5*(v_o(i-1,j+1,Nz)+v_o(i,j+1,Nz))*0.5*(u_o(i,j,Nz)+u_o(i,j+1,Nz))* &
                         0.25*(theta(i-1,j,Nz)+theta(i-1,j+1,Nz)+theta(i,j,Nz)+theta(i,j+1,Nz))/dy - &
                         0.5*(w_o(i-1,j,Nz)+w_o(i,j,Nz))*0.5*(u_o(i,j,Nz-1)+u_o(i,j,Nz))* &
                         0.25*(theta(i-1,j,Nz-1)+theta(i-1,j,Nz)+theta(i,j,Nz-1)+theta(i,j,Nz))/dz + &
                         0.5*(w_o(i-1,j,Nz+1)+w_o(i,j,Nz+1))*0.5*(u_o(i,j,Nz)+tmp_u_o(i,j,2))* &
                         0.25*(theta(i-1,j,Nz)+tmp_theta(i-1,j,2)+theta(i,j,Nz)+tmp_theta(i,j,2))/dz) + &
                         dt/(0.5*(theta(i-1,j,Nz)+theta(i,j,Nz)))*nu*(-(u_o(i,j,Nz)-u_o(i-1,j,Nz))/dx**2*theta(i-1,j,Nz) + &
                         (u_o(i+1,j,Nz)-u_o(i,j,Nz))/dx**2*theta(i,j,Nz) - &
                         (u_o(i,j,Nz)-u_o(i,j-1,Nz))/dy**2*0.25* &
                         (theta(i-1,j-1,Nz)+theta(i-1,j,Nz)+theta(i,j-1,Nz)+theta(i,j,Nz)) + &
                         (u_o(i,j+1,Nz)-u_o(i,j,Nz))/dy**2*0.25* &
                         (theta(i-1,j,Nz)+theta(i-1,j+1,Nz)+theta(i,j,Nz)+theta(i,j+1,Nz)) - &
                         (u_o(i,j,Nz)-u_o(i,j,Nz-1))/dz**2*0.25* &
                         (theta(i-1,j,Nz-1)+theta(i-1,j,Nz)+theta(i,j,Nz-1)+theta(i,j,Nz)) + &
                         (tmp_u_o(i,j,2)-u_o(i,j,Nz))/dz**2*0.25* &
                         (theta(i-1,j,Nz)+tmp_theta(i-1,j,2)+theta(i,j,Nz)+tmp_theta(i,j,2))) + &
                         dt*gvt(1)
                 else                                     
                    if(wall(i,j,k-1)) then
                       u_o(i,j,k-1)=-u_o(i,j,k)
                    endif
                    if(wall(i,j,k+1)) then
                       u_o(i,j,k+1)=-u_o(i,j,k)
                    endif
                    u_star(i,j,k) = u_o(i,j,k) - dt/(0.5*(theta(i-1,j,k)+theta(i,j,k)))* &
                         (-0.5*(u_o(i-1,j,k)+u_o(i,j,k))*0.5*(u_o(i-1,j,k)+u_o(i,j,k))*theta(i-1,j,k)/dx + &
                         0.5*(u_o(i,j,k)+u_o(i+1,j,k))*0.5*(u_o(i,j,k)+u_o(i+1,j,k))*theta(i,j,k)/dx - &
                         0.5*(v_o(i-1,j,k)+v_o(i,j,k))*0.5*(u_o(i,j-1,k)+u_o(i,j,k))* &
                         0.25*(theta(i-1,j-1,k)+theta(i-1,j,k)+theta(i,j-1,k)+theta(i,j,k))/dy + &
                         0.5*(v_o(i-1,j+1,k)+v_o(i,j+1,k))*0.5*(u_o(i,j,k)+u_o(i,j+1,k))* &
                         0.25*(theta(i-1,j,k)+theta(i-1,j+1,k)+theta(i,j,k)+theta(i,j+1,k))/dy - &
                         0.5*(w_o(i-1,j,k)+w_o(i,j,k))*0.5*(u_o(i,j,k-1)+u_o(i,j,k))* &
                         0.25*(theta(i-1,j,k-1)+theta(i-1,j,k)+theta(i,j,k-1)+theta(i,j,k))/dz + &
                         0.5*(w_o(i-1,j,k+1)+w_o(i,j,k+1))*0.5*(u_o(i,j,k)+u_o(i,j,k+1))* &
                         0.25*(theta(i-1,j,k)+theta(i-1,j,k+1)+theta(i,j,k)+theta(i,j,k+1))/dz) + &
                         dt/(0.5*(theta(i-1,j,k)+theta(i,j,k)))*nu*(-(u_o(i,j,k)-u_o(i-1,j,k))/dx**2*theta(i-1,j,k) + &
                         (u_o(i+1,j,k)-u_o(i,j,k))/dx**2*theta(i,j,k) - &
                         (u_o(i,j,k)-u_o(i,j-1,k))/dy**2*0.25* &
                         (theta(i-1,j-1,k)+theta(i-1,j,k)+theta(i,j-1,k)+theta(i,j,k)) + &
                         (u_o(i,j+1,k)-u_o(i,j,k))/dy**2*0.25* &
                         (theta(i-1,j,k)+theta(i-1,j+1,k)+theta(i,j,k)+theta(i,j+1,k)) - &
                         (u_o(i,j,k)-u_o(i,j,k-1))/dz**2*0.25* &
                         (theta(i-1,j,k-1)+theta(i-1,j,k)+theta(i,j,k-1)+theta(i,j,k)) + &
                         (u_o(i,j,k+1)-u_o(i,j,k))/dz**2*0.25* &
                         (theta(i-1,j,k)+theta(i-1,j,k+1)+theta(i,j,k)+theta(i,j,k+1))) + &
                         dt*gvt(1)

                 endif
              endif
           Enddo
        Enddo
     Enddo







     if(mod(istep,Print_step).eq.0) then
        Print*,'u_star---------'
        Print*,u_star(nxx,nyy,:)
     endif




!!! for v_star
     Do i=2,Nx-1  
        Do j=2,Ny
           Do k=1,Nz
              if(wall(i,j-1,k).or.wall(i,j,k)) then
                 v_star(i,j,k) = 0.
              else
                 if(wall(i-1,j,k)) then
                    v_o(i-1,j,k)=-v_o(i,j,k)
                 endif
                 if(wall(i+1,j,k)) then
                    v_o(i+1,j,k)=-v_o(i,j,k)
                 endif

                 if(k.eq.1) then
                    if(tmp_wall(i,j,1)) then
                       tmp_v_o(i,j,1)=-v_o(i,j,1)
                    endif
                    if(wall(i,j,1+1)) then
                       v_o(i,j,1+1)=-v_o(i,j,1)
                    endif
                    v_star(i,j,1) = v_o(i,j,1) - dt/(0.5*(theta(i,j-1,1)+theta(i,j,1)))* &
                         (-0.5*(u_o(i,j-1,1)+u_o(i,j,1))*0.5*(v_o(i-1,j,1)+v_o(i,j,1))* &
                         0.25*(theta(i-1,j-1,1)+theta(i-1,j,1)+theta(i,j-1,1)+theta(i,j-1,1))/dx + &
                         0.5*(u_o(i+1,j-1,1)+u_o(i+1,j,1))*0.5*(v_o(i,j,1)+v_o(i+1,j,1))* &
                         0.25*(theta(i,j-1,1)+theta(i,j,1)+theta(i+1,j-1,1)+theta(i+1,j,1))/dx - &
                         0.5*(v_o(i,j-1,1)+v_o(i,j,1))*0.5*(v_o(i,j-1,1)+v_o(i,j,1))*theta(i,j-1,1)/dy + &
                         0.5*(v_o(i,j,1)+v_o(i,j+1,1))*0.5*(v_o(i,j,1)+v_o(i,j+1,1))*theta(i,j,1)/dy - &
                         0.5*(w_o(i,j-1,1)+w_o(i,j,1))*0.5*(tmp_v_o(i,j,1)+v_o(i,j,1))* &
                         0.25*(tmp_theta(i,j-1,1)+theta(i,j-1,1)+tmp_theta(i,j,1)+theta(i,j,1))/dz + &
                         0.5*(w_o(i,j-1,1+1)+w_o(i,j,1+1))*0.5*(v_o(i,j,1)+v_o(i,j,1+1))* &
                         0.25*(theta(i,j-1,1)+theta(i,j-1,1+1)+theta(i,j,1)+theta(i,j,1+1))/dz) + &
                         dt/(0.5*(theta(i,j-1,1)+theta(i,j,1)))*nu*(-(v_o(i,j,1)-v_o(i-1,j,1))/dx**2* &
                         0.25*(theta(i-1,j-1,1)+theta(i-1,j,1)+theta(i,j-1,1)+theta(i,j,1)) + &
                         (v_o(i+1,j,1)-v_o(i,j,1))/dx**2*0.25*(theta(i,j,1)+theta(i,j-1,1)+theta(i+1,j-1,1)+theta(i+1,j,1)) - &
                         (v_o(i,j,1)-v_o(i,j-1,1))/dy**2*theta(i,j-1,1) + &
                         (v_o(i,j+1,1)-v_o(i,j,1))/dy**2*theta(i,j,1) - &
                         (v_o(i,j,1)-tmp_v_o(i,j,1))/dz**2*0.25*(tmp_theta(i,j-1,1)+theta(i,j-1,1)+tmp_theta(i,j,1)+theta(i,j,1))+&
                         (v_o(i,j,1+1)-v_o(i,j,1))/dz**2*0.25*(theta(i,j-1,1)+theta(i,j-1,1+1)+theta(i,j,1)+theta(i,j,1+1))) + &
                         dt*gvt(2)
                 else if(k.eq.Nz) then
                    if(wall(i,j,Nz-1)) then
                       v_o(i,j,Nz-1)=-v_o(i,j,Nz)
                    endif
                    if(tmp_wall(i,j,2)) then
                       tmp_v_o(i,j,2)=-v_o(i,j,Nz)
                    endif
                    v_star(i,j,Nz) = v_o(i,j,Nz) - dt/(0.5*(theta(i,j-1,Nz)+theta(i,j,Nz)))* &
                         (-0.5*(u_o(i,j-1,Nz)+u_o(i,j,Nz))*0.5*(v_o(i-1,j,Nz)+v_o(i,j,Nz))* &
                         0.25*(theta(i-1,j-1,Nz)+theta(i-1,j,Nz)+theta(i,j-1,Nz)+theta(i,j-1,Nz))/dx + &
                         0.5*(u_o(i+1,j-1,Nz)+u_o(i+1,j,Nz))*0.5*(v_o(i,j,Nz)+v_o(i+1,j,Nz))* &
                         0.25*(theta(i,j-1,Nz)+theta(i,j,Nz)+theta(i+1,j-1,Nz)+theta(i+1,j,Nz))/dx - &
                         0.5*(v_o(i,j-1,Nz)+v_o(i,j,Nz))*0.5*(v_o(i,j-1,Nz)+v_o(i,j,Nz))*theta(i,j-1,Nz)/dy + &
                         0.5*(v_o(i,j,Nz)+v_o(i,j+1,Nz))*0.5*(v_o(i,j,Nz)+v_o(i,j+1,Nz))*theta(i,j,Nz)/dy - &
                         0.5*(w_o(i,j-1,Nz)+w_o(i,j,Nz))*0.5*(v_o(i,j,Nz-1)+v_o(i,j,Nz))* &
                         0.25*(theta(i,j-1,Nz-1)+theta(i,j-1,Nz)+theta(i,j,Nz-1)+theta(i,j,Nz))/dz + &
                         0.5*(w_o(i,j-1,Nz+1)+w_o(i,j,Nz+1))*0.5*(v_o(i,j,Nz)+tmp_v_o(i,j,2))* &
                         0.25*(theta(i,j-1,Nz)+tmp_theta(i,j-1,2)+theta(i,j,Nz)+tmp_theta(i,j,2))/dz) + &
                         dt/(0.5*(theta(i,j-1,Nz)+theta(i,j,Nz)))*nu*(-(v_o(i,j,Nz)-v_o(i-1,j,Nz))/dx**2* &
                         0.25*(theta(i-1,j-1,Nz)+theta(i-1,j,Nz)+theta(i,j-1,Nz)+theta(i,j,Nz)) + &
                         (v_o(i+1,j,Nz)-v_o(i,j,Nz))/dx**2*0.25* &
                         (theta(i,j,Nz)+theta(i,j-1,Nz)+theta(i+1,j-1,Nz)+theta(i+1,j,Nz)) - &
                         (v_o(i,j,Nz)-v_o(i,j-1,Nz))/dy**2*theta(i,j-1,Nz) + &
                         (v_o(i,j+1,Nz)-v_o(i,j,Nz))/dy**2*theta(i,j,Nz) - &
                         (v_o(i,j,Nz)-v_o(i,j,Nz-1))/dz**2*0.25* &
                         (theta(i,j-1,Nz-1)+theta(i,j-1,Nz)+theta(i,j,Nz-1)+theta(i,j,Nz)) + &
                         (tmp_v_o(i,j,2)-v_o(i,j,Nz))/dz**2*0.25* &
                         (theta(i,j-1,Nz)+tmp_theta(i,j-1,2)+theta(i,j,Nz)+tmp_theta(i,j,2))) + &
                         dt*gvt(2)
                 else
                    if(wall(i,j,k-1)) then
                       v_o(i,j,k-1)=-v_o(i,j,k)
                    endif
                    if(wall(i,j,k+1)) then
                       v_o(i,j,k+1)=-v_o(i,j,k)
                    endif
                    v_star(i,j,k) = v_o(i,j,k) - dt/(0.5*(theta(i,j-1,k)+theta(i,j,k)))* &
                         (-0.5*(u_o(i,j-1,k)+u_o(i,j,k))*0.5*(v_o(i-1,j,k)+v_o(i,j,k))* &
                         0.25*(theta(i-1,j-1,k)+theta(i-1,j,k)+theta(i,j-1,k)+theta(i,j-1,k))/dx + &
                         0.5*(u_o(i+1,j-1,k)+u_o(i+1,j,k))*0.5*(v_o(i,j,k)+v_o(i+1,j,k))* &
                         0.25*(theta(i,j-1,k)+theta(i,j,k)+theta(i+1,j-1,k)+theta(i+1,j,k))/dx - &
                         0.5*(v_o(i,j-1,k)+v_o(i,j,k))*0.5*(v_o(i,j-1,k)+v_o(i,j,k))*theta(i,j-1,k)/dy + &
                         0.5*(v_o(i,j,k)+v_o(i,j+1,k))*0.5*(v_o(i,j,k)+v_o(i,j+1,k))*theta(i,j,k)/dy - &
                         0.5*(w_o(i,j-1,k)+w_o(i,j,k))*0.5*(v_o(i,j,k-1)+v_o(i,j,k))* &
                         0.25*(theta(i,j-1,k-1)+theta(i,j-1,k)+theta(i,j,k-1)+theta(i,j,k))/dz + &
                         0.5*(w_o(i,j-1,k+1)+w_o(i,j,k+1))*0.5*(v_o(i,j,k)+v_o(i,j,k+1))* &
                         0.25*(theta(i,j-1,k)+theta(i,j-1,k+1)+theta(i,j,k)+theta(i,j,k+1))/dz) + &
                         dt/(0.5*(theta(i,j-1,k)+theta(i,j,k)))*nu*(-(v_o(i,j,k)-v_o(i-1,j,k))/dx**2* &
                         0.25*(theta(i-1,j-1,k)+theta(i-1,j,k)+theta(i,j-1,k)+theta(i,j,k)) + &
                         (v_o(i+1,j,k)-v_o(i,j,k))/dx**2*0.25* &
                         (theta(i,j,k)+theta(i,j-1,k)+theta(i+1,j-1,k)+theta(i+1,j,k)) - &
                         (v_o(i,j,k)-v_o(i,j-1,k))/dy**2*theta(i,j-1,k) + &
                         (v_o(i,j+1,k)-v_o(i,j,k))/dy**2*theta(i,j,k) - &
                         (v_o(i,j,k)-v_o(i,j,k-1))/dz**2*0.25* &
                         (theta(i,j-1,k-1)+theta(i,j-1,k)+theta(i,j,k-1)+theta(i,j,k)) + &
                         (v_o(i,j,k+1)-v_o(i,j,k))/dz**2*0.25* &
                         (theta(i,j-1,k)+theta(i,j-1,k+1)+theta(i,j,k)+theta(i,j,k+1))) + &
                         dt*gvt(2)
                 endif
              endif
           Enddo
        Enddo
     Enddo





     if(mod(istep,Print_step).eq.0) then
        Print*,'v_star---------'
        Print*,v_star(nxx,nyy,:)
     endif


!!$     if(mod(istep,Print_step).eq.0) then
!!$        Print*,'gvt'
!!$        Print*,gvt(:)
!!$     endif




!!! for w_star
     Do i=2,Nx-1  
        Do j=2,Ny-1
           Do k=1,Nz+1 
              if(k.eq.1) then
                 if(tmp_wall(i,j,1).or.wall(i,j,1)) then
                    w_star(i,j,1) = 0.
                 else
                    if(wall(i-1,j,1)) then
                       w_o(i-1,j,1)=-w_o(i,j,1)
                    endif
                    if(wall(i+1,j,1)) then
                       w_o(i+1,j,1)=-w_o(i,j,1)
                    endif
                    if(wall(i,j-1,1)) then
                       w_o(i,j-1,1)=-w_o(i,j,1)
                    endif
                    if(wall(i,j+1,1)) then
                       w_o(i,j+1,1)=-w_o(i,j,1)
                    endif
                    w_star(i,j,1) = w_o(i,j,1) - dt/(0.5*(tmp_theta(i,j,1)+theta(i,j,1)))* &
                         (-0.5*(tmp_u_o(i,j,1)+u_o(i,j,1))*0.5*(w_o(i-1,j,1)+w_o(i,j,1))* &
                         0.25*(tmp_theta(i-1,j,1)+theta(i-1,j,1)+tmp_theta(i,j,1)+theta(i,j,1))/dx + &
                         0.5*(w_o(i,j,1)+w_o(i,j,1+1))*0.5*(w_o(i,j,1)+w_o(i,j,1+1))*theta(i,j,1)/dz) + &
                         dt/(0.5*(tmp_theta(i,j,1)+theta(i,j,1)))*nu*(-(w_o(i,j,1)-w_o(i-1,j,1))/dx**2* &
                         0.25*(tmp_theta(i-1,j,1)+theta(i-1,j,1)+tmp_theta(i,j,1)+theta(i,j,1)) + &
                         (w_o(i+1,j,1)-w_o(i,j,1))/dx**2*0.25* &
                         (tmp_theta(i,j,1)+theta(i,j,1)+tmp_theta(i+1,j,1)+theta(i+1,j,1)) - &
                         (w_o(i,j,1)-w_o(i,j-1,1))/dy**2*0.25* &
                         (tmp_theta(i,j-1,1)+theta(i,j-1,1)+tmp_theta(i,j,1)+theta(i,j,1)) + &
                         (w_o(i,j+1,1)-w_o(i,j,1))/dy**2*0.25* &
                         (tmp_theta(i,j,1)+theta(i,j,1)+tmp_theta(i,j+1,1)+theta(i,j+1,1)) - &
                         (w_o(i,j,1)-tmp_w_o(i,j,1))/dz**2*tmp_theta(i,j,1) + &
                         (w_o(i,j,1+1)-w_o(i,j,1))/dz**2*theta(i,j,1)) + &
                         dt*gvt(3)
                 endif
              else if(k.eq.Nz+1) then
                 if(wall(i,j,Nz).or.tmp_wall(i,j,2)) then
                    w_star(i,j,Nz+1) = 0.
                 else
                    if(tmp_wall(i-1,j,2)) then
                       w_o(i-1,j,Nz+1)=-w_o(i,j,Nz+1)
                    endif
                    if(tmp_wall(i+1,j,2)) then
                       w_o(i+1,j,Nz+1)=-w_o(i,j,Nz+1)
                    endif
                    if(tmp_wall(i,j-1,2)) then
                       w_o(i,j-1,Nz+1)=-w_o(i,j,Nz+1)
                    endif
                    if(tmp_wall(i,j+1,2)) then
                       w_o(i,j+1,Nz+1)=-w_o(i,j,Nz+1)
                    endif
                    w_star(i,j,Nz+1) = w_o(i,j,Nz+1) - dt/(0.5*(theta(i,j,Nz)+tmp_theta(i,j,2)))* &
                         (-0.5*(u_o(i,j,Nz)+tmp_u_o(i,j,2))*0.5*(w_o(i-1,j,Nz+1)+w_o(i,j,Nz+1))* &
                         0.25*(theta(i-1,j,Nz)+tmp_theta(i-1,j,2)+theta(i,j,Nz)+tmp_theta(i,j,2))/dx + &
                         0.5*(u_o(i+1,j,Nz)+tmp_u_o(i+1,j,2))*0.5*(w_o(i,j,Nz+1)+w_o(i+1,j,Nz+1))* &
                         0.25*(theta(i,j,Nz)+tmp_theta(i,j,2)+theta(i+1,j,Nz)+tmp_theta(i+1,j,2))/dx - &
                         0.5*(v_o(i,j,Nz)+tmp_v_o(i,j,2))*0.5*(w_o(i,j-1,Nz+1)+w_o(i,j,Nz+1))* &
                         0.25*(theta(i,j-1,Nz)+tmp_theta(i,j-1,2)+theta(i,j,Nz)+tmp_theta(i,j,2))/dy + &
                         0.5*(v_o(i,j+1,Nz)+tmp_v_o(i,j+1,2))*0.5*(w_o(i,j,Nz+1)+w_o(i,j+1,Nz+1))* &
                         0.25*(theta(i,j,Nz)+tmp_theta(i,j,2)+theta(i,j+1,Nz)+tmp_theta(i,j+1,2))/dy - &
                         0.5*(w_o(i,j,Nz)+w_o(i,j,Nz+1))*0.5*(w_o(i,j,Nz)+w_o(i,j,Nz+1))*theta(i,j,Nz)/dz + &
                         0.5*(w_o(i,j,Nz+1)+tmp_w_o(i,j,2))*0.5*(w_o(i,j,Nz+1)+tmp_w_o(i,j,2))*tmp_theta(i,j,2)/dz) + &
                         dt/(0.5*(theta(i,j,Nz)+tmp_theta(i,j,2)))*nu*(-(w_o(i,j,Nz+1)-w_o(i-1,j,Nz+1))/dx**2* &
                         0.25*(theta(i-1,j,Nz)+tmp_theta(i-1,j,2)+theta(i,j,Nz)+tmp_theta(i,j,2)) + &
                         (w_o(i+1,j,Nz+1)-w_o(i,j,Nz+1))/dx**2*0.25* &
                         (theta(i,j,Nz)+tmp_theta(i,j,2)+theta(i+1,j,Nz)+tmp_theta(i+1,j,2)) - &
                         (w_o(i,j,Nz+1)-w_o(i,j-1,Nz+1))/dy**2*0.25* &
                         (theta(i,j-1,Nz)+tmp_theta(i,j-1,2)+theta(i,j,Nz)+tmp_theta(i,j,2)) + &
                         (w_o(i,j+1,Nz+1)-w_o(i,j,Nz+1))/dy**2*0.25* &
                         (theta(i,j,Nz)+tmp_theta(i,j,2)+theta(i,j+1,Nz)+tmp_theta(i,j+1,2)) - &
                         (w_o(i,j,Nz+1)-w_o(i,j,Nz))/dz**2*theta(i,j,Nz) + &
                         (tmp_w_o(i,j,2)-w_o(i,j,Nz+1))/dz**2*tmp_theta(i,j,2)) + &
                         dt*gvt(3)
                 endif
              else
                 if(wall(i,j,k-1).or.wall(i,j,k)) then
                    w_star(i,j,k) = 0.
                 else
                    if(wall(i-1,j,k)) then
                       w_o(i-1,j,k)=-w_o(i,j,k)
                    endif
                    if(wall(i+1,j,k)) then
                       w_o(i+1,j,k)=-w_o(i,j,k)
                    endif
                    if(wall(i,j-1,k)) then
                       w_o(i,j-1,k)=-w_o(i,j,k)
                    endif
                    if(wall(i,j+1,k)) then
                       w_o(i,j+1,k)=-w_o(i,j,k)
                    endif
                    w_star(i,j,k) = w_o(i,j,k) - dt/(0.5*(theta(i,j,k-1)+theta(i,j,k)))* &
                         (-0.5*(u_o(i,j,k-1)+u_o(i,j,k))*0.5*(w_o(i-1,j,k)+w_o(i,j,k))* &
                         0.25*(theta(i-1,j,k-1)+theta(i-1,j,k)+theta(i,j,k-1)+theta(i,j,k))/dx + &
                         0.5*(u_o(i+1,j,k-1)+u_o(i+1,j,k))*0.5*(w_o(i,j,k)+w_o(i+1,j,k))* &
                         0.25*(theta(i,j,k-1)+theta(i,j,k)+theta(i+1,j,k-1)+theta(i+1,j,k))/dx - &
                         0.5*(v_o(i,j,k-1)+v_o(i,j,k))*0.5*(w_o(i,j-1,k)+w_o(i,j,k))* &
                         0.25*(theta(i,j-1,k-1)+theta(i,j-1,k)+theta(i,j,k-1)+theta(i,j,k))/dy + &
                         0.5*(v_o(i,j+1,k-1)+v_o(i,j+1,k))*0.5*(w_o(i,j,k)+w_o(i,j+1,k))* &
                         0.25*(theta(i,j,k-1)+theta(i,j,k)+theta(i,j+1,k-1)+theta(i,j+1,k))/dy - &
                         0.5*(w_o(i,j,k-1)+w_o(i,j,k))*0.5*(w_o(i,j,k-1)+w_o(i,j,k))*theta(i,j,k-1)/dz + &
                         0.5*(w_o(i,j,k)+w_o(i,j,k+1))*0.5*(w_o(i,j,k)+w_o(i,j,k+1))*theta(i,j,k)/dz) + &
                         dt/(0.5*(theta(i,j,k-1)+theta(i,j,k)))*nu*(-(w_o(i,j,k)-w_o(i-1,j,k))/dx**2* &
                         0.25*(theta(i-1,j,k-1)+theta(i-1,j,k)+theta(i,j,k-1)+theta(i,j,k)) + &
                         (w_o(i+1,j,k)-w_o(i,j,k))/dx**2*0.25* &
                         (theta(i,j,k-1)+theta(i,j,k)+theta(i+1,j,k-1)+theta(i+1,j,k)) - &
                         (w_o(i,j,k)-w_o(i,j-1,k))/dy**2*0.25* &
                         (theta(i,j-1,k-1)+theta(i,j-1,k)+theta(i,j,k-1)+theta(i,j,k)) + &
                         (w_o(i,j+1,k)-w_o(i,j,k))/dy**2*0.25* &
                         (theta(i,j,k-1)+theta(i,j,k)+theta(i,j+1,k-1)+theta(i,j+1,k)) - &
                         (w_o(i,j,k)-w_o(i,j,k-1))/dz**2*theta(i,j,k-1) + &
                         (w_o(i,j,k+1)-w_o(i,j,k))/dz**2*theta(i,j,k)) + &
                         dt*gvt(3)
                 endif
              endif
           Enddo
        Enddo
     Enddo


!!$     if(mod(istep,Print_step).eq.0) then
!!$        Print*,'w_star---------'
!!$        Print*,w_star(nxx,nyy,:)
!!$     endif






!!! velocity and pressure boundary condition setting

     if(mype.eq.0) then
        if (vbc_d.eq.1) then
           u_star(:,:,1) = u_d
           v_star(:,:,1) = v_d
        endif
        if (vbc_d.eq.2) then
           u_star(:,:,1) = 0.
           v_star(:,:,1) = 0.
        endif
        Do i=1,Nx
           Do j=1,Ny
              if(wall(i,j,1)) then
                 w_star(i,j,1) = 0.
              else
                 if (vbc_d.eq.1) then
                    w_star(i,j,1) = w_d
                    w_star(i,j,2) = w_d
                 endif
                 if (vbc_d.eq.2) then
                    w_star(i,j,1) = w_star(i,j,2)*(theta(i,j,1)+theta(i,j,2))/2/theta(i,j,1)
                 endif
              endif
           enddo
        enddo

     endif




     if(mype.eq.npes-1) then
        if (vbc_t.eq.1) then
           u_star(:,:,Nz) = u_t
           v_star(:,:,Nz) = v_t
        endif
        if (vbc_t.eq.2) then
           u_star(:,:,Nz) = 0.
           v_star(:,:,Nz) = 0.
        endif

        Do i=1,Nx
           Do j=1,Ny
              if(wall(i,j,Nz)) then
                 w_star(i,j,Nz+1) = 0.
              else
                 if (vbc_t.eq.1) then
!notice   this setting make sure no water enter the soil when w_t=0.
                    w_star(i,j,Nz+1) = w_t
                    w_star(i,j,Nz) = w_t
                 endif
                 if (vbc_t.eq.2) then
                    w_star(i,j,Nz+1) = w_star(i,j,Nz)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./theta(i,j,Nz)
                 endif
              endif
           enddo
        enddo

     endif







     !! wall location and pore area
     Do i=1,Nx
        Do j=1,Ny
           Do k=1,Nz
              if(wall(i,j,k)) then
                 u_star(i,j,k) = 0.
                 u_star(i+1,j,k) = 0.
                 v_star(i,j,k) = 0.
                 v_star(i,j+1,k) = 0.
                 w_star(i,j,k) = 0.
                 w_star(i,j,k+1) = 0.
              endif
           Enddo
        Enddo
     Enddo




     if(mod(istep,Print_step).eq.0) then
        Print*,'mype=',mype,'w_star---------'
        Print*,w_star(nxx,nyy,:)
     endif








!!!********* solve h **********

     h_star = h_o
     hijk = 1

     !! outer iteration

     DO !while(hijk.lt.1e4)


        !! calculate Se,S
        Do i=1,Nx
           Do j=1,Ny
              Do k=1,Nz
                 if(.not.wall(i,j,k)) then
                    if(h_star(i,j,k).gt.0.0) then
                       Se(i,j,k) = 1.0
                    else
                       Se(i,j,k) = (1.+(lambda(i,j,k)*Abs(h_star(i,j,k)))**tau(i,j,k))**(-m(i,j,k))
                    endif
                    S(i,j,k) = Se(i,j,k)*(1.-Sm(i,j,k)) + Sm(i,j,k)
                    theta(i,j,k) = ne(i,j,k)*S(i,j,k)
                 endif
              Enddo
           Enddo
        Enddo


        !! calculate KK,KK_inv



        Do i=2,Nx
           Do j=1,Ny
              Do k=1,Nz
                 halfx_ne(i,j,k) = 0.5*(ne(i-1,j,k)+ne(i,j,k))
                 halfx_Se(i,j,k) = 0.5*(Se(i-1,j,k)+Se(i,j,k))
                 halfx_theta(i,j,k) = 0.5*(theta(i-1,j,k)+theta(i,j,k))
                 halfx_m(i,j,k) = 0.5*(m(i-1,j,k)+m(i,j,k))
                    if(halfx_ne(i,j,k).eq.0) then
                       cycle
                    else if(halfx_ne(i,j,k).eq.1) then
                       rx(i,j,k) = 1
                    else                       
                       halfx_KKs(i,j,k) = halfx_ne(i,j,k)**3/((1-halfx_ne(i,j,k))**2) * g0/(5.*Ms**2*nu)
                       halfx_KK(i,j,k) = halfx_KKs(i,j,k) * halfx_Se(i,j,k)**0.5 * &
                            (1.-(1.-halfx_Se(i,j,k)**(1./halfx_m(i,j,k)))**halfx_m(i,j,k))**2
                       halfx_KK_inv(i,j,k) = 1./(halfx_KK(i,j,k))     !!!need to recalculate if KK is a matrix
                       rx(i,j,k) = 1+g0*dt*halfx_KK_inv(i,j,k)*halfx_theta(i,j,k)
                    endif
              Enddo
           Enddo
        Enddo



        Do i=1,Nx
           Do j=2,Ny
              Do k=1,Nz
                 halfy_ne(i,j,k) = 0.5*(ne(i,j-1,k)+ne(i,j,k))
                 halfy_Se(i,j,k) = 0.5*(Se(i,j-1,k)+Se(i,j,k))
                 halfy_theta(i,j,k) = 0.5*(theta(i,j-1,k)+theta(i,j,k))
                 halfy_m(i,j,k) = 0.5*(m(i,j-1,k)+m(i,j,k))
                    if(halfy_ne(i,j,k).eq.0) then
                       cycle
                    else if(halfy_ne(i,j,k).eq.1) then
                       ry(i,j,k) = 1
                    else      
                       halfy_KKs(i,j,k) = halfy_ne(i,j,k)**3/((1-halfy_ne(i,j,k))**2) * g0/(5.*Ms**2*nu)
                       halfy_KK(i,j,k) = halfy_KKs(i,j,k) * halfy_Se(i,j,k)**0.5 * &
                            (1.-(1.-halfy_Se(i,j,k)**(1./halfy_m(i,j,k)))**halfy_m(i,j,k))**2
                       halfy_KK_inv(i,j,k) = 1./(halfy_KK(i,j,k))     !!!need to recalculate if KK is a matrix
                       ry(i,j,k) = 1+g0*dt*halfy_KK_inv(i,j,k)*halfy_theta(i,j,k) 
                    endif
              Enddo
           Enddo
        Enddo




!!! exchange water content at interface

        call mpi_sendrecv(Se(:,:,Nz),ilenz,MPI_Double_precision,inext,910, &
             tmp_Se(:,:,1),ilenz,MPI_Double_precision,iprev,910,nallgrp,stat,error)
        call mpi_sendrecv(Se(:,:,1),ilenz,MPI_Double_precision,iprev,920, &
             tmp_Se(:,:,2),ilenz,MPI_Double_precision,inext,920,nallgrp,stat,error)

        call mpi_sendrecv(theta(:,:,Nz),ilenz,MPI_Double_precision,inext,930, &
             tmp_theta(:,:,1),ilenz,MPI_Double_precision,iprev,930,nallgrp,stat,error)
        call mpi_sendrecv(theta(:,:,1),ilenz,MPI_Double_precision,iprev,940, &
             tmp_theta(:,:,2),ilenz,MPI_Double_precision,inext,940,nallgrp,stat,error)


!!$        CALL MPI_BARRIER(MPI_COMM_WORLD,error)




        Do i=1,Nx
           Do j=1,Ny
              Do k=1,Nz+1
                 if(k.eq.1) then
                    halfz_ne(i,j,1) = 0.5*(tmp_ne(i,j,1)+ne(i,j,1))
                    halfz_Se(i,j,1) = 0.5*(tmp_Se(i,j,1)+Se(i,j,1))
                    halfz_theta(i,j,1) = 0.5*(tmp_theta(i,j,1)+theta(i,j,1))
                    halfz_m(i,j,1) = 0.5*(tmp_m(i,j,1)+m(i,j,1))
                    if(halfz_ne(i,j,1).eq.0) then
                       cycle
                    else if(halfz_ne(i,j,1).eq.1) then
                       rz(i,j,1) = 1
                    else   
                       halfz_KKs(i,j,1) = halfz_ne(i,j,1)**3/((1-halfz_ne(i,j,1))**2) * g0/(5.*Ms**2*nu)
                       halfz_KK(i,j,1) = halfz_KKs(i,j,1) * halfz_Se(i,j,1)**0.5 * &
                            (1.-(1.-halfz_Se(i,j,1)**(1./halfz_m(i,j,1)))**halfz_m(i,j,1))**2
                       halfz_KK_inv(i,j,1) = 1./(halfz_KK(i,j,1))     !!!need to recalculate if KK is a matrix
                       rz(i,j,1) = 1+g0*dt*halfz_KK_inv(i,j,1)*halfz_theta(i,j,1)
                    endif
                 else if(k.eq.Nz+1) then
                    halfz_ne(i,j,Nz+1) = 0.5*(ne(i,j,Nz)+tmp_ne(i,j,2))
                    halfz_Se(i,j,Nz+1) = 0.5*(Se(i,j,Nz)+tmp_Se(i,j,2))
                    halfz_theta(i,j,Nz+1) = 0.5*(theta(i,j,Nz)+tmp_theta(i,j,2))
                    halfz_m(i,j,Nz+1) = 0.5*(m(i,j,Nz)+tmp_m(i,j,2))
                    if(halfz_ne(i,j,Nz+1).eq.0) then
                       cycle
                    else if(halfz_ne(i,j,Nz+1).eq.1) then
                       rz(i,j,Nz+1) = 1
                    else
                       halfz_KKs(i,j,Nz+1) = halfz_ne(i,j,Nz+1)**3/((1-halfz_ne(i,j,Nz+1))**2) * g0/(5.*Ms**2*nu)
                       halfz_KK(i,j,Nz+1) = halfz_KKs(i,j,Nz+1) * halfz_Se(i,j,Nz+1)**0.5 * &
                            (1.-(1.-halfz_Se(i,j,Nz+1)**(1./halfz_m(i,j,Nz+1)))**halfz_m(i,j,Nz+1))**2
                       halfz_KK_inv(i,j,Nz+1) = 1./(halfz_KK(i,j,Nz+1))     !!!need to recalculate if KK is a matrix
                       rz(i,j,Nz+1) = 1+g0*dt*halfz_KK_inv(i,j,Nz+1)*halfz_theta(i,j,Nz+1)
                    endif
                 else
                    halfz_ne(i,j,k) = 0.5*(ne(i,j,k-1)+ne(i,j,k))
                    halfz_Se(i,j,k) = 0.5*(Se(i,j,k-1)+Se(i,j,k))
                    halfz_theta(i,j,k) = 0.5*(theta(i,j,k-1)+theta(i,j,k))
                    halfz_m(i,j,k) = 0.5*(m(i,j,k-1)+m(i,j,k))
                    if(halfz_ne(i,j,k).eq.0) then
                       cycle
                    else if(halfz_ne(i,j,k).eq.1) then
                       rz(i,j,k) = 1
                    else      
                       halfz_KKs(i,j,k) = halfz_ne(i,j,k)**3/((1-halfz_ne(i,j,k))**2) * g0/(5.*Ms**2*nu)
                       halfz_KK(i,j,k) = halfz_KKs(i,j,k) * halfz_Se(i,j,k)**0.5 * &
                            (1.-(1.-halfz_Se(i,j,k)**(1./halfz_m(i,j,k)))**halfz_m(i,j,k))**2
                       halfz_KK_inv(i,j,k) = 1./(halfz_KK(i,j,k))     !!!need to recalculate if KK is a matrix
                       rz(i,j,k) = 1+g0*dt*halfz_KK_inv(i,j,k)*halfz_theta(i,j,k)
                    endif
                 endif
              Enddo
           Enddo
        Enddo








        !! calculate FF
        Do i=1,Nx
           Do j=1,Ny
              Do k=1,Nz
                 if(wall(i,j,k)) then
                    FF(i,j,k) = 0.
                 else 
                    if(h_star(i,j,k).gt.0.0) then
                       dSdh(i,j,k) = 0.0
                    else
                       dSdh(i,j,k) = (1.-Sm(i,j,k))*(-m(i,j,k)) * &
                            (1.+lambda(i,j,k)*Abs(h_star(i,j,k))**tau(i,j,k))**(-m(i,j,k)-1)* &
                            (tau(i,j,k)*(lambda(i,j,k)*Abs(h_star(i,j,k)))**(tau(i,j,k)-1.)*lambda(i,j,k)*h_star(i,j,k)/ &
                            (Abs(h_star(i,j,k))+SMALL))
                    endif
                    FF(i,j,k) = alpha_prime(i,j,k)*theta(i,j,k)/ne(i,j,k) + beta_prime*theta(i,j,k) + ne(i,j,k)*dSdh(i,j,k)
                 endif
              Enddo
           Enddo
        Enddo







        !! calculate QQ
        Do i=2,Nx-1  
           Do j=2,Ny-1
              Do k=1,Nz
                 if(wall(i,j,k)) then
                    q0(i,j,k) = 0.
                 else
                    q0(i,j,k) = 1./rx(i,j,k)*halfx_theta(i,j,k)*u_star(i,j,k)/dx - &
                         1./rx(i+1,j,k)*halfx_theta(i+1,j,k)*u_star(i+1,j,k)/dx + &
                         1./ry(i,j,k)*halfy_theta(i,j,k)*v_star(i,j,k)/dy - &
                         1./ry(i,j+1,k)*halfy_theta(i,j+1,k)*v_star(i,j+1,k)/dy + &
                         1./rz(i,j,k)*halfz_theta(i,j,k)*w_star(i,j,k)/dz - &
                         1./rz(i,j,k+1)*halfz_theta(i,j,k+1)*w_star(i,j,k+1)/dz + &
                         FF(i,j,k)*theta(i,j,k)/dt*h_o(i,j,k)
                 endif
              Enddo
           Enddo
        Enddo




        !! calculate seven-diagonal elements
        Do i=2,Nx-1 
           Do j=2,Ny-1
              Do k=1,Nz
                 if(wall(i,j,k)) then
                    a(i,j,k) = 0.
                    b(i,j,k) = 0.
                    c(i,j,k) = 0.
                    d(i,j,k) = 0.
                    e(i,j,k) = 0.
                    f(i,j,k) = 0.
                    g(i,j,k) = 0.
                 else
                    a(i,j,k) = (halfx_theta(i,j,k)/rx(i,j,k)+halfx_theta(i+1,j,k)/rx(i+1,j,k))*g0*dt/dx/dx + &
                         (halfy_theta(i,j,k)/ry(i,j,k)+halfy_theta(i,j+1,k)/ry(i,j+1,k))*g0*dt/dy/dy 

                    b(i,j,k) = -halfx_theta(i,j,k)/rx(i,j,k)*g0*dt/dx/dx
                    c(i,j,k) = -halfx_theta(i+1,j,k)/rx(i+1,j,k)*g0*dt/dx/dx
                    d(i,j,k) = -halfy_theta(i,j,k)/ry(i,j,k)*g0*dt/dy/dy
                    e(i,j,k) = -halfy_theta(i,j+1,k)/ry(i,j+1,k)*g0*dt/dy/dy
                    f(i,j,k) = -halfz_theta(i,j,k)/rz(i,j,k)*g0*dt/dz/dz
                    g(i,j,k) = -halfz_theta(i,j,k+1)/rz(i,j,k+1)*g0*dt/dz/dz
                 endif
              Enddo
           Enddo
        Enddo





        !! calculate seven-diagonal elements near wall


        Do i=2,Nx-1 
           Do j=2,Ny-1
              Do k=1,Nz
                 if(wall(i-1,j,k)) then
                    a(i,j,k) = a(i,j,k) + b(i,j,k)
                    b(i,j,k) = 0.
                 endif
                 if(wall(i+1,j,k)) then
                    a(i,j,k) = a(i,j,k) + c(i,j,k)
                    c(i,j,k) = 0.
                 endif
                 if(wall(i,j-1,k)) then
                    a(i,j,k) = a(i,j,k) + d(i,j,k)
                    d(i,j,k) = 0.
                 endif
                 if(wall(i,j+1,k)) then
                    a(i,j,k) = a(i,j,k) + e(i,j,k)
                    e(i,j,k) = 0.
                 endif

                 if(k.eq.1) then
                    if(tmp_wall(i,j,1)) then
                       a(i,j,1) = a(i,j,1) + f(i,j,1)
                       f(i,j,1) = 0.
                    endif
                 else if(wall(i,j,k-1)) then
                    a(i,j,k) = a(i,j,k) + f(i,j,k)
                    f(i,j,k) = 0.
                 endif

                 if(k.eq.Nz) then
                    if(tmp_wall(i,j,2)) then
                       a(i,j,Nz) = a(i,j,Nz) + g(i,j,Nz) 
                       g(i,j,Nz) = 0.
                    endif
                 else if(wall(i,j,k+1)) then
                    a(i,j,k) = a(i,j,k) + g(i,j,k) 
                    g(i,j,k) = 0.
                 endif

              Enddo
           Enddo
        Enddo





!!! calculate matrix elements for all side planes 

!!! exchange h_star in interface


        call mpi_sendrecv(h_star(:,:,Nz),ilenz,MPI_double_precision,inext,410, &
             tmp_h_star(:,:,1),ilenz,MPI_double_precision,iprev,410,nallgrp,stat,error)

        call mpi_sendrecv(h_star(:,:,1),ilenz,MPI_double_precision,iprev,420, &
             tmp_h_star(:,:,2),ilenz,MPI_double_precision,inext,420,nallgrp,stat,error)


!!$        CALL MPI_BARRIER(MPI_COMM_WORLD,error)



        !notice
        q = q0

        if(npes.eq.1) then
           f(:,:,2) = 0. 
           q(:,:,2) = q0(:,:,2) + 1./rz(:,:,2)*g0*dt/dz*halfz_theta(:,:,2)/dz*h_star(:,:,1)
           g(:,:,Nz-1) = 0. 
           q(:,:,Nz-1) = q0(:,:,Nz-1) + 1./rz(:,:,Nz)*g0*dt/dz*halfz_theta(:,:,Nz)/dz*h_star(:,:,Nz)
        else
           if(mype.eq.0) then
              f(:,:,2) = 0. 
              q(:,:,2) = q0(:,:,2) + 1./rz(:,:,2)*g0*dt/dz*halfz_theta(:,:,2)/dz*h_star(:,:,1)
              g(:,:,Nz) = 0. 
              q(:,:,Nz) = q0(:,:,Nz) + 1./rz(:,:,Nz+1)*g0*dt/dz*halfz_theta(:,:,Nz+1)/dz*tmp_h_star(:,:,2)
           else if(mype.eq.npes-1) then
              f(:,:,1) = 0. 
              q(:,:,1) = q0(:,:,1) + 1./rz(:,:,1)*g0*dt/dz*halfz_theta(:,:,1)/dz*tmp_h_star(:,:,1)
              g(:,:,Nz-1) = 0. 
              q(:,:,Nz-1) = q0(:,:,Nz-1) + 1./rz(:,:,Nz)*g0*dt/dz*halfz_theta(:,:,Nz)/dz*h_star(:,:,Nz)
           else 
              f(:,:,1) = 0. 
              q(:,:,1) = q0(:,:,1) + 1./rz(:,:,1)*g0*dt/dz*halfz_theta(:,:,1)/dz*tmp_h_star(:,:,1)
              g(:,:,Nz) = 0. 
              q(:,:,Nz) = q0(:,:,Nz) + 1./rz(:,:,Nz+1)*g0*dt/dz*halfz_theta(:,:,Nz+1)/dz*tmp_h_star(:,:,2)
           endif
        endif



        !! calculate elements of matrix 

        Do k=1,Nz
           Do j=1,Ny
              Do i=1,Nx
                 ma((k-1)*Nx*Ny+(j-1)*Nx+i) = a(i,j,k)
                 mb((k-1)*Nx*Ny+(j-1)*Nx+i) = b(i,j,k)
                 mc((k-1)*Nx*Ny+(j-1)*Nx+i) = c(i,j,k)
                 md((k-1)*Nx*Ny+(j-1)*Nx+i) = d(i,j,k)
                 me((k-1)*Nx*Ny+(j-1)*Nx+i) = e(i,j,k)
                 mf((k-1)*Nx*Ny+(j-1)*Nx+i) = f(i,j,k)
                 mg((k-1)*Nx*Ny+(j-1)*Nx+i) = g(i,j,k)

                 QQ((k-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,k)
                 HH_star((k-1)*Nx*Ny+(j-1)*Nx+i) = h_star(i,j,k)

              Enddo
           Enddo
        Enddo







!!! solve h by successive overrelaxation (SOR) 



        nn = 1

        Do while(nn.lt.1e5)

           if(ma(1).eq.0) then
              HH(1) = 0.
           else
              HH(1) = omega*(QQ(1)-mc(1)*HH_star(2)-me(1)*HH_star(1+Nx)-mg(1)*HH_star(1+Nx*Ny))/ &
                   (ma(1)) + (1-omega)*HH_star(1)
           endif
           Do pp=2,Nx
              if(ma(pp).eq.0) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx)- &
                      mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)) + (1-omega)*HH_star(pp)
              endif
           Enddo
           Do pp=1+Nx,Nx*Ny
              if(ma(pp).eq.0) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-md(pp)*HH(pp-Nx)-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)- &
                      me(pp)*HH_star(pp+Nx)-mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)) + &
                      (1-omega)*HH_star(pp)
              endif
           Enddo
           Do pp=1+Nx*Ny,Nt-Nx*Ny
              if(ma(pp).eq.0) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx)- &
                      mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)) + (1-omega)*HH_star(pp)
              endif
           Enddo
           Do pp=Nt-Nx*Ny+1,Nt-Nx
              if(ma(pp).eq.0) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx))/(ma(pp)) + &
                      (1-omega)*HH_star(pp)
              endif
           Enddo
           Do pp=Nt-Nx+1,Nt-1
              if(ma(pp).eq.0) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1))/(ma(pp)) + (1-omega)*HH_star(pp)
              endif
           Enddo
           if(ma(Nt).eq.0) then
              HH(Nt) = 0.
           else
              HH(Nt) = omega*(QQ(Nt)-mf(Nt)*HH(Nt-Nx*Ny)-md(Nt)*HH(Nt-Nx)- &
                   mb(Nt)*HH(Nt-1))/(ma(Nt)) +(1-omega)*HH_star(Nt)
           endif


!!! calculate error
           if(npes.eq.1) then
              HH_err_local = Sum(Abs(HH(Nx*Ny+1:Nt-Nx*Ny)-HH_star(Nx*Ny+1:Nt-Nx*Ny)))
           else 
              if(mype.eq.0) then
                 HH_err_local = Sum(Abs(HH(Nx*Ny+1:Nt)-HH_star(Nx*Ny+1:Nt)))
              else if(mype.eq.npes-1) then
                 HH_err_local = Sum(Abs(HH(1:Nt-Nx*Ny)-HH_star(1:Nt-Nx*Ny)))
              else 
                 HH_err_local = Sum(Abs(HH(1:Nt)-HH_star(1:Nt)))
              endif
           endif

           CALL MPI_ALLREDUCE(HH_err_local,HH_err,1,MPI_double_precision,MPI_SUM,MPI_COMM_WORLD,error)
           HH_err = HH_err/LX/LY/(LZ-2)



!!$           HH_err_local = Sum(Abs(HH-HH_star))/Nt
!!$           CALL MPI_ALLREDUCE(HH_err_local,HH_err,1,MPI_double_precision,MPI_SUM,MPI_COMM_WORLD,error)


           If (HH_err.lt.HH_eps) exit



!!! update interface q

!!$        Do k=1,Nz
!!$           Do j=1,Ny
!!$              Do i=1,Nx
!!$                 h(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
!!$              Enddo
!!$           Enddo
!!$        Enddo

           Do j=1,Ny
              Do i=1,Nx
                 h(i,j,1) = HH((j-1)*Nx+i)
              Enddo
           Enddo

           Do j=1,Ny
              Do i=1,Nx
                 h(i,j,Nz) = HH((Nz-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo



!!$!        if(mype.eq.1 .and. mod(istep,Print_step).eq.0) then
!!$        if(mod(istep,Print_step).eq.0) then
!!$           Print*,'mype=',mype,'h insider SOR---------'
!!$           Print*,nn,h(nxx,nyy,:)
!!$        endif





!!! exchange h during iteration

           call mpi_sendrecv(h(:,:,Nz),ilenz,MPI_double_precision,inext,610, &
                tmp_h(:,:,1),ilenz,MPI_double_precision,iprev,610,nallgrp,stat,error)

           call mpi_sendrecv(h(:,:,1),ilenz,MPI_double_precision,iprev,620, &
                tmp_h(:,:,2),ilenz,MPI_double_precision,inext,620,nallgrp,stat,error)

!!$
!!$           CALL MPI_BARRIER(MPI_COMM_WORLD,error)



           if(mype.ne.0) then
              q(:,:,1) = q0(:,:,1) + 1./rz(:,:,1)*g0*dt/dz*halfz_theta(:,:,1)/dz*tmp_h(:,:,1)
              Do j=1,Ny
                 Do i=1,Nx
                    QQ((j-1)*Nx+i) = q(i,j,1)
                 Enddo
              Enddo
           endif

           if(mype.ne.npes-1) then
              q(:,:,Nz) = q0(:,:,Nz) + 1./rz(:,:,Nz+1)*g0*dt/dz*halfz_theta(:,:,Nz+1)/dz*tmp_h(:,:,2)
              Do j=1,Ny
                 Do i=1,Nx
                    QQ((Nz-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,Nz)
                 Enddo
              Enddo
           endif




!!$   if(mype.eq.0) then
!!$      print*,'inner iteration, nn=',nn,'HH_err=',HH_err
!!$      print*,'HH_star=',maxval(HH_star),minval(HH_star)
!!$      print*,'HH=',maxval(HH),minval(HH)
!!$   endif

           HH_star = HH
           nn = nn+1

        ENDDO





        !! update h

        Do k=1,Nz
           Do j=1,Ny
              Do i=1,Nx
                 h(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo
        Enddo


        !        if(mype.eq.1 .and. mod(istep,Print_step).eq.0) then
!!$        if(mod(istep,Print_step).eq.0) then
!!$           Print*,'mype=',mype,'h after SOR---------'
!!$           Print*,hijk,h(nxy,nxy,:)
!!$        endif


        h_err_local = 0.0
        Do i=1,Nx
           Do j=1,Ny
              if(npes.eq.1) then
                 Do k=2,Nz-1
                    h_err_local = h_err_local + abs(h(i,j,k)-h_star(i,j,k))   
                 Enddo
              else 
                 if(mype.eq.0) then
                    Do k=2,Nz
                       h_err_local = h_err_local + abs(h(i,j,k)-h_star(i,j,k))   
                    Enddo
                 else if(mype.eq.npes-1) then
                    Do k=1,Nz-1
                       h_err_local = h_err_local + abs(h(i,j,k)-h_star(i,j,k))   
                    Enddo
                 else 
                    Do k=1,Nz
                       h_err_local = h_err_local + abs(h(i,j,k)-h_star(i,j,k))   
                    Enddo
                 endif
              endif
           Enddo
        Enddo

        CALL MPI_ALLREDUCE(h_err_local,h_err,1,MPI_double_precision,MPI_SUM,MPI_COMM_WORLD,error)

        h_err = h_err/LX/LY/(LZ-2)

!!$        if(mod(istep,Print_step).eq.0) then
!!$           Print*,'outer iteration'
!!$           Print*,'mype=',mype,'h_err=',h_err
!!$        endif



        if(mype.eq.0) then
           if (hbc_d.eq.1) then
              h(:,:,1) = h_d
           endif
           if (hbc_d.eq.2) then
              h(:,:,1) = h(:,:,2)
           endif
        endif

        if(mype.eq.npes-1) then
           if (hbc_t.eq.1) then
              h(:,:,Nz) = h_t
           endif
           if (hbc_t.eq.2) then
              h(:,:,Nz) = h(:,:,Nz-1)
           endif
        endif




        if(h_err.lt.h_eps) exit




        h_star = h
        hijk = hijk+1

     ENDDO




     !notice
     if(mod(istep,Print_step).eq.0) then
        Print*,'outer iteration'
        Print*,hijk,h_err
        Print*,'mype=',mype,'h after sor---------'
        Print*,hijk,h(nxx,nyy,:)
     endif




     !! exchange h

     call mpi_sendrecv(h(:,:,Nz),ilenz,MPI_double_precision,inext,510, &
          tmp_h(:,:,1),ilenz,MPI_double_precision,iprev,510,nallgrp,stat,error)

     call mpi_sendrecv(h(:,:,1),ilenz,MPI_double_precision,iprev,520, &
          tmp_h(:,:,2),ilenz,MPI_double_precision,inext,520,nallgrp,stat,error)


     CALL MPI_BARRIER(MPI_COMM_WORLD,error)








!!!! calculate u(n+1), consider the effects of wall

     Do i=2,Nx 
        Do j=2,Ny-1
           Do k=1,Nz
              if(wall(i,j,k)) then
                 u(i,j,k) = 0.
              else 
                 u(i,j,k)  = 1./rx(i,j,k) * (u_star(i,j,k) - g0*dt*(h(i,j,k)-h(i-1,j,k))/dx)
              endif
           Enddo
        Enddo
     Enddo




     Do i=2,Nx-1 
        Do j=2,Ny
           Do k=1,Nz
              if(wall(i,j,k)) then
                 v(i,j,k) = 0.
              else 
                 v(i,j,k)  = 1./ry(i,j,k) * (v_star(i,j,k) - g0*dt*(h(i,j,k)-h(i,j-1,k))/dy)
              endif
           Enddo
        Enddo
     Enddo






     Do i=2,Nx-1 
        Do j=2,Ny-1
           Do k=1,Nz+1
              if(k.eq.1) then
                 if(wall(i,j,1)) then
                    w(i,j,1) = 0.
                 else 
                    w(i,j,1)  = 1./rz(i,j,1) * (w_star(i,j,1) - g0*dt*(h(i,j,1)-tmp_h(i,j,1))/dz)
                 endif
              else if(k.eq.Nz+1) then
                 if(tmp_wall(i,j,2)) then
                    w(i,j,Nz+1) = 0.
                 else 
                    w(i,j,Nz+1)  = 1./rz(i,j,Nz+1) * (w_star(i,j,Nz+1) - g0*dt*(tmp_h(i,j,2)-h(i,j,Nz))/dz)
                 endif
              else
                 if(wall(i,j,k)) then
                    w(i,j,k) = 0.
                 else 
                    w(i,j,k)  = 1./rz(i,j,k) * (w_star(i,j,k) - g0*dt*(h(i,j,k)-h(i,j,k-1))/dz)
                 endif
              endif
           Enddo
        Enddo
     Enddo







!!! velocity and pressure boundary condition setting

     if(mype.eq.0) then
        if (vbc_d.eq.1) then
           u(:,:,1) = u_d
           v(:,:,1) = v_d
        endif
        if (vbc_d.eq.2) then
           u(:,:,1) = 0.
           v(:,:,1) = 0.
        endif
        Do i=1,Nx
           Do j=1,Ny
              if(wall(i,j,1)) then
                 w(i,j,1) = 0.
              else
                 if (vbc_d.eq.1) then
                    w(i,j,1) = w_d
                    w(i,j,2) = w_d
                 endif
                 if (vbc_d.eq.2) then
                    w(i,j,1) = w(i,j,2)*(theta(i,j,1)+theta(i,j,2))/2/theta(i,j,1)
                 endif
              endif
           enddo
        enddo

     endif






     if(mype.eq.npes-1) then
        if (vbc_t.eq.1) then
           u(:,:,Nz) = u_t
           v(:,:,Nz) = v_t
        endif
        if (vbc_t.eq.2) then
           u(:,:,Nz) = 0.
           v(:,:,Nz) = 0.
        endif

        Do i=1,Nx
           Do j=1,Ny
              if(wall(i,j,Nz)) then
                 w(i,j,Nz+1) = 0.
              else
                 if (vbc_t.eq.1) then
                    w(i,j,Nz+1) = w_t
                    w(i,j,Nz) = w_t
                 endif
                 if (vbc_t.eq.2) then
                    w(i,j,Nz+1) = w(i,j,Nz)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./theta(i,j,Nz)
                 endif
              endif
           enddo
        enddo

     endif






     !! wall location and pore area
     Do i=1,Nx
        Do j=1,Ny
           Do k=1,Nz
              if(wall(i,j,k)) then
                 u(i,j,k) = 0.
                 u(i+1,j,k) = 0.
                 v(i,j,k) = 0.
                 v(i,j+1,k) = 0.
                 w(i,j,k) = 0.
                 w(i,j,k+1) = 0.
                 h(i,j,k) = 0.
              endif
           Enddo
        Enddo
     Enddo





     u_o = u
     v_o = v
     w_o = w
     h_o = h


     if(mod(istep,Print_step).eq.0) then
        Print*,'mype=',mype,'Se---------','nxx=',nxx,'nyy=',nyy
        Print*,Se(nxx,nyy,:)
        Print*,'mype=',mype,'h---------','nxx=',nxx,'nyy=',nyy
        Print*,h(nxx,nyy,:)
        Print*,'mype=',mype,'w---------','nxx=',nxx,'nyy=',nyy
        Print*,w(nxx,nyy,:)
     endif






!!!*********** output data ***********


     IF(mod(istep,Print_step).eq.0) THEN

        ulocalmax = maxval(u)
        CALL MPI_REDUCE(ulocalmax,umax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        vlocalmax = maxval(v)
        CALL MPI_REDUCE(vlocalmax,vmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        wlocalmax = maxval(w)
        CALL MPI_REDUCE(wlocalmax,wmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)

        ulocalmin = minval(u)
        CALL MPI_REDUCE(ulocalmin,umin,1,MPI_double_precision,MPI_MIN,0,MPI_COMM_WORLD,error)
        vlocalmin = minval(v)
        CALL MPI_REDUCE(vlocalmin,vmin,1,MPI_double_precision,MPI_MIN,0,MPI_COMM_WORLD,error)
        wlocalmin = minval(w)
        CALL MPI_REDUCE(wlocalmin,wmin,1,MPI_double_precision,MPI_MIN,0,MPI_COMM_WORLD,error)


        !! calculate average vel

        frlocalsum = 0.
        azlocalsum = 0.


        Do i=1,Nx
           Do j=1,Ny
              Do k=1,Nz
                 if(.not.wall(i,j,k)) then
                    frlocalsum = frlocalsum + w(i,j,k)*halfz_theta(i,j,k)*dx*dy
                    azlocalsum = azlocalsum + halfz_theta(i,j,k)*dx*dy
                 endif
              Enddo
           Enddo
        Enddo


        CALL MPI_REDUCE(frlocalsum,frsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        CALL MPI_REDUCE(azlocalsum,azsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)



        IF(mype.eq.0) THEN
           print*,'Maximum velocity in x-, y-, z-direction: ',umax,vmax,wmax
           print*,'Minimum velocity in x-, y-, z-direction: ',umin,vmin,wmin
           print*,'Average velocity in vertical direction: ',frsum/azsum
           print*,'Average Darcy velocity in vertical direction: ',frsum/azsum*mean_por
        endif

     endif





!!$!!!*********** write output files ***********

     IF(mod(istep,File_step).eq.0) THEN
        call OutputData !(u,v,w,h,Nz,Nx,Ny,Nstart,MEMO,OUT_DIR,PMname,istep)
     endif


  EndDo



  endtime = MPI_WTIME()


  print*,'End computation'
  print '("Computation time = ",f10.4," seconds.")',endtime-starttime


  CALL MPI_FINALIZE (error)


END PROGRAM USMSSolver



!--------------------------------------------------------------
SUBROUTINE input
  !--------------------------------------------------------------    

  Use GLOBAL
  Implicit None
  include 'mpif.h'

  call MPI_INIT(error)
  if (error.ne.0) stop "ERROR: can't initialize mpi" 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,error)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,error)
  nallgrp=MPI_COMM_WORLD

  if (mype.eq.0) then

     Read(*,*),LX,LY,LZ
     Write(*,*)'## INPUT'
     Write(*,*)'## Spatial resolution LX, LY, LZ = ',LX,LY,LZ
     Read(*,*),dx,dy,dz
     Write(*,*)'## Spatial length dx,dy,dz = ',dx,dy,dz
     Read(*,*),dt
     Write(*,*)'## Time length dt = ',dt

     Read(*,*),hbc_d,hbc_t
     Write(*,*)'## Pressure boundary type on down and top planes = ',hbc_d,hbc_t
     Read(*,*),h_d,h_t
     Write(*,*)'## Pressure on down and top planes h_d,h_t = ',h_d,h_t

     Read(*,*),vbc_d,vbc_t
     Write(*,*)'## Velocity boundary type on down and top planes = ',vbc_d,vbc_t
     Read(*,*),u_d,v_d,w_d
     Write(*,*)'## Velocity on down plane u_d,v_d,w_d = ',u_d,v_d,w_d
     Read(*,*),u_t,v_t,w_t
     Write(*,*)'## Velocity on top plane u_t,v_t,w_t = ',u_t,v_t,w_t

     Read(*,*),rho
     Write(*,*)'## Density rho = ',rho
     Read(*,*),nu
     Write(*,*)'## Kinematic viscosity nu = ',nu

     Read(*,*),g0
     Write(*,*)'## Gravity constant = ',g0

     Read(*,*),gvt(1),gvt(2),gvt(3)
     Write(*,*) '## Gravity acceleration: ',gvt(1),gvt(2),gvt(3)

     Read(*,*),Cm
     Write(*,*)'## Constant Cm = ',Cm
     Read(*,*),m0
     Write(*,*)'## parameter for V-N model m0 = ',m0
     Read(*,*),d10
     Write(*,*)'## particle size distribution d10 = ',d10
     Read(*,*),hcc_s
     Write(*,*)'## parameter for modified lambda = ',hcc_s


     Read(*,*),alpha
     Write(*,*)'## Parameter alpha = ',alpha
     Read(*,*),beta
     Write(*,*)'## Parameter beta = ',beta

     Read(*,*),Ms
     Write(*,*)'## Parameter Ms = ',Ms
     Read(*,*),omega
     Write(*,*)'## Over-relaxation factor = ',omega


     Read(*,*),h_eps
     Write(*,*)'## outer iteration criteria = ',h_eps
     Read(*,*),HH_eps
     Write(*,*)'## inner iteration criteria = ',HH_eps


     Read(*,*), Max_step,Print_step,File_step
     Write(*,*)'## Maximum # of time steps =',Max_step
     Write(*,*)'## # of time steps between print =',Print_step
     Write(*,*)'## # of time steps between file output =',File_step
     Read(*,*), PMname
     Write(*,*)'## Porous medium file: ',Trim(PMname)
     Read(*,*), IN_DIR
     Write(*,*)'## Input directory:  ',Trim(IN_DIR)
     Read(*,*), OUT_DIR
     Write(*,*)'## Output directory: ',Trim(OUT_DIR)
     Read(*,*),mean_por
     Write(*,*)'## mean porosity: ',mean_por
     Read(*,*),t0
     Write(*,*)'## initial time steps = ',t0

  endif

  inext=mod(mype+npes+1,npes)
  iprev=mod(mype+npes-1,npes)


  CALL MPI_BCAST(LX,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(LY,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(LZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(dx,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dy,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dz,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dt,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(hbc_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(hbc_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(h_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(h_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(vbc_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(vbc_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(u_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(v_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(w_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(u_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(v_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(w_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(rho,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(nu,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(g0,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(gvt,3,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(Cm,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(m0,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(d10,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(hcc_s,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(alpha,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(beta,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Ms,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(omega,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(h_eps,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(HH_eps,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(Max_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Print_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(File_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(PMname,64,MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(IN_DIR,64,MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(OUT_DIR,64,MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(mean_por,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(t0,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)


  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  Nstart = mype * LZ / npes + 1
  Nx = LX
  Ny = LY
  Nz = (mype+1) * LZ / npes - mype * LZ / npes

  print *,'Processor =',mype,'Domain start =',Nstart,'Domain size =',Nz


END SUBROUTINE input






!--------------------------------------------------------------
SUBROUTINE SOR(ma,mb,mc,md,me,mf,mg,HH,QQ,omega,HH_star,Nx,Ny,Nz,Nt,HH_eps,WW)
  !--------------------------------------------------------------    

  !  Use GLOBAL
  Implicit None
  Integer::pp,nn,Nx,Ny,Nz,Nt
  Double precision::omega,HH_err,HH_eps
  Double precision::SMALL_NUM
  Double precision,Dimension(Nt)::ma,mb,mc,md,me,mf,mg,HH,QQ,HH_star
  Logical,Dimension(Nt)::WW


  nn = 1
  SMALL_NUM = 1e-15

  Do While(nn.lt.1e6)


     if(WW(1)) then
        HH(1) = 0.
     else
        HH(1) = omega*(QQ(1)-mc(1)*HH_star(2)-me(1)*HH_star(1+(Nx-2))-mg(1)*HH_star(1+(Nx-2)*(Ny-2)))/ &
             (ma(1)+SMALL_NUM) + (1-omega)*HH_star(1)
     endif
     Do pp=2,(Nx-2)
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+(Nx-2))- &
                mg(pp)*HH_star(pp+(Nx-2)*(Ny-2)))/(ma(pp)+SMALL_NUM) + (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=1+(Nx-2),(Nx-2)*(Ny-2)
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-md(pp)*HH(pp-(Nx-2))-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)- &
                me(pp)*HH_star(pp+(Nx-2))-mg(pp)*HH_star(pp+(Nx-2)*(Ny-2)))/(ma(pp)+SMALL_NUM) + &
                (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=1+(Nx-2)*(Ny-2),Nt-(Nx-2)*(Ny-2)
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-(Nx-2)*(Ny-2))-md(pp)*HH(pp-(Nx-2))- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+(Nx-2))- &
                mg(pp)*HH_star(pp+(Nx-2)*(Ny-2)))/(ma(pp)+SMALL_NUM) + (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=Nt-(Nx-2)*(Ny-2)+1,Nt-(Nx-2)
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-(Nx-2)*(Ny-2))-md(pp)*HH(pp-(Nx-2))- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+(Nx-2)))/(ma(pp)+SMALL_NUM) + &
                (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=Nt-(Nx-2)+1,Nt-1
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-(Nx-2)*(Ny-2))-md(pp)*HH(pp-(Nx-2))- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1))/(ma(pp)+SMALL_NUM) + (1-omega)*HH_star(pp)
        endif
     Enddo
     if(WW(Nt)) then
        HH(Nt) = 0.
     else
        HH(Nt) = omega*(QQ(Nt)-mf(Nt)*HH(Nt-(Nx-2)*(Ny-2))-md(Nt)*HH(Nt-(Nx-2))- &
             mb(Nt)*HH(Nt-1))/(ma(Nt)+SMALL_NUM) +(1-omega)*HH_star(Nt)
     endif




     HH_err = Sum(Abs(HH-HH_star))/Nt
!!$
!!$     print*,'inner iteration, nn=',nn,'HH_err=',HH_err
!!$     print*,'HH_star=',maxval(HH_star),minval(HH_star)
!!$     print*,'HH=',maxval(HH),minval(HH)

!!$     do pp=1,Nt
!!$        print*,pp+121,HH(pp)
!!$     enddo
!!$ 
     If (HH_err.lt.HH_eps) exit

     HH_star = HH
     nn = nn + 1

  Enddo
  RETURN
END SUBROUTINE SOR





!--------------------------------------------------------------
SUBROUTINE GS(ma,mb,mc,md,me,mf,mg,HH,QQ,alpha_h,HH_star,Nx,Ny,Nz,Nt,HH_eps,WW)
  !--------------------------------------------------------------    

  !  Use GLOBAL
  Implicit None
  Integer::pp,nn,Nx,Ny,Nz,Nt
  Double precision::alpha_h,HH_err,HH_eps,SMALL_NUM
  Double precision,Dimension(Nt)::HH,QQ,QQ_star,HH_star,ma,mb,mc,md,me,mf,mg
  Logical,Dimension(Nt)::WW

  nn = 1
  SMALL_NUM = 1e-15

  ma = ma/alpha_h


  Do While(nn.lt.1e5)

!!$     print*,'nn=',nn,'HH_err=',HH_err
     !!     print*,'HH_star=',maxval(HH_star),minval(HH_star)
!!$     print*,'HH=',maxval(HH),minval(HH)
!!$     print*,'****'


     Do pp=1,Nt
        QQ_star(pp) = QQ(pp) + (1.0-alpha_h)*ma(pp)*HH_star(pp)
     Enddo


     if(WW(1)) then
        HH(1) = 0.
     else
        HH(1) = (QQ_star(1)-mc(1)*HH_star(2)-me(1)*HH_star(1+Nx)-mg(1)*HH_star(1+Nx*Ny))/ &
             (ma(1)+SMALL_NUM)
     endif
     Do pp=2,Nx
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = (QQ_star(pp)-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx)- &
                mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)+SMALL_NUM) 
        endif
     Enddo
     Do pp=1+Nx,Nx*Ny
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = (QQ_star(pp)-md(pp)*HH(pp-Nx)-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)- &
                me(pp)*HH_star(pp+Nx)-mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)+SMALL_NUM)
        endif
     Enddo
     Do pp=1+Nx*Ny,Nt-Nx*Ny
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = (QQ_star(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx)- &
                mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)+SMALL_NUM) 
        endif
     Enddo
     Do pp=Nt-Nx*Ny+1,Nt-Nx
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = (QQ_star(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx))/(ma(pp)+SMALL_NUM) 
        endif
     Enddo
     Do pp=Nt-Nx+1,Nt-1
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = (QQ_star(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1))/(ma(pp)+SMALL_NUM)
        endif
     Enddo
     if(WW(Nt)) then
        HH(Nt) = 0.
     else
        HH(Nt) = (QQ_star(Nt)-mf(Nt)*HH(Nt-Nx*Ny)-md(Nt)*HH(Nt-Nx)- &
             mb(Nt)*HH(Nt-1))/(ma(Nt)+SMALL_NUM)
     endif


     HH_err = Sum(Abs(HH-HH_star))/Nt

     print*,'inner iteration, nn=',nn,'alpha_h=',alpha_h,'HH_err=',HH_err
!!$     print*,'HH_star=',maxval(HH_star),minval(HH_star)
     print*,'HH=',maxval(HH),minval(HH)
!!$ 

     !!modife over-relaxation factor

!!$     If (HH_err.gt.1000*HH_eps) then
!!$        alpha_h=0.1
!!$     else If (HH_err.gt.100*HH_eps) then
!!$        alpha_h=0.3
!!$     else If (HH_err.gt.10*HH_eps) then
!!$        alpha_h=0.5
!!$     else If (HH_err.gt.5*HH_eps) then
!!$        alpha_h=0.7
!!$     else If (HH_err.gt.2*HH_eps) then
!!$        alpha_h=0.9
!!$     else If (HH_err.gt.HH_eps) then
!!$        alpha_h=1.0
!!$     else
!!$        exit
!!$     endif
!!$
!!$     HH_star = HH
!!$     nn = nn + 1
!!$     


     If (HH_err.lt.HH_eps)  exit
     HH_star = HH
     nn = nn + 1

  Enddo
  RETURN
END SUBROUTINE GS




!--------------------------------------------------------------
Subroutine ReadDensity(den,myblock,LX,LY,sta,MEMO,FileName)
  !--------------------------------------------------------------      
  Implicit None
  Integer:: LX,LY,myblock,MEMO,RecLength,IO,idz,sta
  Double precision,Intent(OUT):: den(LX,LY,myblock)
  Character*(*),Intent(IN):: FileName

  RecLength=MEMO*LX*LY
  Open(30,ACCESS='DIRECT',STATUS='OLD',RECL=RecLength, &
       FILE=FileName,IOSTAT=IO)
  If (IO.Ne.0) Then
     Write(*,*)'Error opening density file, IO = ',IO
     Stop
  End If
  Do idz = 1,myblock
     Read(30,REC=sta+idz-1)den(:,:,idz)
  End Do
  Close(30)
End Subroutine ReadDensity






!--------------------------------------------------------------
Subroutine OutputData
  !--------------------------------------------------------------   

  use GLOBAL
  implicit none
  include 'mpif.h'
!!$  Double precision,Dimension(Nx+1,Ny,Nz)::u
!!$  Double precision,Dimension(Nx,Ny+1,Nz)::v
!!$  Double precision,Dimension(Nx,Ny,Nz+1)::w
!!$  Double precision,Dimension(Nx,Ny,Nz)::h
!!$  Integer:: Nx,Ny,Nz,Nstart,MEMO,RecLength,istep
!!$  Character*(*),Intent(IN):: OUT_DIR,FileName
!!$  Character(64)::Rootname

  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_u_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(u(:,:,:),Nz,Nx+1,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif



  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_v_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(v(:,:,:),Nz,Nx,Ny+1,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_w_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)

        if(npes.eq.1) then
           CALL WriteDEN(w(:,:,:),Nz+1,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
                CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
                CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))
        else
           if(mype.eq.npes-1) then
              CALL WriteDEN(w(:,:,:),Nz+1,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
                   CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
                   CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))
           else 
              CALL WriteDEN(w(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
                   CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
                   CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))
           endif
        endif

!!$  CALL WriteDEN(w(:,:,:),Nz+1,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
!!$       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
!!$       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_h_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(h(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_Se_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(Se(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif





end Subroutine OutputData




!--------------------------------------------------------------
SUBROUTINE WriteDEN(DEN,myblock,LX,LY,Nstart,MEMO,FileName)
  !--------------------------------------------------------------

  IMPLICIT none
  INTEGER:: LX,LY,MEMO,RecLength,IO,idz,myblock,Nstart
  Double precision:: DEN(LX,LY,myblock)
  CHARACTER*(*),INTENT(IN):: FileName

  RecLength=MEMO*LX*LY
  OPEN(20,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=RecLength, &
       FILE=FileName,IOSTAT=IO)
  DO idz=1,myblock
     WRITE(20,REC=idz)DEN(:,:,idz)
  ENDDO
  CLOSE(20)
  RETURN
END SUBROUTINE WriteDEN





!--------------------------------------------------------------
SUBROUTINE WriteMatlab()
  !--------------------------------------------------------------
  USE GLOBAL
  IMPLICIT none
  INTEGER :: idx,time

  OPEN(20,FILE=trim(OUT_DIR)//'/'//'LoadGeo.m')
  WRITE(20,'(A,I4,A,I4,A,I4,A,A,A,A,A)')'por = ReadReal(' ,LX,',', LY, ',' ,LZ,&
       ',''../',trim(IN_DIR),'/',trim(PMname),''',''float64'');'
  CLOSE(20)


  OPEN(20,FILE=trim(OUT_DIR)//'/'//'LoadVel.m')

  WRITE(20,'(A,I8,A)')'DT = ', File_step,';'

  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A)')'u = zeros( ',Max_step/File_step+1,&
       ',', LX+1,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A)')'v = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY+1,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A)')'w = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ+1,');'


  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_u'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'u(', 1,',:,:,:) = ReadReal( ',&
       LX+1,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'     

  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_v'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'v(', 1,',:,:,:) = ReadReal( ', &
       LX,',', LY+1,',', LZ,',''', trim(Rootname),''',''float64'');'

  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_w'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'w(', 1,',:,:,:) = ReadReal( ', &
       LX,',', LY,',', LZ+1,',''', trim(Rootname),''',''float64'');'




  DO idx = 1, Max_step/File_step

     time=idx*File_step
     Rootname = trim(PMname)//'_u_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'u(', idx+1,',:,:,:) = ReadReal( ',&
          LX+1,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'

     Rootname = trim(PMname)//'_v_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'v(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY+1,',', LZ,',''', trim(Rootname),''',''float64'');'

     Rootname = trim(PMname)//'_w_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'w(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ+1,',''', trim(Rootname),''',''float64'');'



  enddo
  close(20)



  OPEN(20,FILE=trim(OUT_DIR)//'/'//'LoadPrs.m')

  WRITE(20,'(A,I8,A)')'DT = ', File_step,';'

  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'h = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'


  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_h'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'h(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'



  DO idx = 1, Max_step/File_step

     time=idx*File_step
     Rootname = trim(PMname)//'_h_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'h(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  enddo

  close(20)






  OPEN(20,FILE=trim(OUT_DIR)//'/'//'LoadSe.m')

  WRITE(20,'(A,I8,A)')'DT = ', File_step,';'

  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'Se = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'

  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_Se'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'Se(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'




  DO idx = 1, Max_step/File_step

     time=idx*File_step
     Rootname = trim(PMname)//'_Se_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'Se(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'


  enddo
  close(20)









  RETURN
END SUBROUTINE WriteMatlab


