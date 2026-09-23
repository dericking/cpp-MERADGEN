      program driver_parity
      implicit real*8(a-h,k,m,o-z)
      real vpgen(4),elab
      real*8 r1,r2,r3,r4
      integer out_unit,iev,ibits,nrad,maxrad,maxcall
      character*1024 qfile,ofile
      character*16 mode
      character*32 carg
      include 'output.inc'
      include 'test.inc'

      if (command_argument_count() .lt. 2) then
         write(*,*) 'usage: driver_parity quads.txt output.txt',
     &              ' [es14|full|hex] [max_rad] [max_calls]',
     &              ' [elab thetacm phi pl]'
         stop 1
      endif
      call get_command_argument(1, qfile)
      call get_command_argument(2, ofile)
      mode='es14'
      if (command_argument_count() .ge. 3) then
         call get_command_argument(3, mode)
      endif
      mode=adjustl(mode)
      maxrad=0
      maxcall=0
      nrad=0
      if (command_argument_count() .ge. 4) then
         call get_command_argument(4, carg)
         read(carg,*) maxrad
      endif
      if (command_argument_count() .ge. 5) then
         call get_command_argument(5, carg)
         read(carg,*) maxcall
      endif

      pi=atan(1d0)*4d0
      data m/.511000d-3/,m2/.261112d-6/
      elab=45.
      thetacm=90.
      phi=10.
      pl=-1d0
      if (command_argument_count() .ge. 6) then
         call get_command_argument(6, carg)
         read(carg,*) elab
      endif
      if (command_argument_count() .ge. 7) then
         call get_command_argument(7, carg)
         read(carg,*) thetacm
      endif
      if (command_argument_count() .ge. 8) then
         call get_command_argument(8, carg)
         read(carg,*) phi
      endif
      if (command_argument_count() .ge. 9) then
         call get_command_argument(9, carg)
         read(carg,*) pl
      endif
      itest=0
      s=2.d0*(elab*m+m2)
      ecm=dsqrt(2.*m*(elab+m))/2
      pcm=dsqrt(ecm**2-m2)

      call merad_init(elab)

      vpgen(4)=0.
      vpgen(1)=-pcm*sin(thetacm*pi/180.)*cos(phi*pi/180.)
      vpgen(2)=-pcm*sin(thetacm*pi/180.)*sin(phi*pi/180.)
      vpgen(3)=pcm*(1.-cos(thetacm*pi/180.))

      open(11,file=qfile,status='old')
      out_unit=20
      open(out_unit,file=ofile,status='unknown')

      write(out_unit,'(A)') 'HEADER START'
      write(out_unit,'(A,ES24.16)') 'elab=',dble(elab)
      write(out_unit,'(A,ES24.16)') 'thetacm=',thetacm
      write(out_unit,'(A,ES24.16)') 'phi=',phi
      write(out_unit,'(A,ES24.16)') 'pl=',pl
      write(out_unit,'(A,ES24.16)') 'm=',m
      write(out_unit,'(A,ES24.16)') 'm2=',m2
      write(out_unit,'(A,ES24.16)') 'ecm=',ecm
      write(out_unit,'(A,ES24.16)') 'pcm=',pcm
      write(out_unit,'(A,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)')
     &  'vpgen=',vpgen(1),vpgen(2),vpgen(3),vpgen(4)
      write(out_unit,'(A,A)') 'mode=',trim(mode)
      write(out_unit,'(A,I0)') 'max_rad=',maxrad
      write(out_unit,'(A,I0)') 'max_calls=',maxcall
      if (maxrad.gt.0) then
         write(out_unit,'(A)') 'dump_policy=radiative_only'
      else
         write(out_unit,'(A)') 'dump_policy=all'
      endif
      write(out_unit,'(A)') 'HEADER END'

      iev=0
 10   read(11,*,end=20) r1,r2,r3,r4
      iev=iev+1
      call meradgen(pl,vpgen,r1,r2,r3,r4)
      if (ich.eq.1) nrad=nrad+1
      if (maxrad.eq.0 .or. ich.eq.1) then
      write(out_unit,'(A,I0,A)') 'EVENT ',iev,' START'
      write(out_unit,'(A,F0.10,1X,F0.10,1X,F0.10,1X,F0.10,A)')
     &  '  RANDOM{',r1,r2,r3,r4,'}'
      if (trim(mode) .eq. 'hex') then
         write(out_unit,'(A)',advance='no') '  VPRAD{'
         ibits=transfer(vprad(1),ibits)
         write(out_unit,'(A,Z8.8)',advance='no') '0x',ibits
         ibits=transfer(vprad(2),ibits)
         write(out_unit,'(A,Z8.8)',advance='no') ' 0x',ibits
         ibits=transfer(vprad(3),ibits)
         write(out_unit,'(A,Z8.8)',advance='no') ' 0x',ibits
         ibits=transfer(vprad(4),ibits)
         write(out_unit,'(A,Z8.8,A)') ' 0x',ibits,'}'
         write(out_unit,'(A)',advance='no') '  PHIRAD{'
         ibits=transfer(phirad(1),ibits)
         write(out_unit,'(A,Z8.8)',advance='no') '0x',ibits
         ibits=transfer(phirad(2),ibits)
         write(out_unit,'(A,Z8.8)',advance='no') ' 0x',ibits
         ibits=transfer(phirad(3),ibits)
         write(out_unit,'(A,Z8.8)',advance='no') ' 0x',ibits
         ibits=transfer(phirad(4),ibits)
         write(out_unit,'(A,Z8.8,A)') ' 0x',ibits,'}'
      elseif (trim(mode) .eq. 'full') then
         write(out_unit,'(A,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,A)')
     &     '  VPRAD{',vprad(1),vprad(2),vprad(3),vprad(4),'}'
         write(out_unit,'(A,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,A)')
     &     '  PHIRAD{',phirad(1),phirad(2),phirad(3),phirad(4),'}'
      else
         write(out_unit,'(A,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,A)')
     &     '  VPRAD{',vprad(1),vprad(2),vprad(3),vprad(4),'}'
         write(out_unit,'(A,ES14.6,1X,ES14.6,1X,ES14.6,1X,ES14.6,A)')
     &     '  PHIRAD{',phirad(1),phirad(2),phirad(3),phirad(4),'}'
      endif
      if (trim(mode) .eq. 'full' .or. trim(mode) .eq. 'hex') then
         write(out_unit,'(A,I0,1X,ES24.16,1X,ES24.16,1X,ES24.16,A)')
     &     '  KIN{',ich,vgen,t1gen,zgen,'}'
      endif
      endif
      if (maxrad.gt.0 .and. nrad.ge.maxrad) goto 20
      if (maxcall.gt.0 .and. iev.ge.maxcall) goto 20
      goto 10
 20   continue
      write(out_unit,'(A)') 'FOOTER START'
      write(out_unit,'(A,I0)') 'calls=',iev
      write(out_unit,'(A,I0)') 'radiative=',nrad
      write(out_unit,'(A)') 'FOOTER END'
      write(0,'(A,I0,A,I0)') 'calls=',iev,' radiative=',nrad
      write(0,'(A,ES16.8,1X,ES16.8,1X,ES16.8,1X,ES16.8)')
     &  'kinematics elab,thetacm,phi,pl=',
     &  dble(elab),thetacm,phi,pl
      close(11)
      close(out_unit)
      end
