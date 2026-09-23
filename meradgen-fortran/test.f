*Here is the test of photonic variable simulation
	 program main
	 
	implicit real*8(a-h,k,m,o-z)
	real vpgen(4),elab
	character gen*7
	include 'output.inc'
	include 'test.inc'
	pi=atan(1d0)*4d0
	data m/.511000d-3/,m2/.261112d-6/
	elab=45.
	itest=1! test v distribution
c	itest=2! test t1 distribution
c	itest=3! test z distribution
c	itest=4! test vectors
	s=2.d0*(elab*m+m2)
	ecm=dsqrt(2.*m*(elab+m))/2
	pcm=dsqrt(ecm**2-m2)
	
	 call merad_init(elab)

	thetacm=90.
	phi=10.
	 pl=1d0

*Input	 
*vpgen - 4-momentum of the virtual photon (vpgen=k1-k2) in CM system


*Output (output.inc):
*vprad=p2-p1 in CM system
*phirad -  4-momentum of the real photon
	
	 vpgen(4)=0.
	 vpgen(1)=-pcm*sin(thetacm*pi/180.)*cos(phi*pi/180.)	 
	 vpgen(2)=-pcm*sin(thetacm*pi/180.)*sin(phi*pi/180.)	 
	 vpgen(3)=pcm*(1.-cos(thetacm*pi/180.))	 
	 if (itest.eq.4)then
	 nev=5	
	 else
	 nev=1e8
	 endif
       t=dble(vpgen(4)**2-vpgen(1)**2-vpgen(2)**2-vpgen(3)**2)				
	 if(itest.ge.2)vgen=(s+t)/4.
	 if(itest.ge.3)then
      t1min=(2d0*m2*t+vgen*(t-vgen-sqrt((t-vgen)**2-4d0*m2*t)))/
     .(2d0*(m2+vgen))
      t1max=m2*t**2/(m2+vgen)/t1min	
	 t1gen=0.5d0*(t1max+t1min)
	 endif
	 do j=1,nbin
	 sigbin(j)=0d0
	 argbin(j)=0d0
	 enddo
	 open(8,file='rnd.dat')
	 read(8,*)iy 	
	 close(8)	

         if(itest.eq.1)gen='v'
         if(itest.eq.2)gen='t1'
         if(itest.eq.3)gen='z'
	 open(12,file='test.dat')
	 if(itest.lt.4)write(12,13)itest,gen,elab,thetacm,
     &    pl,nbin,nev,iy 
13	 format	(
     . 'itest=',i1,/
     . ,a3,'generation',/
     .'rgen is generated probability',/
     .'rcalc is calculated probability',/
     .'Ebeam=',g8.3,'GeV',/
     .'theta=',g8.3,'degrees in CM system',/
     .'P=pb*pt=',g8.3,'beam polarization times target polarization',/
     .'number of bins',i3, /
     .'number of radiative events ',i11,/
     .'initial random number',i3)
	 if(itest.eq.4)write(12,14)itest,elab,thetacm,pl,nev,iy 
14	 format	(
     . 'itest=',i1,/
     . ,'variable reconstruction',/
     .'Ebeam=',g8.3,'GeV',/
     .'theta=',g8.3,'degrees in CM system',/
     .'P=pb*pt=',g8.3,'beam polarization times target polarization',/
     .'number of radiative events ',i2,/
     .'initial random number',i3)
         if(itest.gt.1.and.itest.le.3)write(12,'(a3,g11.4)')'v=',vgen
         if(itest.eq.3)write(12,'(a3,g11.4)')'t1=',t1gen
         if(itest.le.3)write(12,'(a3,4a13)'),'bin',gen,'rgen     ','rcalc     ','rgen/rcalc'

	 DO I=1,nev
232	 continue	
	 call meradgen(pl,vpgen)
	 if(ich.eq.0)goto 232 	
	if(itest.eq.4)then 
      write(12,*),'-------------------------------------------'
      write(12,*),'event=',i
      write(12,*),'test v reconstruction'
      write(12,41),'v=',(vprad(4)+phirad(4)+sqrt(s)/2.)**2
     .-(vprad(1)+phirad(1))**2
     .-(vprad(2)+phirad(2))**2
     .-(vprad(3)+phirad(3)-sqrt(s*(s-4.*m2))/2./sqrt(s))**2-m2
     .,'    reconstructed v from 4-vectors' 	
      write(12,42),'v=',vgen,' generated v ' 	
      write(12,*)
      write(12,*),'test t1 reconstruction'
      write(12,43),'t1=',vprad(4)**2-vprad(1)**2-vprad(2)**2-vprad(3)**2
     .,'   reconstructed t1 from 4-vectors' 	
      write(12,44),'t1=',t1gen,'   generated t1 ' 	
      write(12,*)
      write(12,*),'test z reconstruction'
      write(12,41),'z=',2*((sqrt(s)/2-vprad(4))*phirad(4)
     . +vprad(1)*phirad(1)+vprad(2)*phirad(2)
     .-(sqrt(s*(s-4.*m2))/2./sqrt(s)-vprad(3))*phirad(3))
     .,'    reconstructed z from 4-vectors' 		 
      write(12,42),'z=',zgen,'   generated z ' 	
      write(12,*)
      write(12,45),'m2gamma=',phirad(4)**2-phirad(1)**2-phirad(2)**2-phirad(3)**2
     .,'   real photon mass square' 		 
      write(12,*)
41      format	(a2,g12.6,a35)	
42      format	(a2,g12.6,a17)	
43      format	(a3,g13.6,a35)	
44      format	(a3,g13.6,a17)	
45      format	(a8,g13.6,a26)	
	 goto 233
	endif
c	 print *,itest,ich,t1gen
	 if(itest.eq.1)vargen=vgen
	 if(itest.eq.2)vargen=t1gen
	 if(itest.eq.3)vargen=zgen
c	 print*,i,weight
	 do j=1,nbin
c	 print*,j,vargen,argbin(j) 		
	 if(vargen.le.bin(j))then
	 sigbin(j)=sigbin(j)+1d0/dble(nev)/step
	 argbin(j)=argbin(j)+vargen
	 goto 233
	 endif
	 enddo
233	 continue	
	 enddo
	 if(itest.eq.4)return	
	 do j=1,nbin
	 vargen=argbin(j)/dble(nev)/step/sigbin(j)	
c	 vargen=0.5d0*(bin(j)+bin(j-1))
	 if(itest.eq.1)s=fsir(t,0d0,vargen,0d0,pl,0,2)/dble(weight)	
	 if(itest.eq.2)then
       	call zd(t,vargen,vgen)
	s=fsir(t,vargen,vgen,0d0,pl,0,1)/dble(weight)
	 endif
	 if(itest.eq.3)s=fsir(t,t1gen,vgen,vargen,pl,0,0)/dble(weight)	
c	 print *,'s',s		
	 write(12,'(i3,4g13.4)'),j,vargen,sigbin(j),s,sigbin(j)/s
	 enddo

	 end
