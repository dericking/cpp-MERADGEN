*Here is the sample of generation
	 program main
	 
	implicit real*8(a-h,k,m,o-z)
	real vpgen(4),elab
	include 'output.inc'
	include 'test.inc'
	pi=atan(1d0)*4d0
	data m/.511000d-3/,m2/.261112d-6/
	elab=45.
	itest=0
	s=2.d0*(elab*m+m2)
	ecm=dsqrt(2.*m*(elab+m))/2
	pcm=dsqrt(ecm**2-m2)
	
	 call merad_init(elab)

	thetacm=90.
	phi=10.
	 pl=-1d0

*Input	 
*vpgen - 4-momentum of the virtual photon (vpgen=k1-k2) in CM system


*Output (output.inc):
*vprad=p2-p1 in CM system
*phirad -  4-momentum of the real photon
	
	 vpgen(4)=0.
	 vpgen(1)=-pcm*sin(thetacm*pi/180.)*cos(phi*pi/180.)	 
	 vpgen(2)=-pcm*sin(thetacm*pi/180.)*sin(phi*pi/180.)	 
	 vpgen(3)=pcm*(1.-cos(thetacm*pi/180.))	 
	 n=100

	 DO I=1,n
	 call meradgen(pl,vpgen)
*ich=0 - non-radiative channel	 
*ich=1 - radiative channel
*weight is the ratio of 1-loop corrected cross section to the born one	 
c	print *,i,y,pl,v,t1,z,weight,ich
	print *,vprad
	print *,phirad
c	 stop 
	 enddo

	 end
