*Here is the sample of generation
*input- vpgen (virtual photon four momentum), polarization, lab energy
*output- generated v, z, t1 and real photon four momentum phirad
      subroutine run(pol,enlab,virtphoton, virtphotonradplus
     .,realphotonplus,vvarplus,zvarplus,t1varplus,sigmabornplus
     .,sigmaradplus,sigmanoradplus,virtphotonradminus
     .,realphotonminus,vvarminus,zvarminus,t1varminus
     .,sigmabornminus,sigmaradminus,sigmanoradminus)
	implicit real*8(a-h,k,m,o-z)
	real vpgen(4),elab
	integer it
	real*8 virtphoton(4), realphotonplus(4), virtphotonradplus(4)
	real*8 realphotonminus(4), virtphotonradminus(4), r(4)
	real*8 vvarplus, zvarplus, t1varplus, enlab, pol
	real*8 sigmabornplus, sigmaradplus,sigmanoradplus
	real*8 sigmabornminus, sigmaradminus, sigmanoradminus
	real*8 vvarminus, zvarminus, t1varminus
	include 'output.inc'
	include 'test.inc'
	pi=atan(1d0)*4d0
	data m/.511000d-3/,m2/.261112d-6/
	elab=enlab
	itest=0
	s=2.d0*(elab*m+m2)
	pl = 1.
	do it=1,3 !in c arrays start from 0 so the time component is the first one instead of the last one
		vpgen(it)=virtphoton(it+1)
	enddo
	vpgen(4) = virtphoton(1)

*   Pre generation of random numbers to use them in both p=1 and p=-1
      	call random_seed()
		call random_number(r)
c	print *,vpgen
	 call merad_init(elab)
*Output (output.inc):
*vprad=p2-p1 in CM system
*phirad -  4-momentum of the real photon
	 call meradgen(pl,vpgen, r)
*ich=0 - non-radiative channel	 
*ich=1 - radiative channel
*weight is the ratio of 1-loop corrected cross section to the born one	 
c	print *,i,y,pl,v,t1,z,weight,ich
      do it=1,3 
        virtphotonradplus(it+1) = vprad(it)
		realphotonplus(it+1)=phirad(it)
	  enddo
	  virtphotonradplus(1)=vprad(4)
	  realphotonplus(1)=phirad(4)
	  vvarplus = vgen
	  zvarplus = zgen
	  t1varplus = t1gen
	  sigmabornplus = (sigmar+sinonr)/weight
	  sigmaradplus = sigmar
	  sigmanoradplus = sinonr
c	  print *, elab, '',vvarplus, '', sigmaradplus, '',sigmanoradplus
c     .,'',weight

	  call meradgen(-pl,vpgen, r)
*ich=0 - non-radiative channel	 
*ich=1 - radiative channel
*weight is the ratio of 1-loop corrected cross section to the born one	 
c	print *,i,y,pl,v,t1,z,weight,ich
      do it=1,3 
        virtphotonradminus(it+1) = vprad(it)
		realphotonminus(it+1)=phirad(it)
	  enddo
	  virtphotonradminus(1)=vprad(4)
	  realphotonminus(1)=phirad(4)
	  vvarminus = vgen
	  zvarminus = zgen
	  t1varminus = t1gen
	  sigmabornminus = (sigmar+sinonr)/weight
	  sigmaradminus = sigmar
	  sigmanoradminus = sinonr
c      print *,vvarminus,'',weight,'',sigmaradminus,'',sigmanoradminus
c	print *,phirad
c	 stop 
	end
