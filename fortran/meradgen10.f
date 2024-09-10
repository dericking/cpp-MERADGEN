  subroutine meradgen(ppl,vpgen, r)
implicit  none
  real*8 y,ppl,u0,u,vmax,vmin
  real*8 fsirv,sig,xsvt,xsbt,dcanc
  real*8 distsiv,distarv,fsir
  real*8 distsit1,distart1
  real*8 distsiz,distarz
  real*8 vvn,vvo,sin,sio
  real*8 tt1n,tt1o,sirad1
  real*8 zmin,zmax,zzn,zzo,det
  real*8 sitot,xs0,xsvr,xsV,xsB,xsF,xsR
  real*8 xsadd,sirand,urand
  real*8 fs,v,t1min,t1max
  real*8 t1z,t1z1,t1z2,stept1,t1p
  real vpgen(4)
  real*8 r(4)
  integer nv,nnv,nn,iv,ikey
  integer nt1,it1,itt,iz,iz1,iz2
  integer nz,i
  parameter(nv=60)
  parameter(nt1=30)
  parameter(nz=60)
include 'const.inc'
include 'tv.inc'
include 'gr.inc'
include 'grid.inc'
include 'output.inc'
include 'test.inc'
c	dimension vpgen(4),vprad(4),phirad(4)
  dimension distsiv(0:nv),distarv(0:nv)
  dimension distsit1(0:4*nt1),distart1(0:4*nt1)
  dimension distsiz(0:nz),distarz(0:nz)
  dimension t1p(5)
  save xs0,distsiv,distarv,distsit1,distart1
 .,distsiz,distarz
 external fsirv,fsirv1
data ikey/0/
 pl=ppl
   t=dble(vpgen(4)**2-vpgen(1)**2-vpgen(2)**2-vpgen(3)**2)
  vmax=(s+t)/2d0
 vmin=2.*Egmin*m*1d-1
if(itest.eq.2)goto 22	
if(itest.eq.3)goto 44	
if(ikey.eq.0)then
if(itest.eq.1)then
step=(vmax-vmin)/dble(nbin)
do i=0,nbin
bin(i)=vmin+dble(i)*step
c	print *,i,argbin(i),vmin,vmax
enddo
c	stop
endif	
   vvn=vmin
   sin=0d0
   nn=0
  distsiv(0)=0d0
  distarv(0)=vvn
  xs0=sig(t,pl,0)
 do iv=1,nv
   vvo=vvn
   sio=sin
 vvn=vmin+(vmax-vmin)*grv(iv)
 sin=fsir(t,0d0,vvn,0d0,pl,nn,2)
  distsiv(iv)=distsiv(iv-1)+(sin+sio)*(vvn-vvo)/2d0
  distarv(iv)=vvn
enddo
if(itest.eq.1)ikey=1
 u0=-s-t
*      Born cross sectiom
  xs0=sig(t,pl,0)
*      Self-Energy of Boson contribution
*      Vertex contribution
  xsvr=sig(t,pl,1)
*      Box contribution
  xsB=xsBt(pl,s,t,u0)+xsBt(pl,s,u0,t)
*      Result of IRD cancellation
  xsF=xs0*dcanc(vmin,s,t,u0)
*      Additional soft Bremsstrahlung
   call simpsx(1d-22,vmin,10000,1d-3,fsirv,xsadd)
  sinonr=xs0+xsvr+xsB+xsf+xsadd
endif
  sirad=distsiv(nv)
  sitot=sirad+sinonr
  if(itest.eq.0)weight=sngl(sitot/xs0)
  if(itest.eq.1)weight=sngl(sirad)
  sigmar = sirad

  sirand=r(1)*sitot!urand(iy)*sitot
c	print *,sirand,sitot
c      print *,vmax, vmin, sirand, sitot

  if(sirand.le.sinonr)then
  vgen=0d0
  t1gen=t
  zgen=0d0
ich=0
  do i=1,4
vprad(i)=vpgen(i)
phirad(i)=0.
  enddo
return
  endif
ich=1
  do iv=1,nv
c       print*,iv	
  if(distsiv(iv).gt.sirand-sinonr)then
  vgen=distarv(iv-1)+(distarv(iv)-distarv(iv-1))*
 .(sirand-sinonr-distsiv(iv-1))/
 .(distsiv(iv)-distsiv(iv-1))
if(itest.eq.1)return	
  goto 22
  endif
  enddo

22    continue
if(ikey.eq.0)then
  u=4d0*m2+vgen-s-t
  t1min=(2d0*m2*t+vgen*(t-vgen-sqrt((t-vgen)**2-4d0*m2*t)))/
 .(2d0*(m2+vgen))
  t1max=m2*t**2/(m2+vgen)/t1min	
c	print *,ikey
c	pause
c	ich=1
if(itest.eq.2)ikey=1
step=(t1max-t1min)/dble(nbin)
do i=0,nbin
bin(i)=t1min+dble(i)*step
c	print *,i,argbin(i),t1min,t1max
enddo
i=0
c	stop
  t1z=(-s*t*(t+u)+2d0*m2*vgen**2)/((s-vgen)**2-4d0*m2*s)
  t1z1=(-u*t*(t+s)+2d0*m2*vgen**2)/((u-vgen)**2-4d0*m2*u)
  t1z2=(-s*vgen*(t+s)
 .-2d0*m2*(2d0*t*u+(u-2d0*(s+vgen))*vgen))/((u-vgen)**2-4d0*m2*u)
   t1p(1)=t1min
   t1p(2)=t1z
   t1p(3)=min(t1z1,t1z2)
   t1p(4)=max(t1z1,t1z2)
   t1p(5)=t1max
   tt1n=t1min
   sin=0d0
  distsit1(0)=0d0
  distart1(0)=tt1n
 do i=1,4
 do it1=1,nt1
   tt1o=tt1n
   sio=sin
 tt1n=t1p(i)+(t1p(i+1)-t1p(i))*grt1(it1)
call zd(t,tt1n,vgen)
 sin=fsir(t,tt1n,vgen,0d0,pl,nn,1)
  distsit1((i-1)*nt1+it1)=distsit1((i-1)*nt1+it1-1)
 .+(sin+sio)*(tt1n-tt1o)/2d0
  distart1((i-1)*nt1+it1)=tt1n
enddo
enddo
sirad=distsit1(4*nt1)
  if(itest.eq.2)weight=sngl(sirad)
endif
i=i+1	
  sirand=r(2)*sirad!urand(iy)*sirad
c	print *,i,sirad,sirand
ich=1
  do it1=1,4*nt1
  if(distsit1(it1).gt.sirand)then
  t1gen=distart1(it1-1)+(distart1(it1)-distart1(it1-1))*
 .(sirand-distsit1(it1-1))/
 .(distsit1(it1)-distsit1(it1-1))
c	print *,'2'
if(itest.eq.2)return	
  goto 44
  endif
  enddo
44    continue
if(ikey.eq.0)then
ich=1	 
call zd(t,t1gen,vgen)
det=(bz**2-az*cz)	
zmax=(-bz+dsqrt(det))/az
c	zmin=(-bz-dsqrt(det))/az
zmin=cz/az/zmax
if(itest.eq.3)ikey=1
step=(zmax-zmin)/dble(nbin)
do i=0,nbin
bin(i)=zmin+dble(i)*step
c	print *,i,argbin(i),t1min,t1max
enddo
i=0
c	stop
   zzn=zmin
   sin=0d0
   nn=0
  distsiz(0)=0d0
  distarz(0)=zzn
 do iz=1,nz
   zzo=zzn
   sio=sin
 zzn=zmin+(zmax-zmin)*grz(iz)
 sin=fsir(t,t1gen,vgen,zzn,pl,nn,0)
  distsiz(iz)=distsiz(iz-1)+(sin+sio)*(zzn-zzo)/2d0
  distarz(iz)=zzn
enddo
sirad=distsiz(nz)
  if(itest.eq.3)weight=sngl(sirad)
endif
  sirand=r(3)*sirad!urand(iy)*sirad
  do iz=1,nz
  if(distsiz(iz).gt.sirand)then
  zgen=distarz(iz-1)+(distarz(iz)-distarz(iz-1))*
 .(sirand-distsiz(iz-1))/
 .(distsiz(iz)-distsiz(iz-1))
if(itest.eq.3)return	
  goto 55
  endif
  enddo
55    continue

if(itest.ne.0.and.itest.ne.4)return
call vectrec(vpgen)
return
  
end

  subroutine vectrec(vpgen)
implicit  none
include 'const.inc'
include 'tv.inc'
include 'output.inc'
include 'test.inc'
  real*8 al1,al2,al3,al4,al5,al6,al7,al8
  real*8 sls,sl1,sl3,sl8,phi,sp,cp
real*8 urand
real vpgen(4)
  t=dble(vpgen(4)**2-vpgen(1)**2-vpgen(2)**2-vpgen(3)**2)
phi=atan2(vpgen(2),vpgen(1))
  al1=(s-vgen)**2-4d0*m2*s
  al2=s+2d0*t-vgen-4d0*m2
  al3=-s*t*(s+t-vgen-4d0*m2)-m2*vgen**2
  al4=s*(s-vgen-4d0*m2)-(s+vgen)*zgen
  al5=vgen*zgen*(s-vgen-zgen)-m2*(vgen+zgen)**2
  al6=s*(vgen-zgen)-vgen*(vgen+zgen)
  al7=(s+2d0*t1gen-zgen-4d0*m2)*al1-al2*al4
  al8=16d0*al3*al5-al7**2
  sls=sqrt(als)
  sl1=sqrt(al1)
  sl3=sqrt(al3)
  sl8=sqrt(al8)
  if(urand(iy).gt.0.5)sl8=-sl8
  sp=dsin(phi)
  cp=dcos(phi)
  vprad(1)=-sngl((sls*sl1*sl8*sp+(4d0*al3*al4-s*al2*al7)*cp)
 ./(4d0*al1*sls*sl3))	 
  vprad(2)=sngl((sls*sl1*sl8*cp-(4d0*al3*al4-s*al2*al7)*sp)
 ./(4d0*al1*sls*sl3))	 
  vprad(3)=sngl((als*al1-s*(al7+al2*al4))
 ./(2d0*sqrt(s)*al1*sls))
  vprad(4)=-sngl(zgen/2d0/sqrt(s))	 
  phirad(1)=sngl((sls*sl1*sl8*sp+(4d0*al3*al6-s*al2*al7)*cp)
 ./(4d0*al1*sls*sl3))	 
  phirad(2)=sngl((-sls*sl1*sl8*cp+(4d0*al3*al6-s*al2*al7)*sp)
 ./(4d0*al1*sls*sl3))	 
  phirad(3)=sngl(sqrt(s)*(al7+al2*al6)
 ./(2d0*al1*sls))
  phirad(4)=sngl((vgen+zgen)/2d0/sqrt(s))
return
end

  subroutine merad_init(elab)
implicit none
real elab
include 'const.inc'
pi=atan(1d0)*4d0
  data alfa/.729735d-2/,m/.511000d-3/,m2/.261112d-6/
 .barn/.389379d6/
En=elab
s=2.d0*(En*m+m2)
als=s*(s-4.*m2)
 coeb=4d0*pi*alfa**2*barn*(s-2d0*m2)/als
 coer=alfa**3*barn*(s-2d0*m2)/als/pi/4d0
c Minimum of the photonic enegry which is able to detect by
c calorimeter.
 Egmin=En*1d-2
c Minimum inelasticity
c	  vmin=2.*Egmin*m	     ! ??
 open(12,file='rnd.dat')
c	 read(12,*)iy 	
 close(12)	
call grid_init
 end

  subroutine grid_init
implicit none
  integer i,nv,nt1,nz
  parameter(nv=60)
  parameter(nt1=30)
  parameter(nz=60)
include 'grid.inc'
  do i=1,30
  grv(i)=dble(i-1)/dble(29)/4d0
  enddo
  do i=31,45
  grv(i)=0.25d0+dble(i-30)/dble(45-30)/4d0
  enddo
  do i=46,nv
  grv(i)=0.5d0+dble(i-45)/dble(nv-45)/2d0
  enddo
  do i=1,7
grt1(i)=0.1d0*dble(i)**2/49d0/2d0
grt1(31-i)=1d0+grt1(1)-grt1(i)		
  enddo
  do i=8,15
 grt1(i)=(0.1d0+0.9d0*(dble(i-7))/8d0)/2d0
grt1(31-i)=1d0+grt1(1)-grt1(i)		
  enddo
  do i=1,30
 grz(i)=0.5d0*dble(i**2)/30d0**2
  enddo
  do i=1,30
 grz(61-i)=1d0-0.49d0*dble(i**2)/30d0**2
  enddo
   end

 double precision function sig(t,pl,i)
 implicit none
 real*8 t,pl,u,ss,u1,u2,u3,pl1,pl2,pl3
 real*8 dsvt,dsvu,du1,du2,du3,dp1,dp2,dp3,vacpol,l1f
 integer i
include 'const.inc'
c	    t=-s*y
  u=4*m2-s-t
 ss=s-2d0*m2
 u1=(ss**2+u**2)/2d0+2d0*m2*(s+2d0*t-3d0*m2)
 u2=(ss**2+t**2)/2d0+2d0*m2*(s+2d0*u-3d0*m2)
 u3=(s-2d0*m2)*(s-6*m2)
c      pl1=((16.0*m2**2-12.0*m2*s+2.0*s**2+s*t)*s*t)/(2.0*als)
c      pl2=((16.0*m2**2-8.0*m2*s+s**2-s*t)*(4.0*m2-s-t)*s)/(2.0*als)
c      pl3=(-((2.0*m2-s)**4-4.0*(2.0*m2-s)**2*m2**2+4.0*(2.0*m2-s)**2*
c     . m2*t-8.0*(2.0*m2-t)*m2**2*t))/als
 pl1=-t*(-t*s**2/2d0/als-ss)
 pl2=-t*(-t*s**2/2d0/als-2d0*m2)+als/2d0-ss*s
 pl3=-(ss**4+4d0*m2*(ss**2*t+m2*(-ss**2+2d0*t**2-4d0*m2*t)))
 .	  /als
if(i.eq.0)then
sig=coeb*(u1/t**2+u2/u**2+u3/u/t
 .	    +pl*(pl1/t**2+pl2/u**2+pl3/u/t))
return
else
dsvt=vacpol(-t)+L1f(t,s,m2)
dsvu=vacpol(-u)+L1f(u,s,m2)
du1=dsvt
du2=dsvu
du3=dsvt+dsvu
dp1=dsvt
dp2=dsvu
dp3=dsvt+dsvu
sig=alfa/pi*coeb*(du1*u1/t**2+du2*u2/u**2+du3*u3/u/t/2d0
 .	    +pl*(dp1*pl1/t**2+dp2*pl2/u**2+dp3*pl3/u/t/2d0))
endif
end

   double precision function vacpol(t)
implicit none
real*8 t,am2,a2,sqlmi,allmi,suml
real*8 aaa,bbb,ccc,sumh
integer i
include 'const.inc'
  dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
  data am2/.26110d-6,.111637d-1,3.18301d0/
  suml=0.
  do 10 i=1,3
 a2=2.*am2(i)
 sqlmi=dsqrt(t*t+2.*a2*t)
allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
  if(t.lt.1.d0)then
aaa = -1.345d-9
bbb = -2.302d-3
ccc = 4.091
  elseif(t.lt.64d0)then
aaa = -1.512d-3
bbb =  -2.822d-3
ccc = 1.218
  else
aaa = -1.1344d-3
bbb = -3.0680d-3
ccc = 9.9992d-1
  endif
  sumh = -(aaa+bbb*dlog(1d0+ccc*t)) *2d0*pi/alfa
  vacpol=suml+sumh
  end


  double precision function L1f(xt,xs,xm2)
  implicit none
  real*8 xt,xs,xm2
include 'const.inc'
  L1f= - 2d0*dlog(abs(xt)/xs)*( dlog(abs(xt)/xm2)-1d0 )
 &	   + dlog(abs(xt)/xm2) + (dlog(abs(xt)/xm2))**2
 &	   + 4d0*( pi**2/12d0-1d0 )
  return
  end

c      double precision function xsVt(pl,xs,xt,xu)
c      implicit none
c      real*8 xs,xt,xu,l1f,pl
c	 include 'const.inc'
c      xsVt=2d0*alfa**3/xt**2* barn*
c     & 	 ((1d0+pl)*xu**2/xs-(1d0-pl)*xs**2/xu)*
c     & 	 L1f(xt,xs,m2)
c      return
c      end

  double precision function xsBt(pl,xs,xt,xu)
  implicit none
   real*8 xs,xt,xu,l1f,dgg1,dgg2,pl
include 'const.inc'
  xsBt=2d0*alfa**3/xt**2* barn*
 &	((1d0+pl)*xu**2/xs*dgg1(xs,xt,xu)
 &	 -(1d0-pl)*xs**2/xu*dgg2(xs,xt,xu))
  return
  end

  double precision function dgg1(xs,xt,xu)
  implicit none
  real*8 xs,xt,xu,ls,lx,dgg
include 'const.inc'
LS=dlog(xs/abs(xt))
LX=dlog(xu/xt)
  dgg = LS**2*(xS**2+xU**2)/2./xt-LS*xU-(LX**2+PI**2)*xU**2/xt
  dgg1=2d0*dlog(xs/abs(xu))*dlog(dsqrt(abs(xu/xs)))-xt/xu**2*dgg
  return
  end

  double precision function dgg2(xs,xt,xu)
  implicit none
  real*8 xs,xt,xu,ls,lx,dgg
include 'const.inc'
LS=dlog(xs/abs(xt))
LX=dlog(xu/xt)
  dgg = LS**2*xS**2/xt+LX*xS-(LX**2+PI**2)*(xS**2+xU**2)/2./xt
  dgg2=2d0*dlog(xs/abs(xu))*dlog(dsqrt(abs(xu/xs)))-xt/xs**2*dgg
  return
  end

  double precision function dcanc(xvmin,xs,xt,xu)
  implicit none
  real*8 xvmin,xs,xt,xu,lm,lr,del1s,eps,del1h,del1h1,fspen
include 'const.inc'
include 'xx.inc'
xxs=xs
xxt=xt
xxu=xu
   lm=dlog(-xt/m2)
   lr=dlog(-xu/xs)
  del1s=-2.5d0*lm**2+(3d0-2d0*lr)*lm-lr**2/2d0
 .      -(lm-1)*dlog(xs*(xs+xt)/xt**2)-pi**2/3d0+1d0

  del1h=-lm**2/2d0+(dlog(xt**2*(xs+xt)**2*(xs-xvmin)
 ./xs/(xvmin-xt)/xvmin/(xs+xt-xvmin)**2)+1d0)*lm
 .-dlog(-xvmin/xt)**2/2d0-dlog(1d0-xvmin/xt)**2
 .+dlog((xs+xt)/(xs+xt-xvmin))*dlog((xs+xt)*(xs+xt-xvmin)/xt**2)
 .+dlog((xvmin-xs)/xt)*dlog(1d0-xvmin/xs)+dlog(-xvmin/xt) 
 .+fspen((xs-xvmin)/xs)-fspen((xt-xvmin)/xt)
 .+2d0*(fspen(xvmin/xs)-fspen(xvmin/xt)-fspen(xvmin/(xt+xs)))
 .-pi**2/6d0
  dcanc=alfa/pi*
 &	 (  4d0*dlog(xvmin/dsqrt(m2*xs))*(dlog(xt*xu/m2/xs)-1d0)
 &	   + del1s + del1h  )
  return
  end

  double precision function fspens(x)
c
c    spence function
c
  implicit real*8(a-h,o-z)
  f=0.d0
  a=1.d0
  an=0.d0
  tch=1.d-16
1   an=an+1.d0
  a=a*x
  b=a/an**2
  f=f+b
  if(b-tch)2,2,1
2   fspens=f
  return
  end

  double precision function fspen(x)
  implicit real*8(a-h,o-z)
  data f1/1.644934d0/
  if(x)8,1,1
1   if(x-.5d0)2,2,3
2 fspen=fspens(x)
  return
3 if(x-1d0)4,4,5
4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
  return
5 if(x-2d0)6,6,7
6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
  return
7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
  return
8 if(x+1d0)10,9,9
9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
  return
10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
  return
  end

  subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
  implicit real*8(a-h,o-z)
  dimension f(7),p(5)
  h=dsign(h1,b1-a1)
  s=dsign(1.d0,h)
  a=a1
  b=b1
  ai=0.d0
  aih=0.d0
  aiabs=0.d0
  p(2)=4.d0
  p(4)=4.d0
  p(3)=2.d0
  p(5)=1.d0
  if(b-a) 1,2,1
1 reps=dabs(reps1)
  aeps=dabs(aeps1)
  do 3 k=1,7
3   f(k)=10.d16
  x=a
  c=0.d0
  f(1)=funct(x)/3.
4 x0=x
  if((x0+4.*h-b)*s) 5,5,6
6 h=(b-x0)/4.
  if(h) 7,2,7
7 do 8 k=2,7
8   f(k)=10.d16
  c=1.d0
5 di2=f(1)
  di3=dabs(f(1))
  do 9 k=2,5
  x=x+h
  if((x-b)*s) 23,24,24
24 x=b
23 if(f(k)-10.d16) 10,11,10
11 f(k)=funct(x)/3.
10 di2=di2+p(k)*f(k)
9 di3=di3+p(k)*abs(f(k))
  di1=(f(1)+4.*f(3)+f(5))*2.*h
  di2=di2*h
  di3=di3*h
  if(reps) 12,13,12
13 if(aeps) 12,14,12
12 eps=dabs((aiabs+di3)*reps)
  if(eps-aeps) 15,16,16
15 eps=aeps
16 delta=dabs(di2-di1)
  if(delta-eps) 20,21,21
20 if(delta-eps/8.) 17,14,14
17 h=2.*h
  f(1)=f(5)
  f(2)=f(6)
  f(3)=f(7)
  do 19 k=4,7
19  f(k)=10.d16
  go to 18
14 f(1)=f(5)
  f(3)=f(6)
  f(5)=f(7)
  f(2)=10.d16
  f(4)=10.d16
  f(6)=10.d16
  f(7)=10.d16
18 di1=di2+(di2-di1)/15.
  ai=ai+di1
  aih=aih+di2
  aiabs=aiabs+di3
  go to 22
21 h=h/2.
  f(7)=f(5)
  f(6)=f(4)
  f(5)=f(3)
  f(3)=f(2)
  f(2)=10.d16
  f(4)=10.d16
  x=x0
  c=0.d0
  go to 5
22 if(c) 2,4,2
2 return
  end

  subroutine simpsx(a,b,np,ep,func,res)
  implicit real*8 (a-h,o-z)
  external func
  step=(b-a)/np
  call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
  end

  subroutine simptx(a,b,np,ep,func,res)
  implicit real*8 (a-h,o-z)
  external func
  step=(b-a)/np
  call simpt(a,b,step,ep,1d-18,func,ra,res,r2,r3)
  end

  subroutine simpux(a,b,np,ep,func,res)
  implicit real*8 (a-h,o-z)
  external func
  step=(b-a)/np
  call simpu(a,b,step,ep,1d-18,func,ra,res,r2,r3)
  end

  subroutine simpt(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
  implicit real*8(a-h,o-z)
  dimension f(7),p(5)
  h=dsign(h1,b1-a1)
  s=dsign(1.d0,h)
  a=a1
  b=b1
  ai=0.d0
  aih=0.d0
  aiabs=0.d0
  p(2)=4.d0
  p(4)=4.d0
  p(3)=2.d0
  p(5)=1.d0
  if(b-a) 1,2,1
1 reps=dabs(reps1)
  aeps=dabs(aeps1)
  do 3 k=1,7
3   f(k)=10.d16
  x=a
  c=0.d0
  f(1)=funct(x)/3.
4 x0=x
  if((x0+4.*h-b)*s) 5,5,6
6 h=(b-x0)/4.
  if(h) 7,2,7
7 do 8 k=2,7
8   f(k)=10.d16
  c=1.d0
5 di2=f(1)
  di3=dabs(f(1))
  do 9 k=2,5
  x=x+h
  if((x-b)*s) 23,24,24
24 x=b
23 if(f(k)-10.d16) 10,11,10
11 f(k)=funct(x)/3.
10 di2=di2+p(k)*f(k)
9 di3=di3+p(k)*abs(f(k))
  di1=(f(1)+4.*f(3)+f(5))*2.*h
  di2=di2*h
  di3=di3*h
  if(reps) 12,13,12
13 if(aeps) 12,14,12
12 eps=dabs((aiabs+di3)*reps)
  if(eps-aeps) 15,16,16
15 eps=aeps
16 delta=dabs(di2-di1)
  if(delta-eps) 20,21,21
20 if(delta-eps/8.) 17,14,14
17 h=2.*h
  f(1)=f(5)
  f(2)=f(6)
  f(3)=f(7)
  do 19 k=4,7
19  f(k)=10.d16
  go to 18
14 f(1)=f(5)
  f(3)=f(6)
  f(5)=f(7)
  f(2)=10.d16
  f(4)=10.d16
  f(6)=10.d16
  f(7)=10.d16
18 di1=di2+(di2-di1)/15.
  ai=ai+di1
  aih=aih+di2
  aiabs=aiabs+di3
  go to 22
21 h=h/2.
  f(7)=f(5)
  f(6)=f(4)
  f(5)=f(3)
  f(3)=f(2)
  f(2)=10.d16
  f(4)=10.d16
  x=x0
  c=0.d0
  go to 5
22 if(c) 2,4,2
2 return
  end

  subroutine simpu(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
  implicit real*8(a-h,o-z)
  dimension f(7),p(5)
  h=dsign(h1,b1-a1)
  s=dsign(1.d0,h)
  a=a1
  b=b1
  ai=0.d0
  aih=0.d0
  aiabs=0.d0
  p(2)=4.d0
  p(4)=4.d0
  p(3)=2.d0
  p(5)=1.d0
  if(b-a) 1,2,1
1 reps=dabs(reps1)
  aeps=dabs(aeps1)
  do 3 k=1,7
3   f(k)=10.d16
  x=a
  c=0.d0
  f(1)=funct(x)/3.
4 x0=x
  if((x0+4.*h-b)*s) 5,5,6
6 h=(b-x0)/4.
  if(h) 7,2,7
7 do 8 k=2,7
8   f(k)=10.d16
  c=1.d0
5 di2=f(1)
  di3=dabs(f(1))
  do 9 k=2,5
  x=x+h
  if((x-b)*s) 23,24,24
24 x=b
23 if(f(k)-10.d16) 10,11,10
11 f(k)=funct(x)/3.
10 di2=di2+p(k)*f(k)
9 di3=di3+p(k)*abs(f(k))
  di1=(f(1)+4.*f(3)+f(5))*2.*h
  di2=di2*h
  di3=di3*h
  if(reps) 12,13,12
13 if(aeps) 12,14,12
12 eps=dabs((aiabs+di3)*reps)
  if(eps-aeps) 15,16,16
15 eps=aeps
16 delta=dabs(di2-di1)
  if(delta-eps) 20,21,21
20 if(delta-eps/8.) 17,14,14
17 h=2.*h
  f(1)=f(5)
  f(2)=f(6)
  f(3)=f(7)
  do 19 k=4,7
19  f(k)=10.d16
  go to 18
14 f(1)=f(5)
  f(3)=f(6)
  f(5)=f(7)
  f(2)=10.d16
  f(4)=10.d16
  f(6)=10.d16
  f(7)=10.d16
18 di1=di2+(di2-di1)/15.
  ai=ai+di1
  aih=aih+di2
  aiabs=aiabs+di3
  go to 22
21 h=h/2.
  f(7)=f(5)
  f(6)=f(4)
  f(5)=f(3)
  f(3)=f(2)
  f(2)=10.d16
  f(4)=10.d16
  x=x0
  c=0.d0
  go to 5
22 if(c) 2,4,2
2 return
  end

 double precision function fsirv(v)
 implicit none
 real*8 v,fsir
 integer nn
include 'tv.inc'
 fsirv=fsir(t,0d0,v,0d0,pl,nn,-1)
 end
 double precision function fsirv1(v)
 implicit none
 real*8 v,fsir
 integer nn
include 'tv.inc'
 fsirv1=fsir(t,0d0,v,0d0,pl,nn,2)
 end

subroutine zd(t,t1,v)
  implicit none
real*8 t,t1,v,u
include 'const.inc'
include 'gr.inc'
 u=v-s-t+4d0*m2
az=(v-t)**2-4d0*m2*t
bz=-(v*(2d0*m2*(t+t1)+t1*(v-t)))
 .+s*(-t**2+t1*v+t*(t1+v))
   cz=(s*(t-t1)+t1*v)**2
 .-4d0*m2*(s*(t-t1)**2+t1*v**2)
az1=az
bz1=-(t*(s+t-4d0*m2)*(t-t1))+
 . (t*(2*t-t1)-2d0*m2*(t+t1)+s*(t+t1))*v - t*v**2
cz1=((s+t)*(t-t1)-t*v)**2+4d0*m2*(-((s+t)*(t-t1)**2)+
 .	   (t-t1)*(t+t1)*v-t1*v**2)
az2=az
bz2=(4d0*m2-s-t)*t*(4d0*m2-t1) +
 .	(6d0*m2*t-s*t-2d0*m2*t1+s*t1-t*t1)*v+(-4d0*m2 + s)*v**2
cz2=u*(-4d0*m2*(s+t1-4d0*m2)**2+
 .	   (16d0*m2**2+t1**2-4d0*m2*(s+2d0*t1))*u) -
 .	2d0*(2d0*m2-t1)*(4d0*m2-s-t1)*u*v+(s+t1-4d0*m2)**2*v**2
return
end

c$nodebug
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
  FUNCTION URAND(IY)
C *   This is a standard pseudo-random generator      *
C *   that work on IBM-370 and IBM-PC. We don't       *
C *   know does it work on SUN? 		      *
C *						      *
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
  REAL*8  URAND,S
  INTEGER*4 IY
  INTEGER*4 A,C,MASK
  PARAMETER (A  = 843314861)
  PARAMETER (C  = 453816693)
  PARAMETER (S  = 4.6566128E-10)
  REAL *8 ra
  call random_seed
  call random_number(ra)

  IY=IAND(A*IY+C,Z'7FFFFFFF')
  URAND = ra
c      URAND=FLOAT(IY)*S
c	 if(urand.eq.1.) urand=0.4
  END

  
