Program srlm_quadpack_ozaki
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! **********************************************************************
! This program aims to calculate the Observables 
! for SRLM coupled to phononic bath from  
! the components of Green Function in complex plane
! based on the idea proposed by T. Ozaki namely,
! continued fraction representation of the Fermi-Dirac function
! taking care of convergence factor analytically!
! Using the ODEPACK library to solve the flow equations.
!######################################################################
!**************the calculation of effective hybridization and
!                Conductance has been included 
!######################################################################
! *************using the intlib.f90 to the do the integration for
!              the calculation of occupation number!
!              in particular avint
! *********************************************************************
!---------------------------------------------------------------------------------
! note depending on the value of gama you may play with the number of ploes!
!---------------------------------------------------------------------------------
! compared to the previous version: comparison to the zero temperature calculation
! using the self-energy in matsubara space!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
implicit none
!------------------------------------------------------------------
!------------------------------------------------------------------
!  Parameter declaration:
!------------------------------------------------------------------
common w0,lambda,gama,pi,e0,SHart,SHartout,w1,w2,Yt0FRG

DOUBLE PRECISION::w0,lambda,gama,pi,zero,g1realana,g1ana,g1num,g1FRG,ep,e0,SHart,JAC,SHartout
INTEGER,parameter::NEQ1=1,Npoles=4000,ITOL=1,nrowpd=2,e0inttotal=200 !,NEQ=6000
INTEGER:: Flag, ISTATE, IOPT, LRW, LIW, MF,MU, ML,k,klinear,kcounter&
          ,ilambda,e0int,NEQ,ineq,polecount,pcounter,ixout
DOUBLE PRECISION,DIMENSION(NEQ1)::Y1stno,YD1stno,Y1st,YD1st,Yana,YDana&
                                 ,Yt0ana,YDt0ana,Rt0frg,RDt0frg
DOUBLE PRECISION::x, xout, RTOL, ATOL
DOUBLE PRECISION,allocatable,dimension(:)::RWORK
INTEGER,allocatable,dimension(:)::IWORK
DOUBLE PRECISION,dimension(nrowpd,NEQ1)::pd1
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(2*Npoles)::Ynum,YDnum,YFRG,YDFRG,Yt0FRG,YDt0FRG
!DOUBLE PRECISION,allocatable,DIMENSION(:)::Ynum,YDnum,YFRG,YDFRG,Yt0FRG,YDt0FRG
DOUBLE PRECISION,DIMENSION(Npoles)::w1
DOUBLE PRECISION,DIMENSION(Npoles)::w2
DOUBLE PRECISION,dimension(nrowpd,NEQ1)::pd
!complex*16,dimension(NEQ/2)::ana
complex*16,allocatable,dimension(:)::ana
!------------------------------------------------------------------
complex*16::dp,dm,dtp,dtm,ic,mu0
integer*8::sgn,sgnn
integer::openstatus,NEQw1
double precision,dimension(60000)::SigmaH
!------------------------------------------------------------------
! For Ozaki:
!------------------------------------------------------------------
double precision,dimension(Npoles)::zp,Rp
complex*16,dimension(Npoles)::Gpper,Gnper,Gpnum,Gnnum,GpFRG,GnFRG,G0pDot,G0nDot
double precision,dimension(2*Npoles-1)::E
double precision,dimension(2*Npoles)::Dmatrix
double precision,dimension(max(1,2*(2*Npoles)-2))::WORK
double precision,dimension(2*Npoles,2*Npoles)::Bmatrix,Z
integer*8::INFO,i,ii,j,ki
CHARACTER*1::JOBZ 
double precision::beta,rhoper,rhonumozaki,rhoFRGozaki,rho0Dot&
                 ,ap,am,width,ktilda,dw1,slope,Kper,KFRGozaki&
                 ,geanaana,geanaozaki,genumozaki,gefrgozaki,geananum,get0ana,get0FRG
double precision, PARAMETER::Rlarge=10.0d+00**(10.00d+00)
double precision,dimension(e0inttotal)::ndper,ndnumozaki,ndFRGozaki&
                                       ,ndanarout,nd0Dot,ndananum,e0array&
                                       ,ndt0ana,ndt0FRG
!-----------------------------------------------------------------------
integer*4,parameter::NTAB=20000
double precision,dimension(2*NTAB)::XTAB,rananum
complex*16,dimension(2*NTAB)::sana

double precision,dimension(2*Npoles)::ycheck
CHARACTER*6::readpoles 
double precision::ratio1,ratio2,ratio3,e0eff
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=1,file="ozaki.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
open(unit=40,file="ozaki4000.txt",status="old",action="read",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
 open(unit=2,file="odepackozakicl0.5e0m1gw0scan.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!open(unit=3,file="odepackozakigeflscanw010.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!open(unit=4,file="odepackozakispec.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!------------------------------------------------------------------
!  Parameter of the model:
!------------------------------------------------------------------
!w0=0.0025d+00
zero=0.0d+00
pi=atan2(zero,-gama)
ic=cmplx(0.0d+0,1.0d+0, kind = 16)
kcounter=0
JAC=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
! calculating the poles and Residues:
!------------------------------------------------------------------
!------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
Bmatrix=0.0d+00
Dmatrix=0.0d+00

Do i=1,2*(Npoles)
if(i+1<=2*(Npoles))then
Bmatrix(i,i+1)=1.0d+00/(2.0d+00*sqrt((2.0d+00*i-1)*(2.0d+00*i+1)))
endif
Bmatrix(i+1,i)=Bmatrix(i,i+1)
if(i<2.0d+00*(Npoles))then
E(i)=1.0d+00/(2.0d+00*sqrt((2.0d+00*i-1)*(2.0d+00*i+1)))
endif
end do
!-----------------------------------------------------------------------
JOBZ = 'V'
!-----------------------------------------------------------------------
!call DSTEV( JOBZ, 2*(NEQ/2), Dmatrix, E, Z, 2*(NEQ/2), WORK, INFO )
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
ii=0
Do i=1,2*Npoles
if(Dmatrix(i)>0)then
ii=ii+1
zp(ii)=1.0d+00/Dmatrix(i)
Rp(ii)=(Z(1,i))**2/(4.0d+00*(Dmatrix(i)**2))
endif

end do

!------------------------------------------------------------------
Do ilambda=1,21
!------------------------------------------------------------------
print*,'ilambda=',ilambda

!w0=0.0025d+00

!lambda=0.035d+00*w0
if(ilambda==1)then
lambda=0.0d+00
else
!lambda=w0*10.0d+00**(-6.0d+00+(ilambda-1)*3.0d+00/19.0d+00)!(0.0d+00+(ilambda-1)*(0.01d+00-0.0d+00)/20.0d+00)*w0
lambda=w0*(ilambda)*(0.002d+00)/20.0d+00
endif

!gama=w0/100.0d+00

!gama=w0*(10d+00**(-4.0d+00+7.0d+00*(ilambda-1)*0.02d+00))

!gama=w0*10.0!(10d+00**(-3.0d+00+5.0d+00*(ilambda-1)/20.0d+00))
gama=0.0001d+00
!w0=20.0d+00*gama
w0=1.0d+00*(ilambda)*gama

!lambda=(ilambda-1)*1.0d+00*w0/(20.0d+00)!(ilambda*sqrt(2.0d+00)/5.0d+00)*w0
lambda=0.5d+00*w0

ep=(lambda**2)/w0

print*,'ep/gama',2.0d+00*ep/(pi*gama)

beta=1.0d+00/(0.0001d+00*gama)
!------------------------------------------------------------------
!------------------------------------------------------------------
!Method flag Help:
!------------------------------------------------------------------
! 10  Nonstiff (Adams) method, no Jacobian used.
! 21  Stiff (BDF) method, user-supplied full Jacobian.
! 22  Stiff method, internally generated full Jacobian.
! 24  Stiff method, user-supplied banded Jacobian.
! 25  Stiff method, internally generated banded Jacobian
!------------------------------------------------------------------
!************** Linearization:*************************************
!------------------------------------------------------------------
MF=10
!------------------------------------------------------------------
! Dimension Declaration of WORK:
!------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ1
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ1+NEQ1**2 
LIW=20+NEQ1
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ1+(2*ML+MU)*NEQ1
LIW=20+NEQ1
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
! Linearilized (no feed-back)
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0
IWORK(5:10)=0
IWORK(6)=100000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)
ATOL=w0*10.0d+00**(-13.0d+00)
ISTATE=1
!------------------------------------------------------------------
x=(1000.0d+00)*w0
xout=0.0d+00
ml=1
mu=1

Y1stno(1:NEQ1)=0.0d+00

!CALL Fno (NEQ1, x, Y1stno, YD1stno)
!CALL JAC (NEQ1, x, Y1stno, ml, mu, pd, nrowpd)

!CALL DLSODE (Fno, NEQ1, Y1stno, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"LIN-Num",ISTATE,xout/w0
!------------------------------------------------------------------
!------------------------------------------------------------------
! Linearilized 
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0
IWORK(5:10)=0
IWORK(6)=100000
!-----------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)
ATOL=w0*10.0d+00**(-13.0d+00)
ISTATE=1
x=(1000.0d+00)*w0
xout=0.0d+00
ml=1
mu=1
Y1st(1:NEQ1)=0.0d+00
!-----------------------------------------------------------------------

!CALL F (NEQ1, x, Y1st, YD1st)
!CALL JAC (NEQ1, x, Y1st, ml, mu, pd, nrowpd)

!CALL DLSODE (F, NEQ1, Y1st, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)


print*,"LIN    ",ISTATE,xout/w0
!-----------------------------------------------------------------------
deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! To check the convergence of Ozaki:
!-----------------------------------------------------------------------
Do ineq=40,40
!-----------------------------------------------------------------------
NEQ=(2*Npoles)  !ineq*100
!polecount=(ineq-1)*100+ineq
!-----------------------------------------------------------------------
! Storing Poles and Residues:
!-----------------------------------------------------------------------
print*,NEQ
Do i=1,(NEQ/2)!+polecount
!write(2,*)i,zp(i),Rp(i)
! if(i==40)then
! read(1,*)readpoles,pcounter
! elseif(i>polecount)then
! read(1,*)j,zp(j),Rp(j)
!write(3,*)j,zp(j),Rp(j)
! endif
if(ilambda==1)then
read(40,*)j,zp(j),Rp(j)
endif
end do
!-----------------------------------------------------------------------
!allocate(Ynum(NEQ),YDnum(NEQ),YFRG(NEQ),YDFRG(NEQ),Yt0FRG(NEQ),YDt0FRG(NEQ),ana(NEQ/2))
allocate(ana(NEQ/2))
!------------------------------------------------------------------
! Dot Energy:
!------------------------------------------------------------------
Do e0int=51,51
!------------------------------------------------------------------
slope=0.1d+00!0.5d+00!0.5d+00!(ilambda)*0.02d+00
!e0=ep+300.0d+00*0.1d+00*gama*(((e0int-51.0)*exp((abs(51.0- e0int)-50.0)/23.0d+00)/50.0d+00)) 
!e0=ep+slope*0.1d+00*gama*(((e0int-51)*exp((abs(51-e0int)-50.0d+00)/23.0d+00)/50.0d+00))
e0=ep-1.0d+00*gama
e0array(e0int)=e0
!e0=(e0int-51)*ep+ep*10.0d+00**(-6.0d+00*(e0int-51))
!e0=w0*2.0d+00
print*,'e0',e0int,e0/gama,ep/gama
!------------------------------------------------------------------
ndananum(e0int)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
! Restricted self-consistent Hartree Term:
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
SigmaH(1)=ep
Do k=1,60000
SigmaH(k+1)=-ep*(1.0d+00-(2.0d+00/pi)*atan((e0+SigmaH(k))/gama))
if(abs((SigmaH(k+1)-SigmaH(k))/ep)<10.0d+00**(-40.0d+00))exit
end do
if(abs((SigmaH(k+1)-SigmaH(k))/ep)>10.0d+00**(-40.0d+00))then
print*,"Self-Consistent Hartree not completed"
endif
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
SHart=-ep!SigmaH(k)! -ep*(1.0d+00-(2.0d+00*e0/(pi*gama)))/(1.0d+00-(2.0d+00*ep/(pi*gama))) !
SHartout=SHart
! if(e0int==51)then
! SHart=-ep*(1.0d+00-(2.0d+00/pi)*atan(e0/gama))
! else
! SHart=(-ep*(1.0d+00-(2.0d+00/pi)*atan(ep/gama)))+(e0-ep)*(2.0d+00*ep/(pi*gama))/(1.0d+00+(ep/gama))
! endif
!SHartout=-ep*(1.0d+00-(2.0d+00/pi)*atan(e0/gama))
SHart=-ep!0.0d+00
print*,"SH",2.0d+00*ep/(pi*gama),SHart/ep,-(1.0d+00-(2.0d+00/pi)*atan((e0+SHart)/gama))&
           ,-(1.0d+00-(2.0d+00*e0/(pi*gama)))/(1.0d+00-(2.0d+00*ep/(pi*gama)))

if(e0int==52)then           
!write(2,*)ilambda,(e0-ep)/gama,SHart/ep&
!          ,(SHart/ep)+(1.0d+00-(2.0d+00*e0/(pi*gama)))/(1.0d+00-(2.0d+00*ep/(pi*gama)))&
!          ,(SHart/ep)+(1.0d+00-(2.0d+00/pi)*atan((e0+SHart)/gama))
endif

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Dot Occupancy in Perturbation Theory at zero Temperature
!        using the self-energy at matsubara space ;D
!-----------------------------------------------------------------------
MF=10
!-----------------------------------------------------------------------
! Dimension Declaration of WORK:
!-----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ1
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ1+NEQ1**2 
LIW=20+NEQ1
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ1+(2*ML+MU)*NEQ1
LIW=20+NEQ1
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Do ixout=1,1!00

Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=100000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-15.0d+00)!0.0d+00*w0*(10.0d+00**(-10.0d+00-(4.0d+00*(ixout-1)/99.0d+00)))   !
ATOL=w0*10.0d+00**(-15.0d+00)!10.0d+00*w0*(10.0d+00**(-10.0d+00-(4.0d+00*(ixout-1)/99.0d+00)))  !

ISTATE=1
!------------------------------------------------------------------
!x=-w0*10.0d+00**(11.0d+00)!(10000.0d+00)*w0
x=w0*10.0d+00**(8.0d+00)!w0*((10.0d+00)**(4.0d+00+(ixout-1)*7.0d+00/99.0d+00))!
xout=0.0d+00!w0*((10.0d+00)**(-10.0d+00))!0.0d+00!

ml=1
mu=1

Yt0ana(1:NEQ1)=0.0d+00

!CALL Ft0ana (NEQ1, x, Yt0ana, YDt0ana)
!CALL JAC (NEQ1, x, Yana, ml, mu, pd, nrowpd)

!CALL DLSODE (Ft0ana, NEQ1, Yt0ana, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

if(ISTATE==2)then
print*,"T0-Ana",ISTATE,ATOL/w0,'integration was successfull'
else
print*,"T0-Ana",ISTATE,ATOL/w0,'NOTE:integration was not successfull'
endif

!write(3,*)ixout,xout/w0,0.5d+00-Yt0ana(1)

end do

deallocate(RWORK,IWORK)

!-----------------------------------------------------------------------
ndt0ana(e0int)=0.5d+00+Yt0ana(1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Dot Occupancy in Perturbation Theory from performing the Analytic
!                continuation analytically ;D
!-----------------------------------------------------------------------
MF=10
!------------------------------------------------------------------
! Dimension Declaration of WORK:
!------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ1
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ1+NEQ1**2 
LIW=20+NEQ1
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ1+(2*ML+MU)*NEQ1
LIW=20+NEQ1
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Do ixout=1,1!00

Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=100000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)!0.0d+00*w0*(10.0d+00**(-10.0d+00-(4.0d+00*(ixout-1)/99.0d+00)))   !
ATOL=w0*10.0d+00**(-13.0d+00)!10.0d+00*w0*(10.0d+00**(-10.0d+00-(4.0d+00*(ixout-1)/99.0d+00)))   !

ISTATE=1
!------------------------------------------------------------------
!x=-w0*10.0d+00**(11.0d+00)!(10000.0d+00)*w0
xout=w0*10.0d+00**(11.0d+00)!w0*((10.0d+00)**(4.0d+00+(ixout-1)*7.0d+00/99.0d+00))!
x=-xout

ml=1
mu=1

Yana(1:NEQ1)=0.0d+00

!CALL Fana (NEQ1, x, Yana, YDana)
!!CALL JAC (NEQ1, x, Yana, ml, mu, pd, nrowpd)

!CALL DLSODE (Fana, NEQ1, Yana, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

if(ISTATE==2)then
print*,"Ana-Ana",ISTATE,ATOL/w0,'integration was successfull'
else
print*,"Ana-Ana",ISTATE,ATOL/w0,'NOTE:integration was not successfull'
endif

!write(2,*)ixout,RTOL/w0,Yana(1)

end do

deallocate(RWORK,IWORK)

!------------------------------------------------------------------
ndanarout(e0int)=Yana(1)
!------------------------------------------------------------------
print*,"Checking the Hartree-Fock-Term",ndanarout(e0int),e0/ep 
!------------------------------------------------------------------
!------------------------------------------------------------------
! Incoming Energy Grid:
!------------------------------------------------------------------
Do k=1,(NEQ/2)
!Do k=1,NEQ

! if(k<(NEQ/2)+1)then
! w2(k)=zp(k)/beta

w1(k)=zp(k)/beta!zp(k+Npoles-(NEQ/2))/beta

if(w1(k)<(0.001d+00)*w0 .and. kcounter==0)then
klinear=k
kcounter=kcounter+1
endif

w2(k)=w1(k)!(exp(-(k-1)*0.01d+00))*10000.0d+00*w0*0.90d+00

! else
! w2(k)=-(zp(8001-k)/beta)
! endif
!------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------
! Analytical results for perturbation theory in matsubara space: (+self-consistent Hartree term)
!-----------------------------------------------------------------------------------------------------------------------------

dp=-1.0d+00/(e0+SHart+w0-ic*(w2(k)+gama))
dm=-1.0d+00/(e0+SHart-w0-ic*(w2(k)+gama))
dtp=-1.0d+00/(e0+SHart+w0-ic*(w2(k)-gama))
dtm=-1.0d+00/(e0+SHart-w0-ic*(w2(k)-gama))        

ana(k)=((ic*(lambda**2)/(2.0d+00*pi))*(0.5d+00*(dtp-dtm-dp+dm)*log(gama**2+(e0+SHart)**2)&
                             -0.5d+00*(dtp-dtm-dp+dm)*log(w0**2+w2(k)**2)&
                             +ic*(dp-dtp)*atan2(w0,w2(k))+ic*(dtm-dm)*atan2(-w0,w2(k))&
                             -ic*(dp+dm)*pi*sgnn(w0)+ic*(dtp-dtm)*atan2(e0+SHart,-gama)&
                             -ic*(dp-dm)*atan2(e0+SHart,gama)&
                             -ic*pi*(dtp-dtm)*sgnn(e0+SHart)))&
                             +e0+SHartout
                             !-2.0d+00*0.0d+00*((lambda**2)/w0)*(0.5d+00-(1.0d+00/pi)*atan2(e0,gama))  

!Write(3,*)NEQ/2,k,w2(k)/w0,real(ana(k)),aimag(ana(k))                     
!---------------------------------------------------------------------- 
end do  !k
!------------------------------------------------------------------
! Flow equations:
!---------------------------------------------------------------
MF=10
!------------------------------------------------------------------
! Dimension Declaration of WORK:
!------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ+NEQ**2 
LIW=20+NEQ
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ+(2*ML+MU)*NEQ
LIW=20+NEQ
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! No feed-back (Perturbation Theory):
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0
IWORK(5:10)=0
IWORK(6)=100000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)
ATOL=w0*10.0d+00**(-13.0d+00)
ISTATE=1
x=(1000.0d+00)*w0!w1(1)!
xout=0.0d+00!w1(NEQ/2)
ml=1
mu=1
Ynum(1:NEQ/2)=0.0d+00                !Imaginary Part
Ynum((NEQ/2)+1:NEQ)=e0+SHart         !Real Part self-consistent
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL Fnum (NEQ, x, Ynum, YDnum)
!CALL JAC (NEQ, x, Ynum, ml, mu, pd, nrowpd)
CALL DLSODE (Fnum, NEQ, Ynum, x, xout, ITOL, RTOL, ATOL, Flag, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,ISTATE,xout/w0,RTOL/gama

!------------------------------------------------------------------
!------------------------------------------------------------------
! FRG ozaki:
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0
IWORK(5:10)=0
IWORK(6)=100000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-13.0d+00)!(-16.0d+00)   !
ISTATE=1
x=(1000.0d+00)*w0!w1(1)!
xout=0.0d+00!w1(NEQ/2)
ml=1
mu=1
!Yfrg(1:NEQ)=0.0d+00
Yfrg(1:NEQ/2)=0.0d+00       !Imaginary Part
Yfrg((NEQ/2)+1:NEQ)=e0-ep   !Real Part
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL Ffrg (NEQ, x, Yfrg, YDfrg)

!CALL JAC (NEQ, x, Yfrg, ml, mu, pd, nrowpd)

CALL DLSODE (Ffrg, NEQ, Yfrg, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FFRG",ISTATE,xout/w0
!------------------------------------------------------------------
!------------------------------------------------------------------
! FRG zero Temperrature:
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0
IWORK(5:10)=0
IWORK(6)=100000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)
ATOL=w0*10.0d+00**(-13.0d+00)
ISTATE=1
x=((10.0d+00)**4.0d+00)*w0!w1(1)!
xout=0.0d+00!w1(NEQ/2)
ml=1
mu=1
!Yfrg(1:NEQ)=0.0d+00
Yt0frg(1:NEQ/2)=0.0d+00       !Imaginary Part
Yt0frg((NEQ/2)+1:NEQ)=e0-ep   !Real Part (-Ep->convergence!)
!------------------------------------------------------------------
!------------------------------------------------------------------
!CALL Ft0frg (NEQ, x, Yt0frg, YDt0frg)

!!CALL JAC (NEQ, x, Yfrg, ml, mu, pd, nrowpd)

!CALL DLSODE (Ft0frg, NEQ, Yt0frg, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"Ft0FRG",ISTATE,xout/w0

Yt0FRG((NEQ/2)+1:NEQ)=Yt0FRG((NEQ/2)+1:NEQ)-e0!-ep
!------------------------------------------------------------------
Do k=1,NEQ/2
!write(3,*)k,lambda/w0,(e0-ep)/gama,(-e0+YFRG(k+(NEQ/2)))/ep,(((e0-ep)*(exp(2.0d+00*ep/(pi*gama))))-e0)/ep&
!            ,SHartout/ep,-(1.0d+00-(2.0d+00/pi)*atan(e0/gama))
!write(2,*)k,w1(k)/w0,YFRG(k)/gama,YFRG(k+(NEQ/2))/ep,aimag(ana(k))/gama,real(ana(k))/ep,Ynum(k)/gama,Ynum(k+(NEQ/2))/ep
!write(2,*)k,w2(k)/gama,Yt0FRG(k)/gama,Yt0FRG(k+(NEQ/2)),aimag(ana(k))/gama
end do
!------------------------------------------------------------------------
deallocate(RWORK,IWORK)
!------------------------------------------------------------------------
!------------------------------------------------------------------------
ratio2=(w0**4)/((4.0d+00*(gama**2)*(w0**2))+((gama**2+(e0-ep)**2-w0**2)**2))
ratio3=ratio2*(gama**2-(e0-ep)**2+w0**2)/(w0**2)
ratio1=ratio2*(gama**2+(e0-ep)**2-w0**2)/(w0**2)

e0eff=(e0-ep)*(exp(-((lambda/w0)**2)*ratio1))&
     *(exp((ep/(e0-ep))*(1-ratio3)*(SIGN(1.0d+00,(e0-ep))-((2.0d+00/pi)*atan(gama/(e0-ep))))))&
     *(((w0**2)/((e0-ep)**2+gama**2))**((-2.0d+00/pi)*(gama/w0)*((lambda/w0)**2)*ratio2))
!------------------------------------------------------------------------
write(2,*)ilambda,lambda/w0,w0/gama,w1(NEQ/2)/w0,YFRG(NEQ/2)/gama,YFRG((NEQ/2)+(NEQ/2))/gama&
          ,e0eff/gama&
         ,aimag(ana(NEQ/2))/gama,real(ana(NEQ/2))/gama,Ynum(NEQ/2)/gama,Ynum((NEQ/2)+(NEQ/2))/gama
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!        FRG  Dot Occupancy  at zero Temperature
!        using the self-energy at matsubara space ;D
!------------------------------------------------------------------------
MF=10
!------------------------------------------------------------------
! Dimension Declaration of WORK:
!------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ1
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ1+NEQ1**2 
LIW=20+NEQ1
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ1+(2*ML+MU)*NEQ1
LIW=20+NEQ1
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Do ixout=1,1!00

Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=100000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)!0.0d+00*w0*(10.0d+00**(-10.0d+00-(4.0d+00*(ixout-1)/99.0d+00)))   !
ATOL=w0*10.0d+00**(-13.0d+00)!10.0d+00*w0*(10.0d+00**(-10.0d+00-(4.0d+00*(ixout-1)/99.0d+00)))  !

ISTATE=1
!------------------------------------------------------------------
!x=-w0*10.0d+00**(11.0d+00)!(10000.0d+00)*w0
x=w0*10.0d+00**(4.0d+00)!w0*((10.0d+00)**(4.0d+00+(ixout-1)*7.0d+00/99.0d+00))!
xout=0.0d+00!w0*((10.0d+00)**(-10.0d+00))!0.0d+00!

ml=1
mu=1

Rt0FRG(1:NEQ1)=0.0d+00

!CALL Fnt0FRG (NEQ1, x, Rt0FRG, RDt0FRG)
!!CALL JAC (NEQ1, x, Yana, ml, mu, pd, nrowpd)

!CALL DLSODE (Fnt0FRG, NEQ1, Rt0FRG, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

if(ISTATE==2)then
print*,"T0-FRG",ISTATE,ATOL/w0,'integration was successfull'
else
print*,"T0-FRG",ISTATE,ATOL/w0,'NOTE:integration was not successfull'
endif

!write(3,*)ixout,xout/w0,0.5d+00-Yt0ana(1)

end do

!-----------------------------------------------------------------------
ndt0FRG(e0int)=0.5d+00+Rt0FRG(1)
!-----------------------------------------------------------------------
!------------------------------------------------------------------
! gamma^(1)
!------------------------------------------------------------------
!------------------------------------------------------------------
! Analytic Results for gamma^(1)
!------------------------------------------------------------------
g1realana=-((lambda**2)/(2.0d+00*pi))*(8.0d+00*gama*w0*log(gama/w0)&
                  +2.0d+00*pi*(w0**2-gama**2)+4.0d+00*(gama/w0)*(gama**2+w0**2))/((w0**2+gama**2)**2)
!------------------------------------------------------------------
!------------------------------------------------------------------
! Do k=klinear,(NEQ/2)-1
!
! if(abs(aimag(ana(klinear+1))/w1(klinear+1)-aimag(ana(klinear))/w1(klinear))>w0*(10.0d+00**(-11.0d+00))&
!     .or.abs((Ynum(klinear+1))/w1(klinear+1)-(Ynum(klinear))/w1(klinear))>w0*(10.0d+00**(-11.0d+00))&
!     .or.abs((Yfrg(klinear+1))/w1(klinear+1)-(Yfrg(klinear))/w1(klinear))>w0*(10.0d+00**(-11.0d+00)))then
!  print*,"Not converging fast enough for the numerical derivative"&
!        ,abs(aimag(ana(klinear+1))/w1(klinear+1)-aimag(ana(klinear))/w1(klinear))/w0&
!        ,abs((Ynum(klinear+1))/w1(klinear+1)-(Ynum(klinear))/w1(klinear))/w0&
!        ,abs((Yfrg(klinear+1))/w1(klinear+1)-(Yfrg(klinear))/w1(klinear))/w0

! endif

!end do  !k
!------------------------------------------------------------------
!write(1,*),lambda/w0,gama/w0,g1realana,aimag(ana(klinear))/w1(klinear),Ynum(klinear)/w1(klinear)&
!          ,Yfrg(klinear)/w1(klinear),Y1stno,Y1st

print*,lambda/w0,gama/w0,g1realana,aimag(ana((NEQ/2)))/w1((NEQ/2)),Ynum((NEQ/2))/w1((NEQ/2))&
      ,Yfrg((NEQ/2))/w1((NEQ/2)),Y1stno,Y1st
!------------------------------------------------------------------      
!write(1,*)e0int,e0/gama,Ynum(NEQ)/ep,Yfrg(NEQ)/ep,SHart/ep   
!------------------------------------------------------------------
deallocate(RWORK,IWORK)
!------------------------------------------------------------------
! calculating Dot Occupancy with Ozaki method:
!------------------------------------------------------------------
mu0=(ic*Rlarge/(ic*Rlarge+ic*gama-e0))
rhoper=0.0d+00
rhonumozaki=0.0d+00
rhoFRGozaki=0.0d+00
rho0Dot=0.0d+00

Do k=1, NEQ/2

G0pDot(k)=(1.0d+00/(ic*w1(k)+(ic*gama*sgn(w1(k)))-e0+e0array(51)))
G0nDot(k)=(1.0d+00/(-ic*w1(k)+(ic*gama*sgn(-w1(k)))-e0+e0array(51)))

Gpper(k)=(1.0d+00/(ic*w1(k)+(ic*gama*sgn(w1(k)))-ic*aimag(ana(k))-real(ana(k))))
Gnper(k)=(1.0d+00/(-ic*w1(k)+(ic*gama*sgn(-w1(k)))+ic*aimag(ana(k))-real(ana(k))))

Gpnum(k)=(1.0d+00/(ic*w1(k)+(ic*gama*sgn(w1(k)))-YNum((NEQ/2)+k)-ic*YNum(k)))
Gnnum(k)=(1.0d+00/(-ic*w1(k)+(ic*gama*sgn(-w1(k)))-YNum((NEQ/2)+k)+ic*YNum(k)))

Gpfrg(k)=(1.0d+00/(ic*w1(k)+(ic*gama*sgn(w1(k)))-Yfrg((NEQ/2)+k)-ic*Yfrg(k)))
Gnfrg(k)=(1.0d+00/(-ic*w1(k)+(ic*gama*sgn(-w1(k)))-Yfrg((NEQ/2)+k)+ic*Yfrg(k)))

rho0Dot=rho0Dot-(Rp(k)*real(G0pDot(k)+G0nDot(k))/beta)  
rhoper=rhoper-(Rp(k)*real(Gpper(k)+Gnper(k))/beta)  
rhonumozaki=rhonumozaki-(Rp(k)*real(Gpnum(k)+Gnnum(k))/beta) 
rhoFRGozaki=rhoFRGozaki-(Rp(k)*real(GpFRG(k)+GnFRG(k))/beta) 

!-----------------------------------------------------------------------
! Linear Conductance:
!-----------------------------------------------------------------------
if(k==1)then

Kper=0.0d+00
KFRGozaki=0.0d+00

else

Kper=Kper+2.0d+00*(gama/beta)*Rp(k)*(aimag(Gpper(k)-Gpper(k-1)))/(w1(k)-w1(k-1))
KFRGozaki=KFRGozaki+2.0d+00*(gama/beta)*Rp(k)*(aimag(Gpfrg(k)-Gpfrg(k-1)))/(w1(k)-w1(k-1))

endif
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end do
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
nd0Dot(e0int)=(mu0/2.0d+00)-rho0Dot
ndper(e0int)=(mu0/2.0d+00)-rhoper
ndnumozaki(e0int)=(mu0/2.0d+00)-rhonumozaki
ndFRGozaki(e0int)=(mu0/2.0d+00)-rhoFRGozaki
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Perturbative Self-energy from performing the Analytic
! continuation analytically ;D
!-----------------------------------------------------------------------
if(w0>5*gama)then
width=gama*(log10(w0/gama))*6.0d+00

ktilda=0.5d+00*(NTAB+((ki*width)/(2.0d+00*(w0-width))))/(1.0d+00+((width)/(2.0d+00*(w0-width))))

if(20.0d+00*gama-2.0d+00*w0+width>0)then
ki=100!(ktilda/19.4)/(1+(2*(w0-width)/(20*gama-2*w0+width)))
else
ki=1
endif

else 

width=w0*0.9d+00
ki=1000.0d+00*log10(7*gama/w0)

if(w0<gama)then
ktilda=1000.0d+00+ki+2000.0d+00*log(2.0d+00*gama/w0)
else
ktilda=1000.0d+00+ki
endif

endif
!-----------------------------------------------------------------------
Do k=1,2*NTAB

if(k<NTAB+1)then

if(k<ki)then

XTAB(k)=(-2*w0+width)+(20.0d+00*gama-2*w0+width)*((k-((1+2.0d+00*ki)/2.0))/((2.0d+00*ki-1)/2.0d+00))&
   *exp((abs(k-((1+2.0d+00*ki)/2.0d+00))-((2*ki-1)/2.0d+00))/230.0d+00)

elseif(k>=ki.and.k<ktilda)then

XTAB(k)=-w0+(w0-width)*((k-((ki+ktilda)/2.0d+00))/((ktilda-ki)/2.0d+00))&
   *exp((abs(k-((ki+ktilda)/2.0d+00))-((ktilda-ki)/2.0d+00))/100.0d+00)

else if (k>=ktilda)then

XTAB(k)=width*((k-NTAB)/(1.0d+00*NTAB-ktilda))*exp((abs(k-NTAB)-(NTAB-ktilda))/230.0d+00)

endif

else
XTAB(k)=-XTAB(2*NTAB+1-k)
endif

!(for nonzero T calculation!)
!XTAB(k)=(4.0d+00*w0)*(((k-(NTAB/2.0d+00))*exp((abs((NTAB/2.0d+00)-k)-(NTAB/2.0d+00))/23.0d+00)/(NTAB/2.0d+00)))  
!XTAB(k)=0.0d+00
!-----------------------------------------------------------------------

ap=1.0d+00/(gama**2+(e0+SHart-w0-XTAB(k))**2)
am=-1.0d+00/(gama**2+(e0+SHart+w0-XTAB(k))**2)


sana(k)=ic*((lambda**2)/2.0d+00)*((gama*am*tanh(beta*(XTAB(k)-w0)/2.0d+00)+gama*ap*tanh(beta*(XTAB(k)+w0)/2.0d+00)))&
        -((lambda**2)/(2.0d+00*pi))*(-gama*(ap+am)*log((e0+SHart)**2+gama**2)&
                          +2*gama*(ap*log(abs(XTAB(k)+w0))+am*log(abs(XTAB(k)-w0)))&
                          +2*(ap*(XTAB(k)+w0-(e0+SHart))+am*(XTAB(k)-w0-(e0+SHart)))*atan2((e0+SHart),gama))&
                          +((lambda**2)/2.0d+00)*((1.0d+00/(XTAB(k)-w0-(e0+SHart)+ic*gama))&
                                             +(1.0d+00/(XTAB(k)+w0-(e0+SHart)+ic*gama)))&
         +SHartout
!-----------------------------------------------------------------------
rananum(k)=(-1.0d+00/pi)*(-gama+(aimag(sana(k))))/((XTAB(k)-e0-real(sana(k)))**2+(gama-(aimag(sana(k))))**2)
!write(4,*)k,XTAB(k)/w0,pi*gama*rananum(k)
!write(2,*)k,XTAB(k)/w0,pi*gama*rananum(k)

if(k==1)then
dw1=0.0!abs(100*(gama*2)*(((2-NTAB)*exp((abs(NTAB-2)-NTAB)/230.0)/NTAB))-w1)
else
dw1=abs(XTAB(k)-XTAB(k-1))
endif
ndananum(e0int)=ndananum(e0int)+(dw1)*rananum(k)
!-----------------------------------------------------------------------
!write(3,*)k,XTAB(k)/w0,(-SHartout+real(sana(k)))/ep,aimag(sana(k))/gama
!-----------------------------------------------------------------------
end do
!-----------------------------------------------------------------------
! Linear Conductance:
!-----------------------------------------------------------------------
!write(2,*)e0int,(e0array(e0int)-ep)/gama,pi*gama*rananum(1),Kper*pi,KFRGozaki*pi
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
if(e0int==51)then
g1FRG=Yfrg((NEQ/2))/w1((NEQ/2))
g1Num=Ynum((NEQ/2))/w1((NEQ/2))
endif
!-----------------------------------------------------------------------
end do ! e0int
!-----------------------------------------------------------------------
! effective hybridization:
!-----------------------------------------------------------------------
get0ana=(e0array(52)-e0array(51))/(pi*(ndt0ana(51)-ndt0ana(52)))

geanaana=(e0array(52)-e0array(51))/(pi*(ndanarout(51)-ndanarout(52)))

geanaozaki=(e0array(52)-e0array(51))/(pi*(ndper(51)-ndper(52)))

genumozaki=(e0array(52)-e0array(51))/(pi*(ndnumozaki(51)-ndnumozaki(52)))

gefrgozaki=(e0array(52)-e0array(51))/(pi*(ndFRGozaki(51)-ndFRGozaki(52)))

geananum=(e0array(52)-e0array(51))/(pi*(ndananum(51)-ndananum(52)))

get0FRG=(e0array(52)-e0array(51))/(pi*(ndt0FRG(51)-ndt0FRG(52)))

!write(3,*)NEQ/2,lambda/w0,gama/w0,(e0array(e0int)-e0array(51))/gama,get0ana/gama,geanaana/gama,geanaozaki/gama&
!           ,genumozaki/gama,gefrgozaki/gama,get0FRG/gama,(1.0d+00+g1realana),(1.0d+00+g1FRG),(1.0d+00+g1Num)!&
         !,1.0d+00/(1.0d+00-aimag(ana((NEQ/2)))/w1((NEQ/2))),1.0d+00/(1.0d+00-Ynum((NEQ/2))/w1((NEQ/2)))
!write(3,*)NEQ/2,ilambda,gama/w0,(e0array(e0int)-e0array(51))/gama,gefrgozaki/gama&
!          ,(1.0d+00+g1realana),(1.0d+00+g1FRG),1.0d+00/(1.0d+00-g1FRG),lambda**2/(w0*gama)       
!------------------------------------------------------------------
! deallocate(Ynum,YDnum,YFRG,YDFRG,ana)
deallocate(ana)
!------------------------------------------------------------------
end do !ineq         
!------------------------------------------------------------------
end do ! ilambda
!------------------------------------------------------------------
close(1)
!------------------------------------------------------------------
end program srlm_quadpack_ozaki
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
SUBROUTINE  Fno (NEQ1, x, Y1stno, YD1stno)
implicit none

common w0,lambda,gama,pi,e0

INTEGER::NEQ1
DOUBLE PRECISION::x,e0
DOUBLE PRECISION,DIMENSION(NEQ1)::Y1stno, YD1stno

double precision::w0,lambda,gama,pi

YD1stno(1:NEQ1)=(4.0d+00*w0*(lambda**2)/pi)*(1.0d+00/(gama+x))*(x/((x**2+w0**2)**2))

End SUBROUTINE Fno
!------------------------------------------------------------------
!------------------------------------------------------------------
SUBROUTINE  F (NEQ1, x, Y1st, YD1st)
implicit none

common w0,lambda,gama,pi,e0

INTEGER::NEQ1
DOUBLE PRECISION::x,e0
DOUBLE PRECISION,DIMENSION(NEQ1)::Y1st, YD1st

double precision::w0,lambda,gama,pi

YD1st(1:NEQ1)=(4.0d+00*w0*(lambda**2)/pi)*(1.0d+00/(gama+x*(1.0d+00-Y1st(1:NEQ1))))*(x/((x**2+w0**2)**2))

End SUBROUTINE F
!------------------------------------------------------------------
!------------------------------------------------------------------
!----------------------------------------------------------------------
!*****************NUM:***************************
!----------------------------------------------------------------------
subroutine Fnum(NEQ, x, Ynum, YDnum)

implicit none

common w0,lambda,gama,pi,e0,SHart,SHartout,w1

INTEGER::k,NEQ
DOUBLE PRECISION::x
DOUBLE PRECISION,DIMENSION(NEQ)::Ynum,YDnum
DOUBLE PRECISION,DIMENSION(4000)::w1
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,SHart,SHartout

!----------------------------------------------------------------------
ep=(lambda**2)/w0

Do k=1,NEQ/2
!w1(k)=(exp(-(k-1)*0.012d+00))*10000.0d+00*w0
!w1(NEQ/2)=0.0d+00
end do

!YDnum(1:NEQ)=(w0*(lambda**2)/pi)*(4.0*w1(1:NEQ)*x)/((x+gama)&
!                 *((w1(1:NEQ)+x)**2+w0**2)*((w1(1:NEQ)-x)**2+w0**2))  


YDnum(1:NEQ/2)=(-w0*(lambda**2)/pi)*((-1.0d+00/(w0**2+(x-w1(1:NEQ/2))**2))&
                                     +(1.0d+00/(w0**2+(x+w1(1:NEQ/2))**2)))&
                                   *(x+gama)/((x+gama)**2+(e0+SHart)**2)

YDnum(NEQ/2+1:NEQ)=(w0*(lambda**2)/pi)*((1.0d+00/(w0**2+(x-w1(1:NEQ/2))**2))&
                                   +(1.0d+00/(w0**2+(x+w1(1:NEQ/2))**2)))&
                                   *(e0+SHart)/((x+gama)**2+(e0+SHart)**2)

!write(1,*)x/w0,Ynum(768)/gama,w1(768)/w0
!write(1,*)x/w0,Ynum(NEQ/2)/gama,Ynum(NEQ)/ep,SHart/ep
                 
end subroutine Fnum
!----------------------------------------------------------------------
!*****************FRG:***************************
!----------------------------------------------------------------------
subroutine Ffrg(NEQ, x, Yfrg, YDfrg)

implicit none

common w0,lambda,gama,pi,e0,SHart,SHartout,w1

INTEGER::k,NEQ,founder,ki,xcounter
DOUBLE PRECISION::x,SHart,SHartout
DOUBLE PRECISION,DIMENSION(NEQ)::Yfrg,YDfrg
DOUBLE PRECISION,DIMENSION(4000)::w1
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,refeedback,e0,imfeedback

!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
integer*4, parameter::DIM_P=1          !the spatial dimension.
!integer*4, parameter::DATA_NUM=3000    !the number of data points
integer*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

real*8::T_DATA(NEQ/2)               !the value of the independent variable at the sample points
real*8::P_DATA1(DIM_P,NEQ/2)         !the value of the dependent variables at the sample points
real*8::P_DATA2(DIM_P,NEQ/2)         !the value of the dependent variables at the sample points
real*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
real*8:: P_INTERP1(DIM_P,NEQ/2)      !the interpolated values of the dependent variables
real*8:: P_INTERP2(DIM_P,NEQ/2)      !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
founder=0
ep=(lambda**2)/w0

!Do k=1,NEQ
Do k=NEQ/2,1,-1
!w1(k)=(exp(-(k-1)*0.012d+00))*10000.0d+00*w0
!w1(NEQ/2)=0.0d+00

T_DATA((NEQ/2+1)-k)=w1(k)
P_DATA1(DIM_P,(NEQ/2+1)-k)=Yfrg(k)        !Imaginary Part
P_DATA2(DIM_P,(NEQ/2+1)-k)=Yfrg(NEQ/2+k)   !Real Part

!----------------------------------------------------------------------
if(w1(k)<x.and.founder==0)then
ki=k
founder=founder+1
endif

!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=x
call interp_linear(DIM_P,NEQ/2,T_DATA,P_DATA1,INTERP_NUM,T_INTERP,P_INTERP1)
imfeedback=P_INTERP1(DIM_P,INTERP_NUM)
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=x
call interp_linear(DIM_P,NEQ/2,T_DATA,P_DATA2,INTERP_NUM,T_INTERP,P_INTERP2)
refeedback=P_INTERP2(DIM_P,INTERP_NUM)
!----------------------------------------------------------------------

!YDfrg(1:NEQ/2)=(w0*(lambda**2)/pi)*(4.0*w1(1:NEQ)*x)/((x+gama-feedback)&
!                 *((w1(1:NEQ)+x)**2+w0**2)*((w1(1:NEQ)-x)**2+w0**2))                

YDfrg(1:NEQ/2)=(-w0*(lambda**2)/pi)*((-1.0d+00/(w0**2+(x-w1(1:NEQ/2))**2))&
                                     +(1.0d+00/(w0**2+(x+w1(1:NEQ/2))**2)))&
                                   *(x+gama-imfeedback)/((x+gama-imfeedback)**2+(refeedback)**2)
                                   
                                   
YDfrg(NEQ/2+1:NEQ)=(w0*(lambda**2)/pi)*((1.0d+00/(w0**2+(x-w1(1:NEQ/2))**2))&
                                   +(1.0d+00/(w0**2+(x+w1(1:NEQ/2))**2))-(2.0d+00/(w0**2)))&
                                   *(refeedback)/((x+gama-imfeedback)**2+(refeedback)**2)
                 
if(x==w1(1))then
xcounter=0
print*,xcounter
else
xcounter=xcounter+1
endif                 
                 
!write(2,*)x/w0,Yfrg(768)/gama,w1(768)/w0  
!Do k=1, NEQ
!write(1,*)xcounter,x/w0,feedback

!if(xcounter==35256)then
!Do k=1, NEQ
!write(1,*)k,w1(k)/w0,Yfrg(k)/gama,x/w0,xcounter
!end do
!endif

!end do
!write(2,*)x/w0,Yfrg(NEQ/2)/gama,Yfrg(NEQ)/ep

end subroutine Ffrg
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!**********Zeo Temeprature FRG:***************************
!----------------------------------------------------------------------
subroutine Ft0frg(NEQ, x, Yt0frg, YDt0frg)

implicit none

common w0,lambda,gama,pi,e0,SHart,SHartout,w1,w2

INTEGER::k,NEQ,founder,ki,xcounter
DOUBLE PRECISION::x,SHart,SHartout
DOUBLE PRECISION,DIMENSION(NEQ)::Yt0frg,YDt0frg
DOUBLE PRECISION,DIMENSION(4000)::w1,w2
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,refeedback,e0,imfeedback

!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
integer*4, parameter::DIM_P=1        !the spatial dimension.
!integer*4, parameter::DATA_NUM=3000 !the number of data points
integer*4, parameter::INTERP_NUM=1   !the number of points at which interpolation is to be done

real*8::T_DATA(NEQ/2)                !the value of the independent variable at the sample points
real*8::P_DATA1(DIM_P,NEQ/2)         !the value of the dependent variables at the sample points
real*8::P_DATA2(DIM_P,NEQ/2)         !the value of the dependent variables at the sample points
real*8::T_INTERP(INTERP_NUM)         !the value of the independent variable at the interpolation points
real*8::P_INTERP1(DIM_P,NEQ/2)       !the interpolated values of the dependent variables
real*8::P_INTERP2(DIM_P,NEQ/2)       !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
founder=0
ep=(lambda**2)/w0

!Do k=1,NEQ
Do k=NEQ/2,1,-1
!w1(k)=(exp(-(k-1)*0.012d+00))*10000.0d+00*w0
!w1(NEQ/2)=0.0d+00

T_DATA((NEQ/2+1)-k)=w2(k)
P_DATA1(DIM_P,(NEQ/2+1)-k)=Yt0frg(k)         !Imaginary Part
P_DATA2(DIM_P,(NEQ/2+1)-k)=Yt0frg(NEQ/2+k)   !Real Part

!----------------------------------------------------------------------
if(w2(k)<x.and.founder==0)then
ki=k
founder=founder+1
endif

!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=x
call interp_linear(DIM_P,NEQ/2,T_DATA,P_DATA1,INTERP_NUM,T_INTERP,P_INTERP1)
imfeedback=P_INTERP1(DIM_P,INTERP_NUM)
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=x
call interp_linear(DIM_P,NEQ/2,T_DATA,P_DATA2,INTERP_NUM,T_INTERP,P_INTERP2)
refeedback=P_INTERP2(DIM_P,INTERP_NUM)
!----------------------------------------------------------------------

!YDfrg(1:NEQ/2)=(w0*(lambda**2)/pi)*(4.0*w1(1:NEQ)*x)/((x+gama-feedback)&
!                 *((w1(1:NEQ)+x)**2+w0**2)*((w1(1:NEQ)-x)**2+w0**2))                

YDt0frg(1:NEQ/2)=(-w0*(lambda**2)/pi)*((-1.0d+00/(w0**2+(x-w2(1:NEQ/2))**2))&
                                     +(1.0d+00/(w0**2+(x+w2(1:NEQ/2))**2)))&
                                   *(x+gama-imfeedback)/((x+gama-imfeedback)**2+(refeedback)**2)

YDt0frg(NEQ/2+1:NEQ)=(w0*(lambda**2)/pi)*((1.0d+00/(w0**2+(x-w2(1:NEQ/2))**2))&
                                   +(1.0d+00/(w0**2+(x+w2(1:NEQ/2))**2))-(2.0d+00/(w0**2)))&
                                   *(refeedback)/((x+gama-imfeedback)**2+(refeedback)**2)
                 
if(x==w2(1))then
xcounter=0
print*,xcounter
else
xcounter=xcounter+1
endif                 
                 
!write(2,*)x/w0,Yfrg(768)/gama,w1(768)/w0  
!Do k=1, NEQ
!write(1,*)xcounter,x/w0,feedback

!if(xcounter==35256)then
!Do k=1, NEQ
!write(1,*)k,w1(k)/w0,Yfrg(k)/gama,x/w0,xcounter
!end do
!endif

!end do
!write(2,*)x/w0,Yfrg(NEQ/2)/gama,Yfrg(NEQ)/ep

end subroutine Ft0frg
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!****Dot Occupation from zero Temeprature FRG:*************************
!----------------------------------------------------------------------
subroutine Fnt0frg(NEQ1, x, Rt0frg, RDt0frg)

implicit none

common w0,lambda,gama,pi,e0,SHart,SHartout,w1,w2,Yt0FRG

INTEGER,PARAMETER::NEQ=8000
INTEGER::k,NEQ1,founder,ki,xcounter
DOUBLE PRECISION::x,SHart,SHartout
DOUBLE PRECISION,DIMENSION(NEQ1)::Rt0frg,RDt0frg
DOUBLE PRECISION,DIMENSION(4000)::w1,w2
DOUBLE PRECISION,DIMENSION(8000)::Yt0FRG
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,refeedback,e0,imfeedback

!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
integer*4, parameter::DIM_P=1          !the spatial dimension.
!integer*4, parameter::DATA_NUM=3000    !the number of data points
integer*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

real*8::T_DATA(NEQ/2)               !the value of the independent variable at the sample points
real*8::P_DATA1(DIM_P,NEQ/2)         !the value of the dependent variables at the sample points
real*8::P_DATA2(DIM_P,NEQ/2)         !the value of the dependent variables at the sample points
real*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
real*8:: P_INTERP1(DIM_P,NEQ/2)      !the interpolated values of the dependent variables
real*8:: P_INTERP2(DIM_P,NEQ/2)      !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
founder=0
ep=(lambda**2)/w0

!Do k=1,NEQ
Do k=NEQ/2,1,-1
!w1(k)=(exp(-(k-1)*0.012d+00))*10000.0d+00*w0
!w1(NEQ/2)=0.0d+00

T_DATA((NEQ/2+1)-k)=w2(k)
P_DATA1(DIM_P,(NEQ/2+1)-k)=Yt0frg(k)         !Imaginary Part
P_DATA2(DIM_P,(NEQ/2+1)-k)=Yt0frg(NEQ/2+k)   !Real Part

!----------------------------------------------------------------------
if(w2(k)<x.and.founder==0)then
ki=k
founder=founder+1
endif

!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=x
call interp_linear(DIM_P,NEQ/2,T_DATA,P_DATA1,INTERP_NUM,T_INTERP,P_INTERP1)
imfeedback=P_INTERP1(DIM_P,INTERP_NUM)
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=x
call interp_linear(DIM_P,NEQ/2,T_DATA,P_DATA2,INTERP_NUM,T_INTERP,P_INTERP2)
refeedback=e0+P_INTERP2(DIM_P,INTERP_NUM)
!----------------------------------------------------------------------
              
!IM: YDt0frg(1:NEQ/2) 
!RE: YDt0frg(NEQ/2+1:NEQ) 

RDt0frg(1:NEQ1)=(1.0d+00/pi)*(refeedback)/((x+gama-imfeedback)**2+(refeedback)**2)

end subroutine Fnt0frg
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!---Analytic Perturbation theory after Analytic continuation:-----------
!----------------------------------------------------------------------
subroutine Fana(NEQ1, x, Yana, YDana)

implicit none

common w0,lambda,gama,pi,e0,SHart,SHartout

INTEGER::k,NEQ1,sgnp,sgnm
DOUBLE PRECISION::x,fermi
DOUBLE PRECISION,DIMENSION(NEQ1)::Yana,YDana
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,SHart,ap,am,beta,SHartout
COMPLEX*16::sana
COMPLEX::ic
!----------------------------------------------------------------------
ep=(lambda**2)/w0
beta=1.0d+00/(0.0001d+00*gama)
ic=cmplx(0.0d+0,1.0d+0, kind = 16)
!-----------------------------------------------------------------------

ap=1.0d+00/(gama**2+(e0+SHart-w0-x)**2)
am=-1.0d+00/(gama**2+(e0+SHart+w0-x)**2)

if(x>0.0d+00)then
fermi=0.0d+00
elseif(x==0.0d+00)then
fermi=0.50d+00
else
fermi=1.0d+00
endif

if(x-w0>0.0d+00)then
sgnp=1
elseif(x-w0==0.0d+00)then
sgnp=0
else
sgnp=-1
endif

if(x+w0>0.0d+00)then
sgnm=1
elseif(x+w0==0.0d+00)then
sgnm=0
else
sgnm=-1
endif
!-------------------------------------------------------------
!-------------------------------------------------------------        
sana=ic*((lambda**2)/2.0d+00)*((gama*am*tanh(beta*(x-w0)/2.0d+00)+gama*ap*tanh(beta*(x+w0)/2.0d+00)))&
        -((lambda**2)/(2.0d+00*pi))*(-gama*(ap+am)*log((e0+SHart)**2+gama**2)&
                          +2.0d+00*gama*(ap*log(abs(x+w0))+am*log(abs(x-w0)))&
                          +2.0d+00*(ap*(x+w0-(e0+SHart))+am*(x-w0-(e0+SHart)))*atan2((e0+SHart),gama))&
                          +((lambda**2)/2.0d+00)*((1.0d+00/(x-w0-(e0+SHart)+ic*gama))&
                                             +(1.0d+00/(x+w0-(e0+SHart)+ic*gama)))&
         +SHartout
!-----------------------------------------------------------------------         
sana=ic*((lambda**2)/2.0d+00)*((gama*am*sgnp+gama*ap*sgnm))&
        -((lambda**2)/(2.0d+00*pi))*(-gama*(ap+am)*log((e0+SHart)**2+gama**2)&
                          +2.0d+00*gama*(ap*log(abs(x+w0))+am*log(abs(x-w0)))&
                          +2.0d+00*(ap*(x+w0-(e0+SHart))+am*(x-w0-(e0+SHart)))*atan2((e0+SHart),gama))&
                          +((lambda**2)/2.0d+00)*((1.0d+00/(x-w0-(e0+SHart)+ic*gama))&
                                             +(1.0d+00/(x+w0-(e0+SHart)+ic*gama)))&
         +SHartout  
         
         
!fermi=(1.0d+00-tanh(beta*x/2.0d+00))/2.0d+00         
!-----------------------------------------------------------------------
!(1.0d+00-tanh(beta*x/2.0d+00))

YDana(1:NEQ1)=1.0d+00*fermi*(1.0d+00/pi)&
                     *(gama-(aimag(sana)))/((x-e0-real(sana))**2+(gama-(aimag(sana)))**2)

!write(4,*)x,pi*gama*Yana(1:NEQ1)                    
end subroutine Fana
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!---Analytic Perturbation theory before Analytic continuation:-----------
!------------------------------------------------------------------------
subroutine Ft0ana(NEQ1, x, Yt0ana, YDt0ana)

implicit none

common w0,lambda,gama,pi,e0,SHart,SHartout

INTEGER::k,NEQ1,sgnne0
DOUBLE PRECISION::x,fermi
DOUBLE PRECISION,DIMENSION(NEQ1)::Yt0ana,YDt0ana
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,SHart,beta,SHartout
COMPLEX*16::anat0,manat0,dp,dm,dtp,dtm
COMPLEX::ic
!----------------------------------------------------------------------
ep=(lambda**2)/w0
beta=1.0d+00/(0.0001d+00*gama)
ic=cmplx(0.0d+0,1.0d+0, kind = 16)
!-----------------------------------------------------------------------
if(e0+SHart<0)then
sgnne0=-1
else
sgnne0=1
endif
!-------------------------------------------------------------
!-------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------

dp=-1.0d+00/(e0+SHart+w0-ic*(x+gama))
dm=-1.0d+00/(e0+SHart-w0-ic*(x+gama))
dtp=-1.0d+00/(e0+SHart+w0-ic*(x-gama))
dtm=-1.0d+00/(e0+SHart-w0-ic*(x-gama))        

anat0=((ic*(lambda**2)/(2.0d+00*pi))*(0.5d+00*(dtp-dtm-dp+dm)*log(gama**2+(e0+SHart)**2)&
                             -0.5d+00*(dtp-dtm-dp+dm)*log(w0**2+x**2)&
                             +ic*(dp-dtp)*atan2(w0,x)+ic*(dtm-dm)*atan2(-w0,x)&
                             -ic*(dp+dm)*pi+ic*(dtp-dtm)*atan2(e0+SHart,-gama)&
                             -ic*(dp-dm)*atan2(e0+SHart,gama)&
                             -ic*pi*(dtp-dtm)*sgnne0))&
                             +SHartout 
                           
!---------------------------------------------------------------------- 
!-----------------------------------------------------------------------         

YDt0ana(1:NEQ1)=(1.0d+00/pi)*(e0+real(anat0))/((x+gama-aimag(anat0))**2+(e0+real(anat0))**2)
                   
end subroutine Ft0ana
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine JAC (neq1, x, y, ml, mu, pd, nrowpd)
!------------------------------------------------------------------
integer,intent (in) ::neq1, ml, mu, nrowpd
double precision::x
DOUBLE PRECISION,dimension(NEQ1)::y
DOUBLE PRECISION,dimension(nrowpd,NEQ1)::pd
pd(1,1) = 0.0d0
pd(1,2) = 1.0d0
pd(2,1) = -6.0d0*y(1)*y(2) - 1.0d0
pd(2,2) = 3.0d0*(1.0d0 - y(1)*y(1))
return
end
!-------------------------------------------------------------
!-------------------------------------------------------------
!*************************************************************
!   Sgn function:
!*************************************************************
integer*8 function sgn(w)
implicit none
double precision,intent(in)::w
if(w>0)then
sgn=1
elseif(w==0)then
sgn=0
else
sgn=-1
endif
end function sgn
!-------------------------------------------------------------
integer*8 function sgnn(w)
implicit none
double precision,intent(in)::w
if(w>0)then
sgnn=1
elseif(w==0)then
sgnn=1
else
sgnn=-1
endif
end function sgnn
!-------------------------------------------------------------
!-------------------------------------------------------------
