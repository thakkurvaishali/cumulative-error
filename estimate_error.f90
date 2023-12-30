PROGRAM error_estimation
IMPLICIT NONE
REAL*8,ALLOCATABLE::mean(:),var(:),delta_G(:),xcor(:),k(:),sum2(:),xmean(:)
REAL*8::dr,dum2,sum1,dummy,kbt,num2,num1,temp,r0,dum,x
INTEGER::ir,nr,nsteps,i,i_md,j,skip_steps,nblocks,b,bs
CHARACTER(LEN=50),ALLOCATABLE::filename(:)
REAL*8, PARAMETER :: kb=1.98E-3 !kcal-1 mol-1 K-1 
!----------------------------------------------------------

!reading filenames
! dr=0.05
! prefac=(k*dr)*(k*dr)
nr=24
temp=300.d0
kbt=kb*temp
bs=5000 !block_size
!----------------------------------------------------------

ALLOCATE(filename(nr),xcor(nr),k(nr),sum2(nr))
OPEN(2,file='filename.inp',status='old')

!----------------------------------------------------------

 DO ir=1,nr
   READ(2,'(a)')filename(ir)
   READ(2,*)xcor(ir),k(ir)        !reading mean and restraint
   k(ir)=k(ir)*627.5d0           !au to Kcal/mol
   WRITE(*,'(a)')filename(ir)
 END DO

!----------------------------------------------------------

!OPEN(10,file=filename(1))
! CALL GET_STEPS(10,nsteps)
!CLOSE(10)

nsteps=20000
skip_steps=15000
PRINT*,nsteps
nblocks=nsteps/bs
print*,nblocks
OPEN(21,file="variance.dat")
ALLOCATE(var(nr),xmean(nblocks))

!----------------------------------------------------------

 DO ir=1,nr
   OPEN(11,file=filename(ir),status="old")
   dum=0.0d0
   PRINT*,"WORKING FOR ir=",ir
   DO j=1,skip_steps
      READ(11,*) 
   END DO
      b=0;dum=0.0d0;xmean=0.0d0
      do j=1,nblocks
         do i_md=1,bs                          !loop over steps in block-j
            READ(11,*)dummy,dummy,r0,x
            x=r0+x
            xmean(j)=xmean(j)+x
         end do
         xmean(j)=xmean(j)/dfloat(bs)
         b=b+1
         dum=dum+xmean(j)
      end do
      dum=dum/dfloat(nblocks)

!----------------------------------------------------------

      dummy=0.0d0
      do j=1,nblocks
      var(ir)=dummy+(dum-xmean(j))**2
      end do
      var(ir)=var(ir)/dfloat(nblocks-1)
      write(21,*)ir,xcor(ir),var(ir)
    CLOSE(11)  
 END DO

!----------------------------------------------------------
CLOSE(21)

OPEN(21,file="delta_G.dat")
ALLOCATE(delta_G(nr))

!----------------------------------------------------------
 sum2(1)=0.0d0     !nr is the referenece point
 DO ir=nr,2,-1
   sum2(ir)=(var(ir)*k(ir)**2+var(ir-1)*k(ir-1)**2)*0.25d0
   dr=(xcor(ir)-xcor(ir-1))**2
   sum2(ir)=sum2(ir)*dr
 END DO 
 dum2=0.0d0
 DO ir=1,nr
    dum2=dum2+sum2(ir)
    delta_G(ir)=dum2
! END DO
! 
! DO ir=1,nr
    WRITE(21,*)ir,xcor(ir),delta_G(ir),dsqrt(delta_G(ir))
 END DO
!-----------------------------------------------------------
END PROGRAM error_estimation
   
SUBROUTINE GET_STEPS(filenumber,iter)
 INTEGER::filenumber,iter,ios
iter=0
 DO 
     READ(filenumber,*,iostat=ios)
     IF (ios.ne.0) exit
     iter=iter+1   
END DO
END SUBROUTINE 
