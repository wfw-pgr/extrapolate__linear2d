program main
  implicit none
  integer         , parameter   :: lun     = 50
  integer         , parameter   :: cLen    = 300
  character(cLen) , parameter   :: inpFile = "dat/input.dat"
  character(cLen) , parameter   :: outFile = "dat/output.dat"
  integer                       :: i, j, k, LI, LJ, LK, nlines, ncmp
  double precision, allocatable :: mshape(:,:,:,:), data(:,:,:)
  integer         , allocatable :: flag(:,:)
  character(cLen)               :: cmt
  integer         , parameter   :: x_=1, y_=2, z_=3, i_=4, s_=5, f_=6
  
  ! ------------------------------------------------------ !
  ! --- [1] load mshape_svd.dat                        --- !
  ! ------------------------------------------------------ !
  open(lun,file=trim(inpFile),status="old")
  read(lun,*)
  read(lun,*) cmt, nlines, ncmp
  read(lun,*) cmt, LK, LJ, LI, ncmp
  allocate( mshape(ncmp,LI,LJ,LK) )
  do k=1, LK
     do j=1, LJ
        do i=1, LI
           read(lun,*) mshape(:,i,j,k)
        enddo
     enddo
  enddo
  close(lun)

  allocate( flag(LI,LJ), data(3,LI,LJ) )
  k = 1
  do j=1, LJ
     do i=1, LI
        data(:,i,j) =      mshape(x_:z_,i,j,k)
        flag(  i,j) = int( mshape(   f_,i,j,k) )
     enddo
  enddo

  ! ------------------------------------------------------ !
  ! --- [2] call extrapolation                         --- !
  ! ------------------------------------------------------ !
  call extrapolate__linear2d( data, flag, LI,LJ )
  do j=1, LJ
     do i=1, LI
        mshape(z_,i,j,k) = data(z_,i,j)
     enddo
  enddo

  ! ------------------------------------------------------ !
  ! --- [3] save results                               --- !
  ! ------------------------------------------------------ !
  open(lun,file=trim(outFile),status="replace")
  write(lun,"(a)") "# x_ y_ z_ i_ s_ f_"
  write(lun,"(a,2(i10,1x))") "# ", nlines, ncmp
  write(lun,"(a,4(i10,1x))") "# ", LK, LJ, LI, ncmp
  do k=1, LK
     do j=1, LJ
        do i=1, LI
           write(lun,"(6(e15.8,1x))") mshape(:,i,j,k)
        enddo
     enddo
  enddo
  close(lun)
  
contains

  ! ====================================================== !
  ! === extrapolation subroutine                       === !
  ! ====================================================== !
  subroutine extrapolate__linear2d( data, flag, LI,LJ )
    implicit none
    integer         , intent(in)    :: LI, LJ
    integer         , intent(inout) :: flag(  LI,LJ)
    double precision, intent(inout) :: data(3,LI,LJ)
    integer         , parameter     :: x_=1, y_=2, z_=3
    double precision, allocatable   :: data_(:,:,:)
    integer         , allocatable   :: flag_(:,:)
    integer                         :: i, j, igen, genMax, count, surroundflag
    double precision                :: avg, ext

    genMax = int( 1.5 * max( LI, LJ ) )
    allocate( data_(3,0:LI+1,0:LJ+1), flag_(0:LI+1,0:LJ+1) )
    
    do igen=1, genMax
       
       ! -- copy data -- !
       data_(:,:,:)       = 0.d0
       flag_(:,:)         = 0
       data_(:,1:LI,1:LJ) = data(:,:,:)
       flag_(  1:LI,1:LJ) = flag(:,:)
       count              = 0

       ! -- extrapolate neibour -- !
       do j=2, LJ-1
          do i=2, LI-1
             
             surroundflag = ( + flag_(i,j-1) + flag_(i-1,j) &
                  &           + flag_(i,j+1) + flag_(i+1,j) )
             
             if ( ( flag_(i,j).eq.0 ).and.( surroundflag.gt.0 ) ) then

                count = count + 1
                
                if      ( surroundflag.eq.1 ) then
                   ! e.g.) :: data(i,j) = (i,j-1) + ( (i,j-1) - (i,j-2) )
                   !                    = 2 * (i,j-1) - (i,j-2)
                   if      ( flag_(i,j-1).eq.1 ) then
                      data(z_,i,j) = + 2.d0 * data_(z_,i,j-1) - data_(z_,i,j-2)
                   else if ( flag_(i-1,j).eq.1 ) then
                      data(z_,i,j) = + 2.d0 * data_(z_,i-1,j) - data_(z_,i-2,j)
                   else if ( flag_(i,j+1).eq.1 ) then
                      data(z_,i,j) = + 2.d0 * data_(z_,i,j+1) - data_(z_,i,j+2)
                   else if ( flag_(i+1,j).eq.1 ) then 
                      data(z_,i,j) = + 2.d0 * data_(z_,i+1,j) - data_(z_,i+2,j)
                   else
                      stop "error"
                   endif
                   flag(i,j) = 1
                   
                else if ( surroundflag.eq.2 ) then
                   ! e.g.) :: data(i,j) = (i,j-1) + ( (i,j-1) - (i,j-2) )
                   !                    = 2 * (i,j-1) - (i,j-2)
                   if      ( ( flag_(i,j-1).eq.1 ).and.( flag_(i-1,j).eq.1 ) ) then
                      data(z_,i,j) = + 0.5d0 * ( &
                           &             + 2.d0*data_(z_,i,j-1) - data_(z_,i,j-2) &
                           &             + 2.d0*data_(z_,i-1,j) - data_(z_,i-2,j) )
                   else if ( ( flag_(i,j-1).eq.1 ).and.( flag_(i+1,j).eq.1 ) ) then
                      data(z_,i,j) = + 0.5d0 * ( &
                           &             + 2.d0*data_(z_,i,j-1) - data_(z_,i,j-2) &
                           &             + 2.d0*data_(z_,i+1,j) - data_(z_,i+2,j) )
                   else if ( ( flag_(i,j+1).eq.1 ).and.( flag_(i-1,j).eq.1 ) ) then
                      data(z_,i,j) = + 0.5d0 * ( &
                           &             + 2.d0*data_(z_,i,j+1) - data_(z_,i,j+2) &
                           &             + 2.d0*data_(z_,i-1,j) - data_(z_,i-2,j) )
                   else if ( ( flag_(i,j+1).eq.1 ).and.( flag_(i+1,j).eq.1 ) ) then
                      data(z_,i,j) = + 0.5d0 * ( &
                           &             + 2.d0*data_(z_,i,j+1) - data_(z_,i,j+2) &
                           &             + 2.d0*data_(z_,i+1,j) - data_(z_,i+2,j) )
                   else
                      stop "error"
                   endif
                   flag(i,j) = 1

                else if ( surroundflag.eq.3 ) then
                   ! e.g.) :: data(i,j) = 0.5*( avg + ext )
                   !       ::  avg      = 0.5*( (i-1,j)+(i+1,j) )
                   !       ::  ext      = 2.0*(i,j-1)-(i,j-2)
                   if      ( ( flag_(i,j-1).eq.1 ).and.( flag_(i-1,j).eq.1 ) &
                        &                         .and.( flag_(i,j+1).eq.1 ) ) then
                      avg          = 0.5d0 * ( data_(z_,i,j-1) + data_(z_,i,j+1) )
                      ext          = 2.0d0 * data_(z_,i-1,j) - data_(z_,i-2,j)
                      data(z_,i,j) = 0.5d0 * ( avg + ext )
                   else if ( ( flag_(i,j-1).eq.1 ).and.( flag_(i-1,j).eq.1 ) &
                        &                         .and.( flag_(i+1,j).eq.1 ) ) then
                      avg          = 0.5d0 * ( data_(z_,i-1,j) + data_(z_,i+1,j) )
                      ext          = 2.0d0 * data_(z_,i,j-1) - data_(z_,i,j-2)
                      data(z_,i,j) = 0.5d0 * ( avg + ext )
                   else if ( ( flag_(i,j-1).eq.1 ).and.( flag_(i,j+1).eq.1 ) &
                        &                         .and.( flag_(i+1,j).eq.1 ) ) then
                      avg          = 0.5d0 * ( data_(z_,i,j-1) + data_(z_,i,j+1) )
                      ext          = 2.0d0 * data_(z_,i+1,j) - data_(z_,i+2,j)
                      data(z_,i,j) = 0.5d0 * ( avg + ext )
                   else if ( ( flag_(i,j+1).eq.1 ).and.( flag_(i-1,j).eq.1 ) &
                        &                         .and.( flag_(i+1,j).eq.1 ) ) then
                      avg          = 0.5d0 * ( data_(z_,i-1,j) + data_(z_,i+1,j) )
                      ext          = 2.0d0 * data_(z_,i,j+1) - data_(z_,i,j+2)
                      data(z_,i,j) = 0.5d0 * ( avg + ext )
                   else
                      stop "error"
                   endif
                   flag(i,j) = 1

                   
                else if ( surroundflag.eq.4 ) then
                   ! e.g.) :: data(i,j) = average of 4-point
                   data(z_,i,j) = 0.25d0*( + data_(z_,i,j-1) + data_(z_,i,j+1) &
                        &                  + data_(z_,i-1,j) + data_(z_,i+1,j) )
                   flag(i,j) = 1
                   
                endif
                
             endif
             
          enddo
       enddo

       if ( count.eq.0 ) exit
          
    enddo
    
    return
  end subroutine extrapolate__linear2d


  
end program main
