PROGRAM polyfit
  !-----------------------------------------------------------!
  ! This is extrapolate_tau extended to:                      !
  ! - compute errorbars on parameters using a bootstrap Monte !
  !   Carlo method.                                           !
  ! - compute function values and derivatives with errorbars. !
  !                                                           !
  ! PLR 09.2015                                               !
  !-----------------------------------------------------------!
  IMPLICIT NONE

  ! Global variables.
  ! Precision for real-to-integer conversions and exponent comparisons.
  DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-8

  ! Local variables.
  ! (x,y,dy) data.
  CHARACTER(256) fname
  INTEGER nxy,ncolumn,icol_x,icol_y,icol_dx,icol_dy
  LOGICAL have_dx,have_dy,weighted
  DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:),dx(:),dy(:),weight(:)
  DOUBLE PRECISION xrange,xmin,xmax,xcentre,xmean,xmedian
  DOUBLE PRECISION dxrange,dxmin,dxmax,dxcentre,dxmean,dxmedian
  DOUBLE PRECISION yrange,ymin,ymax,ycentre,ymean,ymedian
  DOUBLE PRECISION dyrange,dymin,dymax,dycentre,dymean,dymedian
  ! Fitting function.
  INTEGER npoly
  DOUBLE PRECISION x0
  DOUBLE PRECISION,ALLOCATABLE :: pow(:),a(:),da(:)
  ! Polynomials resulting from operations on fitting polynomial.
  INTEGER op_npoly
  DOUBLE PRECISION,ALLOCATABLE :: op_pow(:),op_a(:)
  ! Automatic assessment variables.
  DOUBLE PRECISION,ALLOCATABLE ::  chi2_vector(:),rmsc_100_vector(:),&
     &rmsc_150_vector(:),rmsc_200_vector(:)
  INTEGER ntest
  ! Random sampling vectors.
  INTEGER nrandom,nx
  DOUBLE PRECISION grid_x1,grid_x2,f95_lo,f1s_lo,f1s_hi,f95_hi
  DOUBLE PRECISION,ALLOCATABLE :: ran_gauss_x(:),ran_gauss_y(:),ran_x(:),&
     &ran_y(:),ran_a(:),f_array(:,:),w_vector(:),x_target(:)
  ! Distribution analysis.
  DOUBLE PRECISION mean,var,skew,kurt
  ! I/O unit.
  INTEGER,PARAMETER :: io=10
  ! Misc local variables.
  INTEGER i,j,idiv,nderiv,irandom,ierr,ialloc,np,ix,ipos
  CHARACTER(2048) char2048
  CHARACTER(40) print_dble1,long_label,set_label,file_label
  DOUBLE PRECISION t1,xplot,dxplot

  write(6,*)'======='
  write(6,*)'POLYFIT'
  write(6,*)'======='
  write(6,*)

  ! Read in (x,y) data from file specified by user.
  call get_file(fname)
  open(unit=io,file=trim(fname),status='old',iostat=ierr)
  if(ierr/=0)call quit('Error opening "'//trim(fname)//'".')
  call check_file(io,fname,nxy,ncolumn)
  if(nxy<1)call quit('No data found in file "'//trim(fname)//'".')
  if(ncolumn<2)call quit('Not enough data columns found in file "'//&
     &trim(fname)//'".')
  write(6,*)'"'//trim(fname)//'" contains '//trim(i2s(nxy))//' lines with '//&
     &trim(i2s(ncolumn))//' columns.'
  write(6,*)
  call get_columns(ncolumn,icol_x,icol_y,icol_dx,icol_dy)
  have_dx=icol_dx>0
  have_dy=icol_dy>0
  allocate(x(nxy),y(nxy),dx(nxy),dy(nxy),weight(nxy),stat=ialloc)
  if(ialloc/=0)call quit('Allocation error.')
  call read_file(io,fname,ncolumn,icol_x,icol_y,icol_dx,icol_dy,nxy,x,y,dx,dy)
  close(io)

  ! Various stats.
  xmin=minval(x)
  xmax=maxval(x)
  xrange=xmax-xmin
  xcentre=.5d0*(xmin+xmax)
  xmean=sum(x)/dble(nxy)
  xmedian=median(nxy,x)
  ymin=minval(y)
  ymax=maxval(y)
  yrange=ymax-ymin
  ycentre=.5d0*(ymin+ymax)
  ymean=sum(y)/dble(nxy)
  ymedian=median(nxy,y)
  if(have_dx)then
    dxmin=minval(dx)
    dxmax=maxval(dx)
    dxrange=dxmax-dxmin
    dxcentre=.5d0*(dxmin+dxmax)
    dxmean=sum(dx)/dble(nxy)
    dxmedian=median(nxy,dx)
  endif ! have_dy
  if(have_dy)then
    dymin=minval(dy)
    dymax=maxval(dy)
    dyrange=dymax-dymin
    dycentre=.5d0*(dymin+dymax)
    dymean=sum(dy)/dble(nxy)
    dymedian=median(nxy,dy)
  endif ! have_dy
  write(6,*)
  write(6,*)'Data statistics:'
  write(6,*)
  write(6,'(t8,a,t21,a,t34,a,t47,a,t60,a)')'min','mean','median','centre','max'
  write(6,'(t2,70("-"))')
  write(6,'(t2,a,t5,5(1x,es12.4))')'X',xmin,xmean,xmedian,xcentre,xmax
  write(6,'(t2,a,t5,5(1x,es12.4))')'Y',ymin,ymean,ymedian,ycentre,ymax
  if(have_dx)write(6,'(t2,a,t5,5(1x,es12.4))')'DY',dxmin,dxmean,dxmedian,&
     &dxcentre,dxmax
  if(have_dy)write(6,'(t2,a,t5,5(1x,es12.4))')'DY',dymin,dymean,dymedian,&
     &dycentre,dymax
  write(6,'(t2,70("-"))')
  write(6,*)
  if(have_dx.and.have_dy)then
    weight=1.d0/(dx*dy)**2
  elseif(have_dx)then
    weight=1.d0/dx**2
  elseif(have_dy)then
    weight=1.d0/dy**2
  else
    weight=1.d0
  endif
  weighted=have_dx.or.have_dy

  restart_loop: do

    ! Get optional offset.
    do
      write(6,*)'X offset [number, left/mean/median/centre/min/max/right, or &
         &blank for 0.0]:'
      read(5,'(a)',iostat=ierr)char2048
      write(6,*)
      x0=0.d0
      if(ierr==0)then
        char2048=adjustl(char2048)
        select case(trim(char2048))
        case('left')
          x0=xmin
        case('right')
          x0=xmax
        case('min')
          do i=1,nxy
            if(abs(y(i)-ymin)<tol_zero)then
              x0=x(i)
              exit
            endif
          enddo
        case('max')
          do i=1,nxy
            if(abs(y(i)-ymax)<tol_zero)then
              x0=x(i)
              exit
            endif
          enddo
        case('centre')
          x0=xcentre
        case('mean')
          x0=xmean
        case('median')
          x0=xmedian
        case('')
          continue
        case default
          read(char2048,*,iostat=ierr)x0
          if(ierr<0)then
            call quit('Quitting.')
          elseif(ierr>0)then
            write(6,*)'Input problem, try again.'
            write(6,*)
            cycle
          endif
        end select
      endif
      exit
    enddo
    x=x-x0

    ! Assess integer expansion orders up to order 8.
    ntest=min(nxy,9)

    ! Allocate work arrays.
    allocate(chi2_vector(ntest),rmsc_100_vector(ntest),&
       &rmsc_150_vector(ntest),rmsc_200_vector(ntest),stat=ialloc)
    if(ialloc/=0)call quit('Allocation error.')

    ! Loop over expansion orders.
    do npoly=1,ntest
      ! Allocate work arrays.
      allocate(pow(npoly),a(npoly),da(npoly),op_pow((npoly*(npoly+1))/2),&
         &op_a((npoly*(npoly+1))/2),stat=ialloc)
      if(ialloc/=0)call quit('Allocation error.')
      ! Fit to polynomial of order npoly-1.
      do i=1,npoly
        pow(i)=dble(i-1)
      enddo ! i
      call calc_parameters(nxy,npoly,x,y,weight,pow,a,da,weighted)
      ! Evaluate chi^2.
      chi2_vector(npoly)=chi_squared(nxy,npoly,x,y,weight,pow,a,weighted)
      ! Evaluate root mean square curvature of fit as measure of smoothness.
      np=npoly
      call deriv_poly(np,pow,a,op_npoly,op_pow,op_a)
      call deriv_poly(op_npoly,op_pow,op_a)
      call square_poly(op_npoly,op_pow,op_a)
      call int_poly(op_npoly,op_pow,op_a)
      rmsc_100_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,x(nxy))-&
         &eval_poly(op_npoly,op_pow,op_a,x(1))&
         &)/xrange))
      rmsc_150_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,x(nxy)+0.25d0*xrange)-&
         &eval_poly(op_npoly,op_pow,op_a,x(1)-0.25d0*xrange)&
         &)/(1.5d0*xrange)))
      rmsc_200_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,x(nxy)+0.5d0*xrange)-&
         &eval_poly(op_npoly,op_pow,op_a,x(1)-0.5d0*xrange)&
         &)/(2.d0*xrange)))
      ! Destroy work arrays.
      deallocate(pow,a,da,op_pow,op_a)
    enddo ! npoly

    ! Report.
    write(6,*)
    write(6,*)'FIT ASSESSMENT'
    write(6,*)'=============='
    write(6,*)'The following table gives chi^2 values and RMS curvatures for &
       &various integer'
    write(6,*)'expansion orders to help you choose a good fitting function.  &
       &Chi^2 is a'
    write(6,*)'measure of fit accuracy (lower is better), and the RMS &
       &curvature is an'
    write(6,*)'indicator of overfitting (lower is better).  The RMS curvature &
       &is given for'
    write(6,*)'the original data interval and for symmetrically extended &
       &intervals, the'
    write(6,*)'latter being of particular relevance if the fit is to be used &
       &for'
    write(6,*)'extrapolations.'
    write(6,*)
    write(6,*)'Generally:'
    write(6,*)'* If chi^2 does not decrease monotonically, try a different x0.'
    write(6,*)'* If chi^2 has a plateau, the onset may be a reasonable &
       &expansion order.'
    write(6,*)'* Jumps in the RMS curvature with expansion order may indicate &
       &overfitting.'
    write(6,*)
    write(6,'(34x,a)')'RMS curvature'
    write(6,'(22x,36("-"))')
    write(6,'(1x,a,t10,a,t23,a,t36,a,t49,a)')'Order','chi^2','Data range',&
       &'+50%','+100%'
    write(6,'(1x,57("-"))')
    do npoly=1,ntest
      write(6,'(1x,i5,4(1x,es12.4),1x,a)')npoly-1,chi2_vector(npoly),&
         &rmsc_100_vector(npoly),rmsc_150_vector(npoly),rmsc_200_vector(npoly)
    enddo
    write(6,'(1x,57("-"))')
    write(6,*)
    write(6,*)'You can plot these polynomials by typing "plot" below.'
    write(6,*)

    ! Destroy work arrays.
    deallocate(chi2_vector,rmsc_100_vector,rmsc_150_vector,rmsc_200_vector)

    ! Ask user to enter exponents for fitting polynomial.
    parse_exp: do

      write(6,*)'Enter exponents in fitting polynomial as a space-separated &
         &list ["restart"'
      write(6,*)'to provide another X offset, "plot" to plot, blank to quit]:'
      read(5,'(a)')char2048
      write(6,*)

      ! Plot all polynomial orders.
      char2048=adjustl(char2048)
      select case(trim(char2048))
      case('plot')
        open(unit=io,file='poly_orders.dat',status='replace')
        do npoly=1,ntest
          allocate(pow(npoly),a(npoly),da(npoly),stat=ialloc)
          if(ialloc/=0)call quit('Allocation error.')
          do i=1,npoly
            pow(i)=dble(i-1)
          enddo ! i
          call calc_parameters(nxy,npoly,x,y,weight,pow,a,da,weighted)
          xplot=x(1)-0.5d0*xrange
          dxplot=(xrange*2.d0)/dble(1000)
          do i=0,1000
            write(io,*)xplot+x0,eval_poly(npoly,pow,a,xplot)
            xplot=xplot+dxplot
          enddo ! i
          deallocate(pow,a,da)
          write(io,'(a)')'&'
        enddo ! npoly
        close(io)
        write(6,*)'Data written to poly_orders.dat.'
        write(6,*)
        cycle parse_exp
      case('restart')
        cycle restart_loop
      end select

      ! Parse input string.
      npoly=0
      do
        read(char2048,*,iostat=ierr)(t1,i=1,npoly+1)
        if(ierr/=0)exit
        npoly=npoly+1
      enddo
      if(npoly<1)then
        write(6,*)'Interpolating polynomial must have at least one term.  &
           &Try again.'
        write(6,*)
        cycle parse_exp
      elseif(npoly>nxy)then
        write(6,*)'Number of terms cannot exceed number of points.  Try again.'
        write(6,*)
        cycle parse_exp
      endif ! npoly<1 or npoly>nxy
      allocate(pow(npoly),stat=ialloc)
      if(ialloc/=0)call quit('Allocation error.')
      read(char2048,*,iostat=ierr)pow(1:npoly)
      ! Make near-{,half,third,quarter}-integers exactly {*}-integers.
      do idiv=2,4
        do i=1,npoly
          if(abs(anint(pow(i)*idiv)-pow(i)*idiv)<tol_zero)&
             &pow(i)=anint(pow(i)*idiv)/dble(idiv)
        enddo ! i
      enddo ! idiv
      ! Forbid negative powers.
      if(any(pow<0.d0))then
        write(6,*)'Found negative exponent.  Try again.'
        write(6,*)
        deallocate(pow)
        cycle parse_exp
      endif
      ! Forbid repeated powers.
      do i=2,npoly
        if(any(abs(pow(1:i-1)-pow(i))<2.d0*tol_zero))then
          write(6,*)'Two exponents appear to be identical.  Try again.'
          write(6,*)
          deallocate(pow)
          cycle parse_exp
        endif ! Identical exponents
      enddo ! i
      exit parse_exp
    enddo parse_exp

    exit restart_loop
  enddo restart_loop

  ! Sort exponents into ascending order
  do i=1,npoly-1
    do j=i+1,npoly
      if(pow(j)<pow(i))then
        t1=pow(i)
        pow(i)=pow(j)
        pow(j)=t1
      endif ! Swap needed
    enddo ! j
  enddo ! i

  ! Allocate coefficient vectors.
  allocate(a(npoly),da(npoly),stat=ialloc)
  if(ialloc/=0)call quit('Allocation error (2).')

  ! Print polynomial form.
  write(6,*)
  write(6,*)'FITTING'
  write(6,*)'======='
  write(6,*)'Form of fitting polynomial:'
  write(6,*)'  '//trim(print_poly_sym(npoly,pow,x0))
  write(6,*)

  ! Compute fit.
  call calc_parameters(nxy,npoly,x,y,weight,pow,a,da,weighted)

  ! Print fit coefficients.
  write(6,*)'Fit parameters:'
  if(weighted)then
    do j=1,npoly
      write(6,'("   ",a," = ",es24.16," +/- ",es24.16)')'k'//trim(i2s(j)),&
         &a(j),da(j)
    enddo ! j
  else
    do j=1,npoly
      write(6,'("   ",a," = ",es24.16)')'k'//trim(i2s(j)),a(j)
    enddo ! j
  endif ! weighted
  ! Evaluate chi-squared.
  t1=chi_squared(nxy,npoly,x,y,weight,pow,a,weighted)
  if(weighted)then
    write(6,*)'Chi-squared value: ',t1
  else
    write(6,*)'Least-squares function: ',t1
  endif ! weighted
  write(6,*)

  ! Write out fitted polynomial.
  ! NB, this is good for pasting into xmgrace, but xmgrace has a string
  ! length limit of 256 characters.
  write(6,*)'Fitted polynomial:'
  write(6,*)'  '//trim(print_poly_num(npoly,pow,a,x0))
  write(6,*)

  ! Write heading for operations part.
  write(6,*)
  write(6,*)'VALUES AND DERIVATIVES'
  write(6,*)'======================'
  write(6,*)'Now you can calculate specific values of the fitted polynomial &
     &and its'
  write(6,*)'first two derivatives, with error bars calculated using a Monte &
     &Carlo method.'
  write(6,*)

  ! Ask for number of random points.
  nrandom=0
  if(have_dx.or.have_dy)then
    do
      write(6,*)'Enter number of random points for Monte Carlo [blank for &
         &10000]:'
      read(5,'(a)',iostat=ierr)char2048
      write(6,*)
      if(ierr<0)then
        call quit()
      elseif(ierr>0)then
        write(6,*)'Input problem, try again.'
        cycle
      endif
      char2048=adjustl(char2048)
      if(trim(char2048)=='')then
        nrandom=10000
        exit
      endif
      read(char2048,*,iostat=ierr)nrandom
      if(ierr/=0)then
        write(6,*)'Bad input, try again.'
        cycle
      endif
      if(nrandom<10)then
        write(6,*)'Number of random points must be => 10.'
        cycle
      endif
      exit
    enddo
  endif

  ! Allocate work arrays.
  allocate(op_pow(npoly),op_a(npoly),stat=ialloc)
  if(ialloc/=0)call quit('Allocation problem (op_*).')
  if(have_dx.or.have_dy)then
    allocate(ran_gauss_x(nxy),ran_gauss_y(nxy),ran_x(nxy),ran_y(nxy),&
       &ran_a(npoly),w_vector(nrandom),stat=ialloc)
    if(ialloc/=0)call quit('Allocation problem (w_vector).')
  endif

  ! Loop over sets of operations.
  do nderiv=0,2

    select case(nderiv)
    case(0)
      long_label='Values'
      set_label="y"
      file_label='values'
    case(1)
      long_label='First derivatives'
      set_label="y'"
      file_label='fderiv'
    case(2)
      long_label='Second derivatives'
      set_label="y''"
      file_label='sderiv'
    end select
    write(6,*)trim(long_label)
    write(6,'(1x,'//trim(i2s(len_trim(long_label)))//'("-"))')
    write(6,*)'Enter space-separated list of X values [blank for X of source &
       &data,'
    write(6,*)'<x1>:<x2>:<n> a grid, "skip" to skip]:'

    do
      read(5,'(a)')char2048
      write(6,*)

      char2048=adjustl(char2048)
      select case(trim(char2048))
      case('')
        write(char2048,*)x(1)+x0
        char2048=adjustl(char2048)
        do i=2,nxy
          write(print_dble1,*)x(i)+x0
          char2048=trim(char2048)//' '//trim(adjustl(print_dble1))
        enddo ! i
        nx=nxy
      case('skip','-','x')
        exit
      case default
        ! Parse input string.
        nx=0
        do
          read(char2048,*,iostat=ierr)(t1,i=1,nx+1)
          if(ierr/=0)exit
          nx=nx+1
        enddo
        if(nx==0)then
          ! Try to parse grid.
          ipos=scan(char2048,':')
          if(ipos<1)then
            write(6,*)'Input problem - try again.'
            write(6,*)
            cycle
          endif
          read(char2048(1:ipos-1),*,iostat=ierr)grid_x1
          if(ierr/=0)then
            write(6,*)'Input problem - try again.'
            write(6,*)
            cycle
          endif
          char2048=adjustl(char2048(ipos+1:))
          ipos=scan(char2048,':')
          if(ipos<1)then
            write(6,*)'Input problem - try again.'
            write(6,*)
            cycle
          endif
          read(char2048(1:ipos-1),*,iostat=ierr)grid_x2
          if(ierr/=0)then
            write(6,*)'Input problem - try again.'
            write(6,*)
            cycle
          endif
          char2048=adjustl(char2048(ipos+1:))
          read(char2048,*,iostat=ierr)nx
          if(ierr/=0)then
            write(6,*)'Input problem - try again.'
            write(6,*)
            cycle
          endif
          if(nx<1)then
            write(6,*)'Need at least one point in grid.  Try again.'
            write(6,*)
            cycle
          endif
          char2048=''
        endif
      end select

      exit
    enddo

    ! Allocate x array.
    allocate(x_target(nx),stat=ialloc)
    if(ialloc/=0)call quit('Allocation problem (f_array).')
    if(len_trim(char2048)==0)then
      if(nx==1)then
        x_target(1)=grid_x1
      else
        do i=1,nx
          x_target(i)=grid_x1+(grid_x2-grid_x1)*dble(i-1)/dble(nx-1)
        enddo ! i
      endif
    else
      read(char2048,*,iostat=ierr)x_target(1:nx)
    endif

    if(have_dx.or.have_dy)then

      ! Allocate work array.
      allocate(f_array(nrandom,nx),stat=ialloc)
      if(ialloc/=0)call quit('Allocation problem (f_array).')

      ! Do random sampling of data space.
      ran_x=x
      ran_y=y
      do irandom=1,nrandom
        if(have_dx)then
          ran_gauss_x=gaussian_random_number(dx)
          ran_x=x+ran_gauss_x
        endif
        if(have_dy)then
          ran_gauss_y=gaussian_random_number(dy)
          ran_y=y+ran_gauss_y
        endif
        call calc_parameters(nxy,npoly,ran_x,ran_y,weight,pow,ran_a,da,.false.)
        !w_vector(irandom)=1.d0/&
        !   &chi_squared(nxy,npoly,x,ran_y,weight,pow,ran_a,.false.)
        w_vector(irandom)=1.d0
        op_npoly=npoly
        op_pow=pow
        op_a=ran_a
        do i=1,nderiv
          call deriv_poly(op_npoly,op_pow,op_a)
        enddo ! i
        do ix=1,nx
          t1=eval_poly(op_npoly,op_pow,op_a,x_target(ix)-x0)
          f_array(irandom,ix)=t1
        enddo ! ix
      enddo ! irandom

    endif ! have_dx.or.have_dy

    ! Open file.
    fname='polyfit_'//trim(file_label)//'.dat'
    open(unit=io,file=trim(fname),status='replace')

    if(.not.(have_dx.or.have_dy))then

      ! Write header.
      write(io,'(a)')'# x  '//trim(set_label)
      ! Construct target polynomial at original parameter set.
      op_npoly=npoly
      op_pow=pow
      op_a=a
      do i=1,nderiv
        call deriv_poly(op_npoly,op_pow,op_a)
      enddo ! i
      ! Loop over values of x.
      do ix=1,nx
        ! Evaluate target polynomial at x.
        t1=eval_poly(op_npoly,op_pow,op_a,x_target(ix)-x0)
        write(io,*)x_target(ix),t1
      enddo ! ix

    else

      ! Write header.
      write(io,'(a)')'# x  '//trim(set_label)//&
         &'  stderr  stderr+  stderr-  95%err+  95%err-  skew  kurt'
      ! Loop over values of x.
      do ix=1,nx
        ! Report various statistics.
        mean=sum(f_array(:,ix)*w_vector)/sum(w_vector)
        call characterize_dist(nrandom,f_array(:,ix),w_vector,var,skew,&
           &kurt)
        f95_lo=find_pth_smallest(nint(dble(nrandom)*0.025d0),nrandom,&
           &f_array(:,ix))
        f1s_lo=find_pth_smallest(nint(dble(nrandom)*0.158655254d0),nrandom,&
           &f_array(:,ix))
        f1s_hi=find_pth_smallest(nint(dble(nrandom)*0.841344746d0),nrandom,&
           &f_array(:,ix))
        f95_hi=find_pth_smallest(nint(dble(nrandom)*0.975d0),nrandom,&
           &f_array(:,ix))
        write(io,*)x_target(ix),mean,sqrt(var),mean-f95_lo,mean-f1s_lo,&
           &f1s_hi-mean,f95_hi-mean,skew,kurt
      enddo ! ix

    endif ! .not.(have_dx.or.have_dy) or have_dx.or.have_dy

    ! Close file.
    close(io)
    write(6,*)'Data written to "'//trim(fname)//'".'
    write(6,*)

    ! Clean up.
    if(have_dx.or.have_dy)deallocate(f_array)
    deallocate(x_target)

  enddo ! loop over y / y' / y''

  ! Clean up.
  if(have_dx.or.have_dy)deallocate(ran_gauss_x,ran_gauss_y,ran_x,ran_y,ran_a,&
     &w_vector)
  deallocate(op_pow,op_a)

  ! Finish off.
  deallocate(x,y,dx,dy,weight,pow,a,da)
  write(6,*)'Program finished.'
  write(6,*)


CONTAINS


  SUBROUTINE calc_parameters(nxy,npoly,x,y,weight,pow,a,da,weighted)
    !-------------------------------------------------------------!
    ! Evaluate the parameters using the method given in my notes. !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),weight(nxy),pow(npoly)
    DOUBLE PRECISION,INTENT(out) :: a(npoly),da(npoly)
    LOGICAL,INTENT(in) :: weighted
    INTEGER info,i,j,ipiv(npoly),lwork,ialloc
    DOUBLE PRECISION M(npoly,npoly),Minv(npoly,npoly),c(npoly),tempr(1)
    DOUBLE PRECISION,ALLOCATABLE :: work(:)
    INTERFACE
      SUBROUTINE dsytrf(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
        CHARACTER(1),INTENT(in) :: UPLO
        INTEGER,INTENT(in) :: N,LDA,LWORK
        DOUBLE PRECISION,INTENT(out) :: A(LDA,*),WORK(*)
        INTEGER,INTENT(out) :: IPIV(*)
        INTEGER,INTENT(out) :: INFO
      END SUBROUTINE dsytrf
      SUBROUTINE dsytri(UPLO,N,A,LDA,IPIV,WORK,INFO)
        CHARACTER(1),INTENT(in) :: UPLO
        INTEGER,INTENT(in) :: N,LDA,IPIV(*)
        DOUBLE PRECISION,INTENT(inout) :: A(LDA,*)
        DOUBLE PRECISION,INTENT(out) :: WORK(*)
        INTEGER,INTENT(out) :: INFO
      END SUBROUTINE dsytri
    END INTERFACE

    ! Construct vector c.
    call construct_c(nxy,npoly,x,y,weight,pow,c,weighted)

    ! Construct matrix M.
    call construct_M(nxy,npoly,x,weight,pow,M,weighted)

    ! Invert matrix M.
    Minv=M
    call dsytrf('L',npoly,Minv(1,1),npoly,ipiv(1),tempr(1),-1,info)
    if(info/=0)call quit('Matrix inversion failed (1).')
    lwork=nint(tempr(1))
    allocate(work(lwork),stat=ialloc)
    if(ialloc/=0)call quit('Allocation error: WORK (1).')
    call dsytrf('L',npoly,Minv(1,1),npoly,ipiv(1),work(1),lwork,info)
    if(info/=0)call quit('Matrix inversion failed (2).')
    deallocate(work)
    allocate(work(npoly),stat=ialloc)
    if(ialloc/=0)call quit('Allocation error: WORK (2).')
    call dsytri('L',npoly,Minv(1,1),npoly,ipiv(1),work(1),info)
    if(info/=0)call quit('Matrix inversion failed (3).')
    deallocate(work)

    ! Complete the upper triangle of Minv.
    do i=1,npoly-1
      do j=i+1,npoly
        Minv(i,j)=Minv(j,i)
      enddo ! j
    enddo ! i

    ! Hence evaluate the coefficients of the terms in the polynomial.
    a=matmul(Minv,c)

    ! Evaluate the standard errors in the coefficients.
    if(weighted)then
      do j=1,npoly
        da(j)=sqrt(Minv(j,j))
      enddo ! j
    endif ! weighted

  END SUBROUTINE calc_parameters


  DOUBLE PRECISION FUNCTION chi_squared(nxy,npoly,x,y,weight,pow,a,weighted)
  !--------------------------------------------------------------------------!
  ! Evaluate the chi-squared value of the fit.  If error bars are not given, !
  ! return the least-squares function.                                       !
  !--------------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),weight(nxy),a(npoly),&
       &pow(npoly)
    LOGICAL,INTENT(in) :: weighted
    INTEGER i,k
    DOUBLE PRECISION e_fit
    chi_squared=0.d0
    do i=1,nxy
      e_fit=0.d0
      do k=1,npoly
        e_fit=e_fit+a(k)*x(i)**pow(k)
      enddo ! k
      if(weighted)then
        chi_squared=chi_squared+(y(i)-e_fit)**2*weight(i)
      else
        chi_squared=chi_squared+(y(i)-e_fit)**2
      endif ! weighted
    enddo ! i
  END FUNCTION chi_squared


  SUBROUTINE construct_c(nxy,npoly,x,y,weight,pow,c,weighted)
    !----------------------------------------!
    ! Construct the vector c (see my notes). !
    !----------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),weight(nxy),pow(npoly)
    DOUBLE PRECISION,INTENT(out) :: c(npoly)
    LOGICAL,INTENT(in) :: weighted
    INTEGER i,j
    if(weighted)then
      do j=1,npoly
        c(j)=0.d0
        do i=1,nxy
          c(j)=c(j)+y(i)*x(i)**pow(j)*weight(i)
        enddo ! i
      enddo ! j
    else
      do j=1,npoly
        c(j)=0.d0
        do i=1,nxy
          c(j)=c(j)+y(i)*x(i)**pow(j)
        enddo ! i
      enddo ! j
    endif ! weighted
  END SUBROUTINE construct_c


  SUBROUTINE construct_M(nxy,npoly,x,weight,pow,M,weighted)
    !----------------------------------------!
    ! Construct the matrix M (see my notes). !
    ! Only need lower triangular part.       !
    !----------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),weight(nxy),pow(npoly)
    DOUBLE PRECISION,INTENT(out) :: M(npoly,npoly)
    LOGICAL,INTENT(in) :: weighted
    INTEGER i,j,k
    M=0.d0
    if(weighted)then
      do k=1,npoly
        do j=k,npoly
          do i=1,nxy
            M(j,k)=M(j,k)+x(i)**(pow(j)+pow(k))*weight(i)
          enddo ! i
        enddo ! j
      enddo ! k
    else
      do k=1,npoly
        do j=k,npoly
          do i=1,nxy
            M(j,k)=M(j,k)+x(i)**(pow(j)+pow(k))
          enddo ! i
        enddo ! j
      enddo ! k
    endif ! weighted
  END SUBROUTINE construct_M


  CHARACTER(12) FUNCTION i2s(n)
    !-----------------------------------------------------------------------!
    ! Convert integers to left justified strings that can be printed in the !
    ! middle of a sentence without introducing large amounts of white space.!
    !-----------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    INTEGER i,j
    INTEGER,PARAMETER :: ichar0=ichar('0')
    i2s=''
    i=abs(n)
    do j=len(i2s),1,-1
      i2s(j:j)=achar(ichar0+mod(i,10))
      i=i/10
      if(i==0)exit
    enddo ! j
    if(n<0)then
      i2s='-'//adjustl(i2s)
    else
      i2s=adjustl(i2s)
    endif ! n<0
  END FUNCTION i2s


  SUBROUTINE get_file(fname)
    !------------------------!
    ! Find out the filename. !
    !------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(out) :: fname
    INTEGER ierr
    LOGICAL file_exists
    do
      write(6,*)'Enter name of data file:'
      read(5,*,iostat=ierr)fname
      if(ierr/=0)then
        write(6,*)'Error reading file name.  Try again.'
        cycle
      endif
      fname=adjustl(fname)
      inquire(file=trim(fname),exist=file_exists)
      if(.not.file_exists)then
        write(6,*)'File does not appear to exist. Please try again.'
        cycle
      endif ! file non-existent
      exit
    enddo
    write(6,*)
  END SUBROUTINE get_file


  SUBROUTINE check_file(io,fname,nline,ncolumn)
    !----------------------------------------------------!
    ! Check file contains data, and return the number of !
    ! data lines and the number of data columns in it.   !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: io
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(out) :: nline,ncolumn
    CHARACTER(1024) char1024
    INTEGER pline,ierr
    DOUBLE PRECISION t1

    ! Initialize.
    nline=0
    ncolumn=0
    pline=0

    ! Look for first line containing data.
    do
      read(io,'(a)',iostat=ierr)char1024
      if(ierr<0)exit
      pline=pline+1
      if(ierr>0)call quit('Error reading '//trim(fname)//' at line '//&
         &trim(i2s(pline))//'.')
      char1024=adjustl(char1024)
      ! Skip comments.
      ipos=scan(char1024,'#!')
      if(ipos==1)cycle
      if(ipos>1)char1024=char1024(1:ipos-1)
      ! Skip empty lines.
      if(len_trim(char1024)==0)cycle
      if(ncolumn==0)then
        ! Find how many elements there are in this line.
        do
          read(char1024,*,iostat=ierr)(t1,i=1,ncolumn+1)
          if(ierr/=0)exit
          ncolumn=ncolumn+1
        enddo
      else
        ! Ensure we can read ncolumn elements in this line.
        do
          read(char1024,*,iostat=ierr)(t1,i=1,ncolumn)
          if(ierr==0)exit
          ncolumn=ncolumn-1
        enddo
      endif
      if(ncolumn<2)call quit('Too few columns found at line '//&
         &trim(i2s(pline))//' of file "'//trim(fname)//'".')
      nline=nline+1
    enddo

    ! Rewind file.
    rewind(io)

  END SUBROUTINE check_file


  SUBROUTINE get_columns(ncolumn,icol_x,icol_y,icol_dx,icol_dy)
    !----------------------------------------!
    ! Ask the user for which columns to use. !
    !----------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ncolumn
    INTEGER, INTENT(inout) :: icol_x, icol_y, icol_dx, icol_dy
    CHARACTER(1024) char1024,cdata,ccol
    INTEGER ierr

    ! Promt and read use input.
    write(6,*)'Enter the indices for x, y, dx, dy:'
    read(5,'(a)',iostat=ierr)char1024
    if(ierr/=0)call quit()
    write(6,*)

    ! Parse input.
    if(len_trim(char1024)==0)then
      ! Use default if empty input.
      select case(ncolumn)
      case(1)
        icol_x=0
        icol_y=1
        icol_dx=0
        icol_dy=1
      case(2)
        icol_x=1
        icol_y=2
        icol_dx=0
        icol_dy=0
      case(3)
        icol_x=1
        icol_y=2
        icol_dx=0
        icol_dy=3
      case default
        icol_x=1
        icol_y=3
        icol_dx=2
        icol_dy=4
      end select
    else
      ! Adapt parsing to number of fields in input string.
      read(char1024,*,iostat=ierr)icol_x,icol_y,icol_dx,icol_dy
      if(ierr/=0)then
        icol_dx=0
        read(char1024,*,iostat=ierr)icol_x,icol_y,icol_dy
        if(ierr/=0)then
          icol_dy=0
          read(char1024,*,iostat=ierr)icol_x,icol_y
          if(ierr/=0)then
            icol_x=0
            read(char1024,*,iostat=ierr)icol_y
            if(ierr/=0)call quit('Could not parse column indices.')
          endif
        endif
      endif
      if(icol_x>ncolumn.or.icol_y>ncolumn.or.icol_dx>ncolumn.or.&
         &icol_dy>ncolumn)call quit('Column indices out of range.')
    endif

    ! Report.
    cdata='(x,y'
    ccol='('
    if(icol_x==0)then
      ccol=trim(ccol)//'INDEX'
    else
      ccol=trim(ccol)//trim(i2s(icol_x))
    endif
    ccol=trim(ccol)//','
    if(icol_y==0)then
      ccol=trim(ccol)//'INDEX'
    else
      ccol=trim(ccol)//trim(i2s(icol_y))
    endif
    if(icol_dx>0)then
      cdata=trim(cdata)//',dx'
      ccol=trim(ccol)//','//trim(i2s(icol_dx))
    endif
    if(icol_dy>0)then
      cdata=trim(cdata)//',dy'
      ccol=trim(ccol)//','//trim(i2s(icol_dy))
    endif
    cdata=trim(cdata)//')'
    ccol=trim(ccol)//')'
    write(6,*)'Parsing '//trim(cdata)//' from columns '//trim(ccol)//'.'

  END SUBROUTINE get_columns


  SUBROUTINE read_file(io,fname,ncolumn,icol_x,icol_y,icol_dx,icol_dy,nxy,&
     &x,y,dx,dy)
    !-------------------------------------!
    ! Read in the data in the input file. !
    !-------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: io,ncolumn,icol_x,icol_y,icol_dx,icol_dy,nxy
    CHARACTER(*),INTENT(in) :: fname
    DOUBLE PRECISION,INTENT(out) :: x(nxy),y(nxy),dx(nxy),dy(nxy)
    DOUBLE PRECISION tvec(ncolumn)
    INTEGER i,ierr,ipos,pline
    CHARACTER(1024) line

    ! Initialize.
    pline=0
    x=(/(dble(i),i=1,nxy)/)
    y=(/(dble(i),i=1,nxy)/)
    dx=0.d0
    dy=0.d0

    ! Loop over data lines.
    do i=1,nxy
      ! Load next line containing data.
      do
        read(io,'(a)',iostat=ierr)line
        if(ierr<0)call quit('Error reading '//trim(fname)//' at line '//&
           &trim(i2s(pline))//': unexpected end of file.')
        pline=pline+1
        if(ierr>0)call quit('Error reading '//trim(fname)//' at line '//&
           &trim(i2s(pline))//'.')
        line=adjustl(line)
        ! Skip comments.
        ipos=scan(line,'#!')
        if(ipos==1)cycle
        if(ipos>1)line=line(1:ipos-1)
        ! Skip empty lines.
        if(len_trim(line)==0)cycle
        exit
      enddo
      ! Read data point from string.
      read(line,*,iostat=ierr)tvec(1:ncolumn)
      if(ierr/=0)call quit('Error reading '//trim(fname)//' at line '//&
         &trim(i2s(pline))//'.')
      ! Set data and check that error bars are positive.
      if(icol_x>0)x(i)=tvec(icol_x)
      if(icol_y>0)y(i)=tvec(icol_y)
      if(icol_dx>0)then
        dx(i)=tvec(icol_dx)
        if(dx(i)<=0.d0)call quit('Error reading '//trim(fname)//' at line '//&
           &trim(i2s(pline))//': non-positive dx.')
      endif ! have_dx
      if(icol_dy>0)then
        dy(i)=tvec(icol_dy)
        if(dy(i)<=0.d0)call quit('Error reading '//trim(fname)//' at line '//&
           &trim(i2s(pline))//': non-positive dy.')
      endif ! have_dy
    enddo ! i

  END SUBROUTINE read_file


  FUNCTION print_poly_sym(npoly,pow,x0) RESULT(polystr)
  !------------------------------------------------------------------------!
  ! This subroutine returns the fitted polynomial in a suitable format for !
  ! pasting into xmgrace.                                                  !
  !------------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: npoly
    DOUBLE PRECISION,INTENT(in) :: pow(npoly),x0
    CHARACTER(3+npoly*54) :: polystr
    INTEGER j,ipow
    CHARACTER(40) pwstr
    CHARACTER(6) xstr
    if(abs(x0)<tol_zero)then
      xstr='x'
    else
      xstr='(x-x0)'
    endif
    polystr='y ='
    do j=1,npoly
      ipow=nint(pow(j))
      if(abs(dble(ipow)-pow(j))<tol_zero)then
        if(ipow==0)then
          pwstr=''
        elseif(ipow==1)then
          pwstr='*'//trim(xstr)
        else
          pwstr='*'//trim(xstr)//'^'//trim(i2s(ipow))
        endif ! pow=0
      else
        write(pwstr,*)pow(j)
        pwstr='*'//trim(xstr)//'^'//trim(adjustl(pwstr))
      endif ! Power is an integer
      if(j>1)then
        polystr=trim(polystr)//' + k'//trim(i2s(j))//trim(pwstr)
      else
        polystr=trim(polystr)//' k'//trim(i2s(j))//trim(pwstr)
      endif ! "+" needed
    enddo ! j
    ! Change "d" to "e" in polystr.
    do j=1,len_trim(polystr)
      if(polystr(j:j)=='d'.or.polystr(j:j)=='D')polystr(j:j)='E'
    enddo ! j
  END FUNCTION print_poly_sym


  FUNCTION print_poly_num(npoly,pow,a,x0) RESULT(polystr)
  !------------------------------------------------------------------------!
  ! This subroutine returns the fitted polynomial in a suitable format for !
  ! pasting into xmgrace.                                                  !
  !------------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: npoly
    DOUBLE PRECISION,INTENT(in) :: pow(npoly),a(npoly),x0
    CHARACTER(2+npoly*105) :: polystr
    INTEGER j,ipow
    CHARACTER(1) plusstr
    CHARACTER(32) coeffstr
    CHARACTER(72) pwstr
    CHARACTER(36) xstr
    if(abs(x0)<tol_zero)then
      xstr='x'
    else
      write(xstr,*)x0
      xstr='(x-'//trim(adjustl(xstr))//')'
    endif
    polystr='y='
    do j=1,npoly
      ipow=nint(pow(j))
      if(abs(dble(ipow)-pow(j))<tol_zero)then
        if(ipow==0)then
          pwstr=''
        elseif(ipow==1)then
          pwstr='*'//trim(xstr)
        else
          pwstr='*'//trim(xstr)//'^'//trim(i2s(ipow))
        endif ! pow=0
      else
        write(pwstr,*)pow(j)
        pwstr='*'//trim(xstr)//'^'//trim(adjustl(pwstr))
      endif ! Power is an integer
      if(j>1.and.a(j)>=0.d0)then
        plusstr='+'
      else
        plusstr=''
      endif ! "+" needed
      write(coeffstr,*)a(j)
      coeffstr=adjustl(coeffstr)
      polystr=trim(polystr)//trim(plusstr)//trim(coeffstr)//trim(pwstr)
    enddo ! j
    ! Change "d" to "e" in polystr.
    do j=1,len_trim(polystr)
      if(polystr(j:j)=='d'.or.polystr(j:j)=='D')polystr(j:j)='E'
    enddo ! j
  END FUNCTION print_poly_num


  SUBROUTINE parse_expression(expr,nderiv,x_target,x1_target,ierr)
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: expr
    INTEGER,INTENT(out) :: nderiv,ierr
    DOUBLE PRECISION,INTENT(out) :: x_target,x1_target
    COMPLEX(kind(1.d0)) ct1
    CHARACTER(len(expr)) tstr
    ierr=1
    nderiv=0
    tstr=adjustl(expr)
    if(tstr(1:1)=='Y')then
      nderiv=-1
      tstr=adjustl(tstr(2:))
      if(tstr(1:1)/='(')return
      if(tstr(len_trim(tstr):len_trim(tstr))/=')')return
      read(tstr,*,iostat=ierr)ct1
      if(ierr/=0)return
      x_target=dble(ct1)
      x1_target=aimag(ct1)
    else
      if(tstr(1:1)/='y')return
      tstr=adjustl(tstr(2:))
      do while(tstr(1:1)=="'")
        nderiv=nderiv+1
        tstr=adjustl(tstr(2:))
      enddo
      if(tstr(1:1)/='(')return
      tstr=adjustl(tstr(2:))
      if(tstr(len_trim(tstr):len_trim(tstr))/=')')return
      tstr=tstr(1:len_trim(tstr)-1)
      read(tstr,*,iostat=ierr)x_target
    endif
  END SUBROUTINE parse_expression


  SUBROUTINE deriv_poly(npoly1,pow1,a1,npoly,pow,a)
    !-------------------------------------------------!
    ! Differentiate the polynomial given by npoly1,   !
    ! pow1 and a1 and return the resulting polynomial !
    ! in npoly2, pow2 and a2 if present, or in-place  !
    ! overwriting the input if not.                   !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: npoly1
    DOUBLE PRECISION,INTENT(inout) :: pow1(:),a1(:)
    INTEGER,INTENT(out),OPTIONAL :: npoly
    DOUBLE PRECISION,INTENT(out),OPTIONAL :: pow(:),a(:)
    INTEGER npoly2
    DOUBLE PRECISION pow2(npoly1),a2(npoly1)
    INTEGER ipoly
    npoly2=0
    do ipoly=1,npoly1
      if(abs(pow1(ipoly))<tol_zero)then
        ! Delete constant term from polynomial.
        cycle
      else
        ! Evaluate derivative of non-constant term.
        npoly2=npoly2+1
        a2(npoly2)=a1(ipoly)*pow1(ipoly)
        pow2(npoly2)=pow1(ipoly)-1.d0
      endif
    enddo ! ipoly
    if(present(npoly))then
      ! Return result separately.
      npoly=npoly2
      pow(1:npoly)=pow2(1:npoly2)
      a(1:npoly)=a2(1:npoly)
    else
      ! Return result in place.
      npoly1=npoly2
      pow1(1:npoly1)=pow2(1:npoly2)
      a1(1:npoly1)=a2(1:npoly2)
    endif
  END SUBROUTINE deriv_poly


  SUBROUTINE int_poly(npoly1,pow1,a1,npoly,pow,a)
    !------------------------------------------------!
    ! Integrate the polynomial given by npoly1, pow1 !
    ! and a1 and return the resulting polynomial in  !
    ! npoly2, pow2 and a2 if present, or in-place    !
    ! overwriting the input if not.                  !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: npoly1
    DOUBLE PRECISION,INTENT(inout) :: pow1(:),a1(:)
    INTEGER,INTENT(out),OPTIONAL :: npoly
    DOUBLE PRECISION,INTENT(out),OPTIONAL :: pow(:),a(:)
    INTEGER npoly2
    DOUBLE PRECISION pow2(npoly1),a2(npoly1)
    INTEGER ipoly
    npoly2=npoly1
    do ipoly=1,npoly1
      a2(ipoly)=a1(ipoly)/(pow1(ipoly)+1.d0)
      pow2(ipoly)=pow1(ipoly)+1.d0
    enddo ! ipoly
    if(present(npoly))then
      ! Return result separately.
      npoly=npoly2
      pow(1:npoly)=pow2(1:npoly2)
      a(1:npoly)=a2(1:npoly)
    else
      ! Return result in place.
      npoly1=npoly2
      pow1(1:npoly1)=pow2(1:npoly2)
      a1(1:npoly1)=a2(1:npoly2)
    endif
  END SUBROUTINE int_poly


  SUBROUTINE square_poly(npoly1,pow1,a1,npoly,pow,a)
    !---------------------------------------------!
    ! Square the polynomial given by npoly1, pow1 !
    ! and a1 and return the resulting polynomial  !
    ! in npoly2, pow2 and a2 if present, or       !
    ! in-place overwriting the input if not.      !
    !---------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: npoly1
    DOUBLE PRECISION,INTENT(inout) :: pow1(:),a1(:)
    INTEGER,INTENT(out),OPTIONAL :: npoly
    DOUBLE PRECISION,INTENT(out),OPTIONAL :: pow(:),a(:)
    INTEGER npoly2
    DOUBLE PRECISION pow2((npoly1*(npoly1+1))/2),a2((npoly1*(npoly1+1))/2)
    INTEGER ipoly,jpoly
    npoly2=0
    do ipoly=1,npoly1
      do jpoly=ipoly,npoly1
        npoly2=npoly2+1
        a2(npoly2)=a1(ipoly)*a1(jpoly)
        if(ipoly/=jpoly)a2(npoly2)=a2(npoly2)*2.d0
        pow2(npoly2)=pow1(ipoly)+pow1(jpoly)
      enddo
    enddo ! ipoly
    if(present(npoly))then
      ! Return result separately.
      npoly=npoly2
      pow(1:npoly)=pow2(1:npoly2)
      a(1:npoly)=a2(1:npoly)
    else
      ! Return result in place.
      npoly1=npoly2
      pow1(1:npoly1)=pow2(1:npoly2)
      a1(1:npoly1)=a2(1:npoly2)
    endif
  END SUBROUTINE square_poly


  DOUBLE PRECISION FUNCTION eval_poly(npoly,pow,a,x_target)
    !-------------------------------------------------!
    ! Evaluate the polynomial given by npoly, pow and !
    ! a at x_target.                                  !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: npoly
    DOUBLE PRECISION,INTENT(in) :: pow(npoly),a(npoly),x_target
    INTEGER ipoly
    eval_poly=0.d0
    do ipoly=1,npoly
      if(abs(pow(ipoly))<tol_zero)then
        ! Constant term.
        eval_poly=eval_poly+a(ipoly)
      elseif(abs(x_target)<tol_zero)then
        ! Skip evaluation if x is zero (avoids issues with negative powers).
        cycle
      elseif(abs(anint(pow(ipoly))-pow(ipoly))<tol_zero)then
        ! Natural powers force to use integer exponents to avoid issues.
        eval_poly=eval_poly+a(ipoly)*x_target**nint(pow(ipoly))
      elseif(x_target>0.d0)then
        ! Fractional powers only evaluated for positive arguments.
        eval_poly=eval_poly+a(ipoly)*x_target**pow(ipoly)
      endif
    enddo ! ipoly
  END FUNCTION eval_poly


  FUNCTION gaussian_random_number(w) RESULT(gauss_vector)
    !--------------------------------------------------------------!
    ! Given w(1:n), returns n random numbers distributed according !
    ! to a Gaussian of variance w**2 centred at zero.              !
    !--------------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: w(:)
    DOUBLE PRECISION :: gauss_vector(size(w))
    DOUBLE PRECISION v(2),rad,fac
    INTEGER i,n
    n=size(w)
    do i=1,n,2
      do
        call random_number(v)
        v(:)=2.d0*v(:)-1.d0
        rad=sum(v**2)
        if(rad<=1.d0.and.rad>0.d0)exit
      enddo
      fac=log(rad)/rad
      fac=sqrt(-fac-fac)
      gauss_vector(i)=v(2)*fac*w(i)
      if(i<n)gauss_vector(i+1)=v(1)*fac*w(i+1)
    enddo ! i
  END FUNCTION gaussian_random_number


  SUBROUTINE characterize_dist(m,E,W,var,skew,kurt)
    !-------------------------------------------------!
    ! Evaluate the variance, skewness and kurtosis of !
    ! a data set.                                     !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: m
    DOUBLE PRECISION, INTENT(in) :: E(m), W(m)
    DOUBLE PRECISION, INTENT(out) :: var, skew, kurt
    DOUBLE PRECISION V1, V2, V3, V4, ave_E, m2, m3, m4, K2, K3, K4, &
       &num1, num2, denom
    V1=sum(W)
    V2=sum(W**2)
    V3=sum(W**3)
    V4=sum(W**4)
    ave_E=sum(E*W)/V1
    m2=sum(W*(E-ave_E)**2)/V1
    m3=sum(W*(E-ave_E)**3)/V1
    m4=sum(W*(E-ave_E)**4)/V1
    K2=m2*(V1*V1/(V1*V1-V2))
    K3=m3*(V1*V1*V1/(V1*V1*V1-3.d0*V1*V2+2.d0*V3))
    num1=V1*V1*(V1**4-4.d0*V1*V3+3.d0*V2*V2)
    num2=3.d0*V1*V1*(V1**4-2.d0*V1*V1*V2+4.d0*V1*V3-3.d0*V2*V2)
    denom=(V1*V1-V2)*(V1**4-6.d0*V1*V1*V2+8.d0*V1*V3+3.d0*V2*V2-6.d0*V4)
    K4=(m4*(num1/denom)-m2*m2*(num2/denom))
    ! Output variables.
    var=K2
    skew=K3/K2**1.5d0
    kurt=K4/(K2*K2)
  END SUBROUTINE characterize_dist


  DOUBLE PRECISION FUNCTION find_pth_smallest(p,n,x)
    !----------------------------------------------------------------!
    ! This function returns the P-th smallest element of the set of  !
    ! N double-precision data in vector X. Uses pivoting to avoid    !
    ! scrambling the input vector.                                   !
    ! The algorithm is based on partitioning: given an index range   !
    ! in the X vector, a random element is chosen (the middle one in !
    ! the partition range, for example), called the partitioning     !
    ! element, 'a', and the array is rearranged so that all elements !
    ! >= a are moved to the right of a, and all elements <= a to its !
    ! left. After partitioning, a will be in its sorted position,    !
    ! and we will know on which side of it to look for the p-th      !
    ! smallest element.                                              !
    ! For more information see:                                      !
    !  http://en.wikipedia.org/wiki/Selection_algorithm              !
    !   #Partition-based_general_selection_algorithm                 !
    ! PLR 05.2009                                                    !
    !----------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: p,n
    DOUBLE PRECISION,INTENT(in) :: x(n)
    INTEGER i,ileft,iright,region_right,region_left,region_centre,map(n)
    DOUBLE PRECISION x_part

    ! Initialize map vector. When we say "move an element", we mean swap
    ! its 'map' value with another; X is left untouched.
    do i=1,n
      map(i)=i
    enddo

    ! Initialize the partitioning region to the whole vector and start main loop.
    region_left=1 ; region_right=n
    do

      ! Check if we've finished.
      if(region_right<=region_left+1)then ! trivial partition, we're done
        if(region_right==region_left+1)then ! last element to check
          if(x(map(region_right))<x(map(region_left)))&
             &call iswap1(map(region_right),map(region_left))
        endif ! region_right==region_left+1
        find_pth_smallest=x(map(p)) ! p-th largest element is at this position
        return
      endif

      ! Get a potential partitioning element and move it to region_left+1.
      region_centre=(region_left+region_right)/2
      call iswap1(map(region_centre),map(region_left+1))
      ! Rearrange region_left, region_left+1 and region_right so that the
      ! elements are sorted ascendingly. The partitioning element 'x_part' ends
      ! up at region_left+1.
      if(x(map(region_left))>x(map(region_right)))&
         &call iswap1(map(region_left),map(region_right))
      if(x(map(region_left+1))>x(map(region_right)))&
         &call iswap1(map(region_left+1),map(region_right))
      if(x(map(region_left))>x(map(region_left+1)))&
         &call iswap1(map(region_left),map(region_left+1))

      ! Now we scan right from region_left+2 and left from region_right-1, and
      ! move elements >= x_part to the right and those <= x_part to the left.
      ileft=region_left+1
      iright=region_right
      x_part=x(map(region_left+1))
      do
        ! Find the next element on the left that is greater that or equal to
        ! x_part.
        do while(ileft<n)
          ileft=ileft+1
          if(x(map(ileft))>=x_part)exit
        enddo ! ileft
        ! Find the next element on the right that is less than or equal to
        ! x_part.
        do while(iright>1)
          iright=iright-1
          if(x(map(iright))<=x_part)exit
        enddo ! iright
        ! If we've crossed the left and right indices, exit.
        if(iright<ileft)exit
        ! Swap the left and right elements.
        call iswap1(map(ileft),map(iright))
      enddo ! ileft,iright
      ! Move x_part to iright.
      call iswap1(map(region_left+1),map(iright))

      ! Select next partitioning region so that the p-th smallest vector is
      ! inside it.
      if(iright>=p)region_right=iright-1
      if(iright<=p)region_left=ileft

    enddo ! main loop

  END FUNCTION find_pth_smallest


  DOUBLE PRECISION FUNCTION median(n,x)
    !-----------------------------------------------------------------!
    ! This function returns the median of a set of n data in array x. !
    !-----------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    DOUBLE PRECISION,INTENT(in) :: x(n)
    INTEGER p
    DOUBLE PRECISION median1,median2
    if(mod(n,2)==0)then
      p=n/2
      median1=find_pth_smallest(p,n,x)
      median2=find_pth_smallest(p+1,n,x)
      median=0.5d0*(median1+median2)
    else
      p=(n+1)/2
      median=find_pth_smallest(p,n,x)
    endif ! n even
  END FUNCTION median


  SUBROUTINE iswap1(i,j)
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: i,j
    INTEGER k
    k=i
    i=j
    j=k
  END SUBROUTINE iswap1


  SUBROUTINE quit (msg)
    !---------------------!
    ! Quit with an error. !
    !---------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in), OPTIONAL :: msg
    if (present(msg)) then
      write(6,*)'ERROR : '//msg
    else
      write(6,*)'Quitting.'
    endif
    stop
  END SUBROUTINE quit


END PROGRAM polyfit
