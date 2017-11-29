PROGRAM polyfit
  !-------------------------------------------------!
  ! POLYFIT                                         !
  ! =======                                         !
  ! Toolbox for performing and manipulating fits to !
  ! polynomials of data with error bars.            !
  !                                                 !
  ! PLR 09.2015                                     !
  !-------------------------------------------------!
  IMPLICIT NONE

  ! Global variables.
  ! Use weights in Monte Carlo fits?
  LOGICAL :: MONTE_CARLO_FIT_WEIGHTS=.false.
  ! Use 1/chi^2 weights in Monte Carlo?
  LOGICAL :: MONTE_CARLO_CHI2_WEIGHTS=.false.
  ! Precision for real-to-integer conversions and exponent comparisons.
  DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-10
  ! Transformations.
  INTEGER,PARAMETER :: NTRANSF=3
  INTEGER,PARAMETER :: ITRANSF_NONE=0,ITRANSF_REC=1,ITRANSF_LOG=2,ITRANSF_EXP=3
  CHARACTER(32),PARAMETER :: TRANSF_NAME(0:NTRANSF)=&
     &(/'linear     ',  'reciprocal ',  'logarithmic',  'exponential'/)
  LOGICAL,PARAMETER :: TRANSF_REQ_NONZERO(0:NTRANSF)=&
     &(/.false.,.true.,.true.,.false./)
  LOGICAL,PARAMETER :: TRANSF_REQ_POSITIVE(0:NTRANSF)=&
     &(/.false.,.false.,.true.,.false./)

  call main()


CONTAINS


  SUBROUTINE main()
    !--------------!
    ! Main driver. !
    !--------------!
    IMPLICIT NONE
    ! (x,y[,dx][,dy]) data.
    INTEGER nxy
    LOGICAL have_dx,have_dy
    DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:),dx(:),dy(:)
    ! Input file information.
    INTEGER ncolumn,icol_x,icol_y,icol_dx,icol_dy
    CHARACTER(256) fname
    ! Misc variables.
    INTEGER ierr

    ! Write header.
    write(6,*)'======='
    write(6,*)'POLYFIT'
    write(6,*)'======='
    write(6,*)

    ! Read in (x,y[,dx][,dy]) data from file specified by user.
    call get_file(fname)
    call check_file(fname,nxy,ncolumn)
    if(nxy<1)call quit('No data found in file "'//trim(fname)//'".')
    if(ncolumn<2)call quit('Not enough data columns found in file "'//&
       &trim(fname)//'".')
    write(6,*)'"'//trim(fname)//'" contains '//trim(i2s(nxy))//' lines with '//&
       &trim(i2s(ncolumn))//' columns.'
    call get_columns(ncolumn,icol_x,icol_y,icol_dx,icol_dy)
    have_dx=icol_dx>0
    have_dy=icol_dy>0
    allocate(x(nxy),y(nxy),dx(nxy),dy(nxy),stat=ierr)
    if(ierr/=0)call quit('Allocation error.')
    call read_file(fname,ncolumn,icol_x,icol_y,icol_dx,icol_dy,nxy,x,y,dx,dy)

    ! Run main menu.
    call user_interaction(nxy,have_dx,have_dy,x,y,dx,dy)

    ! Finish off.
    deallocate(x,y,dx,dy)
    write(6,*)'Program finished.'
    write(6,*)

  END SUBROUTINE main


  SUBROUTINE user_interaction(nxy,have_dx,have_dy,x,y,dx,dy)
    !-------------------------!
    ! Interact with the user. !
    !-------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(inout) :: x(nxy),y(nxy),dx(nxy),dy(nxy)
    ! Transformed variables.
    DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy)
    ! Masked variables.
    INTEGER rnxy
    DOUBLE PRECISION rx(nxy),ry(nxy),rdx(nxy),rdy(nxy),rtx(nxy),rty(nxy),&
       &rdtx(nxy),rdty(nxy)
    ! Fitting function.
    INTEGER npoly,itransfx,itransfy,i
    LOGICAL mask(nxy)
    DOUBLE PRECISION tx0,pow(nxy)
    ! Misc local variables.
    INTEGER ierr
    CHARACTER(2048) char2048

    ! Initialize.
    tx0=0.d0
    itransfx=0
    itransfy=0
    tx=x
    ty=y
    dtx=dx
    dty=dy
    mask=.true.
    rnxy=nxy
    rx=x
    ry=y
    rdx=dx
    rdy=dy
    rtx=tx
    rty=ty
    rdtx=dtx
    rdty=dty
    npoly=2
    pow(1:2)=(/0.d0,1.d0/)

    ! Loop over user actions.
    do

      ! Print list of actions.
      write(6,*)'---------'
      write(6,*)'Main menu'
      write(6,*)'---------'
      write(6,*)'Available actions:'
      write(6,*)'* Setup:'
      write(6,*)'  [a] Set axis scales [X='//&
         &trim(TRANSF_NAME(itransfx))//',Y='//trim(TRANSF_NAME(itransfy))//']'
      write(6,*)'  [m] Set data mask [using '//trim(i2s(count(mask)))//'/'//&
         &trim(i2s(nxy))//' points]'
      write(6,'(1x,a,es11.4,a)')'  [0] Set fit centre [X0=',tx0,']'
      if(have_dx.or.have_dy)then
        write(6,*)'  [w] Toggle use of weighted fits in Monte Carlo [',&
           &MONTE_CARLO_FIT_WEIGHTS,']'
        write(6,*)'  [x] Toggle use of 1/chi^2 as Monte Carlo weights [',&
           &MONTE_CARLO_CHI2_WEIGHTS,']'
      endif
      write(6,*)'* Pre-fit analysis:'
      write(6,*)'  [r] Show basic data-range statistics'
      write(6,*)'  [e] Show expansion-order analysis'
      write(6,*)'  [d] Show dual number of points/expansion-order analysis'
      write(6,*)'* Fitting:'
      write(6,*)'  [f] Set form and perform fit ['//&
         &trim(print_poly_sym(npoly,pow,tx0))//']'
      write(6,*)'* Post-fit analysis:'
      write(6,*)'  [v] Compute values/derivatives of fit'
      write(6,*)'* Quitting:'
      write(6,*)'  [q] Quit'
      write(6,*)

      ! Read user choice.
      write(6,*)'Enter one of the above options:'
      read(5,'(a)',iostat=ierr)char2048
      if(ierr/=0)char2048='q'
      write(6,*)

      ! Perform selected action.
      select case(trim(adjustl(char2048)))
      case('a')
        call ask_scale(nxy,x,y,itransfx,itransfy)
        call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
        call scale_transform(nxy,itransfy,y,ty,have_dy,dy,dty)
        rtx(1:rnxy)=pack(tx,mask)
        rty(1:rnxy)=pack(ty,mask)
        rdtx(1:rnxy)=pack(dtx,mask)
        rdty(1:rnxy)=pack(dty,mask)
      case('m')
        call ask_mask(nxy,x,y,tx,ty,mask)
        rnxy=count(mask)
        rx(1:rnxy)=pack(x,mask)
        ry(1:rnxy)=pack(y,mask)
        rdx(1:rnxy)=pack(dx,mask)
        rdy(1:rnxy)=pack(dy,mask)
        rtx(1:rnxy)=pack(tx,mask)
        rty(1:rnxy)=pack(ty,mask)
        rdtx(1:rnxy)=pack(dtx,mask)
        rdty(1:rnxy)=pack(dty,mask)
        if(rnxy<npoly)then
          npoly=2
          pow(1:2)=(/0.d0,1.d0/)
          write(6,*)'Fit form reset.'
          write(6,*)
        endif
      case('0')
        call ask_centre(rnxy,rtx,rty,tx0)
      case('f')
        call ask_poly(rnxy,i,pow)
        if(i>0.and.i<=rnxy)npoly=i
        call show_poly(rnxy,have_dx,have_dy,rtx,rty,rdtx,rdty,tx0,npoly,pow)
      case('w')
        MONTE_CARLO_FIT_WEIGHTS=.not.MONTE_CARLO_FIT_WEIGHTS
      case('x')
        MONTE_CARLO_CHI2_WEIGHTS=.not.MONTE_CARLO_CHI2_WEIGHTS
      case('r')
        call show_statistics(rnxy,have_dx,have_dy,rtx,rty,rdtx,rdty)
      case('e')
        call show_exporder_assessment(rnxy,have_dx,have_dy,rtx,rty,rdtx,rdty,&
           &tx0,i,pow)
        if(i<1)cycle
        npoly=i
        call show_poly(rnxy,have_dx,have_dy,rtx,rty,rdtx,rdty,tx0,npoly,pow)
      case('d')
        call show_dual_assessment(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,&
           &itransfy,tx0,i,pow,mask)
        if(i>0)then
          npoly=i
          rnxy=count(mask)
          rx(1:rnxy)=pack(x,mask)
          ry(1:rnxy)=pack(y,mask)
          rdx(1:rnxy)=pack(dx,mask)
          rdy(1:rnxy)=pack(dy,mask)
          rtx(1:rnxy)=pack(tx,mask)
          rty(1:rnxy)=pack(ty,mask)
          rdtx(1:rnxy)=pack(dtx,mask)
          rdty(1:rnxy)=pack(dty,mask)
        endif
      case('v')
        call eval_fit_values_derivs(rnxy,have_dx,have_dy,rx,ry,rdx,rdy,&
           &itransfx,itransfy,tx,ty,tx0,npoly,pow)
      case('q','')
        call quit()
      case default
        call quit()
      end select

    enddo

  END SUBROUTINE user_interaction


  SUBROUTINE show_statistics(nxy,have_dx,have_dy,tx,ty,dtx,dty)
    !-----------------------------------------------------!
    ! Show min/centre/mean/median/max of X, Y, dX and dY. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: tx(nxy),ty(nxy),dtx(nxy),dty(nxy)
    DOUBLE PRECISION vmin,vmax,vcentre,vmean,vmedian

    ! Print stats.
    write(6,*)'Data-range statistics:'
    write(6,*)
    write(6,'(t8,a,t21,a,t34,a,t47,a,t60,a)')'min','mean','median','centre',&
       &'max'
    write(6,'(t2,70("-"))')
    ! Transformed X and dX.
    vmin=minval(tx)
    vmax=maxval(tx)
    vcentre=.5d0*(vmin+vmax)
    vmean=sum(tx)/dble(nxy)
    vmedian=median(nxy,tx)
    write(6,'(t2,a,t5,5(1x,es12.4))')'X',vmin,vmean,vmedian,vcentre,vmax
    if(have_dx)then
      vmin=minval(dtx)
      vmax=maxval(dtx)
      vcentre=.5d0*(vmin+vmax)
      vmean=sum(dtx)/dble(nxy)
      vmedian=median(nxy,dtx)
      write(6,'(t2,a,t5,5(1x,es12.4))')'dX',vmin,vmean,vmedian,vcentre,vmax
    endif
    ! Transformed Y and dY.
    vmin=minval(ty)
    vmax=maxval(ty)
    vcentre=.5d0*(vmin+vmax)
    vmean=sum(ty)/dble(nxy)
    vmedian=median(nxy,ty)
    write(6,'(t2,a,t5,5(1x,es12.4))')'Y',vmin,vmean,vmedian,vcentre,vmax
    if(have_dy)then
      vmin=minval(dty)
      vmax=maxval(dty)
      vcentre=.5d0*(vmin+vmax)
      vmean=sum(dty)/dble(nxy)
      vmedian=median(nxy,dty)
      write(6,'(t2,a,t5,5(1x,es12.4))')'dY',vmin,vmean,vmedian,vcentre,vmax
    endif
    write(6,'(t2,70("-"))')
    write(6,*)

  END SUBROUTINE show_statistics


  SUBROUTINE ask_scale(nxy,x,y,itransfx,itransfy)
    !-------------------------------------!
    ! Ask the user to choose axis scales. !
    !-------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy)
    INTEGER,INTENT(inout) :: itransfx,itransfy
    CHARACTER(2048) char2048
    INTEGER itransf,ierr

    ! List available transformations.
    write(6,*)'The following scales are available:'
    do itransf=0,NTRANSF
      write(6,*)'  ['//trim(i2s(itransf))//'] '//&
         &trim(TRANSF_NAME(itransf))
    enddo ! itransf
    write(6,*)

    ! Ask user for X transformation.
    write(6,*)'Enter index of transformed X scale to use ['//&
      &trim(TRANSF_NAME(itransfx))//']:'
    read(5,'(a)',iostat=ierr)char2048
    if(ierr==0)then
      read(char2048,*,iostat=ierr)itransf
      if(ierr==0)then
        if(itransf>=0.and.itransf<=NTRANSF)then
          if(TRANSF_REQ_NONZERO(itransf).and.any(abs(x)<=0.d0))then
            ierr=1
            write(6,*)'Transformation requires all x to be non-zero.'
          elseif(TRANSF_REQ_POSITIVE(itransf).and.any(x<0.d0))then
            ierr=1
            write(6,*)'Transformation requires all x to be positive.'
          else
            itransfx=itransf
          endif
        else
          ierr=1
          write(6,*)'Index out of range.'
        endif
      else
        do itransf=0,NTRANSF
          if(trim(adjustl(char2048))==trim(TRANSF_NAME(itransf)))exit
        enddo ! itransf
        if(itransf<=NTRANSF)then
          ierr=0
          itransfx=itransf
        endif
      endif
    endif
    if(ierr==0)then
      write(6,*)'X is x in '//trim(TRANSF_NAME(itransfx))//' scale.'
    else
      write(6,*)'X scale unchanged.'
    endif
    write(6,*)

    ! Ask user for Y transformation.
    write(6,*)'Enter index of transformed Y scale to use ['//&
      &trim(TRANSF_NAME(itransfy))//']:'
    read(5,'(a)',iostat=ierr)char2048
    if(ierr==0)then
      read(char2048,*,iostat=ierr)itransf
      if(ierr==0)then
        if(itransf>=0.and.itransf<=NTRANSF)then
          if(TRANSF_REQ_NONZERO(itransf).and.any(abs(y)<=0.d0))then
            ierr=1
            write(6,*)'Transformation requires all y to be non-zero.'
          elseif(TRANSF_REQ_POSITIVE(itransf).and.any(y<0.d0))then
            ierr=1
            write(6,*)'Transformation requires all y to be positive.'
          else
            itransfy=itransf
          endif
        else
          ierr=1
          write(6,*)'Index out of range.'
        endif
      else
        do itransf=0,NTRANSF
          if(trim(adjustl(char2048))==trim(TRANSF_NAME(itransf)))exit
        enddo ! itransf
        if(itransf<=NTRANSF)then
          ierr=0
          itransfy=itransf
        endif
      endif
    endif
    if(ierr==0)then
      write(6,*)'Y is y in '//trim(TRANSF_NAME(itransfy))//' scale.'
    else
      write(6,*)'Y scale unchanged.'
    endif
    write(6,*)

  END SUBROUTINE ask_scale


  SUBROUTINE ask_mask(nxy,x,y,tx,ty,mask)
    !----------------------------------------------!
    ! Ask user to select which data to use in fit. !
    !----------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),tx(nxy),ty(nxy)
    LOGICAL,INTENT(inout) :: mask(nxy)
    CHARACTER(2048) char2048
    LOGICAL by_threshold,by_highest,tmask(nxy)
    DOUBLE PRECISION sortvec(nxy),t1
    INTEGER i,indx(nxy),ierr

    ! Print instructions.
    write(6,*)'Use a mask string of the form <variable><criterion>, where:'
    write(6,*)'* <variable> can be x, y (original), X or Y (scale-transformed)'
    write(6,*)'* <criterion> can be <t, >t, #<n, or #>n, meaning:'
    write(6,*)'  <t  : <variable> less than threshold t'
    write(6,*)'  >t  : <variable> greater than threshold t'
    write(6,*)'  #<n : n points with the smallest <variable>'
    write(6,*)'  #>n : n points with the largest <variable>'

    ! Read user choice.
    write(6,*)'Enter a mask string:'
    read(5,'(a)',iostat=ierr)char2048
    if(ierr/=0)return
    write(6,*)

    ! Parse string.
    by_threshold=.true.
    by_highest=.false.
    select case(char2048(1:1))
    case('X')
      sortvec=tx
    case('Y')
      sortvec=ty
    case('x')
      sortvec=x
    case('y')
      sortvec=y
    case default
      return
    end select
    char2048=char2048(2:)
    if(char2048(1:1)=='#')then
      char2048=char2048(2:)
      by_threshold=.false.
    endif
    select case(char2048(1:1))
    case('<')
      continue
    case('>')
      by_highest=.true.
      sortvec=-sortvec
    case default
      return
    end select
    char2048=char2048(2:)

    ! Perform selected action.
    if(by_threshold)then
      read(char2048,*,iostat=ierr)t1
      if(ierr/=0)return
      if(by_highest)t1=-t1
      tmask=sortvec<t1
    else
      read(char2048,*,iostat=ierr)i
      if(ierr/=0)return
      if(i<1.or.i>nxy)return
      call isort(nxy,sortvec,indx)
      tmask=.false.
      tmask(indx(1:i))=.true.
    endif

    ! Make sure we have at least two points.
    if(count(tmask)<2)return
    mask=tmask
    write(6,*)'New mask: ',mask
    write(6,*)

  END SUBROUTINE ask_mask


  SUBROUTINE ask_centre(nxy,tx,ty,tx0)
    !------------------------------------------!
    ! Ask the user to set offsets for the fit. !
    !------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    DOUBLE PRECISION,INTENT(in) :: tx(nxy),ty(nxy)
    DOUBLE PRECISION,INTENT(inout) :: tx0
    CHARACTER(2048) char2048
    DOUBLE PRECISION t1
    INTEGER i,ierr

    ! Get fit centre X0.
    write(6,*)'Current fit centre X0 = ',tx0
    write(6,*)'Enter new centre or left/mean/median/centre/right or min/max:'
    read(5,'(a)',iostat=ierr)char2048
    write(6,*)
    if(ierr==0)then
      char2048=adjustl(char2048)
      select case(trim(char2048))
      case('left')
        tx0=minval(tx)
      case('right')
        tx0=maxval(tx)
      case('min')
        do i=1,nxy
          if(abs(ty(i)-minval(ty))<tol_zero)then
            tx0=tx(i)
            exit
          endif
        enddo
      case('max')
        do i=1,nxy
          if(abs(ty(i)-maxval(ty))<tol_zero)then
            tx0=tx(i)
            exit
          endif
        enddo
      case('centre')
        tx0=.5d0*(minval(tx)+maxval(tx))
      case('mean')
        tx0=sum(tx)/dble(nxy)
      case('median')
        tx0=median(nxy,tx)
      case('')
        ierr=1
      case default
        read(char2048,*,iostat=ierr)t1
        if(ierr==0)then
          tx0=t1
        else
        endif
      end select
    endif

    ! Report.
    if(ierr==0)then
      write(6,*)'Set X0 = ',tx0
    else
      write(6,*)'Keeping X0 = ',tx0
    endif
    write(6,*)

  END SUBROUTINE ask_centre


  SUBROUTINE scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    !---------------------------------------!
    ! Apply a scale transformation x -> tx. !
    !---------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx
    DOUBLE PRECISION,INTENT(in) :: x(nxy)
    DOUBLE PRECISION,INTENT(inout) :: tx(nxy)
    LOGICAL,INTENT(in) :: have_dx
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: dx(nxy)
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: dtx(nxy)
    if(present(dtx))dtx=0.d0
    select case(itransfx)
    case(ITRANSF_NONE)
      tx=x
      if(have_dx.and.present(dx).and.present(dtx))dtx=dx
    case(ITRANSF_REC)
      tx=1.d0/x
      if(have_dx.and.present(dx).and.present(dtx))dtx=dx*tx**2
    case(ITRANSF_LOG)
      tx=log(x)
      if(have_dx.and.present(dx).and.present(dtx))dtx=dx/x
    case(ITRANSF_EXP)
      tx=exp(x)
      if(have_dx.and.present(dx).and.present(dtx))dtx=dx*tx
    end select
  END SUBROUTINE scale_transform


  SUBROUTINE show_exporder_assessment(nxy,have_dx,have_dy,x,y,dx,dy,x0,&
     &npoly_out,pow_out)
    !--------------------------------------------------------!
    ! Perform an assessment of the performance of polynomial !
    ! fits as a function of expansion order.                 !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0
    INTEGER,INTENT(inout) :: npoly_out
    DOUBLE PRECISION,INTENT(inout) :: pow_out(nxy)
    CHARACTER(2048) char2048
    DOUBLE PRECISION,ALLOCATABLE :: chi2_vector(:),rmsc_100_vector(:),&
       &rmsc_150_vector(:),rmsc_200_vector(:),pow(:),a(:),da(:),op_pow(:),&
       &op_a(:)
    DOUBLE PRECISION weight(nxy),txrange
    INTEGER ntest,npoly,op_npoly,i,np,ierr
    LOGICAL weighted

    ! Initialize.
    npoly_out=0

    ! Evaluate transformed variables and fit weights.
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

    ! Assess integer expansion orders up to order 8.
    ntest=min(nxy-1,9)

    ! Allocate work arrays.
    allocate(chi2_vector(ntest),rmsc_100_vector(ntest),&
       &rmsc_150_vector(ntest),rmsc_200_vector(ntest),stat=ierr)
    if(ierr/=0)call quit('Allocation error.')

    ! Loop over expansion orders.
    do npoly=1,ntest
      ! Allocate work arrays.
      allocate(pow(npoly),a(npoly),da(npoly),op_pow((npoly*(npoly+1))/2),&
         &op_a((npoly*(npoly+1))/2),stat=ierr)
      if(ierr/=0)call quit('Allocation error.')
      ! Fit to polynomial of order npoly-1.
      do i=1,npoly
        pow(i)=dble(i-1)
      enddo ! i
      call perform_fit(nxy,npoly,x-x0,y,pow,a,weighted,weight,da)
      ! Evaluate chi^2.
      chi2_vector(npoly)=chi_squared(nxy,npoly,x-x0,y,weight,pow,a,weighted)
      ! Evaluate root mean square curvature of fit as measure of smoothness.
      np=npoly
      txrange=maxval(x)-minval(x)
      call deriv_poly(np,pow,a,op_npoly,op_pow,op_a)
      call deriv_poly(op_npoly,op_pow,op_a)
      call square_poly(op_npoly,op_pow,op_a)
      call int_poly(op_npoly,op_pow,op_a)
      rmsc_100_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,maxval(x))-&
         &eval_poly(op_npoly,op_pow,op_a,minval(x))&
         &)/txrange))
      rmsc_150_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,maxval(x)+0.25d0*txrange)-&
         &eval_poly(op_npoly,op_pow,op_a,minval(x)-0.25d0*txrange)&
         &)/(1.5d0*txrange)))
      rmsc_200_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,maxval(x)+0.5d0*txrange)-&
         &eval_poly(op_npoly,op_pow,op_a,minval(x)-0.5d0*txrange)&
         &)/(2.d0*txrange)))
      ! Destroy work arrays.
      deallocate(pow,a,da,op_pow,op_a)
    enddo ! npoly

    ! Report.
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
    write(6,'(1x,a,t10,a,t23,a,t36,a,t49,a)')'Order','chi^2/Ndf','Data range',&
       &'+50%','+100%'
    write(6,'(1x,57("-"))')
    do npoly=1,ntest
      write(6,'(1x,i5,4(1x,es12.4),1x,a)')npoly-1,&
         &chi2_vector(npoly)/dble(nxy-npoly),&
         &rmsc_100_vector(npoly),rmsc_150_vector(npoly),rmsc_200_vector(npoly)
    enddo
    write(6,'(1x,57("-"))')
    write(6,*)

    ! Allow user to choose one of the above expansion orders.
    do
      write(6,*)'Enter order to choose fit form [''plot'' to plot, empty to &
         &skip]:'
      read(5,'(a)',iostat=ierr)char2048
      if(ierr/=0)exit
      write(6,*)
      select case(trim(char2048))
      case('')
        exit
      case('plot')
        call plot_exporder(nxy,have_dx,have_dy,x,y,dx,dy,x0)
      case default
        read(char2048,*,iostat=ierr)i
        if(ierr==0)then
          if(i>=0.and.i<=ntest)then
            npoly_out=i+1
            pow_out(1:npoly_out)=(/(i,i=0,npoly_out-1)/)
          endif
        endif
        exit
      end select
    enddo

    ! Destroy work arrays.
    deallocate(chi2_vector,rmsc_100_vector,rmsc_150_vector,rmsc_200_vector)

  END SUBROUTINE show_exporder_assessment


  SUBROUTINE plot_exporder(nxy,have_dx,have_dy,x,y,dx,dy,x0)
    !-----------------------------------------------------!
    ! Plot polynomial fits of different expansion orders. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0
    DOUBLE PRECISION,ALLOCATABLE :: pow(:),a(:),da(:)
    DOUBLE PRECISION weight(nxy),xrange,xplot,dxplot
    INTEGER npoly,i,ierr
    LOGICAL weighted
    ! Parameters.
    INTEGER,PARAMETER :: npoint=1000
    INTEGER,PARAMETER :: io=10

    ! Evaluate transformed variables and fit weights.
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

    ! Prepare to write plot.
    xrange=maxval(x)-minval(x)
    dxplot=(xrange*2.d0)/dble(npoint)
    open(unit=io,file='poly_orders.dat',status='replace')

    ! Plot original data.
    if(have_dx.and.have_dy)then
      write(io,'(a)')'@type xydxdy'
      do i=1,nxy
        write(io,*)x(i),y(i),dx(i),dy(i)
      enddo ! i
    elseif(have_dx)then
      write(io,'(a)')'@type xydx'
      do i=1,nxy
        write(io,*)x(i),y(i),dx(i)
      enddo ! i
    elseif(have_dy)then
      write(io,'(a)')'@type xydy'
      do i=1,nxy
        write(io,*)x(i),y(i),dy(i)
      enddo ! i
    else
      write(io,'(a)')'@type xy'
      do i=1,nxy
        write(io,*)x(i),y(i)
      enddo ! i
    endif
    write(io,'(a)')'&'

    ! Loop over expansion orders.
    do npoly=1,min(nxy-1,9)
      write(io,'(a)')'@type xy'
      ! Prepare for fit.
      allocate(pow(npoly),a(npoly),da(npoly),stat=ierr)
      if(ierr/=0)call quit('Allocation error.')
      do i=1,npoly
        pow(i)=dble(i-1)
      enddo ! i
      ! Perform fit.
      call perform_fit(nxy,npoly,x-x0,y,pow,a,weighted,weight,da)
      ! Dump plot.
      xplot=minval(x)-0.5d0*xrange
      do i=0,npoint
        write(io,*)xplot+x0,eval_poly(npoly,pow,a,xplot)
        xplot=xplot+dxplot
      enddo ! i
      ! Clean up after fit.
      deallocate(pow,a,da)
      write(io,'(a)')'&'
    enddo ! npoly

    ! Close file and report.
    close(io)
    write(6,*)'Data written to poly_orders.dat.'
    write(6,*)

  END SUBROUTINE plot_exporder


  SUBROUTINE show_dual_assessment(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,&
     &itransfy,tx0,npoly_out,pow_out,mask_out)
    !-----------------------------------------------------------!
    ! Perform an "extrapolation" assessment in which the number !
    ! of points included in the fit and the expansion order are !
    ! simultaneously optimized.                                 !
    !-----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),tx0
    INTEGER,INTENT(inout) :: npoly_out
    DOUBLE PRECISION,INTENT(inout) :: pow_out(nxy)
    LOGICAL,INTENT(inout) :: mask_out(nxy)
    DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy),weight(nxy),&
       &rtx(nxy),rty(nxy),rweight(nxy),xtarget,min_chi2,t1,&
       &rx(nxy),ry(nxy),rdx(nxy),rdy(nxy),list_f(nxy),list_df(nxy),&
       &list_gof(nxy),list_a(nxy,nxy),list_da(nxy,nxy),&
       &pow(nxy),a(nxy),da(nxy),txrange,txplot,dtxplot
    INTEGER i,ntest,npoly,rnxy,nderiv,ierr,indx(nxy),rnxy_min_chi2,nrandom,&
       &npoly_best,rnxy_best,ipoly,jpoly,list_rnxy(nxy)
    LOGICAL mask(nxy),weighted
    CHARACTER(20) drop_by,drop_criterion
    CHARACTER(2048) char2048
    DOUBLE PRECISION fmean(1),ferr(1),best_f,best_df
    ! Parameters.
    INTEGER,PARAMETER :: DEFAULT_NRANDOM=100000
    INTEGER,PARAMETER :: npoint=1000
    INTEGER,PARAMETER :: io=10

    ! Initialize.
    npoly_out=0

    ! Initialize internal variables.  FIXME - expose to user.
    nderiv=0
    xtarget=0.d0
    nrandom=DEFAULT_NRANDOM
    drop_by='X'
    drop_criterion='largest'

    ! Evaluate transformed variables and fit weights.
    call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    call scale_transform(nxy,itransfy,y,ty,have_dy,dy,dty)
    if(have_dx.and.have_dy)then
      weight=1.d0/(dtx*dty)**2
    elseif(have_dx)then
      weight=1.d0/dtx**2
    elseif(have_dy)then
      weight=1.d0/dty**2
    else
      weight=1.d0
    endif
    weighted=have_dx.or.have_dy

    ! Perform sort.
    select case(trim(drop_by))
    case('X')
      call isort(nxy,tx,indx)
    case('Y')
      call isort(nxy,ty,indx)
    case('x')
      call isort(nxy,x,indx)
    case('y')
      call isort(nxy,y,indx)
    case('index')
      indx=(/(i,i=1,nxy)/)
    end select

    ! Assess integer expansion orders up to order 8.
    ntest=min(nxy-1,9)
    pow(1:nxy)=(/(dble(i-1),i=1,nxy)/)

    ! Loop over expansion orders.
    write(6,'(1x,a,t8,a,t14,a,t27,a,t44,a)')'Order','  Nxy','  chi^2/Ndf',&
       &'    Target f','    Target df'
    write(6,'(1x,59("-"))')
    list_rnxy=0
    do npoly=2,ntest
      ! Loop over number of points in fit.
      rnxy_min_chi2=0
      min_chi2=0.d0
      do rnxy=nxy,npoly+1,-1
        ! Build data mask.
        mask=.false.
        select case(trim(drop_criterion))
        case('smallest')
          mask(indx(nxy-rnxy+1:nxy))=.true.
        case('largest')
          mask(indx(1:rnxy))=.true.
        end select
        ! Apply data mask.
        rtx(1:rnxy)=pack(tx,mask)
        rty(1:rnxy)=pack(ty,mask)
        rweight(1:rnxy)=pack(weight,mask)
        ! Fit to polynomial of order npoly-1.
        call perform_fit(rnxy,npoly,rtx-tx0,rty,pow,a,weighted,rweight,da)
        ! Evaluate chi^2.
        t1=chi_squared(rnxy,npoly,rtx-tx0,rty,rweight,pow,a,weighted)/&
           &dble(rnxy-npoly)
        if(rnxy_min_chi2==0.or.t1<min_chi2)then
          rnxy_min_chi2=rnxy
          min_chi2=t1
        endif
      enddo ! rnxy
      ! Evaluate estimate at number of points that minimizes gof measure.
      rnxy=rnxy_min_chi2
      ! Build data mask.
      mask=.false.
      select case(trim(drop_criterion))
      case('smallest')
        mask(indx(nxy-rnxy+1:nxy))=.true.
      case('largest')
        mask(indx(1:rnxy))=.true.
      end select
      ! Apply data mask.
      rx(1:rnxy)=pack(x,mask)
      ry(1:rnxy)=pack(y,mask)
      rdx(1:rnxy)=pack(dx,mask)
      rdy(1:rnxy)=pack(dy,mask)
      ! Fit to polynomial of order npoly-1.
      call eval_fit_monte_carlo(rnxy,have_dx,have_dy,rx,ry,rdx,rdy,&
         &itransfx,itransfy,tx0,npoly,pow,nrandom,nderiv,.false.,1,&
         &(/xtarget/),fmean,ferr,amean=a,aerr=da)
      ! Report.
      write(6,'(1x,i5,1x,i5,1x,es12.4,2(1x,es16.8))')npoly-1,rnxy_min_chi2,&
         &min_chi2,fmean(1),ferr(1)
      ! Store.
      list_f(npoly)=fmean(1)
      list_df(npoly)=ferr(1)
      list_rnxy(npoly)=rnxy_min_chi2
      list_gof(npoly)=min_chi2
      list_a(1:npoly,npoly)=a(1:npoly)
      list_da(1:npoly,npoly)=da(1:npoly)
    enddo ! npoly
    write(6,'(1x,59("-"))')
    write(6,*)

    ! Find best choice.
    npoly_best=0
    best_f=0.d0
    best_df=0.d0
    rnxy_best=0
    do ipoly=1,ntest
      if(list_rnxy(ipoly)==0)cycle
      do jpoly=ipoly+1,ntest
        if(list_rnxy(jpoly)==0)cycle
        if(list_gof(jpoly)>list_gof(ipoly))cycle
        if(abs(list_f(ipoly)-list_f(jpoly))>&
           &1.0d0*sqrt(list_df(ipoly)**2+list_df(jpoly)**2))exit
      enddo ! jpoly
      if(jpoly>ntest)then
        npoly_best=ipoly
        best_f=list_f(ipoly)
        best_df=list_df(ipoly)
        rnxy_best=list_rnxy(ipoly)
        exit
      endif
    enddo ! ipoly

    if(npoly_best==0)then
      write(6,*)'Test found no good choice.'
    else
      write(6,*)'Best ('//trim(i2s(npoly_best-1))//','//&
         &trim(i2s(rnxy_best))//'):',best_f,'+/-',best_df
    endif
    write(6,*)

    ! Allow user to choose one of the above expansion orders.
    do
      if(npoly_best==0)then
        write(6,*)'Type ''plot'' to plot, empty to skip:'
      else
        write(6,*)'Type ''y'' to accept fit form and mask, ''plot'' to plot, &
           &empty to skip:'
      endif
      read(5,'(a)',iostat=ierr)char2048
      if(ierr/=0)exit
      write(6,*)
      select case(trim(char2048))
      case('')
        exit
      case('plot')
        ! Prepare to write plot.
        txrange=maxval(tx)-minval(tx)
        dtxplot=(txrange*2.d0)/dble(npoint)
        open(unit=io,file='poly_dual.dat',status='replace')
        ! Plot original data.
        if(have_dx.and.have_dy)then
          write(io,'(a)')'@type xydxdy'
          do i=1,nxy
            write(io,*)tx(i),ty(i),dtx(i),dty(i)
          enddo ! i
        elseif(have_dx)then
          write(io,'(a)')'@type xydx'
          do i=1,nxy
            write(io,*)tx(i),ty(i),dtx(i)
          enddo ! i
        elseif(have_dy)then
          write(io,'(a)')'@type xydy'
          do i=1,nxy
            write(io,*)tx(i),ty(i),dty(i)
          enddo ! i
        else
          write(io,'(a)')'@type xy'
          do i=1,nxy
            write(io,*)tx(i),ty(i)
          enddo ! i
        endif
        write(io,'(a)')'&'
        ! Loop over expansion orders.
        do npoly=2,ntest
          write(io,'(a)')'@type xy'
          ! Dump plot.
          txplot=minval(tx)-0.5d0*txrange
          do i=0,npoint
            write(io,*)txplot+tx0,eval_poly(npoly,pow,list_a(1,npoly),txplot)
            txplot=txplot+dtxplot
          enddo ! i
          write(io,'(a)')'&'
        enddo ! npoly
        close(io)
        write(6,*)'Data written to poly_dual.dat.'
        write(6,*)
      case('y')
        if(npoly_best>0)then
          npoly_out=npoly_best
          pow_out(1:npoly_out)=(/(i,i=0,npoly_out-1)/)
          rnxy=rnxy_best
          mask_out=.false.
          select case(trim(drop_criterion))
          case('smallest')
            mask_out(indx(nxy-rnxy+1:nxy))=.true.
          case('largest')
            mask_out(indx(1:rnxy))=.true.
          end select
        endif
        exit
      case default
        exit
      end select
    enddo

  END SUBROUTINE show_dual_assessment


  SUBROUTINE ask_poly(nxy,npoly,pow)
    !----------------------------------------------------------!
    ! Ask user for string containing exponents for polynomial, !
    ! and return the string and the number of terms.           !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    INTEGER,INTENT(inout) :: npoly
    DOUBLE PRECISION,INTENT(inout) :: pow(nxy)
    CHARACTER(2048) char2048
    INTEGER i1,i2,i,j,ipos,idiv,ierr
    DOUBLE PRECISION t1

    ! Ask user for exponents.
    npoly=0
    write(6,*)'Enter exponents in fitting polynomial (space-separated &
       &list, or i:j, or :j,'
    write(6,*)'or empty to skip):'
    read(5,'(a)',iostat=ierr)char2048
    if(ierr/=0)return
    write(6,*)
    char2048=adjustl(char2048)
    do
      read(char2048,*,iostat=ierr)(t1,i=1,npoly+1)
      if(ierr/=0)exit
      npoly=npoly+1
    enddo
    if(npoly==0)then
      ! Try to parse range.
      ipos=scan(char2048,':')
      if(ipos<1)return
      read(char2048(1:ipos-1),*,iostat=ierr)i1
      if(ierr/=0)i1=0
      read(char2048(ipos+1:),*,iostat=ierr)i2
      if(ierr/=0)return
      npoly=i2-i1+1
      if(npoly<1)return
      if(npoly>nxy)return
      pow(1:npoly)=(/(dble(i),i=i1,i2)/)
    elseif(npoly<nxy)then
      read(char2048,*,iostat=ierr)pow(1:npoly)
    else
      npoly=0
    endif
    if(npoly==0)return

    ! Make near-{,half,third,quarter}-integers exactly {*}-integers.
    do idiv=2,4
      do i=1,npoly
        if(abs(anint(pow(i)*idiv)-pow(i)*idiv)<tol_zero)&
           &pow(i)=anint(pow(i)*idiv)/dble(idiv)
      enddo ! i
    enddo ! idiv

    ! Forbid negative powers.
    if(any(pow<0.d0))then
      write(6,*)'Negative exponents not allowed.'
      write(6,*)
      npoly=0
      return
    endif

    ! Forbid repeated powers.
    do i=2,npoly
      if(any(abs(pow(1:i-1)-pow(i))<2.d0*tol_zero))then
        write(6,*)'Two exponents appear to be identical.'
        write(6,*)
        npoly=0
        return
      endif ! Identical exponents
    enddo ! i

    ! Sort exponents in ascending order
    do i=1,npoly-1
      do j=i+1,npoly
        if(pow(j)<pow(i))then
          t1=pow(i)
          pow(i)=pow(j)
          pow(j)=t1
        endif ! Swap needed
      enddo ! j
    enddo ! i

  END SUBROUTINE ask_poly


  SUBROUTINE show_poly(nxy,have_dx,have_dy,x,y,dx,dy,x0,npoly,pow)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0,pow(npoly)
    DOUBLE PRECISION weight(nxy),a(npoly),da(npoly),t1
    LOGICAL weighted
    INTEGER i

    ! Evaluate fit weights.
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

    ! Print polynomial form.
    write(6,*)'Form of fitting polynomial:'
    write(6,*)'  '//trim(print_poly_sym(npoly,pow,x0))
    write(6,*)

    ! Compute fit.
    call perform_fit(nxy,npoly,x-x0,y,pow,a,weighted,weight,da)

    ! Print fit coefficients.
    write(6,*)'Fit parameters:'
    if(weighted)then
      do i=1,npoly
        write(6,'("   ",a," = ",es24.16," +/- ",es24.16)')'k'//trim(i2s(i)),&
           &a(i),da(i)
      enddo ! i
    else
      do i=1,npoly
        write(6,'("   ",a," = ",es24.16)')'k'//trim(i2s(i)),a(i)
      enddo ! i
    endif ! weighted
    ! Evaluate chi-squared.
    t1=chi_squared(nxy,npoly,x-x0,y,weight,pow,a,weighted)
    write(6,*)'chi^2     = ',t1
    if(nxy>npoly)write(6,*)'chi^2/Ndf = ',t1/dble(nxy-npoly)
    write(6,*)

    ! Write out fitted polynomial.
    ! NB, this is good for pasting into xmgrace, but xmgrace has a string
    ! length limit of 256 characters.
    write(6,*)'Fitted polynomial:'
    write(6,*)'  '//trim(print_poly_num(npoly,pow,a,x0))
    write(6,*)

  END SUBROUTINE show_poly


  SUBROUTINE eval_fit_values_derivs(nxy,have_dx,have_dy,x,y,dx,dy,&
     &itransfx,itransfy,tx,ty,tx0,npoly,pow)
    !--------------------------------------------------!
    ! Evaluate values, first and second derivatives of !
    ! fit at user-requested points.                    !
    !--------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy,npoly
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),&
       &tx(nxy),ty(nxy),tx0,pow(npoly)
    ! Fit parameters.
    DOUBLE PRECISION a(npoly)
    ! Polynomials resulting from operations on fitting polynomial.
    INTEGER op_npoly
    DOUBLE PRECISION op_pow(npoly),op_a(npoly)
    ! Random sampling variables.
    INTEGER nrandom,nx
    DOUBLE PRECISION grid_tx1,grid_tx2
    DOUBLE PRECISION, ALLOCATABLE :: fmean(:),ferr(:),fmean_1s(:),ferr_1s(:),&
       &fmean_2s(:),ferr_2s(:),fmed(:),fskew(:),fkurt(:)
    DOUBLE PRECISION,ALLOCATABLE :: tx_target(:)
    ! Misc local variables.
    CHARACTER(40) set_label,xprint
    CHARACTER(256) fname
    CHARACTER(2048) char2048
    INTEGER i,nderiv,ierr,ix,ipos
    LOGICAL untransformed,eval_relative
    DOUBLE PRECISION t1
    ! Parameters.
    INTEGER,PARAMETER :: DEFAULT_NRANDOM=100000
    INTEGER,PARAMETER :: io=10

    ! Initialize.
    nrandom=DEFAULT_NRANDOM
    nderiv=0
    set_label='f'
    eval_relative=.false.

    ! Loop over user actions.
    do

      ! Print list of actions.
      write(6,*)'-------------'
      write(6,*)'Post-fit menu'
      write(6,*)'-------------'
      write(6,*)'Available actions:'
      if(have_dx.or.have_dy)write(6,*)'[m] Set number of Monte Carlo &
         &samples ['//trim(i2s(nrandom))//']'
      write(6,*)'[d] Set derivative order (0 for value) ['//&
         &trim(i2s(nderiv))//']'
      write(6,*)'[r] Toggle evaluation relative to value at X=0 [',&
         &eval_relative,']'
      write(6,*)'[v] Show fit value at one point x'
      write(6,*)'[V] Show fit value at one point X'
      write(6,*)'[p] Plot fit values at several x'
      write(6,*)'[P] Plot fit values at several X'
      write(6,*)'[q] Return to previous menu'
      write(6,*)

      ! Read user choice.
      write(6,*)'Enter one of the above options:'
      read(5,'(a)',iostat=ierr)char2048
      if(ierr/=0)char2048='q'
      write(6,*)

      ! Perform selected action.
      nx=0
      select case(trim(adjustl(char2048)))
      case('m')
        if(.not.(have_dx.or.have_dy))then
          write(6,*)'Action not available for data without error bars.'
          write(6,*)
          cycle
        endif
        write(6,*)'Enter number of Monte Carlo samples ['//&
           &trim(i2s(nrandom))//']:'
        read(5,'(a)',iostat=ierr)char2048
        if(ierr/=0)return
        write(6,*)
        char2048=adjustl(char2048)
        if(trim(char2048)/='')then
          read(char2048,*,iostat=ierr)i
          if(ierr==0)then
            if(i<10)then
              write(6,*)'Number of random points must be => 10.'
              write(6,*)
            else
              nrandom=i
            endif
          endif
        endif
      case('d')
        write(6,*)'Enter derivative order ['//trim(i2s(nderiv))//']:'
        read(5,'(a)',iostat=ierr)char2048
        if(ierr/=0)return
        write(6,*)
        char2048=adjustl(char2048)
        if(trim(char2048)/='')then
          read(char2048,*,iostat=ierr)i
          if(ierr==0)then
            if(i<0)then
              write(6,*)'Derivative order must be non-negative.'
              write(6,*)
            else
              nderiv=i
              if(nderiv==0)then
                set_label='f'
              elseif(nderiv==1)then
                set_label='df/dX'
              else
                set_label='d'//trim(i2s(nderiv))//'f_dX'//trim(i2s(nderiv))
              endif
            endif
          endif
        endif
      case('r')
        eval_relative=.not.eval_relative
      case('v','V')
        untransformed=trim(char2048)=='v'
        if(untransformed)then
          write(6,*)'Enter (untransformed) x value to evaluate the fit at:'
        else
          write(6,*)'Enter (transformed) X value to evaluate the fit at:'
        endif
        read(5,'(a)',iostat=ierr)char2048
        if(ierr/=0)exit
        write(6,*)
        read(char2048,*,iostat=ierr)t1
        if(ierr/=0)cycle
        nx=1
        allocate(tx_target(nx),stat=ierr)
        if(ierr/=0)call quit('Allocation error.')
        tx_target(1)=t1
        if(untransformed)call scale_transform(1,itransfx,tx_target,tx_target,&
           &.false.)
      case('p')
        untransformed=trim(char2048)=='v'
        if(untransformed)then
          write(6,*)'Enter (untransformed) x values to plot at &
             &[space-separated list of x values,'
          write(6,*)'or <x1>:<x2>:<n> for a grid, or blank for x of source &
             &data]:'
        else
          write(6,*)'Enter (transformed) X values to plot at [space-separated &
             &list of X values,'
          write(6,*)'or <X1>:<X2>:<n> for a grid, or blank for X of source &
             &data]:'
        endif
        nx=0
        read(5,'(a)',iostat=ierr)char2048
        if(ierr/=0)exit
        write(6,*)
        char2048=adjustl(char2048)
        if(len_trim(char2048)==0)then
          nx=nxy
          allocate(tx_target(nx),stat=ierr)
          if(ierr/=0)call quit('Allocation error.')
          tx_target=tx
        else
          ! Parse input string.
          do
            read(char2048,*,iostat=ierr)(t1,i=1,nx+1)
            if(ierr/=0)exit
            nx=nx+1
          enddo
          if(nx>1)then
            allocate(tx_target(nx),stat=ierr)
            if(ierr/=0)call quit('Allocation error.')
            read(char2048,*,iostat=ierr)tx_target
          else
            ! Try to parse grid.
            ipos=scan(char2048,':')
            if(ipos<1)cycle
            read(char2048(1:ipos-1),*,iostat=ierr)grid_tx1
            if(ierr/=0)cycle
            char2048=adjustl(char2048(ipos+1:))
            ipos=scan(char2048,':')
            if(ipos<1)cycle
            read(char2048(1:ipos-1),*,iostat=ierr)grid_tx2
            if(ierr/=0)cycle
            char2048=adjustl(char2048(ipos+1:))
            read(char2048,*,iostat=ierr)i
            if(ierr/=0)cycle
            if(i<2)cycle
            nx=i
            allocate(tx_target(nx),stat=ierr)
            if(ierr/=0)call quit('Allocation error.')
            if(nx==1)then
              tx_target(1)=0.5d0*(grid_tx1+grid_tx2)
            elseif(nx>1)then
              tx_target=(/(grid_tx1+(grid_tx2-grid_tx1)*dble(i-1)/dble(nx-1),&
                 &i=1,nx)/)
            endif
          endif
        endif
        if(nx<1)cycle
        if(untransformed)call scale_transform(nx,itransfx,tx_target,tx_target,&
           &.false.)
      case('q','')
        return
      case default
        return
      end select

      ! Evaluate function at selected points.
      if(nx>0)then

        ! Prepare dump file if needed.
        if(nx>1)then
          fname='polyfit_'//trim(set_label)//'.dat'
          open(unit=io,file=trim(fname),status='replace')
        endif

        ! Choose random sampling or deterministic fit.
        if(have_dx.or.have_dy)then

          allocate(fmean(nx),ferr(nx),fmean_1s(nx),ferr_1s(nx),fmean_2s(nx),&
             &ferr_2s(nx),fmed(nx),fskew(nx),fkurt(nx),stat=ierr)
          if(ierr/=0)call quit('Allocation error (fmean,etc).')
          call eval_fit_monte_carlo(nxy,have_dx,have_dy,x,y,dx,dy,&
             &itransfx,itransfy,tx0,npoly,pow,nrandom,nderiv,eval_relative,&
             &nx,tx_target,fmean,ferr,fmean_1s,ferr_1s,fmean_2s,ferr_2s,&
             &fmed,fskew,fkurt)
          ! Dump.
          if(nx==1)then
            ! Single value is printed to stdout.
            write(xprint,*)tx_target(1)
            write(6,*)'Value of '//trim(set_label)//' at X='//&
               &trim(adjustl(xprint))//':'
            write(6,*)'  Mean            : ',fmean(1),' +/- ',ferr(1)
            write(6,*)'  From 1-sigma CI : ',fmean_1s(1),' +/- ',ferr_1s(1)
            write(6,*)'  From 2-sigma CI : ',fmean_2s(1),' +/- ',ferr_2s(1)
            write(6,*)'  Median          : ',fmed(1)
            write(6,*)'  Skewness        : ',fskew(1)
            write(6,*)'  Kurtosis        : ',fkurt(1)
            write(6,*)
          else
            ! Multiple values are dumped to file.
            write(io,'(a)')'# X  '//trim(set_label)//'  stderr  &
               &1sigma-centre  1sigma-width  2sigma-centre  2sigma-width  &
               &median  skew  kurt'
            do ix=1,nx
              write(io,*)tx_target(ix),fmean(ix),ferr(ix),&
                 &fmean_1s(ix),ferr_1s(ix),fmean_2s(ix),ferr_2s(ix),&
                 &fmed(ix),fskew(ix),fkurt(ix)
            enddo ! ix
          endif
          deallocate(fmean,ferr,fmean_1s,ferr_1s,fmean_2s,ferr_2s,fmed,fskew,&
             &fkurt)

        else ! .not.(have_dx.or.have_dy)

          ! Perform single fit.
          call perform_fit(nxy,npoly,tx-tx0,ty,pow,a,.false.)
          if(nx==1)then
            ! Single value is printed to stdout.
            op_npoly=npoly
            op_pow=pow
            op_a=a
            do i=1,nderiv
              call deriv_poly(op_npoly,op_pow,op_a)
            enddo ! i
            t1=eval_poly(op_npoly,op_pow,op_a,tx_target(1)-tx0)
            write(6,*)trim(set_label)//' = ',t1
            write(6,*)
          else
            ! Multiple values are dumped to file.
            write(io,'(a)')'# x  '//trim(set_label)
            op_npoly=npoly
            op_pow=pow
            op_a=a
            do i=1,nderiv
              call deriv_poly(op_npoly,op_pow,op_a)
            enddo ! i
            do ix=1,nx
              t1=eval_poly(op_npoly,op_pow,op_a,tx_target(ix)-tx0)
              write(io,*)tx_target(ix),t1
            enddo ! ix
          endif

        endif ! have_dx.or.have_dy

        ! Clean up.
        if(nx>1)then
          write(6,*)'Data written to "'//trim(fname)//'".'
          write(6,*)
          close(io)
        endif
        deallocate(tx_target)
      endif

    enddo

  END SUBROUTINE eval_fit_values_derivs


  ! COMPUTATION ROUTINES


  SUBROUTINE perform_fit(nxy,npoly,x,y,pow,a,weighted,weight,da)
    !----------------------------------------------------!
    ! Perform least-squares fit of (weighted) xy data to !
    ! polynomial described by pow(1:npoly).              !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),pow(npoly)
    DOUBLE PRECISION,INTENT(out) :: a(npoly)
    LOGICAL,INTENT(in) :: weighted
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: weight(nxy)
    DOUBLE PRECISION,INTENT(out),OPTIONAL :: da(npoly)
    INTEGER i,j,ipiv(npoly),lwork,ierr
    DOUBLE PRECISION M(npoly,npoly),Minv(npoly,npoly),c(npoly)
    DOUBLE PRECISION,ALLOCATABLE :: work(:)

    ! Construct c vector and M matrix.
    if(weighted)then
      do i=1,npoly
        c(i)=sum(weight(:)*y(:)*x(:)**pow(i))
        M(i,i)=sum(weight(:)*x(:)**(pow(i)+pow(i)))
        do j=i+1,npoly
          M(j,i)=sum(weight(:)*x(:)**(pow(j)+pow(i)))
          M(i,j)=M(j,i)
        enddo ! j
      enddo ! i
    else
      do i=1,npoly
        c(i)=sum(y(:)*x(:)**pow(i))
        M(i,i)=sum(x(:)**(pow(i)+pow(i)))
        do j=i+1,npoly
          M(j,i)=sum(x(:)**(pow(j)+pow(i)))
          M(i,j)=M(j,i)
        enddo ! j
      enddo ! i
    endif ! weighted

    ! Invert M.
    Minv=M
    allocate(work(1))
    lwork=-1
    call dsytrf('L',npoly,Minv,npoly,ipiv,work,lwork,ierr)
    if(ierr/=0)call quit('DSYTRF error on workspace query call.')
    lwork=nint(work(1))
    deallocate(work)
    allocate(work(lwork),stat=ierr)
    if(ierr/=0)call quit('Allocation error (work).')
    call dsytrf('L',npoly,Minv,npoly,ipiv,work,lwork,ierr)
    if(ierr/=0)call quit('DSYTRF error.')
    deallocate(work)
    allocate(work(npoly),stat=ierr)
    if(ierr/=0)call quit('Allocation error (work).')
    call dsytri('L',npoly,Minv,npoly,ipiv,work,ierr)
    if(ierr/=0)call quit('DSYTRI error.')
    deallocate(work)

    ! Complete Minv and evaluate coefficients.
    do i=1,npoly-1
      do j=i+1,npoly
        Minv(i,j)=Minv(j,i)
      enddo ! j
    enddo ! i
    a=matmul(Minv,c)

    ! Evaluate standard errors if data are weighted.
    if(weighted)then
      do j=1,npoly
        da(j)=sqrt(Minv(j,j))
      enddo ! j
    elseif(present(da))then
      da=0.d0
    endif ! weighted

  END SUBROUTINE perform_fit


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


  SUBROUTINE eval_fit_monte_carlo(nxy,have_dx,have_dy,x,y,dx,dy,&
     &itransfx,itransfy,tx0,npoly,pow,nrandom,nderiv,eval_relative,&
     &nx,tx_target,fmean,ferr,fmean_1s,ferr_1s,fmean_2s,ferr_2s,fmed,&
     &fskew,fkurt,amean,aerr)
    !------------------------------------------------------!
    ! Perform Monte Carlo sampling of data space to obtain !
    ! fit values or derivatives at specified points with   !
    ! error bars.                                          !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy,npoly,nrandom,nderiv,nx
    LOGICAL,INTENT(in) :: have_dx,have_dy,eval_relative
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),tx0,&
       &pow(npoly),tx_target(nx)
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: fmean(nx),ferr(nx),&
       &fmean_1s(nx),ferr_1s(nx),fmean_2s(nx),ferr_2s(nx),fmed(nx),&
       &fskew(nx),fkurt(nx),amean(npoly),aerr(npoly)
    ! Transformed data.
    DOUBLE PRECISION tx(nxy),ty(nxy)
    ! Polynomials resulting from operations on fitting polynomial.
    INTEGER op_npoly
    DOUBLE PRECISION op_pow(npoly),op_a(npoly)
    ! Random sampling arrays.
    DOUBLE PRECISION f_array(nrandom,nx),w_vector(nrandom),&
       &a_array(nrandom,npoly)
    DOUBLE PRECISION ran_x(nxy),ran_y(nxy),ran_a(npoly)
    ! Distribution analysis.
    DOUBLE PRECISION var,skew,kurt,f2s_lo,f1s_lo,f1s_hi,f2s_hi
    ! Misc.
    INTEGER ipoly,ideriv,irandom,ix
    DOUBLE PRECISION t1,da(npoly),dtx(nxy),dty(nxy),weight(nxy),f0
    LOGICAL weighted

    if(MONTE_CARLO_FIT_WEIGHTS)then
      call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
      call scale_transform(nxy,itransfy,y,ty,have_dy,dy,dty)
      if(have_dx.and.have_dy)then
        weight=1.d0/(dtx*dty)**2
      elseif(have_dx)then
        weight=1.d0/dtx**2
      elseif(have_dy)then
        weight=1.d0/dty**2
      else
        weight=1.d0
      endif
      weighted=have_dx.or.have_dy
    else
      weight=1.d0
      weighted=.false.
    endif

    ! Initialize.
    ran_x=x
    ran_y=y
    f0=0.d0

    ! Loop over random points.
    do irandom=1,nrandom
      if(have_dx)ran_x=x+gaussian_random_number(dx)
      if(have_dy)ran_y=y+gaussian_random_number(dy)
      call scale_transform(nxy,itransfx,ran_x,tx,.false.)
      call scale_transform(nxy,itransfy,ran_y,ty,.false.)
      call perform_fit(nxy,npoly,tx-tx0,ty,pow,ran_a,weighted,weight,da)
      a_array(irandom,1:npoly)=ran_a(1:npoly)
      if(MONTE_CARLO_CHI2_WEIGHTS)then
        w_vector(irandom)=1.d0/&
           &chi_squared(nxy,npoly,tx-tx0,ty,w_vector,pow,ran_a,.false.)
      else
        w_vector(irandom)=1.d0
      endif
      op_npoly=npoly
      op_pow=pow
      op_a=ran_a
      do ideriv=1,nderiv
        call deriv_poly(op_npoly,op_pow,op_a)
      enddo ! ideriv
      if(eval_relative)f0=eval_poly(op_npoly,op_pow,op_a,-tx0)
      do ix=1,nx
        t1=eval_poly(op_npoly,op_pow,op_a,tx_target(ix)-tx0)
        f_array(irandom,ix)=t1-f0
      enddo ! ix
    enddo ! irandom

    ! Return coefficients.
    do ipoly=1,npoly
      call characterize_dist(nrandom,a_array(:,ipoly),w_vector,mean=t1,var=var)
      if(present(amean))amean(ipoly)=t1
      if(present(aerr))aerr(ipoly)=sqrt(var)
    enddo ! ipoly

    ! Evaluate statistics.
    do ix=1,nx
      call characterize_dist(nrandom,f_array(:,ix),w_vector,mean=t1,var=var,&
         &skew=skew,kurt=kurt)
      if(present(fmean))fmean(ix)=t1
      if(present(ferr))ferr(ix)=sqrt(var)
      if(present(fskew))fskew(ix)=skew
      if(present(fkurt))fkurt(ix)=kurt
      if(present(fmed))fmed(ix)=median(nrandom,f_array(:,ix))
      if(present(fmean_1s).or.present(ferr_1s))then
        f1s_lo=find_pth_smallest(nint(dble(nrandom)*0.158655254d0),&
           &nrandom,f_array(:,ix))
        f1s_hi=find_pth_smallest(nint(dble(nrandom)*0.841344746d0),&
           &nrandom,f_array(:,ix))
        if(present(fmean_1s))fmean_1s(ix)=0.5d0*(f1s_lo+f1s_hi)
        if(present(ferr_1s))ferr_1s(ix)=0.5d0*(f1s_hi-f1s_lo)
      endif
      if(present(fmean_2s).or.present(ferr_2s))then
        f2s_lo=find_pth_smallest(nint(dble(nrandom)*0.022750131948d0),&
           &nrandom,f_array(:,ix))
        f2s_hi=find_pth_smallest(nint(dble(nrandom)*0.977249868052d0),&
           &nrandom,f_array(:,ix))
        if(present(fmean_2s))fmean_2s(ix)=0.5d0*(f2s_lo+f2s_hi)
        if(present(ferr_2s))ferr_2s(ix)=0.25d0*(f2s_hi-f2s_lo)
      endif
    enddo ! ix

  END SUBROUTINE eval_fit_monte_carlo


  ! INPUT FILE HANDLING ROUTINES.


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
      read(5,'(a)',iostat=ierr)fname
      if(ierr/=0)call quit()
      fname=adjustl(fname)
      inquire(file=trim(fname),exist=file_exists)
      if(.not.file_exists)call quit('File does not exist.')
      exit
    enddo
    write(6,*)
  END SUBROUTINE get_file


  SUBROUTINE check_file(fname,nline,ncolumn)
    !----------------------------------------------------!
    ! Check file contains data, and return the number of !
    ! data lines and the number of data columns in it.   !
    !----------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(out) :: nline,ncolumn
    CHARACTER(1024) char1024
    INTEGER i,ipos,pline,ierr
    DOUBLE PRECISION t1
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)call quit('Error opening "'//trim(fname)//'".')

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

    ! Close file.
    close(io)

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
    write(6,*)'Enter the column indices for x, y, dx, dy:'
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
    write(6,*)

  END SUBROUTINE get_columns


  SUBROUTINE read_file(fname,ncolumn,icol_x,icol_y,icol_dx,icol_dy,nxy,&
     &x,y,dx,dy)
    !-------------------------------------!
    ! Read in the data in the input file. !
    !-------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ncolumn,icol_x,icol_y,icol_dx,icol_dy,nxy
    CHARACTER(*),INTENT(in) :: fname
    DOUBLE PRECISION,INTENT(out) :: x(nxy),y(nxy),dx(nxy),dy(nxy)
    DOUBLE PRECISION tvec(ncolumn)
    INTEGER i,ierr,ipos,pline
    CHARACTER(1024) line
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)call quit('Error opening "'//trim(fname)//'".')

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

    ! Close file.
    close(io)

  END SUBROUTINE read_file


  ! POLYNOMIAL HANDLING UTILITIES.


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
      xstr='X'
    else
      xstr='(X-X0)'
    endif
    polystr='Y ='
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


  SUBROUTINE deriv_poly(npoly1,pow1,a1,npoly,pow,a)
    !-------------------------------------------------!
    ! Differentiate the polynomial given by npoly1,   !
    ! pow1 and a1 and return the resulting polynomial !
    ! in npoly, pow and a if present, or in-place     !
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
    ! npoly, pow and a if present, or in-place       !
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


  ! GENERIC NUMERICAL UTILITIES.


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


  SUBROUTINE characterize_dist(M,A,w,mean,stderr,err_stderr,var,skew,kurt)
    !-------------------------------------------------!
    ! Evaluate the variance, skewness and kurtosis of !
    ! a data set.                                     !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M),w(M)
    DOUBLE PRECISION, INTENT(out), OPTIONAL :: mean,stderr,err_stderr,var,&
       &skew,kurt
    DOUBLE PRECISION V1,V2,V3,V4,m1,m2,m3,m4,K2,K3,K4,num1,num2,denom
    ! Initialize.
    if(present(mean))mean=0.d0
    if(present(stderr))stderr=0.d0
    if(present(var))var=0.d0
    if(present(skew))skew=0.d0
    if(present(kurt))kurt=0.d0
    ! Compute mean.
    if(present(mean).or.present(stderr).or.present(err_stderr).or.&
       &present(var).or.present(skew).or.present(kurt))then
      V1=sum(w)
      m1=sum(w*A)/V1
      if(present(mean))mean=m1
    endif
    if(M<2)return
    ! Compute variance.
    if(present(stderr).or.present(err_stderr).or.present(var).or.&
       &present(skew).or.present(kurt))then
      V2=sum(w**2)
      m2=sum(w*(A-m1)**2)/V1
      K2=m2*(V1*V1/(V1*V1-V2))
      if(present(stderr))stderr=sqrt(max(0.d0,K2/dble(M)))
      if(present(err_stderr))err_stderr=sqrt(max(0.d0,&
         &0.5d0*K2/(dble(M)*dble(M-1))))
      if(present(var))var=K2
    endif
    if(K2<=0.d0)return
    if(M<3)return
    ! Compute skewness.
    if(present(skew).or.present(kurt))then
      V3=sum(w**3)
      m3=sum(w*(A-m1)**3)/V1
      K3=m3*(V1*V1*V1/(V1*V1*V1-3.d0*V1*V2+2.d0*V3))
      if(present(skew))skew=K3/K2**1.5d0
    endif
    if(M<4)return
    ! Compute kurtosis.
    if(present(kurt))then
      V4=sum(w**4)
      m4=sum(w*(A-m1)**4)/V1
      num1=V1*V1*(V1**4-4.d0*V1*V3+3.d0*V2*V2)
      num2=3.d0*V1*V1*(V1**4-2.d0*V1*V1*V2+4.d0*V1*V3-3.d0*V2*V2)
      denom=(V1*V1-V2)*(V1**4-6.d0*V1*V1*V2+8.d0*V1*V3+3.d0*V2*V2-6.d0*V4)
      K4=(m4*(num1/denom)-m2*m2*(num2/denom))
      kurt=K4/(K2*K2)
    endif
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


  SUBROUTINE isort(n,x,indx)
    !--------------------------------------------------------!
    ! Perform insertion sort on a real vector X(1:N) so that !
    ! X(I)<X(J) if I<J.  This sorting algorithm typically    !
    ! costs ~ N^2, but is stable (preserves the order of     !
    ! entries with same X) and becomes order N when X(:) is  !
    ! nearly sorted.                                         !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    DOUBLE PRECISION,INTENT(in) :: x(n)
    INTEGER,INTENT(inout) :: indx(n)
    INTEGER i,j
    DOUBLE PRECISION xi
    ! Initialize.
    indx(1:n)=(/(i,i=1,n)/)
    ! Loop over elements from 2:n.
    do i=2,n
      xi=x(indx(i))
      ! Move i-th element upwards until it is greater than previous,
      ! at which point the first I-th elements will be sorted.
      do j=i-1,1,-1
        if(xi>=x(indx(j)))exit
        call iswap1(indx(j),indx(j+1))
      enddo !j
    enddo ! i
  END SUBROUTINE isort


  SUBROUTINE iswap1(i,j)
    !--------------------!
    ! Swap two integers. !
    !--------------------!
    IMPLICIT NONE
    INTEGER,INTENT(inout) :: i,j
    INTEGER k
    k=i
    i=j
    j=k
  END SUBROUTINE iswap1


  ! MISC UTILITIES.


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
      i2s='-'//adjustl(i2s(2:12))
    else
      i2s=adjustl(i2s)
    endif ! n<0
  END FUNCTION i2s


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
