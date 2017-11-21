PROGRAM polyfit
  !-----------------------------------------------------------!
  ! POLYFIT                                                   !
  ! =======                                                   !
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
    ! Fitting function.
    INTEGER npoly,itransfx,itransfy,i
    DOUBLE PRECISION x0
    DOUBLE PRECISION,ALLOCATABLE :: pow(:)
    ! Misc local variables.
    INTEGER ierr
    CHARACTER(2048) char2048

    ! Initialize.
    x0=0.d0
    itransfx=0
    itransfy=0
    npoly=0

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
      write(6,'(1x,a,es11.4,a)')'  [0] Set fit centre [X0=',x0,']'
      write(6,*)'* Data analysis and pre-fit assessment:'
      write(6,*)'  [r] Show basic data-range statistics'
      write(6,*)'  [e] Show expansion-order analysis'
      write(6,*)'  [d] Show dual number of points/expansion-order analysis'
      write(6,*)'  [p] Plot fits at multiple expansion orders'
      write(6,*)'* Fitting:'
      if(npoly<1)then
        write(6,*)'  [f] Set form and perform fit [form not yet chosen]'
        write(6,*)'* Post-fit analysis [must set fit form first]:'
      else
        write(6,*)'  [f] Set form and perform fit ['//&
           &trim(print_poly_sym(npoly,pow,x0))//']'
        write(6,*)'* Post-fit analysis:'
      endif
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
      case('0')
        call ask_centre(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy,x0)
      case('f')
        call ask_poly(nxy,char2048,i)
        if(i<1)cycle
        npoly=i
        if(allocated(pow))deallocate(pow)
        allocate(pow(npoly),stat=ierr)
        if(ierr/=0)call quit('Allocation error.')
        call parse_poly(char2048,npoly,pow)
        if(npoly<1)then
          deallocate(pow)
          cycle
        endif
        call show_poly(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy,x0,&
           &npoly,pow)
      case('r')
        call show_statistics(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy)
      case('e')
        call show_exporder_assessment(nxy,have_dx,have_dy,x,y,dx,dy,&
           &itransfx,itransfy,x0,i,char2048)
        if(i>0)then
          npoly=i
          if(allocated(pow))deallocate(pow)
          allocate(pow(npoly),stat=ierr)
          if(ierr/=0)call quit('Allocation error.')
          call parse_poly(char2048,npoly,pow)
          call show_poly(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy,x0,&
             &npoly,pow)
        endif
      case('d')
        call show_dual_assessment(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,&
           &itransfy,x0)
      case('p')
        call plot_exporder(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy,x0)
      case('v')
        if(npoly<1)then
          write(6,*)'Must choose fit form first!'
          write(6,*)
          cycle
        endif
        call eval_fit_values_derivs(nxy,have_dx,have_dy,x,y,dx,dy,&
           &itransfx,itransfy,x0,npoly,pow)
      case('q','')
        call quit()
      case default
        call quit()
      end select

    enddo

    ! Clean up.
    if(allocated(pow))deallocate(pow)

  END SUBROUTINE user_interaction


  SUBROUTINE show_statistics(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy)
    !-----------------------------------------------------!
    ! Show min/centre/mean/median/max of x, y, dx and dy. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy)
    DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy)
    DOUBLE PRECISION vmin,vmax,vcentre,vmean,vmedian

    ! Apply scale transformations.
    call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    call scale_transform(nxy,itransfy,y,ty,have_dy,dy,dty)

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


  SUBROUTINE ask_centre(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy,x0)
    !------------------------------------------!
    ! Ask the user to set offsets for the fit. !
    !------------------------------------------!
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy)
    DOUBLE PRECISION,INTENT(inout) :: x0
    CHARACTER(2048) char2048
    DOUBLE PRECISION t1,tx(nxy),ty(nxy),dtx(nxy),dty(nxy)
    INTEGER i,ierr

    ! Evaluate transformed variables.
    call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    call scale_transform(nxy,itransfy,y,ty,have_dy,dy,dty)

    ! Get fit centre X0.
    write(6,*)'Current fit centre X0 = ',x0
    write(6,*)'Enter new centre or left/mean/median/centre/right or min/max:'
    read(5,'(a)',iostat=ierr)char2048
    write(6,*)
    if(ierr==0)then
      char2048=adjustl(char2048)
      select case(trim(char2048))
      case('left')
        x0=minval(tx)
      case('right')
        x0=maxval(tx)
      case('min')
        do i=1,nxy
          if(abs(ty(i)-minval(ty))<tol_zero)then
            x0=tx(i)
            exit
          endif
        enddo
      case('max')
        do i=1,nxy
          if(abs(ty(i)-maxval(ty))<tol_zero)then
            x0=tx(i)
            exit
          endif
        enddo
      case('centre')
        x0=.5d0*(minval(tx)+maxval(tx))
      case('mean')
        x0=sum(tx)/dble(nxy)
      case('median')
        x0=median(nxy,tx)
      case('')
        ierr=1
      case default
        read(char2048,*,iostat=ierr)t1
        if(ierr==0)then
          x0=t1
        else
        endif
      end select
    endif

    ! Report.
    if(ierr==0)then
      write(6,*)'Set X0 = ',x0
    else
      write(6,*)'Keeping X0 = ',x0
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


  SUBROUTINE show_exporder_assessment(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,&
     &itransfy,x0,npoly_out,char2048)
    !--------------------------------------------------------!
    ! Perform an assessment of the performance of polynomial !
    ! fits as a function of expansion order.                 !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0
    INTEGER,INTENT(inout) :: npoly_out
    CHARACTER(*),INTENT(inout) :: char2048
    DOUBLE PRECISION,ALLOCATABLE :: chi2_vector(:),rmsc_100_vector(:),&
       &rmsc_150_vector(:),rmsc_200_vector(:),pow(:),a(:),da(:),op_pow(:),&
       &op_a(:)
    DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy),weight(nxy),txrange
    INTEGER ntest,npoly,op_npoly,i,np,ierr
    LOGICAL weighted

    ! Initialize.
    npoly_out=0
    char2048=''

    ! Evaluate transformed variables and fit weights.
    call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    tx=tx-x0
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
      call calc_parameters(nxy,npoly,tx,ty,pow,a,weighted,weight,da)
      ! Evaluate chi^2.
      chi2_vector(npoly)=chi_squared(nxy,npoly,tx,ty,weight,pow,a,weighted)
      ! Evaluate root mean square curvature of fit as measure of smoothness.
      np=npoly
      txrange=maxval(tx)-minval(tx)
      call deriv_poly(np,pow,a,op_npoly,op_pow,op_a)
      call deriv_poly(op_npoly,op_pow,op_a)
      call square_poly(op_npoly,op_pow,op_a)
      call int_poly(op_npoly,op_pow,op_a)
      rmsc_100_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,maxval(tx))-&
         &eval_poly(op_npoly,op_pow,op_a,minval(tx))&
         &)/txrange))
      rmsc_150_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,maxval(tx)+0.25d0*txrange)-&
         &eval_poly(op_npoly,op_pow,op_a,minval(tx)-0.25d0*txrange)&
         &)/(1.5d0*txrange)))
      rmsc_200_vector(npoly)=&
         &sqrt(abs((&
         &eval_poly(op_npoly,op_pow,op_a,maxval(tx)+0.5d0*txrange)-&
         &eval_poly(op_npoly,op_pow,op_a,minval(tx)-0.5d0*txrange)&
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
    write(6,*)'Enter order to choose fit form [empty to skip]:'
    read(5,'(a)',iostat=ierr)char2048
    if(ierr==0)then
      write(6,*)
      if(len_trim(char2048)>0)then
        read(char2048,*,iostat=ierr)i
        if(ierr==0)then
          if(i>=0.and.i<=ntest)then
            npoly_out=i+1
            char2048=''
            do i=0,npoly_out-1
              char2048=trim(char2048)//' '//trim(i2s(i))
            enddo ! i
          endif
        endif
      endif
    endif

    ! Destroy work arrays.
    deallocate(chi2_vector,rmsc_100_vector,rmsc_150_vector,rmsc_200_vector)

  END SUBROUTINE show_exporder_assessment


  SUBROUTINE show_dual_assessment(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,&
     &itransfy,x0)
    !-----------------------------------------------------------!
    ! Perform an "extrapolation" assessment in which the number !
    ! of points included in the fit and the expansion order are !
    ! simultaneously optimized.                                 !
    !-----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0
    DOUBLE PRECISION, ALLOCATABLE :: chi2_array(:,:),pow(:),a(:),da(:)
    DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy),weight(nxy),&
       &rtx(nxy),rty(nxy),rweight(nxy),xtarget,min_chi2,t1,&
       &rx(nxy),ry(nxy),rdx(nxy),rdy(nxy)
    INTEGER i,ntest,npoly,rnxy,nderiv,ierr,indx(nxy),imin_chi2,nrandom,ibest_f
    LOGICAL mask(nxy),weighted
    CHARACTER(20) drop_by,drop_criterion
    DOUBLE PRECISION fmean(1),ferr(1),best_f,best_df

    ! Initialize.
    nderiv=0
    xtarget=0.d0
    nrandom=10000
    drop_by='X'
    drop_criterion='largest'

    ! Evaluate transformed variables and fit weights.
    call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    tx=tx-x0
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
    allocate(chi2_array(nxy,ntest),stat=ierr)
    if(ierr/=0)call quit('Allocation problem (chi2_array).')

    ! Loop over expansion orders.
    ibest_f=0
    best_f=0.d0
    best_df=0.d0
    do npoly=1,ntest
      write(6,*)'At expansion order '//trim(i2s(npoly-1))//':'
      ! Allocate work arrays.
      allocate(pow(npoly),a(npoly),da(npoly),stat=ierr)
      if(ierr/=0)call quit('Allocation error.')
      do i=1,npoly
        pow(i)=dble(i-1)
      enddo ! i
      ! Loop over number of points in fit.
      imin_chi2=0
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
        call calc_parameters(rnxy,npoly,rtx,rty,pow,a,weighted,rweight,da)
        ! Evaluate chi^2.
        t1=chi_squared(rnxy,npoly,rtx,rty,rweight,pow,a,weighted)/&
           &dble(rnxy-npoly)
        chi2_array(rnxy,npoly)=t1
        if(imin_chi2==0.or.t1<min_chi2)then
          imin_chi2=rnxy
          min_chi2=t1
        endif
        ! Report.
        write(6,*)'  NXY = '//trim(i2s(rnxy))//' -> chi^2/Ndf = ',t1
      enddo ! rnxy
      if(imin_chi2>npoly+1)then
        rnxy=imin_chi2
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
           &itransfx,itransfy,x0,npoly,pow,nrandom,nderiv,1,(/xtarget/),&
           &fmean,ferr)
        write(6,*)'  At NXY = '//trim(i2s(imin_chi2))//': ',fmean(1),' +/- ',&
           &ferr(1)
        if(ibest_f==0.or.abs(fmean(1)-best_f)>&
           &2.d0*sqrt(ferr(1)**2+best_df**2))then
          ibest_f=npoly
          best_f=fmean(1)
          best_df=ferr(1)
        endif
      endif
      write(6,*)
      deallocate(pow,a,da)
    enddo ! npoly
    write(6,*)'Best expansion order: '//trim(i2s(ibest_f-1))
    write(6,*)'  at which we get ',best_f,' +/- ',best_df
    write(6,*)

    ! Clean up.
    deallocate(chi2_array)

  END SUBROUTINE show_dual_assessment


  SUBROUTINE plot_exporder(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy,x0)
    !-----------------------------------------------------!
    ! Plot polynomial fits of different expansion orders. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0
    DOUBLE PRECISION,ALLOCATABLE :: pow(:),a(:),da(:)
    DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy),weight(nxy),txrange,&
       &txplot,dtxplot
    INTEGER npoly,i,ierr
    LOGICAL weighted
    ! Parameters.
    INTEGER,PARAMETER :: npoint=1000
    INTEGER,PARAMETER :: io=10

    ! Evaluate transformed variables and fit weights.
    call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    tx=tx-x0
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

    ! Prepare to write plot.
    txrange=maxval(tx)-minval(tx)
    open(unit=io,file='poly_orders.dat',status='replace')

    ! Loop over expansion orders.
    do npoly=1,min(nxy-1,9)
      ! Prepare for fit.
      allocate(pow(npoly),a(npoly),da(npoly),stat=ierr)
      if(ierr/=0)call quit('Allocation error.')
      do i=1,npoly
        pow(i)=dble(i-1)
      enddo ! i
      ! Perform fit.
      call calc_parameters(nxy,npoly,tx,ty,pow,a,weighted,weight,da)
      ! Dump plot.
      txplot=minval(tx)-0.5d0*txrange
      dtxplot=(txrange*2.d0)/dble(npoint)
      do i=0,npoint
        write(io,*)txplot+x0,eval_poly(npoly,pow,a,txplot)
        txplot=txplot+dtxplot
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


  SUBROUTINE ask_poly(nxy,char2048,npoly)
    !----------------------------------------------------------!
    ! Ask user for string containing exponents for polynomial, !
    ! and return the string and the number of terms.           !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    CHARACTER(*),INTENT(inout) :: char2048
    INTEGER,INTENT(inout) :: npoly
    INTEGER i1,i2,i,ipos,ierr
    DOUBLE PRECISION t1

    ! Ask user for exponents.
    npoly=0
    write(6,*)'Enter exponents in fitting polynomial (space-separated &
       &list, or i:j, or :j):'
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
      ! Reconstruct string for later parsing.
      char2048=''
      do i=i1,i2
        char2048=trim(char2048)//' '//trim(i2s(i))
      enddo ! i
    endif

    ! Check number exponents is sensible.
    if(npoly<1)then
      write(6,*)'Interpolating polynomial must have at least one term.'
      npoly=0
      return
    endif
    if(npoly>nxy)then
      write(6,*)'Number of terms cannot exceed number of points.'
      npoly=0
      return
    endif

  END SUBROUTINE ask_poly


  SUBROUTINE parse_poly(char2048,npoly,pow)
    !---------------------------------------------------------!
    ! Parse the string contraining the npoly exponents in the !
    ! fitting polynomial.                                     !
    !---------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: char2048
    INTEGER,INTENT(inout) :: npoly
    DOUBLE PRECISION,INTENT(inout) :: pow(npoly)
    INTEGER i,j,idiv
    DOUBLE PRECISION t1

    ! Read string.
    read(char2048,*)pow(1:npoly)

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
      npoly=0
      return
    endif

    ! Forbid repeated powers.
    do i=2,npoly
      if(any(abs(pow(1:i-1)-pow(i))<2.d0*tol_zero))then
        write(6,*)'Two exponents appear to be identical.'
        npoly=0
        return
      endif ! Identical exponents
    enddo ! i

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

  END SUBROUTINE parse_poly


  SUBROUTINE show_poly(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,itransfy,x0,&
     &npoly,pow)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy,npoly
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0,pow(npoly)
    DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy),weight(nxy),a(npoly),&
       &da(npoly),t1
    LOGICAL weighted
    INTEGER i

    ! Evaluate transformed variables and fit weights.
    call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
    tx=tx-x0
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

    ! Print polynomial form.
    write(6,*)'Form of fitting polynomial:'
    write(6,*)'  '//trim(print_poly_sym(npoly,pow,x0))
    write(6,*)

    ! Compute fit.
    call calc_parameters(nxy,npoly,tx,ty,pow,a,weighted,weight,da)

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
    t1=chi_squared(nxy,npoly,tx,ty,weight,pow,a,weighted)
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
     &itransfx,itransfy,x0,npoly,pow)
    !--------------------------------------------------!
    ! Evaluate values, first and second derivatives of !
    ! fit at user-requested points.                    !
    !--------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy,npoly
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0,pow(npoly)
    ! Transformed data.
    DOUBLE PRECISION tx(nxy),ty(nxy)
    ! Fit parameters.
    DOUBLE PRECISION a(npoly)
    ! Polynomials resulting from operations on fitting polynomial.
    INTEGER op_npoly
    DOUBLE PRECISION op_pow(npoly),op_a(npoly)
    ! Random sampling variables.
    INTEGER nrandom,nx
    DOUBLE PRECISION grid_x1,grid_x2
    DOUBLE PRECISION, ALLOCATABLE :: fmean(:),ferr(:),fmean_1s(:),ferr_1s(:),&
       &fmean_2s(:),ferr_2s(:),fmed(:),fskew(:),fkurt(:)
    DOUBLE PRECISION,ALLOCATABLE :: x_target(:)
    ! Misc local variables.
    CHARACTER(40) set_label,xprint
    CHARACTER(256) fname
    CHARACTER(2048) char2048
    INTEGER i,nderiv,ierr,ix,ipos
    DOUBLE PRECISION t1
    ! Parameters.
    INTEGER,PARAMETER :: DEFAULT_NRANDOM=100000
    INTEGER,PARAMETER :: io=10

    ! Initialize.
    nrandom=DEFAULT_NRANDOM
    nderiv=0
    set_label='f'

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
      write(6,*)'[v] Show fit value at one point'
      write(6,*)'[p] Plot fit values at a range of points'
      write(6,*)'[q] Return to previous menu'
      write(6,*)

      ! Read user choice.
      write(6,*)'Enter one of the above options:'
      read(5,'(a)',iostat=ierr)char2048
      if(ierr/=0)char2048='q'
      write(6,*)

      ! Perform selected action.
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
                set_label='d'//trim(i2s(nderiv))//'f/dX'//trim(i2s(nderiv))
              endif
            endif
          endif
        endif
      case('v')
        write(6,*)'Enter (transformed) X value to evaluate the fit at:'
        read(5,'(a)',iostat=ierr)char2048
        if(ierr/=0)exit
        write(6,*)
        read(char2048,*,iostat=ierr)t1
        if(ierr/=0)cycle
        nx=1
        allocate(x_target(nx),stat=ierr)
        if(ierr/=0)call quit('Allocation problem.')
        x_target(1)=t1
      case('p')
        write(6,*)'Enter (transformed) X values to plot at [space-separated &
           &list of X values,'
        write(6,*)'or <x1>:<x2>:<n> for a grid, or blank for X of source data]:'
        nx=0
        read(5,'(a)',iostat=ierr)char2048
        if(ierr/=0)exit
        write(6,*)
        char2048=adjustl(char2048)
        if(len_trim(char2048)==0)then
          nx=nxy
          allocate(x_target(nx),stat=ierr)
          if(ierr/=0)call quit('Allocation problem.')
          x_target=x+x0
        else
          ! Parse input string.
          do
            read(char2048,*,iostat=ierr)(t1,i=1,nx+1)
            if(ierr/=0)exit
            nx=nx+1
          enddo
          if(nx>1)then
            allocate(x_target(nx),stat=ierr)
            if(ierr/=0)call quit('Allocation problem.')
            read(char2048,*,iostat=ierr)x_target
          else
            ! Try to parse grid.
            ipos=scan(char2048,':')
            if(ipos<1)cycle
            read(char2048(1:ipos-1),*,iostat=ierr)grid_x1
            if(ierr/=0)cycle
            char2048=adjustl(char2048(ipos+1:))
            ipos=scan(char2048,':')
            if(ipos<1)cycle
            read(char2048(1:ipos-1),*,iostat=ierr)grid_x2
            if(ierr/=0)cycle
            char2048=adjustl(char2048(ipos+1:))
            read(char2048,*,iostat=ierr)i
            if(ierr/=0)cycle
            if(i<2)cycle
            nx=i
            allocate(x_target(nx),stat=ierr)
            if(ierr/=0)call quit('Allocation problem.')
            if(nx==1)then
              x_target(1)=0.5d0*(grid_x1+grid_x2)
            elseif(nx>1)then
              x_target=(/(grid_x1+(grid_x2-grid_x1)*dble(i-1)/dble(nx-1),&
                 &i=1,nx)/)
            endif
          endif
        endif
        if(nx<1)cycle
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
          if(ierr/=0)call quit('Allocation problem (fmean,etc).')
          call eval_fit_monte_carlo(nxy,have_dx,have_dy,x,y,dx,dy,&
             &itransfx,itransfy,x0,npoly,pow,nrandom,nderiv,nx,x_target,&
             &fmean,ferr,fmean_1s,ferr_1s,fmean_2s,ferr_2s,fmed,fskew,fkurt)
          ! Dump.
          if(nx==1)then
            ! Single value is printed to stdout.
            write(xprint,*)x_target(1)
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
            write(io,'(a)')'# x  '//trim(set_label)//'  stderr  &
               &1sigma-centre  1sigma-width  2sigma-centre  2sigma-width  &
               &median  skew  kurt'
            do ix=1,nx
              write(io,*)x_target(ix),fmean(ix),ferr(ix),&
                 &fmean_1s(ix),ferr_1s(ix),fmean_2s(ix),ferr_2s(ix),&
                 &fmed(ix),fskew(ix),fkurt(ix)
            enddo ! ix
          endif
          deallocate(fmean,ferr,fmean_1s,ferr_1s,fmean_2s,ferr_2s,fmed,fskew,&
             &fkurt)

        else ! .not.(have_dx.or.have_dy)

          ! Perform single fit.
          call scale_transform(nxy,itransfx,x,tx,.false.)
          tx=tx-x0
          call scale_transform(nxy,itransfy,y,ty,.false.)
          call calc_parameters(nxy,npoly,tx,ty,pow,a,.false.)
          if(nx==1)then
            ! Single value is printed to stdout.
            op_npoly=npoly
            op_pow=pow
            op_a=a
            do i=1,nderiv
              call deriv_poly(op_npoly,op_pow,op_a)
            enddo ! i
            t1=eval_poly(op_npoly,op_pow,op_a,x_target(1)-x0)
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
              t1=eval_poly(op_npoly,op_pow,op_a,x_target(ix)-x0)
              write(io,*)x_target(ix),t1
            enddo ! ix
          endif

        endif ! have_dx.or.have_dy

        ! Clean up.
        if(nx>1)then
          write(6,*)'Data written to "'//trim(fname)//'".'
          write(6,*)
          close(io)
        endif
        nx=0
        deallocate(x_target)
      endif

    enddo

  END SUBROUTINE eval_fit_values_derivs


  ! COMPUTATION ROUTINES


  SUBROUTINE calc_parameters(nxy,npoly,x,y,pow,a,weighted,weight,da)
    !-------------------------------------------------------------!
    ! Evaluate the parameters using the method given in my notes. !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),pow(npoly)
    DOUBLE PRECISION,INTENT(out) :: a(npoly)
    LOGICAL,INTENT(in) :: weighted
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: weight(nxy)
    DOUBLE PRECISION,INTENT(out),OPTIONAL :: da(npoly)
    INTEGER info,i,j,ipiv(npoly),lwork,ierr
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
    call construct_c(nxy,npoly,x,y,pow,c,weighted,weight)

    ! Construct matrix M.
    call construct_M(nxy,npoly,x,pow,M,weighted,weight)

    ! Invert matrix M.
    Minv=M
    call dsytrf('L',npoly,Minv(1,1),npoly,ipiv(1),tempr(1),-1,info)
    if(info/=0)call quit('Matrix inversion failed (1).')
    lwork=nint(tempr(1))
    allocate(work(lwork),stat=ierr)
    if(ierr/=0)call quit('Allocation error: WORK (1).')
    call dsytrf('L',npoly,Minv(1,1),npoly,ipiv(1),work(1),lwork,info)
    if(info/=0)call quit('Matrix inversion failed (2).')
    deallocate(work)
    allocate(work(npoly),stat=ierr)
    if(ierr/=0)call quit('Allocation error: WORK (2).')
    call dsytri('L',npoly,Minv(1,1),npoly,ipiv(1),work(1),info)
    if(info/=0)call quit('Matrix inversion failed (3).')
    deallocate(work)

    ! Complete the upper triangle of Minv.
    do i=1,npoly-1
      do j=i+1,npoly
        Minv(i,j)=Minv(j,i)
      enddo ! j
    enddo ! i

    ! Evaluate the coefficients of the terms in the polynomial.
    a=matmul(Minv,c)

    ! Evaluate the standard errors in the coefficients.
    if(weighted)then
      do j=1,npoly
        da(j)=sqrt(Minv(j,j))
      enddo ! j
    elseif(present(da))then
      da=0.d0
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


  SUBROUTINE construct_c(nxy,npoly,x,y,pow,c,weighted,weight)
    !----------------------------------------!
    ! Construct the vector c (see my notes). !
    !----------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),pow(npoly)
    DOUBLE PRECISION,INTENT(out) :: c(npoly)
    LOGICAL,INTENT(in) :: weighted
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: weight(nxy)
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


  SUBROUTINE construct_M(nxy,npoly,x,pow,M,weighted,weight)
    !----------------------------------------!
    ! Construct the matrix M (see my notes). !
    ! Only need lower triangular part.       !
    !----------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,npoly
    DOUBLE PRECISION,INTENT(in) :: x(nxy),pow(npoly)
    DOUBLE PRECISION,INTENT(out) :: M(npoly,npoly)
    LOGICAL,INTENT(in) :: weighted
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: weight(nxy)
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


  SUBROUTINE eval_fit_monte_carlo(nxy,have_dx,have_dy,x,y,dx,dy,&
     &itransfx,itransfy,x0,npoly,pow,nrandom,nderiv,nx,x_target,&
     &fmean,ferr,fmean_1s,ferr_1s,fmean_2s,ferr_2s,fmed,fskew,fkurt)
    !------------------------------------------------------!
    ! Perform Monte Carlo sampling of data space to obtain !
    ! fit values or derivatives at specified points with   !
    ! error bars.                                          !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx,itransfy,npoly,nrandom,nderiv,nx
    LOGICAL,INTENT(in) :: have_dx,have_dy
    DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),x0,&
       &pow(npoly),x_target(nx)
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: fmean(nx),ferr(nx),&
       &fmean_1s(nx),ferr_1s(nx),fmean_2s(nx),ferr_2s(nx),fmed(nx),&
       &fskew(nx),fkurt(nx)
    ! Transformed data.
    DOUBLE PRECISION tx(nxy),ty(nxy)
    ! Polynomials resulting from operations on fitting polynomial.
    INTEGER op_npoly
    DOUBLE PRECISION op_pow(npoly),op_a(npoly)
    ! Random sampling arrays.
    DOUBLE PRECISION f_array(nrandom,nx),w_vector(nrandom)
    DOUBLE PRECISION ran_gauss_x(nxy),ran_gauss_y(nxy),ran_x(nxy),ran_y(nxy),&
       &ran_a(npoly)
    ! Distribution analysis.
    DOUBLE PRECISION var,skew,kurt,f2s_lo,f1s_lo,f1s_hi,f2s_hi
    ! Misc.
    INTEGER i,irandom,ix
    DOUBLE PRECISION t1

    ! Initialize.
    ran_x=x
    ran_y=y

    ! Loop over random points.
    do irandom=1,nrandom
      if(have_dx)then
        ran_gauss_x=gaussian_random_number(dx)
        ran_x=x+ran_gauss_x
      endif
      if(have_dy)then
        ran_gauss_y=gaussian_random_number(dy)
        ran_y=y+ran_gauss_y
      endif
      call scale_transform(nxy,itransfx,ran_x,tx,.false.)
      tx=tx-x0
      call scale_transform(nxy,itransfy,ran_y,ty,.false.)
      call calc_parameters(nxy,npoly,tx,ty,pow,ran_a,.false.)
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

    ! Evaluate statistics.
    do ix=1,nx
      call characterize_dist(nrandom,f_array(:,ix),w_vector,var,skew,kurt)
      if(present(fmean))fmean(ix)=sum(f_array(:,ix)*w_vector)/sum(w_vector)
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
        if(present(ferr_2s))ferr_2s(ix)=0.5d0*(f2s_hi-f2s_lo)
      endif
    enddo ! ix

  END SUBROUTINE eval_fit_monte_carlo


  ! UTILITIES.


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
