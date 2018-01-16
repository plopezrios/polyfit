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
  ! Precision for real-to-integer conversions and exponent comparisons.
  DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-10
  ! Scale transformations.
  INTEGER,PARAMETER :: NTRANSF=3
  INTEGER,PARAMETER :: ITRANSF_NONE=0,ITRANSF_REC=1,ITRANSF_LOG=2,ITRANSF_EXP=3
  CHARACTER(32),PARAMETER :: TRANSF_NAME(0:NTRANSF)=&
     &(/'linear     ','reciprocal ','logarithmic','exponential'/)
  LOGICAL,PARAMETER :: TRANSF_REQ_NONZERO(0:NTRANSF)=&
     &(/.false.,.true.,.true.,.false./)
  LOGICAL,PARAMETER :: TRANSF_REQ_POSITIVE(0:NTRANSF)=&
     &(/.false.,.false.,.true.,.false./)

  ! Types.
  TYPE xydata
    INTEGER :: nxy=0
    LOGICAL :: have_dx=.false.,have_dy=.false.
    DOUBLE PRECISION,POINTER :: x(:)=>null(),y(:)=>null(),dx(:)=>null(),&
       &dy(:)=>null()
  END TYPE xydata
  TYPE dataset
    INTEGER :: itransfx=ITRANSF_NONE,itransfy=ITRANSF_NONE
    TYPE(xydata),POINTER :: xy=>null() ! original data
    TYPE(xydata),POINTER :: txy=>null() ! transformed data
    TYPE(xydata),POINTER :: rtxy=>null() ! range-restricted transformed data
  END TYPE dataset
  TYPE fit_params
    INTEGER :: npoly=0
    DOUBLE PRECISION :: X0=0.d0
    DOUBLE PRECISION,ALLOCATABLE :: pow(:)
    LOGICAL,ALLOCATABLE :: share(:)
    CHARACTER(64) :: X0_string='0'
  END TYPE fit_params
  TYPE range_def
    CHARACTER :: var='X'
    CHARACTER(2) :: op=''
    DOUBLE PRECISION :: thres=0.d0
    INTEGER :: size=0
  END TYPE range_def
  TYPE eval_def
    CHARACTER :: var='X'
    LOGICAL :: rel=.false.
    INTEGER :: nderiv=0,n=0
    DOUBLE PRECISION,ALLOCATABLE :: x(:)
  END TYPE eval_def
  TYPE monte_carlo_params
    LOGICAL :: weighted=.true.
    INTEGER :: nsample=10000
  END TYPE monte_carlo_params

  call main()


CONTAINS


  SUBROUTINE main()
    !--------------!
    ! Main driver. !
    !--------------!
    IMPLICIT NONE
    ! (x,y[,dx][,dy]) datasets.
    INTEGER ndataset,tndataset
    TYPE(dataset),POINTER :: datasets(:)=>null(),tdatasets(:)=>null()
    ! Fit form.
    TYPE(fit_params),POINTER :: fit=>null()
    ! All-set settings.
    INTEGER itransfx_default,itransfy_default
    TYPE(range_def) drange
    ! Function evaluation.
    INTEGER eval_iset
    TYPE(eval_def) deval
    ! Monte Carlo parameters.
    TYPE(monte_carlo_params) mcparams
    ! Input file information.
    INTEGER ncolumn,icol_x,icol_y,icol_dx,icol_dy
    CHARACTER(256) fname
    ! Input file search facility.
    INTEGER,PARAMETER :: search_size=1024
    INTEGER nsearch,ndiscr
    INTEGER,POINTER :: fsearch(:)=>null(),fdiscr(:)=>null()
    CHARACTER(search_size),POINTER :: search(:)=>null()
    ! Misc variables.
    CHARACTER(8192) command,token
    INTEGER i,ifield,ierr,ierr1,ierr2,ierr3,ierr4
    INTEGER nxy,itransf,iset,i1,i2,ipos,npoly
    DOUBLE PRECISION t1

    ! Initialize.
    ndataset=0
    itransfx_default=ITRANSF_NONE
    itransfy_default=ITRANSF_NONE
    allocate(fit)
    fit%npoly=2
    allocate(fit%pow(fit%npoly),fit%share(fit%npoly))
    fit%pow(1:2)=(/0.d0,1.d0/)
    fit%share=.false.
    mcparams%weighted=.false.
    mcparams%nsample=10000

    ! Write header.
    write(6,'(a)')'===================================='
    write(6,'(a)')'POLYFIT - polynomial fitting toolbox'
    write(6,'(a)')'===================================='
    write(6,'()')
    write(6,'(a)')'Type "help" for a list of commands.'
    write(6,'()')

    ! Loop over user actions.
    user_loop: do

      ! Read user command.
      write(6,'(a)',advance='no')'polyfit> '
      read(5,'(a)',iostat=ierr)command
      if(ierr/=0)then
        write(6,'()')
        call quit()
      endif

      ! Execute command.
      select case(trim(field(1,command)))

      case('inspect')
        ! Report number of lines and columns in data file.
        fname=trim(field(2,command))
        call check_file(fname,nxy,ncolumn)
        write(6,'()')
        select case(nxy)
        case(-1)
          write(6,'(a)')'Could not open file "'//trim(fname)//'".'
          write(6,'()')
          cycle user_loop
        case(-2)
          write(6,'(a)')'Problem reading file "'//trim(fname)//'".'
          write(6,'()')
          cycle user_loop
        case(-3)
          write(6,'(a)')'Column count problem in file "'//trim(fname)//'".'
          write(6,'()')
          cycle user_loop
        case(0)
          write(6,'(a)')'File "'//trim(fname)//'" contains no useful data.'
          write(6,'()')
          cycle user_loop
        case default
          write(6,'(a)')'File "'//trim(fname)//'" contains '//&
             &trim(i2s(nxy))//' lines with '//trim(i2s(ncolumn))//' columns.'
          write(6,'()')
        end select

      case('load')
        ! Load data from specified columns of data file.
        fname=trim(field(2,command))
        call check_file(fname,nxy,ncolumn)
        write(6,'()')
        select case(nxy)
        case(-1)
          write(6,'(a)')'Could not open file "'//trim(fname)//'".'
          write(6,'()')
          cycle user_loop
        case(-2)
          write(6,'(a)')'Problem reading file "'//trim(fname)//'".'
          write(6,'()')
          cycle user_loop
        case(-3)
          write(6,'(a)')'Column count problem in file "'//trim(fname)//'".'
          write(6,'()')
          cycle user_loop
        case(0)
          write(6,'(a)')'File "'//trim(fname)//'" contains no useful data.'
          write(6,'()')
          cycle user_loop
        end select

        ! Initialize search facility.
        nsearch=0
        ndiscr=0

        ! Get column selection.
        ! Initialize to default column indices.
        icol_x=0
        icol_y=0
        icol_dx=0
        icol_dy=0
        select case(ncolumn)
        case(1)
          icol_y=1
        case(2)
          icol_x=1
          icol_y=2
        case(3)
          icol_x=1
          icol_y=2
          icol_dy=3
        case default
          icol_x=1
          icol_y=3
          icol_dx=2
          icol_dy=4
        end select

        ! Parse subcommands.
        ifield=2
        do
          ifield=ifield+1
          select case(trim(field(ifield,command)))
          case('type')
            if(trim(field(ifield+2,command))/='using')then
              write(6,'(a)')'Syntax error in load command: "type" without &
                 &"using".'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              cycle user_loop
            endif
            icol_x=0
            icol_y=0
            icol_dx=0
            icol_dy=0
            ierr1=0
            ierr2=0
            ierr3=0
            ierr4=0
            select case(field(ifield+1,command))
            case('xy')
              icol_x=int_field(ifield+3,command,ierr1)
              icol_y=int_field(ifield+4,command,ierr2)
              ifield=ifield+4
            case('xydy')
              icol_x=int_field(ifield+3,command,ierr1)
              icol_y=int_field(ifield+4,command,ierr2)
              icol_dy=int_field(ifield+5,command,ierr3)
              ifield=ifield+5
            case('xdxy')
              icol_x=int_field(ifield+3,command,ierr1)
              icol_dx=int_field(ifield+4,command,ierr2)
              icol_y=int_field(ifield+5,command,ierr3)
              ifield=ifield+5
            case('xdxydy')
              icol_x=int_field(ifield+3,command,ierr1)
              icol_dx=int_field(ifield+4,command,ierr2)
              icol_y=int_field(ifield+5,command,ierr3)
              icol_dy=int_field(ifield+6,command,ierr4)
              ifield=ifield+6
            case default
              write(6,'(a)')'Syntax error in load command: unknown type "'//&
                 &trim(field(ifield+1,command))//'".'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              if(associated(fdiscr))deallocate(fdiscr)
              cycle user_loop
            end select
          case('where')
            ! Add a search clause.
            if(nfield(command)<ifield+2)then
              write(6,'(a)')'Syntax error in load command: too few arguments &
                 &for "where" subcommand"'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              if(associated(fdiscr))deallocate(fdiscr)
              cycle user_loop
            endif
            i=int_field(ifield+1,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error in load command: invalid column &
                 &index for "where".'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              if(associated(fdiscr))deallocate(fdiscr)
              cycle user_loop
            endif
            if(i<1.or.i>ncolumn)then
              write(6,'(a)')'Column index out of range in "where" subcommand.'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              if(associated(fdiscr))deallocate(fdiscr)
              cycle user_loop
            endif
            nsearch=nsearch+1
            call resize_pointer_int1((/nsearch/),fsearch)
            call resize_pointer_char1(search_size,(/nsearch/),search)
            fsearch(nsearch)=i
            search(nsearch)=field(ifield+2,command)
            ifield=ifield+2
          case('by')
            ! Add a discriminator clause.
            if(nfield(command)<ifield+1)then
              write(6,'(a)')'Syntax error in load command: too few arguments &
                 &for "by" subcommand"'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              if(associated(fdiscr))deallocate(fdiscr)
              cycle user_loop
            endif
            i=int_field(ifield+1,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error in load command: invalid column &
                 &index for "by".'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              if(associated(fdiscr))deallocate(fdiscr)
              cycle user_loop
            endif
            if(i<1.or.i>ncolumn)then
              write(6,'(a)')'Column index out of range in "by" subcommand.'
              write(6,'()')
              if(associated(search))deallocate(fsearch,search)
              if(associated(fdiscr))deallocate(fdiscr)
              cycle user_loop
            endif
            ndiscr=ndiscr+1
            call resize_pointer_int1((/ndiscr/),fdiscr)
            fdiscr(ndiscr)=i
            ifield=ifield+1
          case('')
            exit
          case default
            write(6,'(a)')'Syntax error in load command: unknown subcommand "'&
               &//trim(field(ifield,command))//'".'
            write(6,'()')
            if(associated(search))deallocate(fsearch,search)
            if(associated(fdiscr))deallocate(fdiscr)
            cycle user_loop
          end select
        enddo ! ifield

        ! Check.
        if(ierr1/=0.or.ierr2/=0.or.ierr3/=0.or.ierr4/=0)then
          write(6,'(a)')'Problem parsing column indices.'
          write(6,'()')
          if(associated(search))deallocate(fsearch,search)
          if(associated(fdiscr))deallocate(fdiscr)
          cycle user_loop
        endif
        if(icol_x>ncolumn.or.icol_x<0.or.icol_y>ncolumn.or.icol_y<0.or.&
           &icol_dx>ncolumn.or.icol_dx<0.or.icol_dy>ncolumn.or.icol_dy<0)then
          write(6,'(a)')'Column indices out of range.'
          write(6,'()')
          if(associated(search))deallocate(fsearch,search)
          if(associated(fdiscr))deallocate(fdiscr)
          cycle user_loop
        endif

        ! Read data.
        if(.not.associated(search))allocate(fsearch(nsearch),search(nsearch))
        if(.not.associated(fdiscr))allocate(fdiscr(ndiscr))
        call read_file(fname,icol_x,icol_y,icol_dx,icol_dy,nsearch,fsearch,&
           &search,ndiscr,fdiscr,tndataset,tdatasets,ierr)
        deallocate(fsearch,search)
        deallocate(fdiscr)
        if(ierr/=0)cycle user_loop
        if(tndataset<1)then
          write(6,'(a)')'No data loaded.'
          write(6,'()')
          cycle user_loop
        endif

        ! Check data ar compatible with current transformations.
        do iset=1,tndataset
          tdatasets(iset)%itransfx=itransfx_default
          tdatasets(iset)%itransfy=itransfy_default
          if(TRANSF_REQ_NONZERO(tdatasets(iset)%itransfy))then
            if(any(are_equal(tdatasets(iset)%xy%x,0.d0)))then
              tdatasets(iset)%itransfx=ITRANSF_NONE
              write(6,'(a)')'Note: using linear X for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains x=0.'
              write(6,'()')
            endif
          endif
          if(TRANSF_REQ_NONZERO(tdatasets(iset)%itransfy))then
            if(any(are_equal(tdatasets(iset)%xy%y,0.d0)))then
              tdatasets(iset)%itransfy=ITRANSF_NONE
              write(6,'(a)')'Note: using linear Y for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains y=0.'
              write(6,'()')
            endif
          endif
          if(TRANSF_REQ_POSITIVE(tdatasets(iset)%itransfx))then
            if(any(tdatasets(iset)%xy%x<0.d0))then
              tdatasets(iset)%itransfx=ITRANSF_NONE
              write(6,'(a)')'Note: using linear X for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains x<0.'
              write(6,'()')
            endif
          endif
          if(TRANSF_REQ_POSITIVE(tdatasets(iset)%itransfy))then
            if(any(tdatasets(iset)%xy%y<0.d0))then
              tdatasets(iset)%itransfy=ITRANSF_NONE
              write(6,'(a)')'Note: using linear Y for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains y<0.'
              write(6,'()')
            endif
          endif
        enddo ! iset

        ! Add datasets to list.
        call resize_dataset(ndataset+tndataset,datasets)
        do iset=ndataset+1,ndataset+tndataset
          datasets(iset)=tdatasets(iset-ndataset)
          call refresh_dataset(datasets(iset),drange)
          ! Report.
          write(6,'(a)',advance='no')'Loaded data from "'//trim(fname)//&
             &'" as dataset #'//trim(i2s(iset))//', type '
          if(datasets(iset)%xy%have_dx.and.datasets(iset)%xy%have_dy)then
            write(6,'(a)',advance='no')'xdxydy'
          elseif(datasets(iset)%xy%have_dx)then
            write(6,'(a)',advance='no')'xdxy'
          elseif(datasets(iset)%xy%have_dy)then
            write(6,'(a)',advance='no')'xydy'
          else
            write(6,'(a)',advance='no')'xy'
          endif
          write(6,'(a)')', '//trim(i2s(datasets(iset)%xy%nxy))//' data.'
        enddo ! iset
        deallocate(tdatasets)
        nullify(tdatasets)
        ndataset=ndataset+tndataset
        call refresh_fit(ndataset,datasets,fit)
        write(6,'()')

      case('fit')
        write(6,'()')
        if(ndataset<1)then
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        if(nfield(command)>1)then
          write(6,'(a)')'Syntax error: subcommand "'//&
             &trim(field(2,command))//'" not recognized.'
          write(6,'()')
          cycle user_loop
        endif ! nfield>1
        call show_multipoly(ndataset,datasets,drange,fit,mcparams)

      case('plot')
        write(6,'()')
        if(ndataset<1)then
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        if(nfield(command)>1)then
          write(6,'(a)')'Syntax error: subcommand "'//&
             &trim(field(2,command))//'" not recognized.'
          write(6,'()')
          cycle user_loop
        endif ! nfield>1
        call plot_multipoly(ndataset,datasets,drange,fit,mcparams)

      case('assess')
        if(ndataset<1)then
          write(6,'()')
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        ! FIXME - expose this
        deval%nderiv=0
        deval%rel=.false.
        deval%var='X'
        deval%n=1
        eval_iset=0
        allocate(deval%x(deval%n))
        deval%x(1)=0.d0
        select case(trim(field(2,command)))
        case('fit')
          call assess_fit(ndataset,datasets,drange,fit,mcparams,deval,&
             &eval_iset)
        case('range')
          call assess_range(ndataset,datasets,fit,mcparams,deval,eval_iset)
        case('range,fit','fit,range')
          call assess_fit_range(ndataset,datasets,fit,mcparams,deval,eval_iset)
        case default
          write(6,'()')
          write(6,'(a)')'Unknown variable to assess "'//&
             &trim(field(2,command))//'".'
          write(6,'()')
        end select
        deallocate(deval%x)

      case('report')
        if(ndataset<1)then
          write(6,'()')
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        select case(trim(field(2,command)))
        case('range')
          call report_statistics(ndataset,datasets)
        case default
          write(6,'()')
          write(6,'(a)')'Unknown report name "'//trim(field(2,command))//'".'
          write(6,'()')
          cycle user_loop
        end select

      case('evaluate')
        write(6,'()')
        if(ndataset<1)then
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        ! Initialize unused components.
        deval%rel=.false.
        ! Get object to evaluate.
        select case(field(2,command))
        case("f")
          deval%nderiv=0
        case("f'")
          deval%nderiv=1
        case("f''")
          deval%nderiv=2
        case default
          write(6,'(a)')'Syntax error: unknown function "'//&
             &trim(field(2,command))//'".'
          write(6,'()')
          cycle user_loop
        end select
        ! Get where to evaluate it at.
        if(trim(field(3,command))/='at')then
          write(6,'(a)')'Syntax error: missing "at" subcommand.'
          write(6,'()')
          cycle user_loop
        endif
        select case(trim(field(4,command)))
        case('X')
          deval%var='X'
        case default
          write(6,'(a)')'Syntax error: unknown variable "'//&
             &trim(field(4,command))//'".'
          write(6,'()')
          cycle user_loop
        end select
        if(nfield(command)<5)then
          write(6,'(a)')'Syntax error: no values of "'//&
             &trim(field(4,command))//'" specified.'
          write(6,'()')
          cycle user_loop
        endif
        do ifield=5,nfield(command)
          t1=dble_field(ifield,command,ierr)
          if(ierr/=0)then
            write(6,'(a)')'Syntax error: invalid value of "'//&
               &trim(field(4,command))//'" specified.'
            write(6,'()')
            cycle user_loop
          endif
        enddo ! ifield
        deval%n=nfield(command)-4
        allocate(deval%x(deval%n))
        do ifield=5,nfield(command)
          deval%x(ifield-4)=dble_field(ifield,command,ierr)
        enddo ! ifield
        ! Perform evaluation.
        call evaluate_fit(ndataset,datasets,fit,mcparams,drange,deval)
        deallocate(deval%x)

      case('set')

        ! Set variables.
        write(6,'()')
        select case(trim(field(2,command)))

        case('xscale','yscale')
          ! Set scale transformation.
          do itransf=0,NTRANSF
            if(trim(field(3,command))==trim(TRANSF_NAME(itransf)))exit
          enddo ! itransf
          if(itransf>NTRANSF)then
            write(6,'(a)')'Unknown value "'//trim(field(3,command))//'" for &
               &variable "'//trim(field(2,command))//'".'
            write(6,'()')
            cycle user_loop
          endif
          if(nfield(command)>3)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(4,command))/='for')then
              write(6,'(a)')'Syntax error in set command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".'
              write(6,'()')
              cycle user_loop
            endif
            if(nfield(command)<5)then
              write(6,'(a)')'Syntax error in set command: "for" subcommand &
                 &requires arguments.'
              write(6,'()')
              cycle user_loop
            endif
            ! Check transformation is applicable.
            do ifield=5,nfield(command)
              iset=int_field(ifield,command,ierr)
              if(ierr/=0)then
                write(6,'(a)')'Syntax error in set command: could not parse &
                   &set values.'
                write(6,'()')
                cycle user_loop
              endif
              if(iset<1.or.iset>ndataset)then
                write(6,'(a)')'Set index out of range.'
                write(6,'()')
                cycle user_loop
              endif
              if(TRANSF_REQ_NONZERO(itransf))then
                select case(trim(field(2,command)))
                case('xscale')
                  if(any(are_equal(datasets(iset)%xy%x,0.d0)))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x=0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(are_equal(datasets(iset)%xy%y,0.d0)))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains y=0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                end select
              endif
              if(TRANSF_REQ_POSITIVE(itransf))then
                select case(trim(field(2,command)))
                case('xscale')
                  if(any(datasets(iset)%xy%x<0.d0))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x<0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(datasets(iset)%xy%y<0.d0))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains y<0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                end select
              endif
            enddo ! ifield
            ! Set transformation.
            do ifield=5,nfield(command)
              iset=int_field(ifield,command,ierr)
              select case(trim(field(2,command)))
              case('xscale')
                datasets(iset)%itransfx=itransf
              case('yscale')
                datasets(iset)%itransfy=itransf
              end select
            enddo ! ifield
            write(6,'(a)')'Set '//trim(field(2,command))//' to '//&
               &trim(TRANSF_NAME(itransf))//' for set #'//trim(i2s(iset))//'.'
          else
            ! Check transformation is applicable.
            do iset=1,ndataset
              if(TRANSF_REQ_NONZERO(itransf))then
                select case(trim(field(2,command)))
                case('xscale')
                  if(any(are_equal(datasets(iset)%xy%x,0.d0)))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x=0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(are_equal(datasets(iset)%xy%y,0.d0)))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains y=0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                end select
              endif
              if(TRANSF_REQ_POSITIVE(itransf))then
                select case(trim(field(2,command)))
                case('xscale')
                  if(any(datasets(iset)%xy%x<0.d0))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x<0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(datasets(iset)%xy%y<0.d0))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains y<0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                end select
              endif
            enddo ! iset
            ! Store transformation.
            do iset=1,ndataset
              select case(trim(field(2,command)))
              case('xscale')
                datasets(iset)%itransfx=itransf
              case('yscale')
                datasets(iset)%itransfy=itransf
              end select
            enddo ! ifield
            ! Store default transformation.
            select case(trim(field(2,command)))
            case('xscale')
              itransfx_default=itransf
            case('yscale')
              itransfy_default=itransf
            end select
            write(6,'(a)')'Set '//trim(field(2,command))//' to '//&
               &trim(TRANSF_NAME(itransf))//' for all sets.'
          endif
          write(6,'()')
          ! Apply transformations.
          do iset=1,ndataset
            call refresh_dataset(datasets(iset),drange)
          enddo ! iset
          ! Update X0.
          call refresh_fit(ndataset,datasets,fit)

        case('fit')
          ! Set fit exponents.
          ! Figure out syntax used.
          ierr=0
          select case(nfield(command))
          case(2)
            write(6,'(a)')'Syntax error: no value given for "fit".'
            write(6,'()')
            cycle user_loop
          case(3)
            t1=dble_field(3,command,ierr)
          end select
          if(ierr==0)then
            ! Exponents given as list.
            ! Check exponents.
            do ifield=3,nfield(command)
              t1=dble_field(ifield,command,ierr)
              if(ierr/=0)then
                write(6,'(a)')'Syntax error: could not parse exponents.'
                write(6,'()')
                cycle user_loop
              endif
            enddo ! ifield
            npoly=nfield(command)-2
            i1=0
            i2=-1
          else
            ! Exponents given as range.
            token=field(3,command)
            ipos=scan(token,':')
            if(ipos<1)then
              write(6,'(a)')'Syntax error: could not parse range.'
              write(6,'()')
              cycle user_loop
            endif
            i1=0
            if(ipos>1)then
              read(token(1:ipos-1),*,iostat=ierr)i1
              if(ierr/=0)then
                write(6,'(a)')'Syntax error: could not parse range.'
                write(6,'()')
                cycle user_loop
              endif
            endif
            read(token(ipos+1:),*,iostat=ierr)i2
            if(ierr/=0)then
              write(6,'(a)')'Syntax error: could not parse range.'
              write(6,'()')
              cycle user_loop
            endif
            npoly=i2-i1+1
            if(npoly<1)then
              write(6,'(a)')'Exponent range runs backwards.'
              write(6,'()')
              cycle user_loop
            endif
          endif
          ! Reallocate and set exponents.
          deallocate(fit%pow,fit%share)
          fit%npoly=npoly
          allocate(fit%pow(npoly),fit%share(npoly))
          if(i2<i1)then
            do ifield=3,nfield(command)
              fit%pow(ifield-2)=dble_field(ifield,command,ierr)
            enddo ! ifield
          else
            fit%pow(1:npoly)=(/(dble(i),i=i1,i2)/)
          endif
          fit%share(1:npoly)=.false.
          ! Report.
          write(6,'(a)')'Fit form set to:'
          write(6,'(2x,a)')trim(print_poly_sym(fit%npoly,fit%pow,fit%X0))
          write(6,'(2x,a)')'Shared coefficients reset to: none'
          write(6,'()')

        case('range')
          ! Set fit range.
          ! Check sort variable.
          select case(trim(field(3,command)))
          case('x','y','X','Y')
            continue
          case default
            write(6,'(a)')'Syntax error parsing range: unknown sort variable &
               &"'//trim(field(3,command))//'".'
            write(6,'()')
            cycle user_loop
          end select
          ! Check sort operation and set range parameters.
          select case(trim(field(4,command)))
          case('<','<=','>','>=')
            t1=dble_field(5,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error parsing range: could not parse &
                 &threshold.'
              write(6,'()')
              cycle user_loop
            endif
            drange%thres=t1
            drange%op=field(4,command)
          case('first','last')
            i=int_field(5,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error parsing range: could not parse &
                 &range size.'
              write(6,'()')
              cycle user_loop
            endif
            if(i<0)then
              write(6,'(a)')'Syntax error parsing range: range size < 0.'
              write(6,'()')
              cycle user_loop
            endif
            drange%size=i
            select case(trim(field(4,command)))
            case('first')
              drange%op='f'
            case('last')
              drange%op='l'
            end select
          case default
            write(6,'(a)')'Syntax error parsing range: unknown sort criterion &
               &"'//trim(field(4,command))//'".'
            write(6,'()')
            cycle user_loop
          end select
          ! Set remaining range parameters.
          drange%var=field(3,command)
          write(6,'(a)')'Range set.'
          write(6,'()')
          ! Apply transformations.
          do iset=1,ndataset
            call refresh_dataset(datasets(iset),drange)
          enddo ! iset
          ! Update X0.
          call refresh_fit(ndataset,datasets,fit)

        case('shared')
          ! Set shared coefficients.
          if(trim(field(3,command))=='all')then
            fit%share=.true.
          elseif(trim(field(3,command))=='none')then
            fit%share=.false.
          else
            select case(nfield(command))
            case(2)
              write(6,'(a)')'Syntax error: no value given for "shared".'
              write(6,'()')
              cycle user_loop
            case(3)
              t1=int_field(3,command,ierr)
            end select
            if(ierr==0)then
              ! Indices given as list.
              ! Check indices.
              do ifield=3,nfield(command)
                i1=int_field(ifield,command,ierr)
                if(ierr/=0)then
                  write(6,'(a)')'Syntax error: could not parse indices.'
                  write(6,'()')
                  cycle user_loop
                endif
                if(i1<1.or.i1>fit%npoly)then
                  write(6,'(a)')'Indices out of range.'
                  write(6,'()')
                  cycle user_loop
                endif
              enddo ! ifield
              ! Apply.
              fit%share=.false.
              do ifield=3,nfield(command)
                fit%share(int_field(ifield,command,ierr))=.true.
              enddo ! ifield
            else
              ! Exponents given as range.
              ! Get range and check syntax.
              token=field(3,command)
              ipos=scan(token,':')
              if(ipos<1)then
                write(6,'(a)')'Syntax error: could not parse range.'
                write(6,'()')
                cycle user_loop
              endif
              i1=1
              if(ipos>1)then
                read(token(1:ipos-1),*,iostat=ierr)i1
                if(ierr/=0)then
                  write(6,'(a)')'Syntax error: could not parse range.'
                  write(6,'()')
                  cycle user_loop
                endif
              endif
              i2=npoly
              if(ipos<len_trim(token))then
                read(token(ipos+1:),*,iostat=ierr)i2
                if(ierr/=0)then
                  write(6,'(a)')'Syntax error: could not parse range.'
                  write(6,'()')
                  cycle user_loop
                endif
              endif
              if(i2<i1)then
                write(6,'(a)')'Index range runs backwards.'
                write(6,'()')
                cycle user_loop
              endif
              if(i1<1.or.i2>fit%npoly)then
                write(6,'(a)')'Indices out of range.'
                write(6,'()')
                cycle user_loop
              endif
              ! Apply.
              fit%share=.false.
              do i=i1,i2
                fit%share(i)=.true.
              enddo
            endif
          endif
          write(6,'(a)',advance='no')'Shared coefficients set to:'
          if(.not.any(fit%share))then
            write(6,'(a)')' none'
          else
            do i=1,fit%npoly
              if(fit%share(i))write(6,'(a)',advance='no')' k'//trim(i2s(i))
            enddo ! i
            write(6,'()')
          endif
          write(6,'()')

        case('weights')
          select case(trim(field(3,command)))
          case('yes','YES','y','Y','true','TRUE','t','T')
            mcparams%weighted=.true.
          case('no','NO','n','N','false','FALSE','f','F')
            mcparams%weighted=.false.
          case default
            write(6,'(a)')'Invalid value "'//trim(field(3,command))//'" for &
               &variable "'//trim(field(2,command))//'".'
            write(6,'()')
            cycle user_loop
          end select
          write(6,'(a)',advance='no')'Set weights = '
          if(mcparams%weighted)then
            write(6,'(a)')'yes'
          else
            write(6,'(a)')'no'
          endif
          write(6,'()')

        case('centre')
          ! Check value.
          select case(field(3,command))
          case('left','right','max','min','centre','mean','median')
            continue
          case default
            t1=dble_field(3,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Invalid value "'//trim(field(3,command))//'" for &
                 &variable "'//trim(field(2,command))//'".'
              write(6,'()')
              cycle user_loop
            endif
          end select
          ! Set value.
          fit%x0_string=field(3,command)
          ! Update X0.
          call refresh_fit(ndataset,datasets,fit)

        case('nsample')
          ! Check value.
          i=int_field(3,command,ierr)
          if(ierr/=0)then
            write(6,'(a)')'Invalid value "'//trim(field(3,command))//'" for &
               &variable nsample.'
            write(6,'()')
            cycle user_loop
          endif
          if(i<10)then
            write(6,'(a)')'nsample value too small.'
            write(6,'()')
            cycle user_loop
          endif
          mcparams%nsample=i

        case default
          write(6,'(a)')'Unknown variable "'//trim(field(2,command))//'".'
          write(6,'()')
          cycle user_loop
        end select ! variable to set

      case('status')
        write(6,'()')
        write(6,'(a,es11.4)')'Global settings:'
        select case(drange%op)
        case('')
          write(6,'(a,es11.4)')'  Data range: all data'
        case('f')
          write(6,'(a,es11.4)')'  Data range: first '//&
             &trim(i2s(drange%size))//' data by '//trim(drange%var)//' value'
        case('l')
          write(6,'(a,es11.4)')'  Data range: last '//&
             &trim(i2s(drange%size))//' data by '//trim(drange%var)//' value'
        case default
          write(6,'(a,es11.4)')'  Data range: '//trim(drange%var)//' '//&
             &trim(drange%op)//' ',drange%thres
        end select
        write(6,'(a)',advance='no')'  Use 1/dy^2 weights in Monte Carlo: '
        if(mcparams%weighted)then
          write(6,'(a)')'yes'
        else
          write(6,'(a)')'no'
        endif
        write(6,'(a)')'  Number of Monte Carlo samples: '//&
           &trim(i2s(mcparams%nsample))
        write(6,'()')
        select case(ndataset)
        case(0)
          write(6,'(a)')'No datasets loaded.'
        case(1)
          write(6,'(a)')'1 dataset loaded:'
        case default
          write(6,'(a)')trim(i2s(ndataset))//' datasets loaded:'
        end select
        do iset=1,ndataset
          write(6,'(a)',advance='no')'* Set #'//trim(i2s(iset))//': type '
          if(datasets(iset)%xy%have_dx.and.datasets(iset)%xy%have_dy)then
            write(6,'(a)',advance='no')'xdxydy'
          elseif(datasets(iset)%xy%have_dx)then
            write(6,'(a)',advance='no')'xdxy'
          elseif(datasets(iset)%xy%have_dy)then
            write(6,'(a)',advance='no')'xydy'
          else
            write(6,'(a)',advance='no')'xy'
          endif
          if(datasets(iset)%xy%nxy/=datasets(iset)%rtxy%nxy)then
            write(6,'(a)',advance='no')', '//&
               &trim(i2s(datasets(iset)%xy%nxy))//' ('//&
               &trim(i2s(datasets(iset)%rtxy%nxy))//') data, '
          else
            write(6,'(a)',advance='no')', '//&
               &trim(i2s(datasets(iset)%xy%nxy))//' data, '
          endif
          write(6,'(a)',advance='no')&
             &trim(TRANSF_NAME(datasets(iset)%itransfx))//'-'//&
             &trim(TRANSF_NAME(datasets(iset)%itransfy))//' scale'
          write(6,'(a)')'.'
        enddo ! i
        write(6,'()')
        write(6,'(a)')'Fit function:'
        write(6,'(2x,a)')trim(print_poly_sym(fit%npoly,fit%pow,fit%X0))
        write(6,'(a)',advance='no')'  Shared coefficients:'
        if(.not.any(fit%share))then
          write(6,'(a)')' none'
        else
          do i=1,fit%npoly
            if(fit%share(i))write(6,'(a)',advance='no')' k'//trim(i2s(i))
          enddo ! i
          write(6,'()')
        endif
        write(6,'("  X0 = ",a," =",es12.4)')trim(fit%x0_string),fit%X0
        write(6,'()')

      case('help')
        select case(trim(field(2,command)))
        case('')
          write(6,'()')
          write(6,'(a)')'POLYFIT is a toolbox for performing polynomial fits &
             &on one or more sets of'
          write(6,'(a)')'data.  It handles non-integer exponents, a few common &
             &data transformations and'
          write(6,'(a)')'data with statistical errorbars.  Its most useful &
             &feature is the ability to'
          write(6,'(a)')'provide confidence intervals for values and &
             &derivatives of a fit.'
          write(6,'()')
          write(6,'(a)')'POLYFIT uses a rudimentary command-line interface.  &
             &The list of available'
          write(6,'(a)')'commands is:'
          write(6,'()')
          write(6,'(a)')'* inspect <file>'
          write(6,'(a)')'* load <file> [type <type> using <columns>]'
          write(6,'(a)')'* set <variable> <value> [for <set-list>]'
          write(6,'(a)')'* status'
          write(6,'(a)')'* fit'
          write(6,'(a)')'* assess <variable> [using <function> at <X> [for &
             &<set-list>]]'
          write(6,'(a)')'* help [<command> | set <variable>]'
          write(6,'()')
          write(6,'(a)')'Type help <command> for detailed information.'
          write(6,'()')
        case('inspect')
          write(6,'()')
          write(6,'(a)')'Command: inspect <file>'
          write(6,'()')
          write(6,'(a)')'  Reports the number of data lines and columns &
             &detected in <file>.'
          write(6,'()')
        case('load')
          write(6,'()')
          write(6,'(a)')'Command: load <file> [type <type> using <columns>]'
          write(6,'()')
          write(6,'(a)')'  Loads data from <file> into a new dataset.  By &
             &default:'
          write(6,'(a)')'  * 1-column files are of type xy, with &
             &(x,y)=(index,1)'
          write(6,'(a)')'  * 2-column files are of type xy, with &
             &(x,y)=(1,2)'
          write(6,'(a)')'  * 3-column files are of type xydy, with &
             &(x,y,dy)=(1,2,3)'
          write(6,'(a)')'  * 4- or more-column files are of type xdxydy, with &
             &(x,dx,y,dy)=(1,2,3,4)'
          write(6,'(a)')'  Other column selections can be specified with an &
             &explicity type/using clause,'
          write(6,'(a)')'  e.g., "type xydy using 3 5 7".  A column index of &
             &zero refers to the line'
          write(6,'(a)')'  index.'
          write(6,'()')
        case('set')
          if(nfield(command)==2)then
            write(6,'()')
            write(6,'(a)')'Command: set <variable> <value> [for <set-list>]'
            write(6,'()')
            write(6,'(a)')'  Sets <variable> to <value>, either globally or for &
               &selected sets (for'
            write(6,'(a)')'  certain variables).  The list of available &
               &variables is:'
            write(6,'()')
            write(6,'(a)')'  * xscale'
            write(6,'(a)')'  * yscale'
            write(6,'(a)')'  * fit'
            write(6,'(a)')'  * range'
            write(6,'(a)')'  * shared'
            write(6,'()')
            write(6,'(a)')'Type help set <variable> for detailed information.'
            write(6,'()')
          else
            select case(field(3,command))
            case('xscale','yscale')
              write(6,'()')
              write(6,'(a)')'Variable: xscale, yscale'
              write(6,'()')
              write(6,'(a)')'  These set scale transformations for the &
                 &independent x variable and for the'
              write(6,'(a)')'  dependent y variable, respectively.  In &
                 &POLYFIT notation, the original'
              write(6,'(a)')'  variables are called x and y, and the &
                 &transformed variables are called X and'
              write(6,'(a)')'  Y.'
              write(6,'()')
              write(6,'(a)')'  Possible values are:'
              write(6,'(a)')'  * "linear"      : X = x       [default]'
              write(6,'(a)')'  * "reciprocal"  : X = 1/x     [x=0 forbidden]'
              write(6,'(a)')'  * "logarithmic" : X = log(x)  [x<=0 forbidden]'
              write(6,'(a)')'  * "exponential" : X = exp(x)'
              write(6,'()')
              write(6,'(a)')'  These variables can be set in a per-set manner &
                 &or globally.  Note that the'
              write(6,'(a)')'  global value applies to all loaded datasets &
                 &and becomes the default for new'
              write(6,'(a)')'  datasets.'
              write(6,'()')
            case('fit')
              write(6,'()')
              write(6,'(a)')'Variable: fit'
              write(6,'()')
              write(6,'(a)')'  This sets the exponents to be used in the &
                 &fitting function.  The value can be'
              write(6,'(a)')'  specified as a list, e.g., "set fit 0 &
                 &1.5 3", or as an integer range,'
              write(6,'(a)')'  e.g., "set fit 0:2".  Note that the form &
                 &of the fitting function is'
              write(6,'(a)')'  global, so the "for <set-list>" syntax does &
                 &not apply to this variable.'
              write(6,'()')
            case('range')
              write(6,'()')
              write(6,'(a)')'Variable: range'
              write(6,'()')
              write(6,'(a)')'  This sets the data mask to apply to all &
                 &datasets.  The value can be specified'
              write(6,'(a)')'  as <variable> <selector> where:'
              write(6,'(a)')'  * <variable> is one of x, y, X, or Y &
                 &(lowercase for original and uppercase'
              write(6,'(a)')'    for transformed values).'
              write(6,'(a)')'  * <selector> is one of:'
              write(6,'(a)')'    - <  <threshold>'
              write(6,'(a)')'    - <= <threshold>'
              write(6,'(a)')'    - >= <threshold>'
              write(6,'(a)')'    - >  <threshold>'
              write(6,'(a)')'    - first <number>'
              write(6,'(a)')'    - last <number>'
              write(6,'(a)')'  Note that "range" is a global variable, so the &
                 &"for <set-list>" syntax'
              write(6,'(a)')'  does not apply to this variable.'
              write(6,'()')
            case('shared')
              write(6,'()')
              write(6,'(a)')'Variable: shared'
              write(6,'()')
              write(6,'(a)')'  This specifies which coefficients are to be &
                 &shared among datasets.  The value'
              write(6,'(a)')'  is a list of coefficient indices, or "all".  &
                 &Coefficients not flagged'
              write(6,'(a)')'  as shared take different values for each &
                 &dataset.'
              write(6,'()')
            case default
              write(6,'(a)')'No help for variable "'//&
                 &trim(field(2,command))//'".'
            end select
          endif
        case('status')
          write(6,'()')
          write(6,'(a)')'Command: status'
          write(6,'()')
          write(6,'(a)')'  Report currently loaded datasets and values of &
             &internal variables.'
          write(6,'()')
        case('fit')
          write(6,'()')
          write(6,'(a)')'Command: fit'
          write(6,'()')
          write(6,'(a)')'  Perform fit of currently loaded datasets.'
          write(6,'()')
        case('assess')
          write(6,'()')
          write(6,'(a)')'Command: assess <variables> [using <function> at <X> &
             &[for <sets>]]'
          write(6,'()')
          write(6,'(a)')'  Assess the convergence of the fit with the &
             &specified variables.'
          write(6,'(a)')'  The assessment is carried out based on the value &
             &chi^2/Ndf and on the'
          write(6,'(a)')'  value of the value or derivative of the fit at the &
             &specified position.'
          write(6,'(a)')'  The following <variables> can be specified:'
          write(6,'(a)')'  * fit       : assess convergence with &
             &choice of fit form.'
          write(6,'(a)')'  * range     : assess convergence with &
             &data range.'
          write(6,'(a)')'  * fit,range : assess convergence with choice of &
             &fit form and data range'
          write(6,'()')
        case('report')
          write(6,'()')
          write(6,'(a)')'Command: report <report>'
          write(6,'()')
          write(6,'(a)')'  Produce the requested report of the loaded &
             &datasets.'
          write(6,'(a)')'  Available reports are:'
          write(6,'(a)')'  * order : report convergence of f(0) with number &
             &of parameters.'
          write(6,'(a)')'  * size  : report convergence of f(0) with number &
             &of data points.'
          write(6,'(a)')'  * dual  : report convergence of f(0) with number &
             &of parameters and'
          write(6,'(a)')'    data points.'
          write(6,'(a)')'  * range : report basic range statistics of the &
             &data.'
          write(6,'()')
        case default
          write(6,'(a)')'No help for command "'//trim(field(2,command))//'".'
        end select

      case('quit','exit')
        call quit()

      case('')
        continue

      case default
        write(6,'()')
        write(6,'(a)')'Command "'//trim(field(1,command))//'" not recognized.'
        write(6,'()')
      end select

    enddo user_loop

  END SUBROUTINE main


  SUBROUTINE refresh_dataset(dset,drange)
    !-------------------------------------------------------!
    ! Re-apply transformation and mask to dataset following !
    ! change in choice of transformation or mask.           !
    !-------------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset),INTENT(inout) :: dset
    TYPE(range_def),INTENT(in) :: drange
    TYPE(xydata),POINTER :: xy
    LOGICAL,ALLOCATABLE :: mask(:)
    INTEGER,ALLOCATABLE :: indx(:)
    DOUBLE PRECISION,ALLOCATABLE :: sortvec(:)

    ! Delete pre-existing data.
    if(associated(dset%txy))then
      deallocate(dset%txy%x,dset%txy%y,dset%txy%dx,dset%txy%dy)
      deallocate(dset%txy)
    endif
    if(associated(dset%rtxy))then
      deallocate(dset%rtxy%x,dset%rtxy%y,dset%rtxy%dx,dset%rtxy%dy)
      deallocate(dset%rtxy)
    endif

    ! Transform.
    allocate(xy)
    xy%nxy=dset%xy%nxy
    xy%have_dx=dset%xy%have_dx
    xy%have_dy=dset%xy%have_dy
    allocate(xy%x(xy%nxy),xy%y(xy%nxy),xy%dx(xy%nxy),xy%dy(xy%nxy))
    call scale_transform(xy%nxy,dset%itransfx,dset%xy%x,xy%x,dset%xy%have_dx,&
       &dset%xy%dx,xy%dx)
    call scale_transform(xy%nxy,dset%itransfy,dset%xy%y,xy%y,dset%xy%have_dy,&
       &dset%xy%dy,xy%dy)
    dset%txy=>xy
    nullify(xy)

    ! Mask.
    allocate(sortvec(dset%xy%nxy),mask(dset%xy%nxy))
    mask=.true.
    ! Build vector with sort variable.
    select case(drange%var)
    case('x')
      sortvec=dset%xy%x
    case('y')
      sortvec=dset%xy%y
    case('X')
      sortvec=dset%txy%x
    case('Y')
      sortvec=dset%txy%y
    end select
    ! Act on sort operation.
    select case(trim(drange%op))
    case('<')
      mask=sortvec<drange%thres.and..not.are_equal(sortvec,drange%thres)
    case('<=')
      mask=sortvec<drange%thres.or.are_equal(sortvec,drange%thres)
    case('>')
      mask=sortvec>drange%thres.and..not.are_equal(sortvec,drange%thres)
    case('>=')
      mask=sortvec>drange%thres.or.are_equal(sortvec,drange%thres)
    case('f','l')
      if(drange%size>0)then
        allocate(indx(dset%xy%nxy))
        call isort(dset%xy%nxy,sortvec,indx)
        mask=.false.
        if(trim(drange%op)=='f')then
          mask(indx(1:drange%size))=.true.
        elseif(trim(drange%op)=='l')then
          mask(indx(dset%xy%nxy-drange%size+1:dset%xy%nxy))=.true.
        endif
        deallocate(indx)
      endif
    end select

    ! Create masked transformed dataset.
    allocate(xy)
    xy%nxy=count(mask)
    xy%have_dx=dset%txy%have_dx
    xy%have_dy=dset%txy%have_dy
    allocate(xy%x(xy%nxy),xy%y(xy%nxy),xy%dx(xy%nxy),xy%dy(xy%nxy))
    if(xy%nxy>0)then
      xy%x=pack(dset%txy%x,mask)
      xy%dx=pack(dset%txy%dx,mask)
      xy%y=pack(dset%txy%y,mask)
      xy%dy=pack(dset%txy%dy,mask)
    endif
    dset%rtxy=>xy
    nullify(xy)

  END SUBROUTINE refresh_dataset


  SUBROUTINE refresh_fit(ndataset,datasets,fit)
    !------------------------------------------------!
    ! Refresh value of X0 in fit following change to !
    ! X0_string, range or loaded data.               !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(ndataset)
    TYPE(fit_params),POINTER :: fit
    INTEGER iset,tot_nxy,ixy
    DOUBLE PRECISION t1,tx0
    DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:)

    select case(trim(fit%X0_string))

    case('left','right','max','min','centre','mean','median')
      ! Build combined dataset.
      tot_nxy=0
      do iset=1,ndataset
        tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
      enddo ! iset
      allocate(x(tot_nxy),y(tot_nxy))
      tot_nxy=0
      do iset=1,ndataset
        if(datasets(iset)%rtxy%nxy==0)cycle
        x(tot_nxy+1:tot_nxy+datasets(iset)%rtxy%nxy)=&
           &datasets(iset)%rtxy%x(1:datasets(iset)%rtxy%nxy)
        y(tot_nxy+1:tot_nxy+datasets(iset)%rtxy%nxy)=&
           &datasets(iset)%rtxy%y(1:datasets(iset)%rtxy%nxy)
        tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
      enddo ! iset

      ! Locate X0.
      select case(trim(fit%X0_string))
      case('left')
        tx0=minval(x)
      case('right')
        tx0=maxval(x)
      case('min')
        t1=0.d0
        tx0=0.d0
        do ixy=1,tot_nxy
          if(ixy==1.or.y(ixy)<t1)then
            t1=y(ixy)
            tx0=x(ixy)
          endif
        enddo ! ixy
      case('max')
        t1=0.d0
        tx0=0.d0
        do ixy=1,tot_nxy
          if(ixy==1.or.y(ixy)>t1)then
            t1=y(ixy)
            tx0=x(ixy)
          endif
        enddo ! ixy
      case('centre')
        tx0=.5d0*(minval(x)+maxval(x))
      case('mean')
        tx0=sum(x)/dble(tot_nxy)
      case('median')
        tx0=median(tot_nxy,x)
      end select

      ! Clean up.
      deallocate(x,y)

    case default
      ! Parse real number.
      read(fit%x0_string,*)t1
      tx0=t1
    end select

    ! Update X0.
    fit%X0=tx0

  END SUBROUTINE refresh_fit


  SUBROUTINE show_multipoly(ndataset,datasets,drange,fit,mcparams)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(ndataset)
    TYPE(range_def),INTENT(in) :: drange
    TYPE(fit_params),INTENT(in) :: fit
    TYPE(monte_carlo_params),INTENT(in) :: mcparams
    TYPE(eval_def) deval
    INTEGER tot_nparam,tot_nxy,iset,i
    DOUBLE PRECISION chi2,chi2err,a(fit%npoly,ndataset),da(fit%npoly,ndataset)

    ! Initialize.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)
    tot_nxy=0
    do iset=1,ndataset
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset

    ! Perform fit.
    deval%nderiv=0
    deval%rel=.false.
    deval%var='X'
    deval%n=0
    call eval_multifit_monte_carlo(ndataset,datasets,drange,fit,mcparams,&
       &deval,chi2=chi2,chi2err=chi2err,amean=a,aerr=da)

    ! Print fit coefficients.
    do iset=1,ndataset
      write(6,'(a)')'Set #'//trim(i2s(iset))//':'
      write(6,'(a)')'  Fit parameters:'
      do i=1,fit%npoly
        write(6,'(4x,a," = ",es24.16," +/- ",es24.16)')'k'//trim(i2s(i)),&
           &a(i,iset),da(i,iset)
      enddo ! i
      ! Write out fitted polynomial.
      ! NB, this is good for pasting into xmgrace, but xmgrace has a string
      ! length limit of 256 characters.
      write(6,'(a)')'  Fitted polynomial:'
      write(6,'(a)')'    '//trim(print_poly_num(fit%npoly,fit%pow,a(:,iset),&
         &fit%X0))
      write(6,'()')
    enddo ! iset

    ! Report chi-squared.
    write(6,'(a)')'Fit assessment:'
    write(6,'(a,es20.12," +/- ",es20.12)')'  chi^2     = ',chi2,chi2err
    if(tot_nxy>tot_nparam)write(6,'(a,es20.12," +/- ",es20.12)')&
       &'  chi^2/Ndf = ',chi2/dble(tot_nxy-tot_nparam),&
       &chi2err/dble(tot_nxy-tot_nparam)
    write(6,'()')

  END SUBROUTINE show_multipoly


  SUBROUTINE plot_multipoly(ndataset,datasets,drange,fit,mcparams)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(ndataset)
    TYPE(range_def),INTENT(in) :: drange
    TYPE(fit_params),INTENT(in) :: fit
    TYPE(monte_carlo_params),INTENT(in) :: mcparams
    TYPE(eval_def) deval
    INTEGER tot_nparam,tot_nxy,iset,i,nplot,ierr
    DOUBLE PRECISION chi2,chi2err,a(fit%npoly,ndataset),&
       &da(fit%npoly,ndataset),atilde(fit%npoly),tx0,tx1,tx
    INTEGER,PARAMETER :: io=10

    ! Initialize.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)
    tot_nxy=0
    do iset=1,ndataset
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset

    ! Perform fit.
    deval%nderiv=0
    deval%rel=.false.
    deval%var='X'
    deval%n=0
    call eval_multifit_monte_carlo(ndataset,datasets,drange,fit,mcparams,&
       &deval,chi2=chi2,chi2err=chi2err,amean=a,aerr=da)

    ! Plot.
    open(unit=io,file='fit.plot',status='replace',iostat=ierr)
    if(ierr/=0)then
      write(6,'(a)')'Problem opening file.'
      write(6,'()')
      return
    endif

    ! Find plot range.
    nplot=100
    tx0=minval(datasets(1)%rtxy%x)
    tx1=maxval(datasets(1)%rtxy%x)
    do iset=2,ndataset
      tx0=min(tx0,minval(datasets(iset)%rtxy%x))
      tx1=max(tx1,maxval(datasets(iset)%rtxy%x))
    enddo ! iset

    if(all(.not.fit%share))then
      ! Plot data.
      do iset=1,ndataset
        if(datasets(iset)%rtxy%have_dx.and.datasets(iset)%rtxy%have_dy)then
          write(io,'(a)')'@type xydxdy'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i),&
               &datasets(iset)%rtxy%dx(i),datasets(iset)%rtxy%dy(i)
          enddo ! i
        elseif(datasets(iset)%rtxy%have_dx)then
          write(io,'(a)')'@type xydx'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i),&
               &datasets(iset)%rtxy%dx(i)
          enddo ! i
        elseif(datasets(iset)%rtxy%have_dy)then
          write(io,'(a)')'@type xydy'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i),&
               &datasets(iset)%rtxy%dy(i)
          enddo ! i
        else
          write(io,'(a)')'@type xy'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i)
          enddo ! i
        endif
        write(io,'(a)')'&'
      enddo ! iset
      ! Plot fit functions.
      do iset=1,ndataset
        write(io,'(a)')'@type xy'
        do i=0,nplot
          tx=tx0+(tx1-tx0)*dble(i)/dble(nplot)
          write(io,*)tx,eval_poly(fit%npoly,fit%pow,a(1,iset),tx)
        enddo ! i
        write(io,'(a)')'&'
      enddo ! iset
    else
      ! Plot data minus independent bit.
      do iset=1,ndataset
        atilde(1:fit%npoly)=a(1:fit%npoly,iset)
        where(fit%share)atilde=0.d0
        if(datasets(iset)%rtxy%have_dx.and.datasets(iset)%rtxy%have_dy)then
          write(io,'(a)')'@type xydxdy'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i)-&
               &eval_poly(fit%npoly,fit%pow,atilde,datasets(iset)%rtxy%x(i)),&
               &datasets(iset)%rtxy%dx(i),datasets(iset)%rtxy%dy(i)
          enddo ! i
        elseif(datasets(iset)%rtxy%have_dx)then
          write(io,'(a)')'@type xydx'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i)-&
               &eval_poly(fit%npoly,fit%pow,atilde,datasets(iset)%rtxy%x(i)),&
               &datasets(iset)%rtxy%dx(i)
          enddo ! i
        elseif(datasets(iset)%rtxy%have_dy)then
          write(io,'(a)')'@type xydy'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i)-&
               &eval_poly(fit%npoly,fit%pow,atilde,datasets(iset)%rtxy%x(i)),&
               &datasets(iset)%rtxy%dy(i)
          enddo ! i
        else
          write(io,'(a)')'@type xy'
          do i=1,datasets(iset)%rtxy%nxy
            write(io,*)datasets(iset)%rtxy%x(i),datasets(iset)%rtxy%y(i)-&
               &eval_poly(fit%npoly,fit%pow,atilde,datasets(iset)%rtxy%x(i))
          enddo ! i
        endif
        write(io,'(a)')'&'
      enddo ! iset
      ! Plot shared part of fit functions.
      atilde(1:fit%npoly)=a(1:fit%npoly,1)
      where(.not.fit%share)atilde=0.d0
      write(io,'(a)')'@type xy'
      do i=0,nplot
        tx=tx0+(tx1-tx0)*dble(i)/dble(nplot)
        write(io,*)tx,eval_poly(fit%npoly,fit%pow,atilde,tx)
      enddo ! i
    endif

    ! Close file.
    close(io)

    ! Report.
    write(6,'(a)')'Plot saved to fit.plot.'
    write(6,'()')

  END SUBROUTINE plot_multipoly


  SUBROUTINE assess_fit(ndataset,datasets,drange,fit,mcparams,deval,eval_iset)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! expansion order.                               !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),POINTER :: datasets(:)
    TYPE(range_def),INTENT(in) :: drange
    TYPE(fit_params),INTENT(in) :: fit
    TYPE(monte_carlo_params),INTENT(in) :: mcparams
    TYPE(eval_def),INTENT(in) :: deval
    INTEGER,INTENT(in) :: eval_iset
    TYPE(fit_params),POINTER :: tfit=>null()
    INTEGER tot_nxy,tot_nparam,iset,npoly,ix,prev_npoly
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION chi2_all(fit%npoly),chi2err_all(fit%npoly),&
       &fmean_all(deval%n,ndataset,fit%npoly),&
       &ferr_all(deval%n,ndataset,fit%npoly)

    ! Prepare combined dataset arrays.
    tot_nxy=0
    do iset=1,ndataset
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset

    ! Print table header.
    write(6,'()')
    write(6,'(2x,a5,1x,2(1x,a12))',advance='no')'Order','chi^2/Ndf',&
       &'dchi^2/Ndf'
    do iset=1,ndataset
      do ix=1,deval%n
        if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,a12))',advance='no')&
           &'f'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  ',&
           &'df'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  '
      enddo ! ix
    enddo ! iset
    write(6,'()')
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(32))//'("-"))')
    endif

    ! Initialize test results to "invalid".
    chi2_all=-1.d0
    chi2err_all=-1.d0
    fmean_all=0.d0
    ferr_all=-1.d0

    ! Loop over expansion orders.
    do npoly=1,fit%npoly
      ! Allocate work arrays.
      allocate(tfit)
      tfit%npoly=npoly
      allocate(tfit%pow(npoly),tfit%share(npoly))
      tfit%pow=fit%pow(1:npoly)
      tfit%share=fit%share(1:npoly)
      tfit%X0_string=fit%X0_string
      call refresh_fit(ndataset,datasets,tfit)
      ! Adjust counters.
      tot_nparam=count(tfit%share)+ndataset*count(.not.tfit%share)
      if(tot_nparam<tot_nxy)then
        ! Perform fit and report.
        call eval_multifit_monte_carlo(ndataset,datasets,drange,tfit,mcparams,&
           &deval,fmean,ferr,chi2,chi2err)
        chi2_all(npoly)=chi2/dble(tot_nxy-tot_nparam)
        chi2err_all(npoly)=chi2err/dble(tot_nxy-tot_nparam)
        fmean_all(:,:,npoly)=fmean
        ferr_all(:,:,npoly)=ferr
        write(6,'(2x,i5,1x,2(1x,es12.4))',advance='no')npoly-1,&
           &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
        do ix=1,deval%n
          do iset=1,ndataset
            if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,es12.4))',&
               &advance='no')fmean(ix,iset),ferr(ix,iset)
          enddo ! iset
        enddo ! ix
        write(6,'()')
      endif
      ! Destroy work arrays.
      deallocate(tfit%pow,tfit%share)
      deallocate(tfit)
      nullify(tfit)
    enddo ! npoly

    ! Print table footer.
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(32))//'("-"))')
    endif
    write(6,'()')

    ! Report best choice of parameters.
    ! Locate minimum chi2/Ndf.
    min_chi2=-1.d0
    do npoly=1,fit%npoly
      if(chi2_all(npoly)<0.d0)cycle
      t1=chi2_all(npoly)+2.d0*chi2err_all(npoly)
      if(min_chi2<0.d0.or.t1<min_chi2)min_chi2=t1
    enddo ! npoly
    ! Converge function value with respect to npoly.
    prev_npoly=0
    do npoly=1,fit%npoly
      ! Skip invalid points.
      if(chi2_all(npoly)<0.d0)cycle
      ! Skip points whose chi2 is not within uncertainty of minimum.
      if(min_chi2<chi2_all(npoly)-2.d0*chi2err_all(npoly))cycle
      if(prev_npoly>0)then
        ! Check all functions.
        do iset=1,ndataset
          do ix=1,deval%n
            t1=fmean_all(ix,iset,npoly)-fmean_all(ix,iset,prev_npoly)
            t2=sqrt(ferr_all(ix,iset,npoly)**2+&
               &ferr_all(ix,iset,prev_npoly)**2)
            if(abs(t1)>3.d0*t2)exit
          enddo ! ix
          if(ix<=deval%n)exit
        enddo ! iset
        if(iset>ndataset)exit
        if(.true.)exit ! ignore functions
      endif
      prev_npoly=npoly
    enddo ! npoly
    if(npoly<=fit%npoly)then
      write(6,'(a)')'Suggested fit: '//trim(i2s(prev_npoly-1))
    else
      write(6,'(a)')'Could not find optimal fit by criteria.'
    endif
    write(6,'()')

  END SUBROUTINE assess_fit


  SUBROUTINE assess_range(ndataset,datasets,fit,mcparams,deval,eval_iset)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! number of data points.                         !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(:)
    TYPE(fit_params),POINTER :: fit
    TYPE(monte_carlo_params),INTENT(in) :: mcparams
    TYPE(eval_def),INTENT(in) :: deval
    INTEGER,INTENT(in) :: eval_iset
    TYPE(dataset),POINTER :: tdatasets(:)
    TYPE(xydata),POINTER :: xy=>null()
    TYPE(range_def) drange
    INTEGER ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,ix,prev_igrid
    INTEGER,ALLOCATABLE :: indx(:)
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION,ALLOCATABLE :: chi2_all(:),chi2err_all(:),&
       &fmean_all(:,:,:),ferr_all(:,:,:)
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:)

    ! FIXME - expose these parameters.
    drange%var='X'
    drange%op='<='
    drange%size=0

    ! Initialize.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)

    ! Compute grid.
    tot_nxy=0
    do iset=1,ndataset
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset
    allocate(txall(tot_nxy),indx(tot_nxy),txgrid(tot_nxy))
    tot_nxy=0
    do iset=1,ndataset
      txall(tot_nxy+1:tot_nxy+datasets(iset)%rtxy%nxy)=&
         &datasets(iset)%rtxy%x(1:datasets(iset)%rtxy%nxy)
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset
    call isort(tot_nxy,txall,indx)
    ngrid=1
    txgrid(ngrid)=txall(indx(1))
    do ixy=2,tot_nxy
      if(are_equal(txall(indx(ixy)),txgrid(ngrid)))cycle
      ngrid=ngrid+1
      txgrid(ngrid)=txall(indx(ixy))
    enddo ! ixy
    deallocate(txall,indx)

    ! Allocate arrays to store test results and initialize to "invalid".
    allocate(chi2_all(tot_nxy),chi2err_all(tot_nxy),&
       &fmean_all(deval%n,ndataset,tot_nxy),&
       &ferr_all(deval%n,ndataset,tot_nxy))
    chi2_all=-1.d0
    chi2err_all=-1.d0
    fmean_all=0.d0
    ferr_all=-1.d0

    ! Make copy of datasets.
    allocate(tdatasets(ndataset))
    do iset=1,ndataset
      allocate(xy)
      xy%nxy=datasets(iset)%xy%nxy
      allocate(xy%x(xy%nxy),xy%dx(xy%nxy),xy%y(xy%nxy),xy%dy(xy%nxy))
      xy%x=datasets(iset)%xy%x
      xy%dx=datasets(iset)%xy%dx
      xy%y=datasets(iset)%xy%y
      xy%dy=datasets(iset)%xy%dy
      xy%have_dx=datasets(iset)%xy%have_dx
      xy%have_dy=datasets(iset)%xy%have_dy
      tdatasets(iset)%xy=>xy
      tdatasets(iset)%itransfx=datasets(iset)%itransfx
      tdatasets(iset)%itransfy=datasets(iset)%itransfy
      nullify(tdatasets(iset)%txy)
      nullify(tdatasets(iset)%rtxy)
      nullify(xy)
    enddo ! iset

    ! Print table header.
    write(6,'()')
    write(6,'(2x,a5,1x,2(1x,a12))',advance='no')'Ndata','chi^2/Ndf',&
       &'dchi^2/Ndf'
    do ix=1,deval%n
      do iset=1,ndataset
        if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,a12))',advance='no')&
           &'f'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  ',&
           &'df'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  '
      enddo ! iset
    enddo ! ix
    write(6,'()')
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(32))//'("-"))')
    endif

    ! Loop over points in grid.
    do igrid=ngrid,1,-1
      ! Apply range.
      drange%thres=txgrid(igrid)
      tot_nxy=0
      do iset=1,ndataset
        call refresh_dataset(tdatasets(iset),drange)
        tot_nxy=tot_nxy+tdatasets(iset)%rtxy%nxy
      enddo ! iset
      if(tot_nparam<tot_nxy)then
        ! Perform fit and report.
        call eval_multifit_monte_carlo(ndataset,tdatasets,drange,fit,mcparams,&
           &deval,fmean,ferr,chi2,chi2err)
        chi2_all(igrid)=chi2/dble(tot_nxy-tot_nparam)
        chi2err_all(igrid)=chi2err/dble(tot_nxy-tot_nparam)
        fmean_all(:,:,igrid)=fmean
        ferr_all(:,:,igrid)=ferr
        write(6,'(2x,i5,1x,2(1x,es12.4))',advance='no')tot_nxy,&
           &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
        do ix=1,deval%n
          do iset=1,ndataset
            if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,es12.4))',&
               &advance='no')fmean(ix,iset),ferr(ix,iset)
          enddo ! iset
        enddo ! ix
        write(6,'()')
      endif
    enddo ! igrid

    ! Print table footer.
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(32))//'("-"))')
    endif
    write(6,'()')

    ! Report best choice of parameters.
    ! Locate minimum chi2/Ndf.
    min_chi2=-1.d0
    do igrid=1,ngrid
      if(chi2_all(igrid)<0.d0)cycle
      t1=chi2_all(igrid)+2.d0*chi2err_all(igrid)
      if(min_chi2<0.d0.or.t1<min_chi2)min_chi2=t1
    enddo ! igrid
    ! Converge function value with respect to igrid.
    prev_igrid=0
    do igrid=ngrid,1,-1
      ! Skip invalid points.
      if(chi2_all(igrid)<0.d0)cycle
      ! Skip points whose chi2 is not within uncertainty of minimum.
      if(min_chi2<chi2_all(igrid)-2.d0*chi2err_all(igrid))cycle
      if(prev_igrid>0)then
        ! Check all functions.
        do iset=1,ndataset
          do ix=1,deval%n
            t1=fmean_all(ix,iset,igrid)-fmean_all(ix,iset,prev_igrid)
            t2=sqrt(ferr_all(ix,iset,igrid)**2+&
               &ferr_all(ix,iset,prev_igrid)**2)
            if(abs(t1)>3.d0*t2)exit
          enddo ! ix
          if(ix<=deval%n)exit
        enddo ! iset
        if(iset>ndataset)exit
      endif
      prev_igrid=igrid
    enddo ! igrid
    if(igrid>=1)then
      write(6,'(a)')'Suggested grid point: '//trim(i2s(prev_igrid))
    else
      write(6,'(a)')'Could not find optimal range by criteria.'
    endif
    write(6,'()')

    ! Clean up.
    do iset=1,ndataset
      deallocate(tdatasets(iset)%xy%x,tdatasets(iset)%xy%dx,&
         &tdatasets(iset)%xy%y,tdatasets(iset)%xy%dy)
      deallocate(tdatasets(iset)%txy%x,tdatasets(iset)%txy%dx,&
         &tdatasets(iset)%txy%y,tdatasets(iset)%txy%dy)
      deallocate(tdatasets(iset)%rtxy%x,tdatasets(iset)%rtxy%dx,&
         &tdatasets(iset)%rtxy%y,tdatasets(iset)%rtxy%dy)
      deallocate(tdatasets(iset)%xy)
      deallocate(tdatasets(iset)%txy)
      deallocate(tdatasets(iset)%rtxy)
    enddo ! iset
    deallocate(tdatasets)
    nullify(tdatasets)
    deallocate(chi2_all,chi2err_all,fmean_all,ferr_all)

  END SUBROUTINE assess_range


  SUBROUTINE assess_fit_range(ndataset,datasets,fit,mcparams,deval,eval_iset)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! expansion order and number of data points.     !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),POINTER :: datasets(:)
    TYPE(fit_params),POINTER :: fit
    TYPE(monte_carlo_params),INTENT(in) :: mcparams
    TYPE(eval_def),INTENT(in) :: deval
    INTEGER,INTENT(in) :: eval_iset
    TYPE(dataset),POINTER :: tdatasets(:)
    TYPE(xydata),POINTER :: xy=>null()
    TYPE(range_def) drange
    TYPE(fit_params),POINTER :: tfit=>null()
    INTEGER ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,npoly,ix,prev_igrid,&
       &prev_npoly,prev_prev_npoly
    INTEGER,ALLOCATABLE :: indx(:)
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION,ALLOCATABLE :: chi2_all(:,:),chi2err_all(:,:),&
       &fmean_all(:,:,:,:),ferr_all(:,:,:,:)
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:)

    ! FIXME - expose these parameters.
    drange%var='X'
    drange%op='<='
    drange%size=0

    ! Compute grid.
    tot_nxy=0
    do iset=1,ndataset
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset
    allocate(txall(tot_nxy),indx(tot_nxy),txgrid(tot_nxy))
    tot_nxy=0
    do iset=1,ndataset
      txall(tot_nxy+1:tot_nxy+datasets(iset)%rtxy%nxy)=&
         &datasets(iset)%rtxy%x(1:datasets(iset)%rtxy%nxy)
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset
    call isort(tot_nxy,txall,indx)
    ngrid=1
    txgrid(ngrid)=txall(indx(1))
    do ixy=2,tot_nxy
      if(are_equal(txall(indx(ixy)),txgrid(ngrid)))cycle
      ngrid=ngrid+1
      txgrid(ngrid)=txall(indx(ixy))
    enddo ! ixy
    deallocate(txall,indx)

    ! Allocate arrays to store test results and initialize to "invalid".
    allocate(chi2_all(fit%npoly,tot_nxy),chi2err_all(fit%npoly,tot_nxy),&
       &fmean_all(deval%n,ndataset,fit%npoly,tot_nxy),&
       &ferr_all(deval%n,ndataset,fit%npoly,tot_nxy))
    chi2_all=-1.d0
    chi2err_all=-1.d0
    fmean_all=0.d0
    ferr_all=-1.d0

    ! Make copy of datasets.
    allocate(tdatasets(ndataset))
    do iset=1,ndataset
      allocate(xy)
      xy%nxy=datasets(iset)%xy%nxy
      allocate(xy%x(xy%nxy),xy%dx(xy%nxy),xy%y(xy%nxy),xy%dy(xy%nxy))
      xy%x=datasets(iset)%xy%x
      xy%dx=datasets(iset)%xy%dx
      xy%y=datasets(iset)%xy%y
      xy%dy=datasets(iset)%xy%dy
      xy%have_dx=datasets(iset)%xy%have_dx
      xy%have_dy=datasets(iset)%xy%have_dy
      tdatasets(iset)%xy=>xy
      tdatasets(iset)%itransfx=datasets(iset)%itransfx
      tdatasets(iset)%itransfy=datasets(iset)%itransfy
      nullify(tdatasets(iset)%txy)
      nullify(tdatasets(iset)%rtxy)
      nullify(xy)
    enddo ! iset

    ! Print table header.
    write(6,'()')
    write(6,'(2x,a5,2x,a5,1x,2(1x,a12))',advance='no')'Ndata','Order',&
       &'chi^2/Ndf','dchi^2/Ndf'
    do ix=1,deval%n
      do iset=1,ndataset
        if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,a12))',advance='no')&
           &'f'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  ',&
           &'df'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  '
      enddo ! iset
    enddo ! ix
    write(6,'()')
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(39))//'("-"))')
    endif

    ! Loop over points in grid.
    do igrid=ngrid,1,-1
      ! Apply range.
      drange%var='X'
      drange%op='<='
      drange%thres=txgrid(igrid)
      drange%size=0
      tot_nxy=0
      do iset=1,ndataset
        call refresh_dataset(tdatasets(iset),drange)
        tot_nxy=tot_nxy+tdatasets(iset)%rtxy%nxy
      enddo ! iset
      ! Loop over expansion orders.
      do npoly=1,fit%npoly
        ! Allocate work arrays.
        allocate(tfit)
        tfit%npoly=npoly
        allocate(tfit%pow(npoly),tfit%share(npoly))
        tfit%pow=fit%pow(1:npoly)
        tfit%share=fit%share(1:npoly)
        tfit%X0_string=fit%X0_string
        call refresh_fit(ndataset,tdatasets,tfit)
        ! Adjust counters.
        tot_nparam=count(tfit%share)+ndataset*count(.not.tfit%share)
        if(tot_nparam<tot_nxy)then
          ! Perform fit and report.
          call eval_multifit_monte_carlo(ndataset,tdatasets,drange,tfit,&
             &mcparams,deval,fmean,ferr,chi2,chi2err)
          chi2_all(npoly,igrid)=chi2/dble(tot_nxy-tot_nparam)
          chi2err_all(npoly,igrid)=chi2err/dble(tot_nxy-tot_nparam)
          fmean_all(:,:,npoly,igrid)=fmean
          ferr_all(:,:,npoly,igrid)=ferr
          write(6,'(2x,i5,2x,i5,1x,2(1x,es12.4))',advance='no')&
             &tot_nxy,npoly-1,&
             &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
          do ix=1,deval%n
            do iset=1,ndataset
              if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,es12.4))',&
                 &advance='no')fmean(ix,iset),ferr(ix,iset)
            enddo ! iset
          enddo ! ix
          write(6,'()')
        endif
        ! Destroy work arrays.
        deallocate(tfit%pow,tfit%share)
        deallocate(tfit)
        nullify(tfit)
      enddo ! npoly
    enddo ! igrid

    ! Print table footer.
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(39))//'("-"))')
    endif
    write(6,'()')

    ! Report best choice of parameters.
    ! Locate minimum chi2/Ndf.
    min_chi2=-1.d0
    do igrid=1,ngrid
      do npoly=1,fit%npoly
        if(chi2_all(npoly,igrid)<0.d0)cycle
        t1=chi2_all(npoly,igrid)+2.d0*chi2err_all(npoly,igrid)
        if(min_chi2<0.d0.or.t1<min_chi2)min_chi2=t1
      enddo ! npoly
    enddo ! igrid
    ! Converge function value with respect to igrid.
    prev_igrid=0
    do igrid=ngrid,1,-1
      prev_npoly=0
      do npoly=1,fit%npoly
        ! Skip invalid points.
        if(chi2_all(npoly,igrid)<0.d0)cycle
        ! Skip points whose chi2 is not within uncertainty of minimum.
        if(min_chi2<chi2_all(npoly,igrid)-2.d0*chi2err_all(npoly,igrid))cycle
        if(prev_npoly>0)then
          ! Check all functions.
          do iset=1,ndataset
            do ix=1,deval%n
              t1=fmean_all(ix,iset,npoly,igrid)-&
                 &fmean_all(ix,iset,prev_npoly,igrid)
              t2=sqrt(ferr_all(ix,iset,npoly,igrid)**2+&
                 &ferr_all(ix,iset,prev_npoly,igrid)**2)
              if(abs(t1)>3.d0*t2)exit
            enddo ! ix
            if(ix<=deval%n)exit
          enddo ! iset
          if(iset>ndataset)exit
        endif
        prev_npoly=npoly
      enddo ! npoly
      ! Skip points for which test did not converge.
      if(npoly>fit%npoly)cycle
      if(prev_igrid>0)then
        ! Check function value across grid sizes.
        do iset=1,ndataset
          do ix=1,deval%n
            t1=fmean_all(ix,iset,prev_prev_npoly,igrid)-&
               &fmean_all(ix,iset,prev_prev_npoly,prev_igrid)
            t2=sqrt(ferr_all(ix,iset,prev_prev_npoly,igrid)**2+&
               &ferr_all(ix,iset,prev_prev_npoly,prev_igrid)**2)
            if(abs(t1)>3.d0*t2)exit
          enddo ! ix
          if(ix<=deval%n)exit
        enddo ! iset
        if(iset>ndataset)exit
      endif
      prev_igrid=igrid
      prev_prev_npoly=prev_npoly
    enddo ! igrid
    if(igrid>=1)then
      write(6,'(a)')'Suggested fit: '//trim(i2s(prev_prev_npoly-1))
      write(6,'(a)')'Suggested grid point: '//trim(i2s(prev_igrid))
    else
      write(6,'(a)')'Could not find optimal range by criteria.'
    endif
    write(6,'()')

    ! Clean up.
    do iset=1,ndataset
      deallocate(tdatasets(iset)%xy%x,tdatasets(iset)%xy%dx,&
         &tdatasets(iset)%xy%y,tdatasets(iset)%xy%dy)
      deallocate(tdatasets(iset)%txy%x,tdatasets(iset)%txy%dx,&
         &tdatasets(iset)%txy%y,tdatasets(iset)%txy%dy)
      deallocate(tdatasets(iset)%rtxy%x,tdatasets(iset)%rtxy%dx,&
         &tdatasets(iset)%rtxy%y,tdatasets(iset)%rtxy%dy)
      deallocate(tdatasets(iset)%xy)
      deallocate(tdatasets(iset)%txy)
      deallocate(tdatasets(iset)%rtxy)
    enddo ! iset
    deallocate(tdatasets)
    nullify(tdatasets)
    deallocate(chi2_all,chi2err_all,fmean_all,ferr_all)

  END SUBROUTINE assess_fit_range


  SUBROUTINE report_statistics(ndataset,datasets)
    !-----------------------------------------------------!
    ! Show min/centre/mean/median/max of X, Y, dX and dY. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(ndataset)
    INTEGER iset
    DOUBLE PRECISION vmin,vmax,vcentre,vmean,vmedian

    ! Print stats.
    write(6,'()')
    write(6,'(a)')'Data-range statistics (ignoring user-provided range):'
    write(6,'()')
    write(6,'(2x,a3,1x,a4,5(1x,a12))')'var','set','min','mean','median',&
       &'centre','max'
    write(6,'(2x,74("-"))')
    do iset=1,ndataset
      ! Transformed X and dX.
      vmin=minval(datasets(iset)%txy%x)
      vmax=maxval(datasets(iset)%txy%x)
      vcentre=.5d0*(vmin+vmax)
      vmean=sum(datasets(iset)%txy%x)/dble(datasets(iset)%txy%nxy)
      vmedian=median(datasets(iset)%txy%nxy,datasets(iset)%txy%x)
      write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'X',iset,vmin,vmean,vmedian,&
         &vcentre,vmax
    enddo ! iset
    do iset=1,ndataset
      if(datasets(iset)%txy%have_dx)then
        vmin=minval(datasets(iset)%txy%dx)
        vmax=maxval(datasets(iset)%txy%dx)
        vcentre=.5d0*(vmin+vmax)
        vmean=sum(datasets(iset)%txy%dx)/dble(datasets(iset)%txy%nxy)
        vmedian=median(datasets(iset)%txy%nxy,datasets(iset)%txy%dx)
        write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'dX',iset,vmin,vmean,vmedian,&
           &vcentre,vmax
      endif
    enddo ! iset
    do iset=1,ndataset
      ! Transformed Y and dY.
      vmin=minval(datasets(iset)%txy%y)
      vmax=maxval(datasets(iset)%txy%y)
      vcentre=.5d0*(vmin+vmax)
      vmean=sum(datasets(iset)%txy%y)/dble(datasets(iset)%txy%nxy)
      vmedian=median(datasets(iset)%txy%nxy,datasets(iset)%txy%y)
      write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'Y',iset,vmin,vmean,vmedian,&
         &vcentre,vmax
    enddo ! iset
    do iset=1,ndataset
      if(datasets(iset)%txy%have_dy)then
        vmin=minval(datasets(iset)%txy%dy)
        vmax=maxval(datasets(iset)%txy%dy)
        vcentre=.5d0*(vmin+vmax)
        vmean=sum(datasets(iset)%txy%dy)/dble(datasets(iset)%txy%nxy)
        vmedian=median(datasets(iset)%txy%nxy,datasets(iset)%txy%dy)
        write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'dY',iset,vmin,vmean,vmedian,&
           &vcentre,vmax
      endif
    enddo ! iset
    write(6,'(2x,74("-"))')
    write(6,'()')

  END SUBROUTINE report_statistics


  SUBROUTINE evaluate_fit(ndataset,datasets,fit,mcparams,drange,deval)
    !-------------------------------------------------!
    ! Evaluate value or derivative of fit at provided !
    ! points.                                         !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(ndataset)
    TYPE(fit_params),INTENT(in) :: fit
    TYPE(monte_carlo_params),INTENT(in) :: mcparams
    TYPE(range_def),INTENT(in) :: drange
    TYPE(eval_def),INTENT(in) :: deval
    INTEGER iset,ix
    DOUBLE PRECISION,ALLOCATABLE :: fmean(:,:),ferr(:,:)

    ! Evaluate.
    allocate(fmean(deval%n,ndataset),ferr(deval%n,ndataset))
    call eval_multifit_monte_carlo(ndataset,datasets,drange,fit,mcparams,&
       &deval,fmean,ferr)

    ! Report.
    write(6,'()')
    write(6,'(2x,a4,1x,3(1x,a16))')'set','X       ','f       ','df       '
    write(6,'(2x,56("-"))')
    do iset=1,ndataset
      do ix=1,deval%n
        write(6,'(2x,i4,1x,3(1x,f16.12))')iset,deval%x(ix),&
           &fmean(ix,iset),ferr(ix,iset)
      enddo ! ix
    enddo ! iset
    write(6,'(2x,56("-"))')
    write(6,'()')
    deallocate(fmean,ferr)

  END SUBROUTINE evaluate_fit


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


  ! COMPUTATION ROUTINES


  SUBROUTINE perform_multifit(ndataset,datasets,fit,weighted,chi2,a,da)
    !----------------------------------------------------!
    ! Perform least-squares fit of sets of (weighted) xy !
    ! data to the polynomial of exponents pow(1:npoly),  !
    ! with equal/independent coefficients for each set   !
    ! depending on the value of share(1:npoly).          !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(ndataset)
    TYPE(fit_params),INTENT(in) :: fit
    LOGICAL,INTENT(in) :: weighted
    DOUBLE PRECISION,INTENT(inout) :: chi2,a(fit%npoly,ndataset)
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: da(fit%npoly,ndataset)
    TYPE(xydata),POINTER :: xy=>null()
    DOUBLE PRECISION,ALLOCATABLE :: x(:,:),y(:,:),weight(:,:)
    INTEGER tot_nparam,max_nxy,tot_nxy,ieq,ip,jp,iset,jset,&
       &ixy,ipoly,jpoly,i,j,lwork,ierr
    DOUBLE PRECISION e_fit
    DOUBLE PRECISION,ALLOCATABLE :: M(:,:),Minv(:,:),c(:),dc(:),work(:)
    INTEGER,ALLOCATABLE :: ipiv(:)

    ! Extract fit properties.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)

    ! Prepare combined dataset arrays.
    max_nxy=0
    tot_nxy=0
    do iset=1,ndataset
      max_nxy=max(max_nxy,datasets(iset)%rtxy%nxy)
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset
    allocate(x(max_nxy,ndataset),y(max_nxy,ndataset),weight(max_nxy,ndataset))
    x=1.d0
    y=0.d0
    weight=0.d0
    do iset=1,ndataset
      xy=>datasets(iset)%rtxy
      if(xy%nxy==0)cycle
      x(1:xy%nxy,iset)=xy%x(1:xy%nxy)-fit%X0
      y(1:xy%nxy,iset)=xy%y(1:xy%nxy)
      if(.not.weighted)then
        weight(1:xy%nxy,iset)=1.d0
      else
        if(xy%have_dx.and.xy%have_dy)then
          weight(1:xy%nxy,iset)=1.d0/(xy%dx(1:xy%nxy)*xy%dy(1:xy%nxy))**2
        elseif(xy%have_dx)then
          weight(1:xy%nxy,iset)=1.d0/xy%dx(1:xy%nxy)**2
        elseif(xy%have_dy)then
          weight(1:xy%nxy,iset)=1.d0/xy%dy(1:xy%nxy)**2
        else
          weight(1:xy%nxy,iset)=1.d0
        endif
      endif
    enddo ! iset
    nullify(xy)

    ! Construct c vector and M matrix.
    allocate(M(tot_nparam,tot_nparam),Minv(tot_nparam,tot_nparam),&
       &ipiv(tot_nparam),c(tot_nparam),dc(tot_nparam))

    ! Initialize equation counter.
    ieq=0
    ! Loop over shared parameters.
    do ipoly=1,fit%npoly
      if(.not.fit%share(ipoly))cycle
      ! There is one equation for this parameter.
      ieq=ieq+1
      ! Right-hand side.
      c(ieq)=0.d0
      do iset=1,ndataset
        c(ieq)=c(ieq)+sum(weight(:,iset)*y(:,iset)*x(:,iset)**fit%pow(ipoly))
      enddo ! iset
      ! Coefficients for shared parameters.
      jp=0
      do jpoly=1,fit%npoly
        if(.not.fit%share(jpoly))cycle
        jp=jp+1
        M(jp,ieq)=0.d0
        do iset=1,ndataset
          M(jp,ieq)=M(jp,ieq)+sum(weight(:,iset)*x(:,iset)**&
             &(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! iset
      enddo ! jpoly
      ! Coefficients for independent parameters.
      do jpoly=1,fit%npoly
        if(fit%share(jpoly))cycle
        do iset=1,ndataset
          jp=jp+1
          M(jp,ieq)=sum(weight(:,iset)*x(:,iset)**&
             &(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! iset
      enddo ! jpoly
    enddo ! ipoly
    ! Loop over independent parameters.
    do ipoly=1,fit%npoly
      if(fit%share(ipoly))cycle
      ! There are ndataset equations for the ndataset instances of this
      ! parameter.
      do iset=1,ndataset
        ieq=ieq+1
        ! Right-hand side.
        c(ieq)=sum(weight(:,iset)*y(:,iset)*x(:,iset)**fit%pow(ipoly))
        ! Coefficients for shared parameters.
        jp=0
        do jpoly=1,fit%npoly
          if(.not.fit%share(jpoly))cycle
          jp=jp+1
          M(jp,ieq)=M(jp,ieq)+sum(weight(:,iset)*x(:,iset)**&
             &(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! jpoly
        ! Coefficients for independent parameters.
        do jpoly=1,fit%npoly
          if(fit%share(jpoly))cycle
          do jset=1,ndataset
            jp=jp+1
            if(jset==iset)then
              M(jp,ieq)=sum(weight(:,iset)*x(:,iset)**&
                 &(fit%pow(jpoly)+fit%pow(ipoly)))
            else
              M(jp,ieq)=0.d0
            endif
          enddo ! iset
        enddo ! jpoly
      enddo ! iset
    enddo ! ipoly

    ! Invert M.
    Minv=M
    allocate(work(1))
    lwork=-1
    call dsytrf('L',tot_nparam,Minv,tot_nparam,ipiv,work,lwork,ierr)
    if(ierr/=0)call quit('DSYTRF error on workspace query call.')
    lwork=nint(work(1))
    deallocate(work)
    allocate(work(lwork),stat=ierr)
    if(ierr/=0)call quit('Allocation error (work).')
    call dsytrf('L',tot_nparam,Minv,tot_nparam,ipiv,work,lwork,ierr)
    if(ierr/=0)call quit('DSYTRF error.')
    deallocate(work)
    allocate(work(tot_nparam),stat=ierr)
    if(ierr/=0)call quit('Allocation error (work).')
    call dsytri('L',tot_nparam,Minv,tot_nparam,ipiv,work,ierr)
    if(ierr/=0)call quit('DSYTRI error.')
    deallocate(work)

    ! Complete Minv and evaluate coefficients.
    do i=1,tot_nparam-1
      do j=i+1,tot_nparam
        Minv(i,j)=Minv(j,i)
      enddo ! j
    enddo ! i
    c=matmul(Minv,c)

    ! Evaluate standard errors if data are weighted.
    dc=0.d0
    if(weighted)then
      do i=1,tot_nparam
        dc(i)=sqrt(Minv(i,i))
      enddo ! i
    endif ! weighted

    ! Put parameters in output arrays.
    ip=0
    ! Loop over shared parameters.
    do ipoly=1,fit%npoly
      if(.not.fit%share(ipoly))cycle
      ip=ip+1
      a(ipoly,1:ndataset)=c(ip)
      if(present(da))da(ipoly,1:ndataset)=dc(ip)
    enddo ! ipoly
    ! Loop over independent parameters.
    do ipoly=1,fit%npoly
      if(fit%share(ipoly))cycle
      do iset=1,ndataset
        ip=ip+1
        a(ipoly,iset)=c(ip)
        if(present(da))da(ipoly,iset)=dc(ip)
      enddo ! iset
    enddo ! ipoly

    ! Return chi^2 value.
    chi2=0.d0
    do iset=1,ndataset
      do ixy=1,max_nxy
        e_fit=0.d0
        do ipoly=1,fit%npoly
          e_fit=e_fit+a(ipoly,iset)*x(ixy,iset)**fit%pow(ipoly)
        enddo ! ipoly
        chi2=chi2+(y(ixy,iset)-e_fit)**2*weight(ixy,iset)
      enddo ! ixy
    enddo ! iset

    ! Free memory.
    deallocate(M,Minv,ipiv,c,dc)

  END SUBROUTINE perform_multifit


  SUBROUTINE eval_multifit_monte_carlo(ndataset,datasets,drange,fit,mcparams,&
     &deval,fmean,ferr,chi2,chi2err,amean,aerr,fmean_1s,ferr_1s,fmean_2s,&
     &ferr_2s,fmed,fskew,fkurt)
    !------------------------------------------------------!
    ! Perform Monte Carlo sampling of data space to obtain !
    ! fit values or derivatives at specified points with   !
    ! error bars.                                          !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),INTENT(in) :: datasets(ndataset)
    TYPE(range_def),INTENT(in) :: drange
    TYPE(fit_params),INTENT(in) :: fit
    TYPE(eval_def),INTENT(in) :: deval
    TYPE(monte_carlo_params),INTENT(in) :: mcparams
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: &
       &fmean(deval%n,ndataset),ferr(deval%n,ndataset),&
       &chi2,chi2err,&
       &amean(fit%npoly,ndataset),aerr(fit%npoly,ndataset),&
       &fmean_1s(deval%n,ndataset),ferr_1s(deval%n,ndataset),&
       &fmean_2s(deval%n,ndataset),ferr_2s(deval%n,ndataset),&
       &fmed(deval%n,ndataset),fskew(deval%n,ndataset),&
       &fkurt(deval%n,ndataset)
    ! Pointers.
    TYPE(dataset),POINTER :: tdatasets(:)=>null()
    TYPE(xydata),POINTER :: xy=>null()
    ! Polynomials resulting from operations on fitting polynomial.
    INTEGER op_npoly
    DOUBLE PRECISION op_pow(fit%npoly),op_a(fit%npoly)
    ! Random sampling arrays.
    DOUBLE PRECISION f_array(mcparams%nsample,deval%n,ndataset),&
       &w_vector(mcparams%nsample),&
       &a_array(mcparams%nsample,fit%npoly,ndataset),&
       &chi2_array(mcparams%nsample)
    ! Distribution analysis.
    DOUBLE PRECISION var,skew,kurt,f2s_lo,f1s_lo,f1s_hi,f2s_hi
    ! Misc.
    LOGICAL weighted
    INTEGER nsample,ipoly,ideriv,irandom,ix,iset
    DOUBLE PRECISION t1,a(fit%npoly,ndataset),da(fit%npoly,ndataset),f0

    ! Make copy of datasets.
    allocate(tdatasets(ndataset))
    do iset=1,ndataset
      allocate(xy)
      xy%nxy=datasets(iset)%xy%nxy
      allocate(xy%x(xy%nxy),xy%dx(xy%nxy),xy%y(xy%nxy),xy%dy(xy%nxy))
      xy%x=datasets(iset)%xy%x
      xy%dx=datasets(iset)%xy%dx
      xy%y=datasets(iset)%xy%y
      xy%dy=datasets(iset)%xy%dy
      xy%have_dx=datasets(iset)%xy%have_dx
      xy%have_dy=datasets(iset)%xy%have_dy
      tdatasets(iset)%xy=>xy
      tdatasets(iset)%itransfx=datasets(iset)%itransfx
      tdatasets(iset)%itransfy=datasets(iset)%itransfy
      nullify(tdatasets(iset)%txy)
      nullify(tdatasets(iset)%rtxy)
      nullify(xy)
    enddo ! iset

    ! Override weights if not all sets have stderrs.
    weighted=mcparams%weighted
    do iset=1,ndataset
      weighted=weighted.and.(tdatasets(iset)%xy%have_dx.or.&
         &tdatasets(iset)%xy%have_dy)
    enddo ! iset

    ! Initialize.
    f0=0.d0
    nsample=mcparams%nsample

    ! Loop over random points.
    do irandom=1,nsample
      do iset=1,ndataset
        if(datasets(iset)%xy%have_dx)tdatasets(iset)%xy%x=&
           &datasets(iset)%xy%x+gaussian_random_number(datasets(iset)%xy%dx)
        if(datasets(iset)%xy%have_dy)tdatasets(iset)%xy%y=&
           &datasets(iset)%xy%y+gaussian_random_number(datasets(iset)%xy%dy)
        call refresh_dataset(tdatasets(iset),drange)
      enddo ! iset
      call perform_multifit(ndataset,tdatasets,fit,weighted,&
         &chi2_array(irandom),a,da)
      a_array(irandom,1:fit%npoly,1:ndataset)=a(1:fit%npoly,1:ndataset)
      w_vector(irandom)=1.d0
      do iset=1,ndataset
        op_npoly=fit%npoly
        op_pow=fit%pow
        op_a=a(:,iset)
        do ideriv=1,deval%nderiv
          call deriv_poly(op_npoly,op_pow,op_a)
        enddo ! ideriv
        if(deval%rel)f0=eval_poly(op_npoly,op_pow,op_a,-fit%X0)
        do ix=1,deval%n
          t1=eval_poly(op_npoly,op_pow,op_a,deval%x(ix)-fit%X0)
          f_array(irandom,ix,iset)=t1-f0
        enddo ! ix
      enddo ! iset
    enddo ! irandom

    ! Return coefficients and statistical properties of requested function
    ! of fit.
    if(present(chi2).or.present(chi2err))then
      call characterize_dist(nsample,chi2_array,w_vector,mean=t1,var=var)
      if(present(chi2))chi2=t1
      if(present(chi2err))chi2err=sqrt(var)
    endif
    do iset=1,ndataset
      do ipoly=1,fit%npoly
        call characterize_dist(nsample,a_array(:,ipoly,iset),w_vector,mean=t1,&
           &var=var)
        if(present(amean))amean(ipoly,iset)=t1
        if(present(aerr))aerr(ipoly,iset)=sqrt(var)
      enddo ! ipoly
      do ix=1,deval%n
        call characterize_dist(nsample,f_array(:,ix,iset),w_vector,mean=t1,&
           &var=var,skew=skew,kurt=kurt)
        if(present(fmean))fmean(ix,iset)=t1
        if(present(ferr))ferr(ix,iset)=sqrt(var)
        if(present(fskew))fskew(ix,iset)=skew
        if(present(fkurt))fkurt(ix,iset)=kurt
        if(present(fmed))fmed(ix,iset)=median(nsample,f_array(:,ix,iset))
        if(present(fmean_1s).or.present(ferr_1s))then
          f1s_lo=find_pth_smallest(nint(dble(nsample)*0.158655254d0),&
             &nsample,f_array(:,ix,iset))
          f1s_hi=find_pth_smallest(nint(dble(nsample)*0.841344746d0),&
             &nsample,f_array(:,ix,iset))
          if(present(fmean_1s))fmean_1s(ix,iset)=0.5d0*(f1s_lo+f1s_hi)
          if(present(ferr_1s))ferr_1s(ix,iset)=0.5d0*(f1s_hi-f1s_lo)
        endif
        if(present(fmean_2s).or.present(ferr_2s))then
          f2s_lo=find_pth_smallest(nint(dble(nsample)*0.022750131948d0),&
             &nsample,f_array(:,ix,iset))
          f2s_hi=find_pth_smallest(nint(dble(nsample)*0.977249868052d0),&
             &nsample,f_array(:,ix,iset))
          if(present(fmean_2s))fmean_2s(ix,iset)=0.5d0*(f2s_lo+f2s_hi)
          if(present(ferr_2s))ferr_2s(ix,iset)=0.25d0*(f2s_hi-f2s_lo)
        endif
      enddo ! ix
    enddo ! iset

    do iset=1,ndataset
      deallocate(tdatasets(iset)%xy%x,tdatasets(iset)%xy%dx,&
         &tdatasets(iset)%xy%y,tdatasets(iset)%xy%dy)
      deallocate(tdatasets(iset)%txy%x,tdatasets(iset)%txy%dx,&
         &tdatasets(iset)%txy%y,tdatasets(iset)%txy%dy)
      deallocate(tdatasets(iset)%rtxy%x,tdatasets(iset)%rtxy%dx,&
         &tdatasets(iset)%rtxy%y,tdatasets(iset)%rtxy%dy)
      deallocate(tdatasets(iset)%xy)
      deallocate(tdatasets(iset)%txy)
      deallocate(tdatasets(iset)%rtxy)
    enddo ! iset
    deallocate(tdatasets)
    nullify(tdatasets)

  END SUBROUTINE eval_multifit_monte_carlo


  ! INPUT FILE HANDLING ROUTINES.


  SUBROUTINE check_file(fname,nline,ncolumn)
    !----------------------------------------------------!
    ! Check file contains data, and return the number of !
    ! data lines and the number of data columns in it.   !
    !----------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(out) :: nline,ncolumn
    CHARACTER(8192) line
    INTEGER ipos,ierr,ncol
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Initialize.
    nline=-1 ! flag non-existing file
    ncolumn=0

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)return
    nline=0

    ! Loop over lines.
    do
      read(io,'(a)',iostat=ierr)line
      if(ierr<0)exit
      if(ierr>0)then
        nline=-2 ! flag reading error
        exit
      endif
      line=adjustl(line)
      ! Skip comments.
      ipos=scan(line,'#!')
      if(ipos==1)cycle
      if(ipos>1)line=line(1:ipos-1)
      ! Skip empty lines.
      if(len_trim(line)==0)cycle
      ! Find how many elements there are in this line.
      ncol=nfield(line)
      if(ncolumn==0)then
        ncolumn=ncol
      else
        ncolumn=min(ncolumn,ncol)
      endif
      if(ncolumn<1)then
        nline=-3 ! flag column count problem
        exit
      endif
      nline=nline+1
    enddo

    ! Close file.
    close(io)

  END SUBROUTINE check_file


  SUBROUTINE read_file(fname,icol_x,icol_y,icol_dx,icol_dy,nsearch,fsearch,&
     &search,ndiscr,fdiscr,ndataset,datasets,ierr)
    !-------------------------------------!
    ! Read in the data in the input file. !
    !-------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(in) :: icol_x,icol_y,icol_dx,icol_dy,nsearch,&
       &fsearch(nsearch),ndiscr,fdiscr(ndiscr)
    CHARACTER(*),INTENT(in) :: search(nsearch)
    INTEGER,INTENT(inout) :: ndataset,ierr
    TYPE(dataset),POINTER :: datasets(:)
    CHARACTER(8192) line
    CHARACTER(8192),POINTER :: discr(:,:)=>null()
    INTEGER i,ipos,isearch,iset,idiscr
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)then
      write(6,'(a)')'Problem opening "'//trim(fname)//'".'
      write(6,'()')
      return
    endif

    ! Initialize.
    ndataset=0

    ! Loop over lines.
    do
      read(io,'(a)',iostat=ierr)line
      if(ierr<0)then
        ierr=0
        exit
      endif
      if(ierr>0)then
        write(6,'(a)')'Problem getting line from "'//trim(fname)//'".'
        write(6,'()')
        exit
      endif
      line=adjustl(line)
      ! Skip comments.
      ipos=scan(line,'#!')
      if(ipos==1)cycle
      if(ipos>1)line=line(1:ipos-1)
      ! Skip empty lines.
      if(len_trim(line)==0)cycle
      ! Verify that this line contains all search strings.
      do isearch=1,nsearch
        if(.not.are_equal_string(trim(field(fsearch(isearch),line)),&
          &trim(search(isearch))))exit
      enddo ! isearch
      if(isearch<=nsearch)cycle
      ! Decide which dataset this goes in.
      do iset=1,ndataset
        do idiscr=1,ndiscr
          if(.not.are_equal_string(trim(field(fdiscr(idiscr),line)),&
             &trim(discr(idiscr,iset))))exit
        enddo ! idiscr
        if(idiscr>ndiscr)exit
      enddo ! iset
      if(iset>ndataset)then
        ! Create new dataset.
        ndataset=ndataset+1
        call resize_dataset(ndataset,datasets)
        iset=ndataset
        ! Store discriminators.
        call resize_pointer_char2(8192,(/ndiscr,ndataset/),discr)
        do idiscr=1,ndiscr
          discr(idiscr,iset)=trim(field(fdiscr(idiscr),line))
        enddo ! idiscr
      endif
      ! Enlarge arrays.
      if(.not.associated(datasets(iset)%xy))then
        i=1
      else
        i=datasets(iset)%xy%nxy+1
      endif
      call enlarge_xydata(i,icol_dx>0,icol_dy>0,datasets(iset)%xy)
      ! Read data point from string.
      if(icol_x>0)then
        datasets(iset)%xy%x(i)=dble_field(icol_x,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of x in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
      else
        datasets(iset)%xy%x(i)=dble(i)
      endif
      if(icol_y>0)then
        datasets(iset)%xy%y(i)=dble_field(icol_y,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of y in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
      else
        datasets(iset)%xy%y(i)=dble(i)
      endif
      if(icol_dx>0)then
        datasets(iset)%xy%dx(i)=dble_field(icol_dx,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of dx in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
        if(datasets(iset)%xy%dx(i)<=0.d0)then
          write(6,'(a)')'Found non-positive dx in "'//trim(fname)//'".'
          write(6,'()')
          ierr=-5
          exit
        endif
      endif ! have_dx
      if(icol_dy>0)then
        datasets(iset)%xy%dy(i)=dble_field(icol_dy,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of dx in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
        if(datasets(iset)%xy%dy(i)<=0.d0)then
          write(6,'(a)')'Found non-positive dy in "'//trim(fname)//'".'
          write(6,'()')
          ierr=-6
          exit
        endif
      endif ! have_dy
    enddo ! i

    ! Close file.
    close(io)

  END SUBROUTINE read_file


  ! POLYNOMIAL HANDLING UTILITIES.


  FUNCTION print_poly_sym(npoly,pow,X0) RESULT(polystr)
    !---------------------------------------------------!
    ! Returns a string with the symbolic version of the !
    ! fitted polynomial.                                !
    !---------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: npoly
    DOUBLE PRECISION,INTENT(in) :: pow(npoly),X0
    CHARACTER(3+npoly*54) :: polystr
    INTEGER j,ipow
    CHARACTER(40) pwstr
    CHARACTER(6) xstr
    if(abs(X0)<tol_zero)then
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


  FUNCTION print_poly_num(npoly,pow,a,X0) RESULT(polystr)
    !----------------------------------------------------!
    ! Returns a string with the numerical version of the !
    ! fitted polynomial in a suitable format for pasting !
    ! into xmgrace.                                      !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: npoly
    DOUBLE PRECISION,INTENT(in) :: pow(npoly),a(npoly),X0
    CHARACTER(2+npoly*105) :: polystr
    INTEGER j,ipow
    CHARACTER(1) plusstr
    CHARACTER(32) coeffstr
    CHARACTER(72) pwstr
    CHARACTER(36) xstr
    if(abs(X0)<tol_zero)then
      xstr='x'
    else
      write(xstr,*)X0
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


  !SUBROUTINE int_poly(npoly1,pow1,a1,npoly,pow,a)
  !  !------------------------------------------------!
  !  ! Integrate the polynomial given by npoly1, pow1 !
  !  ! and a1 and return the resulting polynomial in  !
  !  ! npoly, pow and a if present, or in-place       !
  !  ! overwriting the input if not.                  !
  !  !------------------------------------------------!
  !  IMPLICIT NONE
  !  INTEGER,INTENT(inout) :: npoly1
  !  DOUBLE PRECISION,INTENT(inout) :: pow1(:),a1(:)
  !  INTEGER,INTENT(out),OPTIONAL :: npoly
  !  DOUBLE PRECISION,INTENT(out),OPTIONAL :: pow(:),a(:)
  !  INTEGER npoly2
  !  DOUBLE PRECISION pow2(npoly1),a2(npoly1)
  !  INTEGER ipoly
  !  npoly2=npoly1
  !  do ipoly=1,npoly1
  !    a2(ipoly)=a1(ipoly)/(pow1(ipoly)+1.d0)
  !    pow2(ipoly)=pow1(ipoly)+1.d0
  !  enddo ! ipoly
  !  if(present(npoly))then
  !    ! Return result separately.
  !    npoly=npoly2
  !    pow(1:npoly)=pow2(1:npoly2)
  !    a(1:npoly)=a2(1:npoly)
  !  else
  !    ! Return result in place.
  !    npoly1=npoly2
  !    pow1(1:npoly1)=pow2(1:npoly2)
  !    a1(1:npoly1)=a2(1:npoly2)
  !  endif
  !END SUBROUTINE int_poly


  !SUBROUTINE square_poly(npoly1,pow1,a1,npoly,pow,a)
  !  !---------------------------------------------!
  !  ! Square the polynomial given by npoly1, pow1 !
  !  ! and a1 and return the resulting polynomial  !
  !  ! in npoly2, pow2 and a2 if present, or       !
  !  ! in-place overwriting the input if not.      !
  !  !---------------------------------------------!
  !  IMPLICIT NONE
  !  INTEGER,INTENT(inout) :: npoly1
  !  DOUBLE PRECISION,INTENT(inout) :: pow1(:),a1(:)
  !  INTEGER,INTENT(out),OPTIONAL :: npoly
  !  DOUBLE PRECISION,INTENT(out),OPTIONAL :: pow(:),a(:)
  !  INTEGER npoly2
  !  DOUBLE PRECISION pow2((npoly1*(npoly1+1))/2),a2((npoly1*(npoly1+1))/2)
  !  INTEGER ipoly,jpoly
  !  npoly2=0
  !  do ipoly=1,npoly1
  !    do jpoly=ipoly,npoly1
  !      npoly2=npoly2+1
  !      a2(npoly2)=a1(ipoly)*a1(jpoly)
  !      if(ipoly/=jpoly)a2(npoly2)=a2(npoly2)*2.d0
  !      pow2(npoly2)=pow1(ipoly)+pow1(jpoly)
  !    enddo
  !  enddo ! ipoly
  !  if(present(npoly))then
  !    ! Return result separately.
  !    npoly=npoly2
  !    pow(1:npoly)=pow2(1:npoly2)
  !    a(1:npoly)=a2(1:npoly)
  !  else
  !    ! Return result in place.
  !    npoly1=npoly2
  !    pow1(1:npoly1)=pow2(1:npoly2)
  !    a1(1:npoly1)=a2(1:npoly2)
  !  endif
  !END SUBROUTINE square_poly


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
        ! For natural powers, force integer exponents to avoid issues.
        eval_poly=eval_poly+a(ipoly)*x_target**nint(pow(ipoly))
      elseif(.not.are_equal(x_target,0.d0))then
        ! Fractional powers only evaluated for non-zero arguments.
        eval_poly=eval_poly+a(ipoly)*x_target**pow(ipoly)
      endif
    enddo ! ipoly
  END FUNCTION eval_poly


  ! DERIVED-TYPE TOOLS


  SUBROUTINE enlarge_xydata(nxy,have_dx,have_dy,xy)
    !--------------------------------------------!
    ! Enlarge an xydata-type pointer to (1:nxy), !
    ! preserving any existing data.              !
    !--------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy
    LOGICAL,INTENT(in) :: have_dx,have_dy
    TYPE(xydata),POINTER :: xy
    INTEGER i
    if(.not.associated(xy))then
      allocate(xy)
      xy%nxy=0
      xy%have_dx=have_dx
      xy%have_dy=have_dy
      nullify(xy%x)
      nullify(xy%y)
      nullify(xy%dx)
      nullify(xy%dy)
    endif
    ! Enlarge arrays.
    call resize_pointer_dble1((/nxy/),xy%x)
    call resize_pointer_dble1((/nxy/),xy%y)
    call resize_pointer_dble1((/nxy/),xy%dx)
    call resize_pointer_dble1((/nxy/),xy%dy)
    ! Initialize x and y.
    do i=xy%nxy+1,nxy
      xy%x(i)=dble(i)
      xy%y(i)=dble(i)
    enddo ! i
    ! Store dataset length.
    xy%nxy=nxy
  END SUBROUTINE enlarge_xydata


  SUBROUTINE kill_xydata(xy)
    !---------------------------------!
    ! Destroy an xydata-type pointer. !
    !---------------------------------!
    IMPLICIT NONE
    TYPE(xydata),POINTER :: xy
    if(associated(xy))then
      if(associated(xy%x))deallocate(xy%x)
      if(associated(xy%y))deallocate(xy%y)
      if(associated(xy%dx))deallocate(xy%dx)
      if(associated(xy%dy))deallocate(xy%dy)
      deallocate(xy)
      nullify(xy)
    endif
  END SUBROUTINE kill_xydata


  SUBROUTINE resize_dataset(ndataset,datasets)
    !------------------------------------------------!
    ! Resize a dataset-type pointer to (1:ndataset), !
    ! preserving any existing data.                  !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset),POINTER :: datasets(:)
    TYPE(dataset),POINTER :: tdatasets(:)
    INTEGER old_ndataset,i

    ! Backup original vector of pointers to datasets and destroy original.
    old_ndataset=0
    if(associated(datasets))then
      old_ndataset=size(datasets)
      if(ndataset>0)then
        allocate(tdatasets(min(old_ndataset,ndataset)))
        do i=1,min(old_ndataset,ndataset)
          tdatasets(i)=datasets(i)
        enddo ! i
      endif
      deallocate(datasets)
      if(ndataset<old_ndataset)then
        do i=ndataset+1,old_ndataset
          call kill_xydata(datasets(i)%xy)
          call kill_xydata(datasets(i)%txy)
          call kill_xydata(datasets(i)%rtxy)
        enddo ! i
      endif
    endif

    ! Create new vector of pointers to datasets.
    if(ndataset>0)allocate(datasets(ndataset))

    ! Copy vector of pointers to datasets from backup and destroy backup.
    if(old_ndataset>0)then
      if(ndataset>0)then
        do i=1,min(old_ndataset,ndataset)
          datasets(i)=tdatasets(i)
        enddo ! i
        deallocate(tdatasets)
        nullify(tdatasets)
      endif
    endif

    ! Nullify pointer if empty.
    if(ndataset==0)nullify(datasets)

  END SUBROUTINE resize_dataset


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


  LOGICAL ELEMENTAL FUNCTION are_equal(x,y,tol)
    !------------------------------------------------------!
    ! Check if two floating-point numbers are equal within !
    ! a reasonable tolerance.                              !
    !------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    DOUBLE PRECISION abs_x,abs_y,big,small
    ! Parameters.
    DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-12
    DOUBLE PRECISION,PARAMETER :: tol_rel=1.d-9
    if(present(tol))then
      are_equal=abs(x-y)<=tol
    else
      abs_x=abs(x)
      abs_y=abs(y)
      if(abs_x<=tol_zero.and.abs_y<=tol_zero)then
        are_equal=.true.
      elseif(x>0.d0.eqv.y>0.d0)then
        big=max(abs_x,abs_y)
        small=min(abs_x,abs_y)
        are_equal=big-small<=big*tol_rel
      else
        are_equal=.false.
      endif
    endif
  END FUNCTION are_equal


  LOGICAL FUNCTION are_equal_string(cx,cy,tol)
    !-------------------------------------------------------!
    ! Check if two strings are equal, or if their numerical !
    ! values are equal within a reasonable tolerance.       !
    !-------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: cx,cy
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    DOUBLE PRECISION x,y
    INTEGER ierr
    are_equal_string=trim(cx)==trim(cy)
    if(are_equal_string)return
    read(cx,*,iostat=ierr)x
    if(ierr/=0)return
    read(cy,*,iostat=ierr)y
    if(ierr/=0)return
    are_equal_string=are_equal(x,y,tol)
  END FUNCTION are_equal_string


  ! POINTER RESIZING TOOLS.


  SUBROUTINE resize_pointer_char1(sz,dims,pt,init)
    !------------------------------------------------------!
    ! Allocate or resize a first-rank character pointer PT !
    ! to size DIMS, keeping any existing data untouched.   !
    ! New elements initialized to empty string.            !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: sz,dims(1)
    CHARACTER(sz),INTENT(in),OPTIONAL :: init
    INTEGER old_dims(1)
    CHARACTER(0),PARAMETER :: init_default=''
    CHARACTER(sz),POINTER :: pt(:),pt_new(:)
    if(.not.associated(pt))then
      allocate(pt(dims(1)))
      if(present(init))then
        pt=init
      else
        pt=init_default
      endif
      return
    endif
    old_dims=shape(pt)
    if(all(old_dims==dims))return
    allocate(pt_new(dims(1)))
    if(any(old_dims<dims))then
      if(present(init))then
        pt_new=init
      else
        pt_new=init_default
      endif
    endif
    pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
    deallocate(pt)
    pt=>pt_new
    nullify(pt_new)
  END SUBROUTINE resize_pointer_char1


  SUBROUTINE resize_pointer_char2(sz,dims,pt,init)
    !-------------------------------------------------------!
    ! Allocate or resize a second-rank character pointer PT !
    ! to size DIMS, keeping any existing data untouched.    !
    ! New elements initialized to empty string.             !
    !----------------------------------------------------- -!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: sz,dims(2)
    CHARACTER(sz),INTENT(in),OPTIONAL :: init
    INTEGER old_dims(2)
    CHARACTER(0),PARAMETER :: init_default=''
    CHARACTER(sz),POINTER :: pt(:,:),pt_new(:,:)
    if(.not.associated(pt))then
      allocate(pt(dims(1),dims(2)))
      if(present(init))then
        pt=init
      else
        pt=init_default
      endif
      return
    endif
    old_dims=shape(pt)
    if(all(old_dims==dims))return
    allocate(pt_new(dims(1),dims(2)))
    if(any(old_dims<dims))then
      if(present(init))then
        pt_new=init
      else
        pt_new=init_default
      endif
    endif
    pt_new(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))=&
       &pt(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))
    deallocate(pt)
    pt=>pt_new
    nullify(pt_new)
  END SUBROUTINE resize_pointer_char2


  SUBROUTINE resize_pointer_int1(dims,pt,init)
    !-----------------------------------------------------!
    ! Allocate or resize a first-rank int pointer PT to   !
    ! size DIMS, keeping any existing data untouched. New !
    ! elements initialized to zero unless INIT specified. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dims(1)
    INTEGER,INTENT(in),OPTIONAL :: init
    INTEGER,POINTER :: pt(:),pt_new(:)=>null()
    INTEGER,PARAMETER :: init_default=0
    INTEGER old_dims(1)
    if(.not.associated(pt))then
      allocate(pt(dims(1)))
      if(present(init))then
        pt=init
      else
        pt=init_default
      endif
      return
    endif
    old_dims=shape(pt)
    if(all(old_dims==dims))return
    allocate(pt_new(dims(1)))
    if(any(old_dims<dims))then
      if(present(init))then
        pt_new=init
      else
        pt_new=init_default
      endif
    endif
    pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
    deallocate(pt)
    pt=>pt_new
    nullify(pt_new)
  END SUBROUTINE resize_pointer_int1


  SUBROUTINE resize_pointer_dble1(dims,pt,init)
    !-----------------------------------------------------!
    ! Allocate or resize a first-rank dble pointer PT to  !
    ! size DIMS, keeping any existing data untouched. New !
    ! elements initialized to zero unless INIT specified. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dims(1)
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: init
    DOUBLE PRECISION,POINTER :: pt(:),pt_new(:)=>null()
    DOUBLE PRECISION,PARAMETER :: init_default=0.d0
    INTEGER old_dims(1)
    if(.not.associated(pt))then
      allocate(pt(dims(1)))
      if(present(init))then
        pt=init
      else
        pt=init_default
      endif
      return
    endif
    old_dims=shape(pt)
    if(all(old_dims==dims))return
    allocate(pt_new(dims(1)))
    if(any(old_dims<dims))then
      if(present(init))then
        pt_new=init
      else
        pt_new=init_default
      endif
    endif
    pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
    deallocate(pt)
    pt=>pt_new
    nullify(pt_new)
  END SUBROUTINE resize_pointer_dble1


  ! STRING UTILITIES.


  FUNCTION field(n,line)
    !---------------------------------------------------------!
    ! Return the N-th field in string LINE, where the fields  !
    ! are separated by one or more spaces.  If N is negative, !
    ! return the |N|-th-to-last field in LINE.                !
    !---------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    CHARACTER(*),INTENT(in) :: line
    CHARACTER(len(line)) :: field
    CHARACTER(len(line)) tline
    INTEGER i,k,absn
    LOGICAL back
    ! FIXME - should honour and strip quotes.
    ! Initialize.
    field=''
    if(n==0)return
    absn=abs(n)
    tline=adjustl(line)
    back=(n<0)
    ! Loop over fields.
    do i=1,absn-1
      k=scan(trim(tline),' ',back)
      if(k==0)return
      if(back)then
        tline=adjustl(tline(1:k-1))
      else
        tline=adjustl(tline(k+1:))
      endif
    enddo
    ! Return nth field.
    k=scan(trim(tline),' ',back)
    if(k==0)then
      field=adjustl(tline)
    elseif(back)then
      field=adjustl(tline(k+1:))
    else
      field=adjustl(tline(1:k-1))
    endif
  END FUNCTION field


  INTEGER FUNCTION nfield(line)
    !--------------------------------------!
    ! Return the number of fields in LINE. !
    !--------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: line
    nfield=0
    do
      if(len_trim(field(nfield+1,line))==0)exit
      nfield=nfield+1
    enddo
  END FUNCTION nfield


  INTEGER FUNCTION int_field(ifield,command,ierr)
    !----------------------------------------------------!
    ! Like field, but returning the value as an integer. !
    ! ierr is set to a non-zero value if the requested   !
    ! field could not be parsed as an integer.           !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ifield
    CHARACTER(*),INTENT(in) :: command
    INTEGER,INTENT(inout) :: ierr
    CHARACTER(len(command)) f
    ierr=0
    f=field(ifield,command)
    read(f,*,iostat=ierr)int_field
  END FUNCTION int_field


  DOUBLE PRECISION FUNCTION dble_field(ifield,command,ierr)
    !--------------------------------------------------------!
    ! Like field, but returning the value as an real number. !
    ! ierr is set to a non-zero value if the requested field !
    ! could not be parsed as a real number.                  !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ifield
    CHARACTER(*),INTENT(in) :: command
    INTEGER,INTENT(inout) :: ierr
    CHARACTER(len(command)) f
    ierr=0
    f=field(ifield,command)
    read(f,*,iostat=ierr)dble_field
  END FUNCTION dble_field


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


  ! QUIT ROUTINE.


  SUBROUTINE quit (msg)
    !---------------------!
    ! Quit with an error. !
    !---------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in), OPTIONAL :: msg
    if (present(msg)) then
      write(6,'(a)')'ERROR : '//msg
    else
      write(6,'(a)')'Quitting.'
    endif
    stop
  END SUBROUTINE quit


END PROGRAM polyfit
