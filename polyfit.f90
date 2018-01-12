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

  ! Types.
  TYPE xydata
    INTEGER nxy
    LOGICAL have_dx,have_dy
    DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:),dx(:),dy(:)
  END TYPE xydata
  TYPE dataset
    INTEGER itransfx,itransfy
    TYPE(xydata),POINTER :: xy=>null() ! original data
    TYPE(xydata),POINTER :: txy=>null() ! transformed data
    TYPE(xydata),POINTER :: rtxy=>null() ! range-restricted transformed data
  END TYPE dataset
  TYPE fit_params
    INTEGER npoly
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

  call main()


CONTAINS


  SUBROUTINE main()
    !--------------!
    ! Main driver. !
    !--------------!
    IMPLICIT NONE
    ! (x,y[,dx][,dy]) datasets.
    INTEGER ndataset
    TYPE(xydata),POINTER :: xy=>null()
    TYPE(dataset),POINTER :: datasets(:)=>null(),temp_datasets(:)=>null()
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
            end select
            ifield=ifield+3
          case('')
            exit
          case default
            write(6,'(a)')'Syntax error in load command: unknown subcommand "'&
               &//trim(field(ifield,command))//'".'
            write(6,'()')
            cycle user_loop
          end select
        enddo ! ifield

        ! Check.
        if(ierr1/=0.or.ierr2/=0.or.ierr3/=0.or.ierr4/=0)then
          write(6,'(a)')'Problem parsing column indices.'
          write(6,'()')
          cycle user_loop
        endif
        if(icol_x>ncolumn.or.icol_x<0.or.icol_y>ncolumn.or.icol_y<0.or.&
           &icol_dx>ncolumn.or.icol_dx<0.or.icol_dy>ncolumn.or.icol_dy<0)then
          write(6,'(a)')'Column indices out of range.'
          write(6,'()')
          cycle user_loop
        endif

        ! Read data.
        allocate(xy)
        xy%nxy=nxy
        xy%have_dx=icol_dx>0
        xy%have_dy=icol_dy>0
        allocate(xy%x(nxy),xy%y(nxy),xy%dx(nxy),xy%dy(nxy))
        call read_file(fname,ncolumn,icol_x,icol_y,icol_dx,icol_dy,xy,ierr)

        ! Check data ar compatible with current transformations.
        if(ierr==0)then
          if(TRANSF_REQ_NONZERO(itransfx_default))then
            if(any(are_equal(xy%x,0.d0)))ierr=-7
          endif
          if(TRANSF_REQ_NONZERO(itransfy_default))then
            if(any(are_equal(xy%y,0.d0)))ierr=-8
          endif
          if(TRANSF_REQ_POSITIVE(itransfx_default))then
            if(any(xy%x<0.d0))ierr=-9
          endif
          if(TRANSF_REQ_POSITIVE(itransfy_default))then
            if(any(xy%y<0.d0))ierr=-10
          endif
        endif

        ! Bail out if there has been an error.
        if(ierr/=0)then
          select case(ierr)
          case(-1)
            write(6,'(a)')'Problem opening "'//trim(fname)//'".'
          case(-2)
            write(6,'(a)')'Unexpected enf of file reading "'//trim(fname)//'".'
          case(-3)
            write(6,'(a)')'Problem reading line in "'//trim(fname)//'".'
          case(-4)
            write(6,'(a)')'Problem parsing data in "'//trim(fname)//'".'
          case(-5)
            write(6,'(a)')'"'//trim(fname)//'" contains non-positive dx.'
          case(-6)
            write(6,'(a)')'"'//trim(fname)//'" contains non-positive dy.'
          case(-7)
            write(6,'(a)')'"'//trim(fname)//'" contains x=0, incompatible &
               &with current xscale.'
          case(-8)
            write(6,'(a)')'"'//trim(fname)//'" contains y=0, incompatible &
               &with current yscale.'
          case(-9)
            write(6,'(a)')'"'//trim(fname)//'" contains x<0, incompatible &
               &with current xscale.'
          case(-10)
            write(6,'(a)')'"'//trim(fname)//'" contains y<0, incompatible &
               &with current yscale.'
          end select
          write(6,'()')
          deallocate(xy%x,xy%y,xy%dx,xy%dy)
          deallocate(xy)
          nullify(xy)
          cycle user_loop
        endif

        if(ndataset>0)then
          ! Copy of vector of pointers to datasets and destroy original.
          allocate(temp_datasets(ndataset+1))
          do i=1,ndataset
            temp_datasets(i)=datasets(i)
          enddo
          deallocate(datasets)
        endif

        ! Create new vector of pointers to datasets.
        ndataset=ndataset+1
        allocate(datasets(ndataset))

        if(ndataset>1)then
          ! Copy vector of pointers to datasets from backup and destroy backup.
          do i=1,ndataset-1
            datasets(i)=temp_datasets(i)
          enddo
          deallocate(temp_datasets)
          nullify(temp_datasets)
        endif

        ! Point last element in vector of pointers to the new dataset.
        datasets(ndataset)%xy=>xy
        nullify(xy)

        ! Set other properties.
        datasets(ndataset)%itransfx=itransfx_default
        datasets(ndataset)%itransfy=itransfy_default

        ! Populate transformed dataset.
        call refresh_dataset(datasets(ndataset),drange)
        ! Update X0.
        call refresh_fit(ndataset,datasets,fit)

        ! Report success.
        write(6,'(a)',advance='no')'Loaded data from "'//trim(fname)//&
           &'" as dataset #'//trim(i2s(ndataset))//', type '
        if(datasets(ndataset)%xy%have_dx.and.datasets(ndataset)%xy%have_dy)then
          write(6,'(a)',advance='no')'xdxydy'
        elseif(datasets(ndataset)%xy%have_dx)then
          write(6,'(a)',advance='no')'xdxy'
        elseif(datasets(ndataset)%xy%have_dy)then
          write(6,'(a)',advance='no')'xydy'
        else
          write(6,'(a)',advance='no')'xy'
        endif
        write(6,'(a)')', '//trim(i2s(datasets(ndataset)%xy%nxy))//' data.'
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
          write(6,'(2x,a)')trim(print_poly_sym(fit%npoly,fit%pow,fit%X0))
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
          write(6,'(a)')'Set shared coefficients.'
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
          write(6,'()')
          write(6,'(a)')'Unknown variable "'//trim(field(2,command))//'".'
          write(6,'()')
          cycle user_loop
        end select ! variable to set

      case('plot')
        if(ndataset<1)then
          write(6,'()')
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif

      case('evaluate')
        if(ndataset<1)then
          write(6,'()')
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        deval%nderiv=0
        deval%rel=.false.
        deval%var='X'
        deval%n=1
        allocate(deval%x(deval%n))
        deval%x(1)=0.d0
        call evaluate_fit(ndataset,datasets,fit,mcparams,drange,deval)
        deallocate(deval%x)

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
    DOUBLE PRECISION chi2,chi2err
    DOUBLE PRECISION,ALLOCATABLE :: a(:,:),da(:,:)

    ! Initialize.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)
    tot_nxy=0
    do iset=1,ndataset
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset

    ! Allocate space for fit parameters.
    allocate(a(fit%npoly,ndataset),da(fit%npoly,ndataset))

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

    ! Clean up.
    deallocate(a,da)

  END SUBROUTINE show_multipoly


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
    INTEGER tot_nxy,tot_nparam,iset,npoly,ix
    DOUBLE PRECISION chi2,chi2err
    DOUBLE PRECISION,ALLOCATABLE :: fmean(:,:),ferr(:,:)

    ! Initialize.
    allocate(fmean(deval%n,ndataset),ferr(deval%n,ndataset))

    ! Prepare combined dataset arrays.
    tot_nxy=0
    do iset=1,ndataset
      tot_nxy=tot_nxy+datasets(iset)%rtxy%nxy
    enddo ! iset

    ! Loop over expansion orders.
    write(6,'()')
    write(6,'(2x,a5,1x,2(1x,a12))',advance='no')'Order','chi^2/Ndf',&
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
        write(6,'(2x,i5,1x,2(1x,es12.4))',advance='no')npoly-1,&
           &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
        do ix=1,deval%n
          do iset=1,ndataset
            if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,es12.4))',&
               &advance='no')fmean(1,iset),ferr(1,iset)
          enddo ! iset
        enddo ! ix
        write(6,'()')
      endif
      ! Destroy work arrays.
      deallocate(tfit%pow,tfit%share)
      deallocate(tfit)
      nullify(tfit)
    enddo ! npoly
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(32))//'("-"))')
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
    INTEGER ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,ix
    INTEGER,ALLOCATABLE :: indx(:)
    DOUBLE PRECISION chi2,chi2err
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:),fmean(:,:),ferr(:,:)

    ! FIXME - expose these parameters.
    drange%var='X'
    drange%op='<='
    drange%size=0

    ! Initialize.
    allocate(fmean(deval%n,ndataset),ferr(deval%n,ndataset))
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

    ! Loop over data to remove.
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
        write(6,'(2x,i5,1x,2(1x,es12.4))',advance='no')tot_nxy,&
           &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
        do ix=1,deval%n
          do iset=1,ndataset
            if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,es12.4))',&
               &advance='no')fmean(1,iset),ferr(1,iset)
          enddo ! iset
        enddo ! ix
        write(6,'()')
      endif
    enddo ! igrid
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(32+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(32))//'("-"))')
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
    INTEGER ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,npoly,ix
    INTEGER,ALLOCATABLE :: indx(:)
    DOUBLE PRECISION chi2,chi2err
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:),fmean(:,:),ferr(:,:)

    ! FIXME - expose these parameters.
    drange%var='X'
    drange%op='<='
    drange%size=0

    ! Initialize.
    allocate(fmean(deval%n,ndataset),ferr(deval%n,ndataset))

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

    ! Loop over data to remove.
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
          write(6,'(2x,i5,2x,i5,1x,2(1x,es12.4))',advance='no')&
             &tot_nxy,npoly-1,&
             &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
          do ix=1,deval%n
            do iset=1,ndataset
              if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,es12.4))',&
                 &advance='no')fmean(1,iset),ferr(1,iset)
            enddo ! iset
          enddo ! ix
          write(6,'()')
        endif
        ! Destroy work arrays.
        deallocate(tfit%pow,tfit%share)
        deallocate(tfit)
        nullify(tfit)
      enddo ! npoly
    enddo ! nxy
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(39))//'("-"))')
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


  !SUBROUTINE show_dual_assessment(nxy,have_dx,have_dy,x,y,dx,dy,itransfx,&
  !   &itransfy,tx0,npoly_out,pow_out,mask_out)
  !  !-----------------------------------------------------------!
  !  ! Perform an "extrapolation" assessment in which the number !
  !  ! of points included in the fit and the expansion order are !
  !  ! simultaneously optimized.                                 !
  !  !-----------------------------------------------------------!
  !  IMPLICIT NONE
  !  INTEGER,INTENT(in) :: nxy,itransfx,itransfy
  !  LOGICAL,INTENT(in) :: have_dx,have_dy
  !  DOUBLE PRECISION,INTENT(in) :: x(nxy),y(nxy),dx(nxy),dy(nxy),tx0
  !  INTEGER,INTENT(inout) :: npoly_out
  !  DOUBLE PRECISION,INTENT(inout) :: pow_out(nxy)
  !  LOGICAL,INTENT(inout) :: mask_out(nxy)
  !  DOUBLE PRECISION tx(nxy),ty(nxy),dtx(nxy),dty(nxy),weight(nxy),&
  !     &rtx(nxy),rty(nxy),rweight(nxy),xtarget,min_chi2,t1,&
  !     &rx(nxy),ry(nxy),rdx(nxy),rdy(nxy),list_f(nxy),list_df(nxy),&
  !     &list_gof(nxy),list_a(nxy,nxy),list_da(nxy,nxy),&
  !     &pow(nxy),a(nxy),da(nxy),txrange,txplot,dtxplot
  !  INTEGER i,ntest,npoly,rnxy,nderiv,ierr,indx(nxy),rnxy_min_chi2,nrandom,&
  !     &npoly_best,rnxy_best,ipoly,jpoly,list_rnxy(nxy)
  !  LOGICAL mask(nxy),weighted
  !  CHARACTER(20) drop_by,drop_criterion
  !  CHARACTER(2048) char2048
  !  DOUBLE PRECISION fmean(1),ferr(1),best_f,best_df
  !  ! Parameters.
  !  INTEGER,PARAMETER :: DEFAULT_NRANDOM=100000
  !  INTEGER,PARAMETER :: npoint=1000
  !  INTEGER,PARAMETER :: io=10

  !  ! Initialize.
  !  npoly_out=0

  !  ! Initialize internal variables.  FIXME - expose to user.
  !  nderiv=0
  !  xtarget=0.d0
  !  nrandom=DEFAULT_NRANDOM
  !  drop_by='X'
  !  drop_criterion='largest'

  !  ! Evaluate transformed variables and fit weights.
  !  call scale_transform(nxy,itransfx,x,tx,have_dx,dx,dtx)
  !  call scale_transform(nxy,itransfy,y,ty,have_dy,dy,dty)
  !  if(have_dx.and.have_dy)then
  !    weight=1.d0/(dtx*dty)**2
  !  elseif(have_dx)then
  !    weight=1.d0/dtx**2
  !  elseif(have_dy)then
  !    weight=1.d0/dty**2
  !  else
  !    weight=1.d0
  !  endif
  !  weighted=have_dx.or.have_dy

  !  ! Perform sort.
  !  select case(trim(drop_by))
  !  case('X')
  !    call isort(nxy,tx,indx)
  !  case('Y')
  !    call isort(nxy,ty,indx)
  !  case('x')
  !    call isort(nxy,x,indx)
  !  case('y')
  !    call isort(nxy,y,indx)
  !  case('index')
  !    indx=(/(i,i=1,nxy)/)
  !  end select

  !  ! Assess integer expansion orders up to order 8.
  !  ntest=min(nxy-1,9)
  !  pow(1:nxy)=(/(dble(i-1),i=1,nxy)/)

  !  ! Loop over expansion orders.
  !  write(6,'(1x,a,t8,a,t14,a,t27,a,t44,a)')'Order','  Nxy','  chi^2/Ndf',&
  !     &'    Target f','    Target df'
  !  write(6,'(1x,59("-"))')
  !  list_rnxy=0
  !  do npoly=2,ntest
  !    ! Loop over number of points in fit.
  !    rnxy_min_chi2=0
  !    min_chi2=0.d0
  !    do rnxy=nxy,npoly+1,-1
  !      ! Build data mask.
  !      mask=.false.
  !      select case(trim(drop_criterion))
  !      case('smallest')
  !        mask(indx(nxy-rnxy+1:nxy))=.true.
  !      case('largest')
  !        mask(indx(1:rnxy))=.true.
  !      end select
  !      ! Apply data mask.
  !      rtx(1:rnxy)=pack(tx,mask)
  !      rty(1:rnxy)=pack(ty,mask)
  !      rweight(1:rnxy)=pack(weight,mask)
  !      ! Fit to polynomial of order npoly-1.
  !      call perform_fit(rnxy,npoly,rtx-tx0,rty,pow,a,weighted,rweight,da)
  !      ! Evaluate chi^2.
  !      t1=chi_squared(rnxy,npoly,rtx-tx0,rty,rweight,pow,a,weighted)/&
  !         &dble(rnxy-npoly)
  !      ! Report.
  !      write(6,'(1x,i5,1x,i5,1x,es12.4)')npoly-1,rnxy,t1
  !      if(rnxy_min_chi2==0.or.t1<min_chi2)then
  !        rnxy_min_chi2=rnxy
  !        min_chi2=t1
  !      endif
  !    enddo ! rnxy
  !    ! Evaluate estimate at number of points that minimizes gof measure.
  !    rnxy=rnxy_min_chi2
  !    ! Build data mask.
  !    mask=.false.
  !    select case(trim(drop_criterion))
  !    case('smallest')
  !      mask(indx(nxy-rnxy+1:nxy))=.true.
  !    case('largest')
  !      mask(indx(1:rnxy))=.true.
  !    end select
  !    ! Apply data mask.
  !    rx(1:rnxy)=pack(x,mask)
  !    ry(1:rnxy)=pack(y,mask)
  !    rdx(1:rnxy)=pack(dx,mask)
  !    rdy(1:rnxy)=pack(dy,mask)
  !    ! Fit to polynomial of order npoly-1.
  !    call eval_fit_monte_carlo(rnxy,have_dx,have_dy,rx,ry,rdx,rdy,&
  !       &itransfx,itransfy,tx0,npoly,pow,nrandom,nderiv,.false.,1,&
  !       &(/xtarget/),fmean,ferr,amean=a,aerr=da)
  !    ! Report.
  !    write(6,'(1x,i5,1x,i5,1x,es12.4,2(1x,es16.8))')npoly-1,rnxy_min_chi2,&
  !       &min_chi2,fmean(1),ferr(1)
  !    ! Store.
  !    list_f(npoly)=fmean(1)
  !    list_df(npoly)=ferr(1)
  !    list_rnxy(npoly)=rnxy_min_chi2
  !    list_gof(npoly)=min_chi2
  !    list_a(1:npoly,npoly)=a(1:npoly)
  !    list_da(1:npoly,npoly)=da(1:npoly)
  !  enddo ! npoly
  !  write(6,'(1x,59("-"))')
  !  write(6,*)

  !  ! Find best choice.
  !  npoly_best=0
  !  best_f=0.d0
  !  best_df=0.d0
  !  rnxy_best=0
  !  do ipoly=1,ntest
  !    if(list_rnxy(ipoly)==0)cycle
  !    do jpoly=ipoly+1,ntest
  !      if(list_rnxy(jpoly)==0)cycle
  !      if(list_gof(jpoly)>list_gof(ipoly))cycle
  !      if(abs(list_f(ipoly)-list_f(jpoly))>&
  !         &1.0d0*sqrt(list_df(ipoly)**2+list_df(jpoly)**2))exit
  !    enddo ! jpoly
  !    if(jpoly>ntest)then
  !      npoly_best=ipoly
  !      best_f=list_f(ipoly)
  !      best_df=list_df(ipoly)
  !      rnxy_best=list_rnxy(ipoly)
  !      exit
  !    endif
  !  enddo ! ipoly

  !  if(npoly_best==0)then
  !    write(6,*)'Test found no good choice.'
  !  else
  !    write(6,*)'Best ('//trim(i2s(npoly_best-1))//','//&
  !       &trim(i2s(rnxy_best))//'):',best_f,'+/-',best_df
  !  endif
  !  write(6,*)

  !  ! Allow user to choose one of the above expansion orders.
  !  do
  !    if(npoly_best==0)then
  !      write(6,*)'Type ''plot'' to plot, empty to skip:'
  !    else
  !      write(6,*)'Type ''y'' to accept fit form and mask, ''plot'' to plot, &
  !         &empty to skip:'
  !    endif
  !    read(5,'(a)',iostat=ierr)char2048
  !    if(ierr/=0)exit
  !    write(6,*)
  !    select case(trim(char2048))
  !    case('')
  !      exit
  !    case('plot')
  !      ! Prepare to write plot.
  !      txrange=maxval(tx)-minval(tx)
  !      dtxplot=(txrange*2.d0)/dble(npoint)
  !      open(unit=io,file='poly_dual.dat',status='replace')
  !      ! Plot original data.
  !      if(have_dx.and.have_dy)then
  !        write(io,'(a)')'@type xydxdy'
  !        do i=1,nxy
  !          write(io,*)tx(i),ty(i),dtx(i),dty(i)
  !        enddo ! i
  !      elseif(have_dx)then
  !        write(io,'(a)')'@type xydx'
  !        do i=1,nxy
  !          write(io,*)tx(i),ty(i),dtx(i)
  !        enddo ! i
  !      elseif(have_dy)then
  !        write(io,'(a)')'@type xydy'
  !        do i=1,nxy
  !          write(io,*)tx(i),ty(i),dty(i)
  !        enddo ! i
  !      else
  !        write(io,'(a)')'@type xy'
  !        do i=1,nxy
  !          write(io,*)tx(i),ty(i)
  !        enddo ! i
  !      endif
  !      write(io,'(a)')'&'
  !      ! Loop over expansion orders.
  !      do npoly=2,ntest
  !        write(io,'(a)')'@type xy'
  !        ! Dump plot.
  !        txplot=minval(tx)-0.5d0*txrange
  !        do i=0,npoint
  !          write(io,*)txplot+tx0,eval_poly(npoly,pow,list_a(1,npoly),txplot)
  !          txplot=txplot+dtxplot
  !        enddo ! i
  !        write(io,'(a)')'&'
  !      enddo ! npoly
  !      close(io)
  !      write(6,*)'Data written to poly_dual.dat.'
  !      write(6,*)
  !    case('y')
  !      if(npoly_best>0)then
  !        npoly_out=npoly_best
  !        pow_out(1:npoly_out)=(/(i,i=0,npoly_out-1)/)
  !        rnxy=rnxy_best
  !        mask_out=.false.
  !        select case(trim(drop_criterion))
  !        case('smallest')
  !          mask_out(indx(nxy-rnxy+1:nxy))=.true.
  !        case('largest')
  !          mask_out(indx(1:rnxy))=.true.
  !        end select
  !      endif
  !      exit
  !    case default
  !      exit
  !    end select
  !  enddo

  !END SUBROUTINE show_dual_assessment


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
    CHARACTER(1024) char1024
    INTEGER i,ipos,ierr
    DOUBLE PRECISION t1
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Initialize.
    nline=-1 ! flag non-existing file
    ncolumn=0

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)return
    nline=0

    ! Look for first line containing data.
    do
      read(io,'(a)',iostat=ierr)char1024
      if(ierr<0)exit
      if(ierr>0)then
        nline=-2 ! flag reading error
        exit
      endif
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
      if(ncolumn<1)then
        nline=-3 ! flag column count problem
        exit
      endif
      nline=nline+1
    enddo

    ! Close file.
    close(io)

  END SUBROUTINE check_file


  SUBROUTINE read_file(fname,ncolumn,icol_x,icol_y,icol_dx,icol_dy,xy,ierr)
    !-------------------------------------!
    ! Read in the data in the input file. !
    !-------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ncolumn,icol_x,icol_y,icol_dx,icol_dy
    TYPE(xydata),POINTER :: xy
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(inout) :: ierr
    DOUBLE PRECISION tvec(ncolumn)
    INTEGER i,ipos,pline
    CHARACTER(1024) line
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)then
      ierr=-1 ! flag problem opening file
      return
    endif

    ! Initialize.
    pline=0
    xy%x=(/(dble(i),i=1,xy%nxy)/)
    xy%y=(/(dble(i),i=1,xy%nxy)/)
    xy%dx=0.d0
    xy%dy=0.d0

    ! Loop over data lines.
    do i=1,xy%nxy
      ! Load next line containing data.
      do
        read(io,'(a)',iostat=ierr)line
        if(ierr<0)then
          ierr=-2 ! flag unexpected EOF
          exit
        endif
        pline=pline+1
        if(ierr>0)then
          ierr=-3 ! flag reading error
          exit
        endif
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
      if(ierr/=0)then
        ierr=-4 ! flag data parsing error
        exit
      endif
      ! Set data and check that error bars are positive.
      if(icol_x>0)xy%x(i)=tvec(icol_x)
      if(icol_y>0)xy%y(i)=tvec(icol_y)
      if(icol_dx>0)then
        xy%dx(i)=tvec(icol_dx)
        if(xy%dx(i)<=0.d0)then
          ierr=-5 ! flag non-positive dx
          exit
        endif
      endif ! have_dx
      if(icol_dy>0)then
        xy%dy(i)=tvec(icol_dy)
        if(xy%dy(i)<=0.d0)then
          ierr=-6 ! flag non-poistive dy
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


  ! MISC UTILITIES.


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
