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

  ! Derived types.
  TYPE xy_type
    INTEGER :: nxy=0
    LOGICAL :: have_dx=.false.,have_dy=.false.,have_w=.false.
    DOUBLE PRECISION,POINTER :: x(:)=>null(),y(:)=>null(),dx(:)=>null(),&
       &dy(:)=>null(),w(:)=>null()
  END TYPE xy_type
  TYPE dataset_type
    INTEGER :: itransfx=ITRANSF_NONE,itransfy=ITRANSF_NONE
    DOUBLE PRECISION :: wexp=1.d0
    DOUBLE PRECISION :: weight=1.d0
    TYPE(xy_type),POINTER :: xy=>null() ! original data
    TYPE(xy_type),POINTER :: txy=>null() ! transformed data
    TYPE(xy_type),POINTER :: rtxy=>null() ! range-restricted transformed data
  END TYPE dataset_type
  TYPE dataset_list_type
    TYPE(dataset_type),POINTER :: dataset=>null()
  END TYPE dataset_list_type
  TYPE fit_form_type
    INTEGER :: npoly=0
    DOUBLE PRECISION :: X0=0.d0
    DOUBLE PRECISION,ALLOCATABLE :: pow(:)
    LOGICAL,ALLOCATABLE :: share(:)
    CHARACTER(64) :: X0_string='0'
  END TYPE fit_form_type
  TYPE range_type
    CHARACTER :: var='X'
    CHARACTER(2) :: op=''
    LOGICAL :: no_rhs=.false.
    DOUBLE PRECISION :: thres=0.d0
    INTEGER :: size=0
  END TYPE range_type
  TYPE eval_type
    CHARACTER :: var='X'
    CHARACTER(10) :: what='function'
    LOGICAL :: rel=.false.
    INTEGER :: nderiv=0,n=0
    DOUBLE PRECISION,POINTER :: x(:)=>null()
  END TYPE eval_type
  TYPE mc_params_type
    ! NB, changing the default here does nothing since this is reset in main().
    INTEGER :: nsample=10000
  END TYPE mc_params_type

  call main()


CONTAINS


  SUBROUTINE main()
    !--------------!
    ! Main driver. !
    !--------------!
    IMPLICIT NONE
    ! (x,y[,dx][,dy]) datasets.
    INTEGER ndataset,file_ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:),tmp_dlist(:),file_dlist(:)
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Fit form.
    TYPE(fit_form_type),POINTER :: fit
    ! All-set settings.
    INTEGER itransfx_default,itransfy_default
    DOUBLE PRECISION wexp_default
    TYPE(range_type) drange,tdrange
    ! Function evaluation.
    INTEGER eval_iset
    TYPE(eval_type) deval
    ! Monte Carlo parameters.
    TYPE(mc_params_type) mcparams
    ! Input file information.
    INTEGER ncolumn,icol_x,icol_y,icol_dx,icol_dy,icol_w
    CHARACTER(256) fname
    ! Input file search facility.
    INTEGER,PARAMETER :: search_size=1024
    INTEGER nsearch,ndiscr
    INTEGER,POINTER :: fsearch(:),fdiscr(:)
    CHARACTER(search_size),POINTER :: search(:)
    ! Misc variables.
    CHARACTER(7) settype
    CHARACTER(8192) command,token
    LOGICAL,ALLOCATABLE :: smask(:)
    INTEGER i,ifield,ierr,ipos_x,ipos_y,ipos_dx,ipos_dy,ipos_w,&
       &ierr1,ierr2,ierr3,ierr4,ierr5
    INTEGER nxy,itransf,iset,i1,i2,ipos,npoly
    DOUBLE PRECISION t1,t2,wexp

    ! Initialize.
    ndataset=0
    nullify(dlist)
    itransfx_default=ITRANSF_NONE
    itransfy_default=ITRANSF_NONE
    wexp_default=1.d0
    allocate(fit)
    fit%npoly=2
    allocate(fit%pow(fit%npoly),fit%share(fit%npoly))
    fit%pow(1:2)=(/0.d0,1.d0/)
    fit%share=.false.
    mcparams%nsample=10000
    nullify(fsearch,fdiscr,search)

    ! Write header.
    write(6,'(a)')'===================================='
    write(6,'(a)')'POLYFIT - polynomial fitting toolbox'
    write(6,'(a)')'===================================='
    write(6,'()')
    write(6,'(a)')'Type "help" for a list of commands.'
    write(6,'()')

    ! Loop over user actions.
    user_loop: do

      ! Clean up any mess from previous loops.
      if(associated(search))deallocate(search)
      if(associated(fsearch))deallocate(fsearch)
      if(associated(fdiscr))deallocate(fdiscr)
      if(associated(deval%x))deallocate(deval%x)
      nullify(search)
      nullify(fsearch)
      nullify(fdiscr)
      nullify(deval%x)

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

      case('load','wload')
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
        if(trim(field(1,command))=='wload')then
          settype='y'
        else
          select case(ncolumn)
          case(1)
            settype='y'
          case(2)
            settype='xy'
          case(3)
            settype='xydy'
          case(4)
            settype='xydyw'
          case default
            settype='xdxydyw'
          end select
        endif
        call parse_type_string(settype,icol_x,icol_dx,icol_y,icol_dy,icol_w,&
           &ierr)

        ! Parse subcommands.
        ifield=2
        do
          ifield=ifield+1
          select case(trim(field(ifield,command)))
          case('using') ! lone "using" for wload only
            if(trim(field(1,command))=='load')then
              write(6,'(a)')'Syntax error in load command: "using" without &
                 &preceeding "type" not allowed.'
              write(6,'()')
              cycle user_loop
            endif
            icol_y=int_field(ifield+1,command,ierr1)
            ifield=ifield+1
            if(ierr1/=0)then
              write(6,'(a)')'Problem parsing "using" index.'
              write(6,'()')
              cycle user_loop
            endif
            ! Check column index.
            if(icol_y>ncolumn.or.icol_y<0)then
              write(6,'(a)')'Column indices out of range.'
              write(6,'()')
              cycle user_loop
            endif
          case('type')
            if(trim(field(1,command))=='wload')then
              write(6,'(a)')'Syntax error in wload command: "type" is not &
                 &an allowed subcommand.'
              write(6,'()')
              cycle user_loop
            endif
            call parse_type_string(field(ifield+1,command),ipos_x,ipos_dx,&
               &ipos_y,ipos_dy,ipos_w,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error in load command: unrecognized &
                 &dataset type.'
              write(6,'()')
              cycle user_loop
            endif
            if(trim(field(ifield+2,command))=='using')then
              ierr1=0
              ierr2=0
              ierr3=0
              ierr4=0
              ierr5=0
              icol_x=0
              icol_dx=0
              icol_y=0
              icol_dy=0
              icol_w=0
              if(ipos_x>0)icol_x=int_field(ifield+2+ipos_x,command,ierr1)
              if(ipos_dx>0)icol_dx=int_field(ifield+2+ipos_dx,command,ierr2)
              if(ipos_y>0)icol_y=int_field(ifield+2+ipos_y,command,ierr3)
              if(ipos_dy>0)icol_dy=int_field(ifield+2+ipos_dy,command,ierr4)
              if(ipos_w>0)icol_w=int_field(ifield+2+ipos_w,command,ierr5)
              if(ierr1/=0.or.ierr2/=0.or.ierr3/=0.or.ierr4/=0.or.ierr5/=0)then
                write(6,'(a)')'Problem parsing "using" indices.'
                write(6,'()')
                cycle user_loop
              endif
              ifield=ifield+2
              if(ipos_x>0)ifield=ifield+1
              if(ipos_dx>0)ifield=ifield+1
              if(ipos_y>0)ifield=ifield+1
              if(ipos_dy>0)ifield=ifield+1
              if(ipos_w>0)ifield=ifield+1
            else
              icol_x=ipos_x
              icol_dx=ipos_dx
              icol_y=ipos_y
              icol_dy=ipos_dy
              icol_w=ipos_w
              ifield=ifield+1
            endif
            ! Check column indices.
            if(icol_x>ncolumn.or.icol_x<0.or.icol_dx>ncolumn.or.icol_dx<0.or.&
               &icol_y>ncolumn.or.icol_y<0.or.icol_dy>ncolumn.or.icol_dy<0.or.&
               &icol_w>ncolumn.or.icol_w<0)then
              write(6,'(a)')'Column indices out of range.'
              write(6,'()')
              cycle user_loop
            endif

          case('where')
            ! Add a search clause.
            if(nfield(command)<ifield+2)then
              write(6,'(a)')'Syntax error in load command: too few arguments &
                 &for "where" subcommand"'
              write(6,'()')
              cycle user_loop
            endif
            i=int_field(ifield+1,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error in load command: invalid column &
                 &index for "where".'
              write(6,'()')
              cycle user_loop
            endif
            if(i<1.or.i>ncolumn)then
              write(6,'(a)')'Column index out of range in "where" subcommand.'
              write(6,'()')
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
            if(trim(field(1,command))=='wload')then
              write(6,'(a)')'Syntax error in wload command: "by" is not &
                 &an allowed subcommand.'
              write(6,'()')
              cycle user_loop
            endif
            if(nfield(command)<ifield+1)then
              write(6,'(a)')'Syntax error in load command: too few arguments &
                 &for "by" subcommand"'
              write(6,'()')
              cycle user_loop
            endif
            i=int_field(ifield+1,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error in load command: invalid column &
                 &index for "by".'
              write(6,'()')
              cycle user_loop
            endif
            if(i<1.or.i>ncolumn)then
              write(6,'(a)')'Column index out of range in "by" subcommand.'
              write(6,'()')
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
            cycle user_loop
          end select
        enddo ! ifield

        ! Read data.
        if(.not.associated(search))allocate(fsearch(nsearch),search(nsearch))
        if(.not.associated(fdiscr))allocate(fdiscr(ndiscr))
        call read_file(fname,icol_x,icol_y,icol_dx,icol_dy,icol_w,nsearch,&
           &fsearch,search,ndiscr,fdiscr,file_ndataset,file_dlist,ierr)
        if(ierr/=0)cycle user_loop
        if(file_ndataset<1)then
          write(6,'(a)')'No data loaded.'
          write(6,'()')
          cycle user_loop
        endif

        if(field(1,command)=='wload')then
          ! We are loading dataset weights, so lets see if we have the
          ! right number of them.
          if(file_ndataset/=1)then
            write(6,'(a)')'Problem with wload: loaded too many datasets.'
            write(6,'()')
            call kill_dlist(file_dlist)
            cycle user_loop
          endif
          if(file_dlist(1)%dataset%xy%nxy/=ndataset)then
            write(6,'(a)')'Problem with wload: loaded '//&
               &trim(i2s(file_dlist(1)%dataset%xy%nxy))//' weights but &
               &expected '//trim(i2s(ndataset))//'.'
            write(6,'()')
            call kill_dlist(file_dlist)
            cycle user_loop
          endif
          ! Check that dataset weights are > 0.
          if(any(file_dlist(1)%dataset%xy%y<=0.d0))then
            write(6,'(a)')'Problem with wload: found negative weights.'
            write(6,'()')
            call kill_dlist(file_dlist)
            cycle user_loop
          endif
          ! Set weights.
          do iset=1,ndataset
            dlist(iset)%dataset%weight=file_dlist(1)%dataset%xy%y(iset)
          enddo ! iset
          ! Clean up and loop here to avoid giant if-block below.
          call kill_dlist(file_dlist)
          write(6,'(a)')'Loaded data from "'//trim(fname)//&
             &'" as dataset weights.'
          write(6,'()')
          cycle user_loop
        endif

        ! Check data ar compatible with current transformations.
        do iset=1,file_ndataset
          dataset=>file_dlist(iset)%dataset
          xy=>dataset%xy
          dataset%itransfx=itransfx_default
          dataset%itransfy=itransfy_default
          dataset%wexp=wexp_default
          if(TRANSF_REQ_NONZERO(dataset%itransfy))then
            if(any(are_equal(xy%x,0.d0)))then
              dataset%itransfx=ITRANSF_NONE
              write(6,'(a)')'Note: using linear X for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains x=0.'
              write(6,'()')
            endif
          endif
          if(TRANSF_REQ_NONZERO(dataset%itransfy))then
            if(any(are_equal(xy%y,0.d0)))then
              dataset%itransfy=ITRANSF_NONE
              write(6,'(a)')'Note: using linear Y for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains y=0.'
              write(6,'()')
            endif
          endif
          if(TRANSF_REQ_POSITIVE(dataset%itransfx))then
            if(any(xy%x<0.d0))then
              dataset%itransfx=ITRANSF_NONE
              write(6,'(a)')'Note: using linear X for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains x<0.'
              write(6,'()')
            endif
          endif
          if(TRANSF_REQ_POSITIVE(dataset%itransfy))then
            if(any(xy%y<0.d0))then
              dataset%itransfy=ITRANSF_NONE
              write(6,'(a)')'Note: using linear Y for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains y<0.'
              write(6,'()')
            endif
          endif
          if(dataset%xy%have_w)then
            if(any(xy%w<0.d0.and..not.are_equal(xy%w,0.d0)))then
              if(.not.are_equal(dataset%wexp,0.d0))then
                dataset%wexp=0.d0
                write(6,'(a)')'Note: zeroing weight exponent for set '//&
                   &trim(i2s(iset+ndataset))//' since it contains w<0.'
                write(6,'()')
              endif
            elseif(any(are_equal(xy%w,0.d0)))then
              if(dataset%wexp<0.d0.and..not.are_equal(dataset%wexp,0.d0))then
                dataset%wexp=1.d0
                write(6,'(a)')'Note: setting weight exponent to 1 for set '//&
                   &trim(i2s(iset+ndataset))//' since it contains w=0.'
                write(6,'()')
              endif
            endif
          endif
        enddo ! iset

        ! Add datasets to list.
        if(ndataset>0)then
          allocate(tmp_dlist(ndataset))
          tmp_dlist=dlist(1:ndataset)
          deallocate(dlist)
        endif
        allocate(dlist(ndataset+file_ndataset))
        if(ndataset>0)then
          dlist(1:ndataset)=tmp_dlist
          deallocate(tmp_dlist)
        endif
        ndataset=ndataset+file_ndataset
        dlist(ndataset-file_ndataset+1:ndataset)=file_dlist(1:file_ndataset)
        deallocate(file_dlist)
        do iset=ndataset-file_ndataset+1,ndataset
          dataset=>dlist(iset)%dataset
          call refresh_dataset(dataset,drange)
          ! Report.
          xy=>dataset%xy
          write(6,'(a)')'Loaded data from "'//trim(fname)//'" as dataset #'//&
             &trim(i2s(iset))//', type '//trim(type_string(xy%have_dx,&
             &xy%have_dy,xy%have_w))//', '//trim(i2s(xy%nxy))//' data.'
        enddo ! iset
        call refresh_fit(ndataset,dlist,fit)
        write(6,'()')

      case('unload')
        write(6,'()')
        if(ndataset<1)then
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        allocate(smask(ndataset))
        if(nfield(command)==1)then
          ! Unload all.
          smask=.true.
        else
          ! Unload specified sets.
          smask=.false.
          ifield=1
          do
            ifield=ifield+1
            if(ifield>nfield(command))exit
            i=parse_int(field(ifield,command),ierr)
            if(ierr/=0)then
              write(6,'(a)')'Invalid dataset index.'
              write(6,'()')
              deallocate(smask)
              cycle user_loop
            endif
            if(i<1.or.i>ndataset)then
              write(6,'(a)')'Dataset index out of range.'
              write(6,'()')
              deallocate(smask)
              cycle user_loop
            endif
            smask(i)=.true.
          enddo
        endif
        ! Delete datasets backwards.
        if(.not.all(smask))then
          allocate(tmp_dlist(count(.not.smask)))
          tmp_dlist=pack(dlist,.not.smask)
        endif
        do i=1,ndataset
          if(.not.smask(i))cycle
          call kill_dataset(dlist(i)%dataset)
        enddo
        deallocate(dlist)
        nullify(dlist)
        if(.not.all(smask))then
          allocate(dlist(count(.not.smask)))
          dlist=tmp_dlist
          deallocate(tmp_dlist)
        endif
        ndataset=count(.not.smask)
        call refresh_fit(ndataset,dlist,fit)
        write(6,'(a)')trim(i2s(count(smask)))//' datasets unloaded.'
        write(6,'()')
        deallocate(smask)

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
        call show_multipoly(ndataset,dlist,drange,fit,mcparams)

      case('plot')
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
          deval%what='function'
          deval%nderiv=0
        case("f'")
          deval%what='function'
          deval%nderiv=1
        case("f''")
          deval%what='function'
          deval%nderiv=2
        case("sharedf")
          deval%what='shared'
          deval%nderiv=0
        case("sharedf'")
          deval%what='shared'
          deval%nderiv=1
        case("sharedf''")
          deval%what='shared'
          deval%nderiv=2
        case default
          write(6,'(a)')'Syntax error: unknown function "'//&
             &trim(field(2,command))//'".'
          write(6,'()')
          cycle user_loop
        end select
        ! Parse sub-commands.
        fname='fit.plot'
        ifield=2
        do
          ifield=ifield+1
          if(ifield>nfield(command))exit
          select case(trim(field(ifield,command)))
          case('to')
            if(nfield(command)<ifield+1)then
              write(6,'(a)')'Syntax error: "to" subcommand must be followed &
                 &by a filename.'
              write(6,'()')
              cycle user_loop
            endif
            fname=field(ifield+1,command)
            ifield=ifield+1
          case('at')
            if(nfield(command)<ifield+1)then
              write(6,'(a)')'Syntax error: "at" subcommand must be followed &
                 &by a data range.'
              write(6,'()')
              cycle user_loop
            endif
            call parse_xeval(trim(field(ifield+1,command)),deval)
            if(.not.associated(deval%x))then
              write(6,'(a)')'Syntax error: could not parse range.'
              write(6,'()')
              cycle user_loop
            endif
            ifield=ifield+1
          case default
            write(6,'(a)')'Syntax error: unknown subcommand "'//&
               &trim(field(ifield,command))//'".'
            write(6,'()')
            cycle user_loop
          end select
        enddo
        call plot_multipoly(ndataset,dlist,drange,fit,deval,mcparams,&
           &trim(fname))

      case('assess')
        if(ndataset<1)then
          write(6,'()')
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        ! Check assessment target.
        select case(trim(field(2,command)))
        case('fit','range','range,fit','fit,range')
        case default
          write(6,'()')
          write(6,'(a)')'Unknown variable to assess "'//&
             &trim(field(2,command))//'".'
          write(6,'()')
          cycle user_loop
        end select
        ! Initialize.
        deval%n=0
        deval%rel=.false.
        eval_iset=0
        tdrange%var='X'
        tdrange%op='<='
        tdrange%thres=0.d0
        tdrange%size=0
        tdrange%no_rhs=.true.
        ! Loop over subcommands.
        ifield=2
        do
          ifield=ifield+1
          if(ifield>nfield(command))exit
          select case(field(ifield,command))
          case('using')
            ! Get object to evaluate.
            select case(field(ifield+1,command))
            case("f")
              deval%what='function'
              deval%nderiv=0
            case("f'")
              deval%what='function'
              deval%nderiv=1
            case("f''")
              deval%what='function'
              deval%nderiv=2
            case("sumf")
              deval%what='sum'
              deval%nderiv=0
            case("sumf'")
              deval%what='sum'
              deval%nderiv=1
            case("sumf''")
              deval%what='sum'
              deval%nderiv=2
            case default
              write(6,'(a)')'Syntax error: unknown function "'//&
                 &trim(field(ifield+1,command))//'".'
              write(6,'()')
              cycle user_loop
            end select
            ! Get where to evaluate it at.
            if(trim(field(ifield+2,command))/='at')then
              write(6,'(a)')'Syntax error: missing "at" subcommand.'
              write(6,'()')
              cycle user_loop
            endif
            call parse_xeval(field(ifield+3,command),deval)
            if(.not.associated(deval%x))then
              write(6,'(a)')'Syntax error: problem parsing list of X values.'
              write(6,'()')
              cycle user_loop
            endif
            ifield=ifield+3
          case('by')
            call parse_range(field(ifield+1,command),tdrange)
            if(tdrange%op=='')then
              write(6,'(a)')'Syntax error: problem parsing "by" argument.'
              write(6,'()')
              cycle user_loop
            endif
            ifield=ifield+1
          case('for')
            i=int_field(ifield+1,command,ierr)
            if(ierr/=0)then
              write(6,'(a)')'Syntax error: invalid set index "'//&
                 &trim(field(ifield+1,command))//'" specified.'
              write(6,'()')
              cycle user_loop
            endif
            if(i<1.or.i>ndataset)then
              write(6,'(a)')'Syntax error: set index out of range.'
              write(6,'()')
              cycle user_loop
            endif
            eval_iset=i
            ifield=ifield+1
          case default
            write(6,'(a)')'Syntax error: unknown subcommand "'//&
                 &trim(field(ifield,command))//'".'
            write(6,'()')
            cycle user_loop
          end select
        enddo ! ifield
        ! Evaluate.
        select case(trim(field(2,command)))
        case('fit')
          call assess_fit(ndataset,dlist,drange,fit,mcparams,deval,&
             &eval_iset)
        case('range')
          call assess_range(ndataset,dlist,tdrange,fit,mcparams,deval,&
             &eval_iset)
        case('range,fit','fit,range')
          call assess_fit_range(ndataset,dlist,tdrange,fit,mcparams,deval,&
             &eval_iset)
        end select

      case('report')
        if(ndataset<1)then
          write(6,'()')
          write(6,'(a)')'No datasets loaded.'
          write(6,'()')
          cycle user_loop
        endif
        select case(trim(field(2,command)))
        case('range')
          call report_statistics(ndataset,dlist)
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
          deval%what='function'
          deval%nderiv=0
        case("f'")
          deval%what='function'
          deval%nderiv=1
        case("f''")
          deval%what='function'
          deval%nderiv=2
        case("sumf")
          deval%what='sum'
          deval%nderiv=0
        case("sumf'")
          deval%what='sum'
          deval%nderiv=1
        case("sumf''")
          deval%what='sum'
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
        call parse_xeval(field(4,command),deval)
        if(.not.associated(deval%x))then
          write(6,'(a)')'Syntax error: could not parse range.'
          write(6,'()')
          cycle user_loop
        endif
        ! Perform evaluation.
        call evaluate_fit(ndataset,dlist,fit,mcparams,drange,deval)

      case('intersect')
        ! This is only useful when two or more datasets are loaded.
        if(ndataset<2)then
          write(6,'(a)')'Need at least two datasets to find intersections.'
          write(6,'()')
          cycle user_loop
        endif
        ! Parse options.
        select case(field(2,command))
        case("between")
          t1=dble_field(3,command,ierr1)
          t2=dble_field(4,command,ierr2)
          if(ierr1/=0.or.ierr2/=0)then
            write(6,'(a)')'Syntax error: could not parse arguments of &
               &"between" subcommand.'
            write(6,'()')
            cycle user_loop
          endif
          if(t1>=t2)then
            write(6,'(a)')'Syntax error: intersection range has non-positive &
               &length.'
            write(6,'()')
            cycle user_loop
          endif
        case default
          write(6,'(a)')'Syntax error: "between" subcommand missing.'
          write(6,'()')
          cycle user_loop
        end select
        ! Perform intersection.
        call intersect_fit(ndataset,dlist,fit,mcparams,drange,t1,t2)

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
              xy=>dlist(iset)%dataset%xy
              if(TRANSF_REQ_NONZERO(itransf))then
                select case(trim(field(2,command)))
                case('xscale')
                  if(any(are_equal(xy%x,0.d0)))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x=0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(are_equal(xy%y,0.d0)))then
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
                  if(any(xy%x<0.d0))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x<0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(xy%y<0.d0))then
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
              dataset=>dlist(iset)%dataset
              select case(trim(field(2,command)))
              case('xscale')
                dataset%itransfx=itransf
              case('yscale')
                dataset%itransfy=itransf
              end select
              write(6,'(a)')'Set '//trim(field(2,command))//' to '//&
                 &trim(TRANSF_NAME(itransf))//' for set #'//trim(i2s(iset))//'.'
            enddo ! ifield
          else
            ! Check transformation is applicable.
            do iset=1,ndataset
              xy=>dlist(iset)%dataset%xy
              if(TRANSF_REQ_NONZERO(itransf))then
                select case(trim(field(2,command)))
                case('xscale')
                  if(any(are_equal(xy%x,0.d0)))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x=0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(are_equal(xy%y,0.d0)))then
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
                  if(any(xy%x<0.d0))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains x<0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                case('yscale')
                  if(any(xy%y<0.d0))then
                    write(6,*)'Cannot apply axis transformation: set #'//&
                       &trim(i2s(iset))//' contains y<0.'
                    write(6,'()')
                    cycle user_loop
                  endif
                end select
              endif
            enddo ! iset
            ! Store transformation.
            select case(trim(field(2,command)))
            case('xscale')
              itransfx_default=itransf
              do iset=1,ndataset
                dlist(iset)%dataset%itransfx=itransf
              enddo ! iset
            case('yscale')
              itransfy_default=itransf
              do iset=1,ndataset
                dlist(iset)%dataset%itransfy=itransf
              enddo ! iset
            end select
            write(6,'(a)')'Set '//trim(field(2,command))//' to '//&
               &trim(TRANSF_NAME(itransf))//' for all sets.'
          endif
          write(6,'()')
          ! Apply transformations.
          do iset=1,ndataset
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! iset
          ! Update X0.
          call refresh_fit(ndataset,dlist,fit)

        case('wexp')
          ! Set weight exponent.
          wexp=parse_dble(field(3,command),ierr)
          if(ierr/=0)then
            write(6,'(a)')'Problem parsing value of weight exponent.'
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
            ! Check exponent is applicable.
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
              if(wexp<0.d0.and..not.are_equal(wexp,0.d0))then
                xy=>dlist(iset)%dataset%xy
                if(any(are_equal(xy%w,0.d0)))then
                  write(6,*)'Cannot apply weight exponent: set #'//&
                     &trim(i2s(iset))//' contains w=0.'
                  write(6,'()')
                  cycle user_loop
                endif
              endif
            enddo ! ifield
            ! Set exponent.
            do ifield=5,nfield(command)
              iset=int_field(ifield,command,ierr)
              dataset=>dlist(iset)%dataset
              dataset%wexp=wexp
              write(6,'(a)')'Set wexp for set #'//trim(i2s(iset))//'.'
            enddo ! ifield
          else
            ! Check exponent is applicable.
            ! FIXME - create functions eq_dble, le_dble, lt_dble, etc.
            if(wexp<0.d0.and..not.are_equal(wexp,0.d0))then
              do iset=1,ndataset
                xy=>dlist(iset)%dataset%xy
                if(any(are_equal(xy%w,0.d0)))then
                  write(6,*)'Cannot apply weight exponent: set #'//&
                     &trim(i2s(iset))//' contains w=0.'
                  write(6,'()')
                  cycle user_loop
                endif
              enddo ! iset
            endif ! wexp<0
            ! Set exponent.
            wexp_default=wexp
            do iset=1,ndataset
              dlist(iset)%dataset%wexp=wexp
            enddo ! iset
            write(6,'(a)')'Set wexp for all sets.'
          endif
          write(6,'()')

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
              i1=parse_int(token(1:ipos-1),ierr)
              if(ierr/=0)then
                write(6,'(a)')'Syntax error: could not parse range.'
                write(6,'()')
                cycle user_loop
              endif
            endif
            i2=parse_int(token(ipos+1:),ierr)
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
          write(6,'(2x,a)')trim(print_poly_sym(fit))
          write(6,'(2x,a)')'Shared coefficients reset to: none'
          write(6,'()')

        case('range')
          ! Set fit range.
          ! FIXME - "for" clause
          ! Check sort variable.
          call parse_range(trim(field(3,command)),drange)
          if(drange%op=='')then
            write(6,'(a)')'Syntax error parsing range string.'
            write(6,'()')
            cycle user_loop
          endif
          if(drange%no_rhs)then
            write(6,'(a)')'Syntax error parsing right-hand side of range.'
            write(6,'()')
            drange%op=''
            drange%no_rhs=.false.
            cycle user_loop
          endif
          ! Apply transformations.
          do iset=1,ndataset
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! iset
          ! Update X0.
          call refresh_fit(ndataset,dlist,fit)
          ! Report.
          write(6,'(a)')'Range set.'
          write(6,'()')

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
              i1=int_field(3,command,ierr)
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
                i1=parse_int(token(1:ipos-1),ierr)
                if(ierr/=0)then
                  write(6,'(a)')'Syntax error: could not parse range.'
                  write(6,'()')
                  cycle user_loop
                endif
              endif
              i2=npoly
              if(ipos<len_trim(token))then
                i2=parse_int(token(ipos+1:),ierr)
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

        case('centre')
          ! Check value.
          ! FIXME - "for" clause
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
          call refresh_fit(ndataset,dlist,fit)

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

      case('unset')
        ! Set variables.
        write(6,'()')
        select case(trim(field(2,command)))
        case('xscale','yscale')
          ! FIXME - "for" clause
          select case(trim(field(2,command)))
          case('xscale')
            itransfx_default=ITRANSF_NONE
            do iset=1,ndataset
              dlist(iset)%dataset%itransfx=ITRANSF_NONE
            enddo ! ifield
          case('yscale')
            itransfy_default=ITRANSF_NONE
            do iset=1,ndataset
              dlist(iset)%dataset%itransfy=ITRANSF_NONE
            enddo ! ifield
          end select
          do iset=1,ndataset
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! iset
          call refresh_fit(ndataset,dlist,fit)
        case('wexp')
          ! FIXME - "for" clause
          wexp_default=1.d0
          do iset=1,ndataset
            dlist(iset)%dataset%wexp=1.d0
          enddo ! ifield
        case('fit')
          deallocate(fit%pow,fit%share)
          fit%npoly=2
          allocate(fit%pow(2),fit%share(2))
          fit%pow=(/0.d0,1.d0/)
          fit%share=.false.
        case('range')
          drange%var='X'
          drange%op=''
          drange%thres=0.d0
          drange%size=0
          do iset=1,ndataset
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! iset
          call refresh_fit(ndataset,dlist,fit)
        case('shared')
          fit%share=.false.
        case('centre')
          fit%X0_string='0'
          call refresh_fit(ndataset,dlist,fit)
        case('nsample')
          mcparams%nsample=10000
        case default
          write(6,'(a)')'Unknown variable "'//trim(field(2,command))//'".'
          write(6,'()')
          cycle user_loop
        end select
        write(6,'(a)')'Variable "'//trim(field(2,command))//'" reset.'
        write(6,'()')

      case('status')
        write(6,'()')
        write(6,'(a,es11.4)')'Global settings:'
        select case(drange%op)
        case('')
          write(6,'(a,es11.4)')'  Data range: all data'
        case('[')
          write(6,'(a,es11.4)')'  Data range: first '//&
             &trim(i2s(drange%size))//' data by '//trim(drange%var)//' value'
        case(']')
          write(6,'(a,es11.4)')'  Data range: last '//&
             &trim(i2s(drange%size))//' data by '//trim(drange%var)//' value'
        case default
          write(6,'(a,es11.4)')'  Data range: '//trim(drange%var)//' '//&
             &trim(drange%op)//' ',drange%thres
        end select
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
          dataset=>dlist(iset)%dataset
          xy=>dataset%xy
          write(6,'(a)',advance='no')'* Set #'//trim(i2s(iset))//': '//&
             &trim(i2s(xy%nxy))
          if(xy%nxy/=dataset%rtxy%nxy)write(6,'(a)',advance='no')' ('//&
             &trim(i2s(dataset%rtxy%nxy))//')'
          write(6,'(a)',advance='no')' '//&
             &trim(type_string(xy%have_dx,xy%have_dy,xy%have_w))//' data, '
          write(6,'(a)',advance='no')&
             &trim(TRANSF_NAME(dataset%itransfx))//'-'//&
             &trim(TRANSF_NAME(dataset%itransfy))//' scale'
          write(6,'(a,es10.4)',advance='no')', wexp=',dataset%wexp
          write(6,'(a,es10.4)',advance='no')', sw=',dataset%weight
          write(6,'(a)')'.'
        enddo ! i
        write(6,'()')
        write(6,'(a)')'Fit function:'
        write(6,'(2x,a)')trim(print_poly_sym(fit))
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
          call pprint('')
          call pprint('POLYFIT is a toolbox for performing polynomial fits on &
             &one or more sets of data.  It handles non-integer exponents, a &
             &few common data transformations and data with statistical &
             &errorbars.  POLYFIT can provide confidence intervals for values &
             &and derivatives of a fit, can perform fits to multiple &
             &datasets simultaneously with shared parameters, and provides &
             &useful tools for assessing quality of fits and automatically &
             &setting the fit form and data ranges.')
          call pprint('')
          call pprint('POLYFIT uses a command-line interface.  The list &
             &of available commands is:')
          call pprint('')
          call pprint('* inspect <file>',0,2)
          call pprint('* load <file> [type <type> using <columns>] &
             &[where <column> <value>] [by <column>]',0,2)
          call pprint('* wload <file> [using <column>] [where <column> &
             &<value>]',0,2)
          call pprint('* unload <set-index>',0,2)
          call pprint('* set <variable> <value> [for <set-list>]',0,2)
          call pprint('* unset <variable>',0,2)
          call pprint('* status',0,2)
          call pprint('* assess <variables> [using <function> at X <X> &
             &[for <set>]]',0,2)
          call pprint('* report <report>',0,2)
          call pprint('* fit',0,2)
          call pprint('* plot <file-name>',0,2)
          call pprint('* evaluate <function> at X <X>',0,2)
          call pprint('* help [<command> | set <variable>]',0,2)
          call pprint('')
          call pprint('Type help <command> for detailed information.')
          call pprint('')
        case('inspect')
          call pprint('')
          call pprint('Command: inspect <file>',0,9)
          call pprint('')
          call pprint('Reports the number of data lines and columns detected &
             &in <file>.',2,2)
          call pprint('')
        case('load')
          call pprint('')
          call pprint('Command: load <file> [type <type> [using <columns>]] &
             &[where <column> <value>] [by <column>]',0,9)
          call pprint('')
          call pprint('Loads data from <file> into a new dataset.  By &
             &default:',2,2)
          call pprint('* 1-column files are of type y',2,4)
          call pprint('* 2-column files are of type xy',2,4)
          call pprint('* 3-column files are of type xydy',2,4)
          call pprint('* 4-column files are of type xydyw',2,4)
          call pprint('* 5- or more-column files are of type xdxydyw',2,4)
          call pprint('Other types can be specified with an "type" clause. &
             &POLYFIT assumes by default that data can be found in the &
             &order given by <type> starting at column 1.  Other columns can &
             &be specified with the "using" clause; column indices are &
             &assigned to the variables in the order given by <type>. &
             &Specifying a column index of zero for x or y is equivalent to &
             &omitting x or y from <type>, and causes the variable to be set &
             &to the line index.',2,2)
          call pprint('')
          call pprint('The "where" subcommand restricts the file parsing to &
             &those lines in <file> where column <column> takes the value &
             &<value>.  <value> can be a string.  If multiple "where" &
             &subcommands are specified, only lines for which ALL &
             &specified columns take the required values are loaded.',2,2)
          call pprint('')
          call pprint('The "by" subcommand allows loading multiple datasets &
             &from <file>, each corresponding to a different value of column &
             &<column>.  If multiple "by" commands are specified, two lines &
             &belong to different datasets if ANY of the specified columns &
             &differs in value.',2,2)
          call pprint('')
          call pprint('Example:',2,2)
          call pprint('')
          call pprint('load "../data.dat" type ywdy using 3 5 4 &
             &where 2 "good" where 7 1/4 by 1',4,6)
          call pprint('')
          call pprint('This will read file "../data.dat", loading y, dy, and &
             &w from colums 3, 4, and 5, setting x to the data-point index. &
             &Lines whose column 2 does not contain "good" or whose column 7 &
             &does not contain 0.25 are skipped, and the loaded data are &
             &split into individual datasets for each distinct value of &
             &column 1.  Note that data-point indices for each dataset run &
             &independently, i.e., the first value of x in each dataset will &
             &be 1.',2,2)
          call pprint('')
        case('wload')
          call pprint('')
          call pprint('Command: wload <file> [using <column>] [where <column> &
             &<value>]',0,9)
          call pprint('')
          call pprint('Loads global dataset weights from column <column> &
             &(column 1 by default) of <file>.  These weights are applied to &
             &all data in a dataset during fitting, and are used in the "sum" &
             &operation as linear coefficients.  E.g., "evaluate sumf at X=0" &
             &on two weighted datasets will compute w1*f1(0) + w2*f2(0).',2,2)
          call pprint('')
          call pprint('The "where" subcommand restricts the file parsing to &
             &those lines in <file> where column <column> takes the value &
             &<value>.  <value> can be a string.  If multiple "where" &
             &subcommands are specified, only lines for which ALL &
             &specified columns take the required values are loaded.',2,2)
          call pprint('')
        case('unload')
          call pprint('')
          call pprint('Command: unload [<set1> [<set2> [...]]]',0,9)
          call pprint('')
          call pprint('Unload datasets from memory.  If no sets are &
             &specified, all datasets are unloaded.',2,2)
          call pprint('')
        case('status')
          call pprint('')
          call pprint('Command: status',0,9)
          call pprint('')
          call pprint('Report currently loaded datasets and values of &
             &internal variables.',2,2)
          call pprint('')
        case('assess')
          call pprint('')
          call pprint('Command: assess <variables> [by <criterion>] [using &
             &<function> at <x-values> [for <set>]]',0,9)
          call pprint('')
          call pprint('Assess the convergence of the fit with the specified &
             &variables.',2,2)
          call pprint('')
          call pprint('The assessment is carried out based on the value of &
             &chi^2/Ndf and, if specified, on the value of a function of the &
             &fit.',2,2)
          call pprint('')
          call pprint('The following <variables> can be specified:',2,2)
          call pprint('* fit : assess convergence with choice of fit form.',&
             &2,4)
          call pprint('* range : assess convergence with data range.',2,4)
          call pprint('* fit,range : assess convergence with choice of fit &
             &form and data range',2,4)
          call pprint('')
          call pprint('The "fit" assessment is carried out by successively &
             &introducing terms of the currently-defined fit in the order &
             &they were given.',2,2)
          call pprint('')
          call pprint('The "range" assessment is performed by successively &
             &restricting the data range according to <criterion>.  &
             &<criterion> has the form "<variable><selector>", where &
             &<variable> is one of x, y, X, or Y (lowercase for original and &
             &uppercase for transformed values), and <selector> is one of:',&
             &2,2)
          call pprint('* <T',2,4)
          call pprint('* <=T',2,4)
          call pprint('* >=T',2,4)
          call pprint('* >T',2,4)
          call pprint('* [N (first N)',2,4)
          call pprint('* ]N (last N)',2,4)
          call pprint('Note that in the above, T and N are a literal "T" &
             &and a literal "N", respectively.',2,2)
          call pprint('')
          call pprint('<function> can be f, f'', or f'''' for the value, &
             &first, and second derivative of the fit, respectively, or &
             &sumf, sumf'', or sumf'''' for the sum over datasets of the &
             &value, first, and second derivative of the fit, respectively.',&
             &2,2)
          call pprint('')
          call pprint('<xvalues> is specified as &
             &"<variable>=<comma-separated-list>" (e.g., "X=0,1,2,3") or as &
             &"<variable>=<first>:<last>:<count>" (e.g., "X=0:3:4").',2,2)
          call pprint('')
        case('report')
          call pprint('')
          call pprint('Command: report <report>',0,9)
          call pprint('')
          call pprint('Produce the requested report of the loaded datasets.  &
             &Available reports are:',2,2)
          call pprint('* range : report basic range statistics of the data.',&
             &2,4)
          call pprint('')
        case('fit')
          call pprint('')
          call pprint('Command: fit',0,9)
          call pprint('')
          call pprint('Perform fit of currently loaded datasets.',2,2)
          call pprint('')
        case('plot')
          call pprint('')
          call pprint('Command: plot <function> [at <xvalues>] &
             &[to <filename>]',0,9)
          call pprint('')
          call pprint('Plot a function of the fit to <filename>.',2,2)
          call pprint('')
          call pprint('<function> can be f, f'', or f'''' for the value, &
             &first, and second derivative, respectively, or sharedf, &
             &sharedf'', or sharedf'''' for the value, first, and second &
             &derivative of the shared part of the fit function (for &
             &multi-dataset fits containing shared parameters). If &
             &<function> is f or sharedf, the original data are also &
             &plotted.',2,2)
          call pprint('')
          call pprint('<xvalues> is specified as &
             &"<variable>=<comma-separated-list>" (e.g., "X=0,1,2,3") or as &
             &"<variable>=<first>:<last>:<count>" (e.g., "X=0:3:4").',2,2)
          call pprint('')
          call pprint('If any fit parameters are shared among datasets, the &
             &fit function is split into a shared part and a set-specific &
             &part -- data points are offset by the value of the set-specific &
             &part, and the shared part of the fit is plotted.',2,2)
          call pprint('')
        case('evaluate')
          call pprint('')
          call pprint('Command: evaluate <function> at <xvalues>',0,9)
          call pprint('')
          call pprint('Evaluate a function of the fit and print the value.',&
             &2,2)
          call pprint('')
          call pprint('<function> can be:',2,2)
          call pprint('* f for the fit value',2,4)
          call pprint('* f'' for the first derivative of the fit function',2,4)
          call pprint('* f'''' for the second derivative of the fit function',&
             &2,4)
          call pprint('* sumf for the sum of fit values over all detasets &
             &(weighted sum if dataset weights loaded)',2,4)
          call pprint('* sumf'' for the sum of the first derivative of the &
             &fit function over all datasets (weighted sum if dataset weights &
             &loaded)',2,4)
          call pprint('* sumf'''' for the sum of the second derivative of the &
             &fit function over all datasets (weighted sum if dataset weights &
             &loaded)',2,4)
          call pprint('')
          call pprint('<xvalues> is specified as &
             &"<variable>=<comma-separated-list>" (e.g., "X=0,1,2,3") or as &
             &"<variable>=<first>:<last>:<count>" (e.g., "X=0:3:4").',2,2)
          call pprint('')
        case('intersect')
          call pprint('')
          call pprint('Command: intersect between <X1> <X2>',0,9)
          call pprint('')
          call pprint('Evaluate the average location of the intersections &
             &between each pair of datasets in the range X1:X2.',2,2)
          call pprint('')
          call pprint('The intersect command finds the intersection between &
             &the fits to each pair of datasets and computes the average &
             &location of the intersection.  This operation requires two or &
             &more datasets to be loaded.  The difference between each pair &
             &of fits is required to be of opposite signs at X1 and X2.  &
             &This command will return an error message if this is &
             &consistently not the case during random sampling.',2,2)
          call pprint('')
        case('unset')
          call pprint('')
          call pprint('Command: unset <variable>',0,9)
          call pprint('')
          call pprint('Sets <variable> to its default value.',2,2)
          call pprint('')
          call pprint('Type "help set <variable>" for detailed information on &
             &variables and their default values.',2,2)
          call pprint('')
        case('set')
          if(nfield(command)==2)then
            call pprint('')
            call pprint('Command: set <variable> <value> [for <set-list>]',0,9)
            call pprint('')
            call pprint('Sets <variable> to <value>, either globally or for &
               &selected sets (for certain variables).  The list of available &
               &variables is:',2,2)
            call pprint('')
            call pprint('* xscale',2,4)
            call pprint('* yscale',2,4)
            call pprint('* fit',2,4)
            call pprint('* range',2,4)
            call pprint('* shared',2,4)
            call pprint('* centre',2,4)
            call pprint('* nsample',2,4)
            call pprint('')
            call pprint('Type "help set <variable>" for detailed &
               &information.',2,2)
            call pprint('')
          else
            select case(field(3,command))
            case('xscale','yscale')
              call pprint('')
              call pprint('Variable: xscale, yscale',0,10)
              call pprint('')
              call pprint('"xscale" and "yscale" set scale transformations &
                 &for the independent x variable and for the dependent y &
                 &variable, respectively.  In POLYFIT notation, the original &
                 &variables are called x and y, and the transformed variables &
                 &are called X and Y.',2,2)
              call pprint('')
              call pprint('The allowed values for "xscale" and "yscale" are:',&
                 &2,2)
              call pprint('* "linear" : X = x',2,4)
              call pprint('* "reciprocal" : X = 1/x [x=0 forbidden]',2,4)
              call pprint('* "logarithmic" : X = log(x) [x<=0 forbidden]',2,4)
              call pprint('* "exponential" : X = exp(x)',2,4)
              call pprint('')
              call pprint('These variables can be set in a per-set manner or &
                 &globally.  Note that the global value applies to all loaded &
                 &datasets and becomes the default for new datasets.  POLYFIT &
                 &will refuse to apply a transformation to datasets &
                 &containing incompatible data, e.g., "set xscale logarithmic &
                 &for 1" is not allowed if dataset #1 contains negative x &
                 &values.',2,2)
              call pprint('')
              call pprint('The default value of "xscale" and "yscale" &
                 &is "linear".',2,2)
              call pprint('')
            case('wexp')
              call pprint('')
              call pprint('Variable: wexp',0,10)
              call pprint('')
              call pprint('"wexp" sets the exponent of the data weights to &
                 &be used in the fits, so that the i-th data point is &
                 &weighted by w_i^wexp. The default value of wexp is 1.',2,2)
              call pprint('')
              call pprint('"wexp" can be set globally or per dataset using &
                 &the "for" clause.  Note that the global value applies to &
                 &all loaded datasets and becomes the default for new &
                 &datasets.  POLYFIT will refuse to apply a negative weight &
                 &exponent transformation to datasets containing zero &
                 &weights.',2,2)
              call pprint('')
              call pprint('Note that one can use wexp=-2 to perform the usual &
                 &dy^-2 weighting used in "chi-squared" fits, e.g.,',2,2)
              call pprint('')
              call pprint('load "file.dat" type xydyw using 1 2 3 3',4,6)
              call pprint('set wexp -2',4,6)
              call pprint('')
              call pprint('The repeated index in the "load" command hints at &
                 &the double-counting problem with this type of fit, which, &
                 &along with its inability to handle data with an uncertainty &
                 &of zero, is the reason why we do not recommend using this &
                 &technique.',2,2)
              call pprint('')
            case('fit')
              call pprint('')
              call pprint('Variable: fit',0,10)
              call pprint('')
              call pprint('"fit" sets the exponents to be used in the fit &
                 &function.  The value can be specified as a list, e.g., "set &
                 &fit 0 1.5 3", or as an integer range, e.g., "set fit 0:2".  &
                 &Note that the form of the fitting function is global, so &
                 &the "for <set-list>" syntax does not apply to this &
                 &variable.',2,2)
              call pprint('')
              call pprint('The default value of "fit" is "0 1", &
                 &corresponding to a linear fit.',2,2)
              call pprint('')
            case('range')
              call pprint('')
              call pprint('Variable: range',0,10)
              call pprint('')
              call pprint('"range" sets the data mask to apply to all &
                 &datasets.  The value can be specified as &
                 &"<variable><selector>", where:',2,2)
              call pprint('* <variable> is one of x, y, X, or Y (lowercase &
                 &for original and uppercase for transformed values).',2,4)
              call pprint('* <selector> is one of:',2,4)
              call pprint('- <T',4,6)
              call pprint('- <=T',4,6)
              call pprint('- >=T',4,6)
              call pprint('- >T',4,6)
              call pprint('- [N (first N)',4,6)
              call pprint('- ]N (last N)',4,6)
              call pprint('Note that in the above, T is a real-valued &
                 &threshold, and N is an integer number.')
              call pprint('')
              call pprint('Note that "range" is a global variable, so the &
                 &"for <set-list>" syntax does not apply to this variable.',&
                 &2,2)
              call pprint('')
              call pprint('The default value of "range" is <unset>, which &
                 &selects all data.',2,2)
              call pprint('')
            case('shared')
              call pprint('')
              call pprint('Variable: shared',0,10)
              call pprint('')
              call pprint('"shared" specifies which coefficients in the fit &
                 &are shared among datasets.  The value can be specified as a &
                 &list of coefficient indices (e.g. "set shared 2 3 4"), or &
                 &as a range (e.g., "set shared 2:4"), or using the special &
                 &values "none" or "all".  Coefficients not flagged as shared &
                 &take different values for each dataset.',2,2)
              call pprint('')
              call pprint('The default value of "shared" is "none".',2,2)
              call pprint('')
            case('centre')
              call pprint('')
              call pprint('Variable: centre',0,10)
              call pprint('')
              call pprint('"centre" determines the X offset to be used in &
                 &fitting.  While polynomial fits are analytically identical &
                 &irrespective of the offset, numerical fits might benefit &
                 &from choosing, e.g., an X offset in the middle of the data &
                 &range, close to the data minimum, etc.',2,2)
              call pprint('')
              call pprint('Allowed values are:',2,2)
              call pprint('* left : smallest X in datasets',2,4)
              call pprint('* right : largest X in datasets',2,4)
              call pprint('* centre : midpoint between "left" and &
                 &"right"',2,4)
              call pprint('* mean : mean value of X in datasets',2,4)
              call pprint('* median : median value of X in datasets',2,4)
              call pprint('* min : X of the smallest Y in datasets',2,4)
              call pprint('* max : X of the largest Y in datasets',2,4)
              call pprint('* <X> : user-specified value',2,4)
              call pprint('')
              call pprint('Note that if "centre" is set to a named value &
                 &(i.e., all of the above except "<X>"), the X offset &
                 &dynamically adapts to changes in other variables, e.g., if &
                 &data are loaded/unloaded, if the x/y scales are redefined, &
                 &if data ranges are modified, etc.',2,2)
              call pprint('')
              call pprint('The default value of "centre" is "0".',2,2)
              call pprint('')
            case('nsample')
              call pprint('')
              call pprint('Variable: nsample',0,10)
              call pprint('')
              call pprint('"nsample" is an integer which determines the &
                 &number of Monte Carlo samples to use in the evaluation of &
                 &uncertainties. Note that the uncertainty ddY in an &
                 &uncertainty dY is ddY = dY/sqrt(nsample).',2,2)
              call pprint('')
              call pprint('The default value of "nsample" is 10000, which &
                 &yields an uncertainty in the estimated uncertainty of 1% of &
                 &its value.',2,2)
              call pprint('')
            case default
              call pprint('No help for variable "'//trim(field(3,command))//&
                 &'".')
            end select
          endif
        case default
          call pprint('No help for command "'//trim(field(2,command))//'".')
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


  SUBROUTINE refresh_dataset(dataset,drange)
    !-------------------------------------------------------!
    ! Re-apply transformation and mask to dataset following !
    ! change in choice of transformation or mask.           !
    !-------------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),INTENT(inout) :: dataset
    TYPE(range_type),INTENT(in) :: drange
    TYPE(xy_type),POINTER :: xy,txy,rtxy
    LOGICAL,ALLOCATABLE :: mask(:)
    INTEGER,ALLOCATABLE :: indx(:)
    DOUBLE PRECISION,POINTER :: sortvec(:)
    INTEGER n

    ! Transform.
    xy=>dataset%xy
    call kill_xy(dataset%txy)
    call clone_xy(xy,txy)
    dataset%txy=>txy
    call scale_transform(xy%nxy,dataset%itransfx,xy%x,txy%x,xy%have_dx,xy%dx,&
       &txy%dx)
    call scale_transform(xy%nxy,dataset%itransfy,xy%y,txy%y,xy%have_dy,xy%dy,&
       &txy%dy)

    ! Mask.
    allocate(mask(xy%nxy))
    mask=.true.
    ! Point at sort variable.
    select case(drange%var)
    case('x')
      sortvec=>xy%x
    case('y')
      sortvec=>xy%y
    case('X')
      sortvec=>txy%x
    case('Y')
      sortvec=>txy%y
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
    case('[',']')
      if(drange%size>0)then
        n=min(drange%size,xy%nxy)
        allocate(indx(xy%nxy))
        call isort(xy%nxy,sortvec,indx)
        mask=.false.
        if(trim(drange%op)=='[')then
          mask(indx(1:n))=.true.
        elseif(trim(drange%op)==']')then
          mask(indx(xy%nxy-n+1:xy%nxy))=.true.
        endif
        deallocate(indx)
      endif
    end select
    nullify(sortvec)
    ! Apply weight exponent.
    if(.not.are_equal(dataset%wexp,1.d0))then
      if(are_equal(dataset%wexp,0.d0))then
        txy%w=1.d0
      else
        txy%w=txy%w**dataset%wexp
      endif
    endif

    ! Create masked transformed dataset.
    call kill_xy(dataset%rtxy)
    allocate(rtxy)
    dataset%rtxy=>rtxy
    rtxy%nxy=count(mask)
    rtxy%have_dx=txy%have_dx
    rtxy%have_dy=txy%have_dy
    rtxy%have_w=txy%have_w
    allocate(rtxy%x(rtxy%nxy),rtxy%y(rtxy%nxy),rtxy%dx(rtxy%nxy),&
       &rtxy%dy(rtxy%nxy),rtxy%w(rtxy%nxy))
    if(rtxy%nxy>0)then
      rtxy%x=pack(txy%x,mask)
      rtxy%dx=pack(txy%dx,mask)
      rtxy%y=pack(txy%y,mask)
      rtxy%dy=pack(txy%dy,mask)
      rtxy%w=pack(txy%w,mask)
    endif

  END SUBROUTINE refresh_dataset


  SUBROUTINE refresh_fit(ndataset,dlist,fit)
    !------------------------------------------------!
    ! Refresh value of X0 in fit following change to !
    ! X0_string, range or loaded data.               !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),INTENT(in) :: dlist(ndataset)
    TYPE(fit_form_type),POINTER :: fit
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    INTEGER iset,tot_nxy,ixy,ierr
    DOUBLE PRECISION t1,tx0
    DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:)

    select case(trim(fit%X0_string))

    case('left','right','max','min','centre','mean','median')
      ! Build combined dataset.
      tot_nxy=0
      do iset=1,ndataset
        dataset=>dlist(iset)%dataset
        tot_nxy=tot_nxy+dataset%rtxy%nxy
      enddo ! iset
      allocate(x(tot_nxy),y(tot_nxy))
      tot_nxy=0
      do iset=1,ndataset
        xy=>dlist(iset)%dataset%rtxy
        if(xy%nxy==0)cycle
        x(tot_nxy+1:tot_nxy+xy%nxy)=xy%x(1:xy%nxy)
        y(tot_nxy+1:tot_nxy+xy%nxy)=xy%y(1:xy%nxy)
        tot_nxy=tot_nxy+xy%nxy
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
      t1=parse_dble(fit%x0_string,ierr)
      tx0=t1
    end select

    ! Update X0.
    fit%X0=tx0

  END SUBROUTINE refresh_fit


  SUBROUTINE show_multipoly(ndataset,dlist,drange,fit,mcparams)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_type),INTENT(in) :: fit
    TYPE(mc_params_type),INTENT(in) :: mcparams
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    ! Local variables.
    TYPE(eval_type) deval ! dummy arg.  FIXME - make optional in eval_...
    INTEGER tot_nparam,tot_nxy,iset,i,ierr
    DOUBLE PRECISION chi2,chi2err,rmsy,rmsyerr,a(fit%npoly,ndataset),&
       &da(fit%npoly,ndataset)

    ! Initialize.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)
    tot_nxy=0
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      tot_nxy=tot_nxy+dataset%rtxy%nxy
    enddo ! iset

    ! Perform fit.
    call eval_multifit_monte_carlo(ndataset,dlist,drange,fit,mcparams,&
       &deval,ierr,chi2mean=chi2,chi2err=chi2err,amean=a,aerr=da,&
       &rmsymean=rmsy,rmsyerr=rmsyerr)
    if(ierr/=0)then
      write(6,'(a)')'Could not perform fit.'
      write(6,'()')
      return
    endif

    ! Print table header.
    write(6,'()')
    write(6,'(a)')'Fit parameters:'
    write(6,'()')
    write(6,'(a)',advance='no')'   '
    write(6,'(2x,a5,2x,a5,1x,2(1x,a20))')'Set','i','ki       ','dki       '
    write(6,'(a)',advance='no')'   '
    write(6,'(2x,'//trim(i2s(55))//'("-"))')
    ! Print table.
    do iset=1,ndataset
      do i=1,fit%npoly
        write(6,'(a)',advance='no')'FIT'
        write(6,'(2x,i5,2x,i5,1x,2(1x,es20.12))')iset,i,a(i,iset),da(i,iset)
      enddo ! i
    enddo ! iset
    ! Print table footer.
    write(6,'(a)',advance='no')'   '
    write(6,'(2x,'//trim(i2s(55))//'("-"))')
    write(6,'()')

    ! Report chi-squared.
    write(6,'(a)')'Fit assessment:'
    write(6,'()')
    write(6,'(a)',advance='no')'   '
    write(6,'(2x,a12,1x,2(1x,a20))')'Measure  ','Value       ','Stderr      '
    write(6,'(a)',advance='no')'   '
    write(6,'(2x,'//trim(i2s(55))//'("-"))')
    write(6,'(a)',advance='no')'FIT'
    write(6,'(2x,a12,1x,2(1x,es20.12))')'chi^2    ',chi2,chi2err
    if(tot_nxy>tot_nparam)then
      write(6,'(a)',advance='no')'FIT'
      write(6,'(2x,a12,1x,2(1x,es20.12))')'chi^2/Ndf',&
         &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
    endif
    write(6,'(a)',advance='no')'FIT'
    write(6,'(2x,a12,1x,2(1x,es20.12))')'RMS(Y-f) ',rmsy,rmsyerr
    write(6,'(a)',advance='no')'   '
    write(6,'(2x,'//trim(i2s(55))//'("-"))')
    write(6,'()')

    ! Write out fit in XMGRACE format.
    ! NB, xmgrace has a string length limit of 256 characters, so this may not
    ! be very useful in some cases.
    write(6,'(a)')'Fit in XMGRACE format:'
    do iset=1,ndataset
      write(6,'(a)')'  Set #'//trim(i2s(iset))//': '//&
         &trim(print_poly_num(fit,a(:,iset)))
    enddo ! iset
    write(6,'()')

  END SUBROUTINE show_multipoly


  SUBROUTINE evaluate_fit(ndataset,dlist,fit,mcparams,drange,deval)
    !-------------------------------------------------!
    ! Evaluate value or derivative of fit at provided !
    ! points.                                         !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_type),INTENT(in) :: fit
    TYPE(mc_params_type),INTENT(in) :: mcparams
    TYPE(range_type),INTENT(in) :: drange
    TYPE(eval_type),INTENT(in) :: deval
    ! Local variables.
    INTEGER iset,ix,ierr
    DOUBLE PRECISION,ALLOCATABLE :: fmean(:,:),ferr(:,:)

    ! Evaluate.
    allocate(fmean(deval%n,ndataset),ferr(deval%n,ndataset))
    call eval_multifit_monte_carlo(ndataset,dlist,drange,fit,mcparams,&
       &deval,ierr,fmean,ferr)
    if(ierr/=0)then
      write(6,'(a)')'Could not perform fit.'
      write(6,'()')
      return
    endif

    ! Report.
    write(6,'()')
    write(6,'(4x)',advance='no')
    write(6,'(2x,a4,1x,3(1x,a16))')'set','X       ','f       ','df       '
    write(6,'(4x)',advance='no')
    write(6,'(2x,56("-"))')
    do iset=1,ndataset
      do ix=1,deval%n
        write(6,'(a4)',advance='no')'EVAL'
        write(6,'(2x,i4,1x,3(1x,es20.12))')iset,deval%x(ix),&
           &fmean(ix,iset),ferr(ix,iset)
      enddo ! ix
      select case(trim(deval%what))
      case('shared','sum')
        exit
      end select
    enddo ! iset
    write(6,'(4x)',advance='no')
    write(6,'(2x,56("-"))')
    write(6,'()')
    deallocate(fmean,ferr)

  END SUBROUTINE evaluate_fit


  SUBROUTINE intersect_fit(ndataset,dlist,fit,mcparams,drange,x1,x2)
    !-----------------------------------------------------!
    ! Find the (average) intersection of all datasets and !
    ! report to stdout.                                   !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_type),INTENT(in) :: fit
    TYPE(mc_params_type),INTENT(in) :: mcparams
    TYPE(range_type),INTENT(in) :: drange
    DOUBLE PRECISION,INTENT(in) :: x1,x2
    ! Monte Carlo sample storage.
    DOUBLE PRECISION,ALLOCATABLE :: x0_array(:,:),y0_array(:,:),w_vector(:)
    LOGICAL,ALLOCATABLE :: valid(:)
    ! Parameter vector.
    DOUBLE PRECISION a(fit%npoly,ndataset)
    ! Pointers.
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy,xy_orig
    ! Local variables.
    INTEGER iset,jset,iintersect,nintersect,irandom,nsample,ierr
    DOUBLE PRECISION x0,y0,dx0,dy0,sum_weight,var,chi2

    ! Make copy of datasets.
    call clone_dlist(dlist,tmp_dlist)

    ! Allocate storage for location of intersection.
    nintersect=(ndataset*(ndataset-1))/2
    allocate(x0_array(mcparams%nsample,nintersect),&
       &y0_array(mcparams%nsample,nintersect),w_vector(mcparams%nsample),&
       &valid(nintersect))
    valid=.true.

    ! Compute total dataset weight.
    sum_weight=0.d0
    do iset=1,ndataset
      sum_weight=sum_weight+tmp_dlist(iset)%dataset%weight
    enddo ! iset

    ! Initialize.
    nsample=mcparams%nsample

    ! Loop over random points.
    irandom=0
    do irandom=1,nsample
      do iset=1,ndataset
        dataset=>tmp_dlist(iset)%dataset
        xy=>dataset%xy
        xy_orig=>dlist(iset)%dataset%xy
        if(xy%have_dx)xy%x=xy_orig%x+gaussian_random_number(xy_orig%dx)
        if(xy%have_dy)xy%y=xy_orig%y+gaussian_random_number(xy_orig%dy)
        call refresh_dataset(dataset,drange)
      enddo ! iset
      call perform_multifit(ndataset,tmp_dlist,fit,chi2,a,ierr)
      if(ierr/=0)then
        call kill_dlist(tmp_dlist)
        return
      endif
      w_vector(irandom)=1.d0
      ! Loop over pairs of datasets.
      iintersect=0
      do iset=1,ndataset-1
        do jset=iset+1,ndataset
          iintersect=iintersect+1
          if(.not.valid(iintersect))cycle
          ! Find intersection between sets ISET and JSET.
          call intersect(ndataset,fit,a,iset,jset,x1,x2,x0,y0,ierr)
          if(ierr/=0)valid(iintersect)=.false.
          x0_array(irandom,iintersect)=x0
          y0_array(irandom,iintersect)=y0
        enddo ! jset
      enddo ! iset
    enddo ! irandom

    ! Compute and report intersection(s).
    write(6,'(a)')'Intersections:'
    write(6,'(4x)',advance='no')
    write(6,'(1x,a7,4(1x,a20))')'Sets ','X0          ','DX0         ',&
       &'Y0          ','DY0         '
    iintersect=0
    do iset=1,ndataset-1
      do jset=iset+1,ndataset
        iintersect=iintersect+1
        if(valid(iintersect))then
          call characterize_dist(nsample,x0_array(1,iintersect),w_vector,&
             &mean=x0,var=var)
          dx0=sqrt(var)
          call characterize_dist(nsample,y0_array(1,iintersect),w_vector,&
             &mean=y0,var=var)
          dy0=sqrt(var)
          write(6,'(a4)',advance='no')'INTR'
          write(6,'(2(1x,i3),4(1x,es20.12))')iset,jset,x0,dx0,y0,dy0
        else
          write(6,'(a4)',advance='no')'INTR'
          write(6,'(2(1x,i3),a84)')iset,jset,&
             &'NO RELIABLE INTERSECTION                            '
        endif
      enddo ! jset
    enddo ! iset

    ! Clean up.
    call kill_dlist(tmp_dlist)

  END SUBROUTINE intersect_fit


  SUBROUTINE intersect(ndataset,fit,a,iset,jset,x1,x2,x0,y0,ierr)
    !---------------------------------------------------!
    ! Find intersection between datasets ISET and JSET. !
    !---------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(fit_form_type),INTENT(in) :: fit
    DOUBLE PRECISION,INTENT(in) :: a(fit%npoly,ndataset)
    INTEGER,INTENT(in) :: iset,jset
    DOUBLE PRECISION,INTENT(in) :: x1,x2
    DOUBLE PRECISION,INTENT(inout) :: x0,y0
    INTEGER,INTENT(inout) :: ierr
    LOGICAL lplus,rplus
    DOUBLE PRECISION xl,xr,yl,yr,a_diff(fit%npoly)

    ! Initialize output variables.
    ierr=1
    x0=0.d0
    y0=0.d0

    ! Get coefficients of difference between the two relevant sets.
    a_diff(1:fit%npoly)=a(1:fit%npoly,jset)-a(1:fit%npoly,iset)

    ! Initialize bisection.
    xl=x1
    yl=eval_poly(fit%npoly,fit%pow,a_diff,xl-fit%x0)
    xr=x2
    yr=eval_poly(fit%npoly,fit%pow,a_diff,xr-fit%x0)
    if(abs(yl)<=0.d0.or.abs(yr)<=0.d0)return
    lplus=yl>0.d0
    rplus=yr>0.d0
    if(lplus.eqv.rplus)return

    ! Loop over bisection iterations.
    ierr=0
    do
      x0=0.5d0*(xl+xr)
      y0=eval_poly(fit%npoly,fit%pow,a_diff,x0-fit%x0)
      if(y0>0.d0.eqv.lplus)then
        ! Replace left bracket.
        xl=x0
      else
        ! Replace right bracket.
        xr=x0
      endif
      if(abs(xl-xr)<1.d-10)exit
    enddo

    ! Evaluate final intersection point.
    x0=0.5d0*(xl+xr)
    y0=eval_poly(fit%npoly,fit%pow,a(:,iset),x0-fit%x0)

  END SUBROUTINE intersect


  SUBROUTINE plot_multipoly(ndataset,dlist,drange,fit,deval,mcparams,fname)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_type),INTENT(in) :: fit
    TYPE(eval_type),INTENT(inout) :: deval
    TYPE(mc_params_type),INTENT(in) :: mcparams
    CHARACTER(*),INTENT(in) :: fname
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    INTEGER iset,i,ierr,op_npoly
    DOUBLE PRECISION op_a(fit%npoly),op_pow(fit%npoly),tx0,tx1,t1
    DOUBLE PRECISION,ALLOCATABLE :: fmean(:,:),ferr(:,:),a(:,:)
    INTEGER,PARAMETER :: io=10

    ! Open plot file.
    open(unit=io,file=fname,status='replace',iostat=ierr)
    if(ierr/=0)then
      write(6,'(a)')'Problem opening file.'
      write(6,'()')
      return
    endif

    ! Define plot range if not provided.
    if(.not.associated(deval%x))then
      deval%var='X'
      deval%n=101
      allocate(deval%x(deval%n))
      ! Get data range.
      dataset=>dlist(1)%dataset
      tx0=minval(dataset%txy%x)
      tx1=maxval(dataset%txy%x)
      do iset=2,ndataset
        dataset=>dlist(iset)%dataset
        tx0=min(tx0,minval(dataset%txy%x))
        tx1=max(tx1,maxval(dataset%txy%x))
      enddo ! iset
      ! Extend range past far end and ensure we get to zero on near end.
      t1=tx1-tx0
      if(tx0>=0.d0)then
        tx0=0.d0
      else
        tx0=tx0-0.25d0*t1
      endif
      if(tx1<=0.d0)then
        tx1=0.d0
      else
        tx1=tx1+0.25d0*t1
      endif
      ! Generate grid.
      do i=1,deval%n
        deval%x(i)=tx0+(tx1-tx0)*dble(i-1)/dble(deval%n-1)
      enddo ! i
    endif

    ! Perform fit and plot
    allocate(fmean(deval%n,ndataset),ferr(deval%n,ndataset),&
       &a(fit%npoly,ndataset))
    if(trim(deval%what)/='shared')then
      call eval_multifit_monte_carlo(ndataset,dlist,drange,fit,mcparams,&
         &deval,ierr,fmean=fmean,ferr=ferr,amean=a)
      if(ierr/=0)then
        write(6,'(a)')'Could not perform fit.'
        write(6,'()')
        close(io,status='delete')
        return
      endif
      if(deval%nderiv==0)then
        ! Plot data.
        do iset=1,ndataset
          xy=>dlist(iset)%dataset%txy
          if(xy%have_dx.and.xy%have_dy)then
            !write(io,'(a)')'@type xydxdy'
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i),xy%dx(i),xy%dy(i)
            enddo ! i
          elseif(xy%have_dx)then
            !write(io,'(a)')'@type xydx'
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i),xy%dx(i)
            enddo ! i
          elseif(xy%have_dy)then
            !write(io,'(a)')'@type xydy'
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i),xy%dy(i)
            enddo ! i
          else
            !write(io,'(a)')'@type xy'
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i)
            enddo ! i
          endif
          !write(io,'(a)')'&'
          write(io,'()')
          write(io,'()')
        enddo ! iset
      endif
      ! Plot fit functions.
      do iset=1,ndataset
        !write(io,'(a)')'@type xydy'
        do i=1,deval%n
          write(io,*)deval%x(i),fmean(i,iset),ferr(i,iset)
        enddo ! i
        if(iset<ndataset)then
          !write(io,'(a)')'&'
          write(io,'()')
          write(io,'()')
        endif
      enddo ! iset
    else
      call eval_multifit_monte_carlo(ndataset,dlist,drange,fit,mcparams,&
         &deval,ierr,fmean=fmean,ferr=ferr,amean=a)
      if(ierr/=0)then
        write(6,'(a)')'Could not perform fit.'
        write(6,'()')
        close(io,status='delete')
        return
      endif
      ! Plot data minus independent bit.
      if(deval%nderiv==0)then
        op_npoly=fit%npoly
        op_pow=fit%pow
        do iset=1,ndataset
          xy=>dlist(iset)%dataset%txy
          op_a(1:op_npoly)=a(1:fit%npoly,iset)
          where(fit%share)op_a=0.d0
          if(xy%have_dx.and.xy%have_dy)then
            !write(io,'(a)')'@type xydxdy'
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i)),xy%dx(i),&
                 &xy%dy(i)
            enddo ! i
          elseif(xy%have_dx)then
            !write(io,'(a)')'@type xydx'
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i)),xy%dx(i)
            enddo ! i
          elseif(xy%have_dy)then
            !write(io,'(a)')'@type xydy'
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i)),xy%dy(i)
            enddo ! i
          else
            !write(io,'(a)')'@type xy'
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i))
            enddo ! i
          endif
          !write(io,'(a)')'&'
          write(io,'()')
          write(io,'()')
        enddo ! iset
      endif
      ! Plot shared part of fit functions.
      !write(io,'(a)')'@type xydy'
      do i=1,deval%n
        write(io,*)deval%x(i),fmean(i,1),ferr(i,1)
      enddo ! i
    endif
    nullify(xy)

    ! Close file.
    close(io)

    ! Report.
    write(6,'(a)')'Plot saved to "'//trim(fname)//'".'
    write(6,'()')

  END SUBROUTINE plot_multipoly


  SUBROUTINE assess_fit(ndataset,dlist,drange,fit,mcparams,deval,eval_iset)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! expansion order.                               !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_type),INTENT(in) :: fit
    TYPE(mc_params_type),INTENT(in) :: mcparams
    TYPE(eval_type),INTENT(in) :: deval
    INTEGER,INTENT(in) :: eval_iset
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    ! Local variables.
    TYPE(fit_form_type),POINTER :: tfit
    INTEGER tot_nxy,tot_nparam,iset,npoly,ix,prev_npoly,ierr
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION chi2_all(fit%npoly),chi2err_all(fit%npoly),&
       &fmean_all(deval%n,ndataset,fit%npoly),&
       &ferr_all(deval%n,ndataset,fit%npoly)

    ! Prepare combined dataset arrays.
    tot_nxy=0
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      tot_nxy=tot_nxy+dataset%rtxy%nxy
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
      call refresh_fit(ndataset,dlist,tfit)
      ! Adjust counters.
      tot_nparam=count(tfit%share)+ndataset*count(.not.tfit%share)
      if(tot_nparam<tot_nxy)then
        ! Perform fit and report.
        call eval_multifit_monte_carlo(ndataset,dlist,drange,tfit,mcparams,&
           &deval,ierr,fmean,ferr,chi2,chi2err)
        if(ierr/=0)cycle
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


  SUBROUTINE assess_range(ndataset,dlist,drange,fit,mcparams,deval,eval_iset)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! number of data points.                         !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_type),POINTER :: fit
    TYPE(range_type),INTENT(inout) :: drange
    TYPE(mc_params_type),INTENT(in) :: mcparams
    TYPE(eval_type),INTENT(in) :: deval
    INTEGER,INTENT(in) :: eval_iset
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Pointer-resizing pointers.
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    ! Local variables.
    INTEGER ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,ix,prev_igrid,ierr
    INTEGER,ALLOCATABLE :: indx(:),tngrid(:),tot_nxy_grid(:)
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION,ALLOCATABLE :: chi2_all(:),chi2err_all(:),&
       &fmean_all(:,:,:),ferr_all(:,:,:)
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:)

    ! Initialize.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)

    ! Compute grid.
    tot_nxy=0
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      tot_nxy=tot_nxy+dataset%txy%nxy
    enddo ! iset
    allocate(tot_nxy_grid(tot_nxy))
    select case(drange%op(1:1))
    case('<','>')
      allocate(txall(tot_nxy),indx(tot_nxy),txgrid(tot_nxy),tngrid(tot_nxy))
      tngrid=0
      tot_nxy=0
      do iset=1,ndataset
        xy=>dlist(iset)%dataset%txy
        txall(tot_nxy+1:tot_nxy+xy%nxy)=xy%x(1:xy%nxy)
        tot_nxy=tot_nxy+xy%nxy
      enddo ! iset
      call isort(tot_nxy,txall,indx)
      if(drange%op(1:1)=='>')then
        do ixy=1,tot_nxy
          if(ixy>=tot_nxy-ixy+1)exit
          call iswap1(indx(ixy),indx(tot_nxy-ixy+1))
        enddo ! ixy
      endif
      ngrid=1
      txgrid(1)=txall(indx(1))
      do ixy=2,tot_nxy
        if(are_equal(txall(indx(ixy)),txgrid(ngrid)))cycle
        ngrid=ngrid+1
        txgrid(ngrid)=txall(indx(ixy))
      enddo ! ixy
      deallocate(txall,indx)
    case('[',']')
      ngrid=0
      do iset=1,ndataset
        xy=>dlist(iset)%dataset%txy
        if(iset==1.or.xy%nxy<ngrid)ngrid=xy%nxy
      enddo ! iset
      allocate(txgrid(ngrid),tngrid(ngrid))
      txgrid=0.d0
      do igrid=1,ngrid
        tngrid(igrid)=igrid
      enddo ! igrid
    end select

    ! Allocate arrays to store test results and initialize to "invalid".
    allocate(chi2_all(tot_nxy),chi2err_all(tot_nxy),&
       &fmean_all(deval%n,ndataset,tot_nxy),&
       &ferr_all(deval%n,ndataset,tot_nxy))
    chi2_all=-1.d0
    chi2err_all=-1.d0
    fmean_all=0.d0
    ferr_all=-1.d0

    ! Make copy of datasets.
    call clone_dlist(dlist,tmp_dlist)

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
      drange%size=tngrid(igrid)
      tot_nxy=0
      do iset=1,ndataset
        call refresh_dataset(tmp_dlist(iset)%dataset,drange)
        tot_nxy=tot_nxy+tmp_dlist(iset)%dataset%rtxy%nxy
      enddo ! iset
      tot_nxy_grid(igrid)=tot_nxy
      if(tot_nparam<tot_nxy)then
        ! Perform fit and report.
        call eval_multifit_monte_carlo(ndataset,tmp_dlist,drange,fit,mcparams,&
           &deval,ierr,fmean,ferr,chi2,chi2err)
        if(ierr/=0)cycle
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
      write(6,'(a)')'Suggested grid point: '//&
         &trim(i2s(tot_nxy_grid(prev_igrid)))
    else
      write(6,'(a)')'Could not find optimal range by criteria.'
    endif
    write(6,'()')

    ! Clean up.
    call kill_dlist(tmp_dlist)
    deallocate(chi2_all,chi2err_all,fmean_all,ferr_all)

  END SUBROUTINE assess_range


  SUBROUTINE assess_fit_range(ndataset,dlist,drange,fit,mcparams,deval,&
     &eval_iset)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! expansion order and number of data points.     !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(inout) :: drange
    TYPE(fit_form_type),POINTER :: fit
    TYPE(mc_params_type),INTENT(in) :: mcparams
    TYPE(eval_type),INTENT(in) :: deval
    INTEGER,INTENT(in) :: eval_iset
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    TYPE(fit_form_type),POINTER :: tfit
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    INTEGER ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,npoly,ix,prev_igrid,&
       &prev_npoly,prev_prev_npoly,ierr
    INTEGER,ALLOCATABLE :: indx(:),tngrid(:),tot_nxy_grid(:)
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION,ALLOCATABLE :: chi2_all(:,:),chi2err_all(:,:),&
       &fmean_all(:,:,:,:),ferr_all(:,:,:,:)
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:)

    ! Compute grid.
    tot_nxy=0
    do iset=1,ndataset
      xy=>dlist(iset)%dataset%txy
      tot_nxy=tot_nxy+xy%nxy
    enddo ! iset
    allocate(tot_nxy_grid(tot_nxy))
    select case(drange%op(1:1))
    case('<','>')
      allocate(txall(tot_nxy),indx(tot_nxy),txgrid(tot_nxy),tngrid(tot_nxy))
      tngrid=0
      tot_nxy=0
      do iset=1,ndataset
        xy=>dlist(iset)%dataset%txy
        txall(tot_nxy+1:tot_nxy+xy%nxy)=xy%x(1:xy%nxy)
        tot_nxy=tot_nxy+xy%nxy
      enddo ! iset
      call isort(tot_nxy,txall,indx)
      if(drange%op(1:1)=='>')then
        do ixy=1,tot_nxy
          if(ixy>=tot_nxy-ixy+1)exit
          call iswap1(indx(ixy),indx(tot_nxy-ixy+1))
        enddo ! ixy
      endif
      ngrid=1
      txgrid(1)=txall(indx(1))
      do ixy=2,tot_nxy
        if(are_equal(txall(indx(ixy)),txgrid(ngrid)))cycle
        ngrid=ngrid+1
        txgrid(ngrid)=txall(indx(ixy))
      enddo ! ixy
      deallocate(txall,indx)
    case('[',']')
      ngrid=0
      do iset=1,ndataset
        xy=>dlist(iset)%dataset%txy
        if(iset==1.or.xy%nxy<ngrid)ngrid=xy%nxy
      enddo ! iset
      allocate(txgrid(ngrid),tngrid(ngrid))
      txgrid=0.d0
      do igrid=1,ngrid
        tngrid(igrid)=igrid
      enddo ! igrid
    end select

    ! Allocate arrays to store test results and initialize to "invalid".
    allocate(chi2_all(fit%npoly,tot_nxy),chi2err_all(fit%npoly,tot_nxy),&
       &fmean_all(deval%n,ndataset,fit%npoly,tot_nxy),&
       &ferr_all(deval%n,ndataset,fit%npoly,tot_nxy))
    chi2_all=-1.d0
    chi2err_all=-1.d0
    fmean_all=0.d0
    ferr_all=-1.d0

    ! Make copy of datasets.
    call clone_dlist(dlist,tmp_dlist)

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
      drange%thres=txgrid(igrid)
      drange%size=tngrid(igrid)
      tot_nxy=0
      do iset=1,ndataset
        dataset=>dlist(iset)%dataset
        call refresh_dataset(dataset,drange)
        tot_nxy=tot_nxy+dataset%rtxy%nxy
      enddo ! iset
      tot_nxy_grid(igrid)=tot_nxy
      ! Loop over expansion orders.
      do npoly=1,fit%npoly
        ! Allocate work arrays.
        allocate(tfit)
        tfit%npoly=npoly
        allocate(tfit%pow(npoly),tfit%share(npoly))
        tfit%pow=fit%pow(1:npoly)
        tfit%share=fit%share(1:npoly)
        tfit%X0_string=fit%X0_string
        call refresh_fit(ndataset,tmp_dlist,tfit)
        ! Adjust counters.
        tot_nparam=count(tfit%share)+ndataset*count(.not.tfit%share)
        if(tot_nparam<tot_nxy)then
          ! Perform fit and report.
          call eval_multifit_monte_carlo(ndataset,tmp_dlist,drange,tfit,&
             &mcparams,deval,ierr,fmean,ferr,chi2,chi2err)
          if(ierr/=0)cycle
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

    ! Print second table header.
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
        ! Store npoly.
        prev_npoly=npoly
      enddo ! npoly
      ! Skip points for which test did not converge.
      if(npoly>fit%npoly)cycle
      ! Write table entry.
      write(6,'(2x,i5,2x,i5,1x,2(1x,es12.4))',advance='no')&
         &tot_nxy_grid(igrid),prev_npoly-1,&
         &chi2_all(prev_npoly,igrid),chi2err_all(prev_npoly,igrid)
      do ix=1,deval%n
        do iset=1,ndataset
          if(eval_iset==0.or.eval_iset==iset)write(6,'(2(1x,es12.4))',&
             &advance='no')fmean_all(ix,iset,prev_npoly,igrid),&
             &ferr_all(ix,iset,prev_npoly,igrid)
        enddo ! iset
      enddo ! ix
      write(6,'()')
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

    ! Print table footer.
    if(eval_iset==0)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n*ndataset))//'("-"))')
    elseif(eval_iset>0.and.eval_iset<ndataset)then
      write(6,'(2x,'//trim(i2s(39+26*deval%n))//'("-"))')
    else
      write(6,'(2x,'//trim(i2s(39))//'("-"))')
    endif
    write(6,'()')

    ! Report suggestion.
    if(igrid>=1)then
      write(6,'(a)')'Suggested fit: '//trim(i2s(prev_prev_npoly-1))
      write(6,'(a)')'Suggested grid point: '//&
         &trim(i2s(tot_nxy_grid(prev_igrid)))
    else
      write(6,'(a)')'Could not find optimal range by criteria.'
    endif
    write(6,'()')

    ! Clean up.
    call kill_dlist(tmp_dlist)
    deallocate(chi2_all,chi2err_all,fmean_all,ferr_all)

  END SUBROUTINE assess_fit_range


  SUBROUTINE report_statistics(ndataset,dlist)
    !-----------------------------------------------------!
    ! Show min/centre/mean/median/max of X, Y, dX and dY. !
    !-----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),INTENT(in) :: dlist(ndataset)
    ! Quick-access pointers.
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
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
      ! Transformed X.
      xy=>dlist(iset)%dataset%txy
      vmin=minval(xy%x)
      vmax=maxval(xy%x)
      vcentre=.5d0*(vmin+vmax)
      vmean=sum(xy%x)/dble(xy%nxy)
      vmedian=median(xy%nxy,xy%x)
      write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'X',iset,vmin,vmean,vmedian,&
         &vcentre,vmax
    enddo ! iset
    do iset=1,ndataset
      ! Transformed dX.
      xy=>dlist(iset)%dataset%txy
      if(xy%have_dx)then
        vmin=minval(xy%dx)
        vmax=maxval(xy%dx)
        vcentre=.5d0*(vmin+vmax)
        vmean=sum(xy%dx)/dble(xy%nxy)
        vmedian=median(xy%nxy,xy%dx)
        write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'dX',iset,vmin,vmean,vmedian,&
           &vcentre,vmax
      endif
    enddo ! iset
    do iset=1,ndataset
      ! Transformed Y.
      xy=>dlist(iset)%dataset%txy
      vmin=minval(xy%y)
      vmax=maxval(xy%y)
      vcentre=.5d0*(vmin+vmax)
      vmean=sum(xy%y)/dble(xy%nxy)
      vmedian=median(xy%nxy,xy%y)
      write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'Y',iset,vmin,vmean,vmedian,&
         &vcentre,vmax
    enddo ! iset
    do iset=1,ndataset
      ! Transformed dY.
      xy=>dlist(iset)%dataset%txy
      if(xy%have_dy)then
        vmin=minval(xy%dy)
        vmax=maxval(xy%dy)
        vcentre=.5d0*(vmin+vmax)
        vmean=sum(xy%dy)/dble(xy%nxy)
        vmedian=median(xy%nxy,xy%dy)
        write(6,'(2x,a3,1x,i4,5(1x,es12.4))')'dY',iset,vmin,vmean,vmedian,&
           &vcentre,vmax
      endif
    enddo ! iset
    write(6,'(2x,74("-"))')
    write(6,'()')

  END SUBROUTINE report_statistics


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


  SUBROUTINE perform_multifit(ndataset,dlist,fit,chi2,a,ierr)
    !----------------------------------------------------!
    ! Perform least-squares fit of sets of (weighted) xy !
    ! data to the polynomial of exponents pow(1:npoly),  !
    ! with equal/independent coefficients for each set   !
    ! depending on the value of share(1:npoly).          !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),INTENT(in) :: dlist(ndataset)
    TYPE(fit_form_type),INTENT(in) :: fit
    DOUBLE PRECISION,INTENT(inout) :: chi2,a(fit%npoly,ndataset)
    INTEGER,INTENT(inout) :: ierr
    ! Quick-access pointers.
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    DOUBLE PRECISION,ALLOCATABLE :: x(:,:),y(:,:),w(:,:)
    INTEGER tot_nparam,max_nxy,tot_nxy,ieq,ip,jp,iset,jset,&
       &ixy,ipoly,jpoly,lwork,i,j
    DOUBLE PRECISION e_fit,set_weight
    DOUBLE PRECISION,ALLOCATABLE :: M(:,:),Minv(:,:),c(:),work(:)
    INTEGER,ALLOCATABLE :: ipiv(:)

    ! Initialize.
    ierr=0

    ! Extract fit properties.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)
    max_nxy=0
    tot_nxy=0
    do iset=1,ndataset
      xy=>dlist(iset)%dataset%rtxy
      max_nxy=max(max_nxy,xy%nxy)
      tot_nxy=tot_nxy+xy%nxy
    enddo ! iset

    ! Explicitly prevent fitting with insufficient data.
    if(tot_nxy<tot_nparam)then
      ierr=1
      return
    endif

    ! Prepare combined dataset arrays.
    allocate(x(max_nxy,ndataset),y(max_nxy,ndataset),w(max_nxy,ndataset))
    x=1.d0
    y=0.d0
    w=0.d0
    do iset=1,ndataset
      set_weight=dlist(iset)%dataset%weight
      xy=>dlist(iset)%dataset%rtxy
      if(xy%nxy==0)cycle
      x(1:xy%nxy,iset)=xy%x(1:xy%nxy)-fit%X0
      y(1:xy%nxy,iset)=xy%y(1:xy%nxy)
      w(1:xy%nxy,iset)=set_weight*xy%w(1:xy%nxy)
    enddo ! iset
    nullify(xy)

    ! Construct c vector and M matrix.
    allocate(M(tot_nparam,tot_nparam),Minv(tot_nparam,tot_nparam),&
       &ipiv(tot_nparam),c(tot_nparam))

    ! Initialize equation counter.
    ieq=0
    ! Loop over shared parameters.
    do ipoly=1,fit%npoly
      if(.not.fit%share(ipoly))cycle
      ! There is one equation associated with this parameter.
      ieq=ieq+1
      ! Right-hand side.
      c(ieq)=0.d0
      do iset=1,ndataset
        c(ieq)=c(ieq)+sum(w(:,iset)*y(:,iset)*x(:,iset)**fit%pow(ipoly))
      enddo ! iset
      ! Coefficients of shared parameters.
      jp=0
      do jpoly=1,fit%npoly
        if(.not.fit%share(jpoly))cycle
        jp=jp+1
        M(jp,ieq)=0.d0
        do iset=1,ndataset
          M(jp,ieq)=M(jp,ieq)+sum(w(:,iset)*x(:,iset)**&
             &(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! iset
      enddo ! jpoly
      ! Coefficients of independent parameters.
      do jpoly=1,fit%npoly
        if(fit%share(jpoly))cycle
        do iset=1,ndataset
          jp=jp+1
          M(jp,ieq)=sum(w(:,iset)*x(:,iset)**(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! iset
      enddo ! jpoly
    enddo ! ipoly
    ! Loop over independent parameters.
    do ipoly=1,fit%npoly
      if(fit%share(ipoly))cycle
      ! There are ndataset equations associated with this parameter -- each
      ! associated with each of the ndataset instances of the parameter.
      do iset=1,ndataset
        ieq=ieq+1
        ! Right-hand side.
        c(ieq)=sum(w(:,iset)*y(:,iset)*x(:,iset)**fit%pow(ipoly))
        ! Coefficients of shared parameters.
        jp=0
        do jpoly=1,fit%npoly
          if(.not.fit%share(jpoly))cycle
          jp=jp+1
          M(jp,ieq)=M(jp,ieq)+sum(w(:,iset)*x(:,iset)**&
             &(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! jpoly
        ! Coefficients of independent parameters.
        do jpoly=1,fit%npoly
          if(fit%share(jpoly))cycle
          do jset=1,ndataset
            jp=jp+1
            if(jset==iset)then
              M(jp,ieq)=sum(w(:,iset)*x(:,iset)**&
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
    if(ierr/=0)return
    lwork=nint(work(1))
    deallocate(work)
    allocate(work(lwork),stat=ierr)
    if(ierr/=0)return
    call dsytrf('L',tot_nparam,Minv,tot_nparam,ipiv,work,lwork,ierr)
    if(ierr/=0)return
    deallocate(work)
    allocate(work(tot_nparam),stat=ierr)
    if(ierr/=0)return
    call dsytri('L',tot_nparam,Minv,tot_nparam,ipiv,work,ierr)
    if(ierr/=0)return
    deallocate(work)
    ! Complete Minv.
    do i=1,tot_nparam-1
      do j=i+1,tot_nparam
        Minv(i,j)=Minv(j,i)
      enddo ! j
    enddo ! i

    ! FIXME - for some reason the following does not work in the presence
    ! of shared parameters.  I would have expected the opposite, if anything.
    !! Invert M.
    !Minv=M
    !call dgetrf(tot_nparam,tot_nparam,Minv,tot_nparam,ipiv,ierr)
    !if(ierr/=0)return
    !allocate(work(1))
    !lwork=-1
    !call dgetri(tot_nparam,Minv,tot_nparam,ipiv,work,lwork,ierr)
    !if(ierr/=0)return
    !lwork=nint(work(1))
    !deallocate(work)
    !allocate(work(lwork),stat=ierr)
    !if(ierr/=0)return
    !call dgetri(tot_nparam,Minv,tot_nparam,ipiv,work,lwork,ierr)
    !if(ierr/=0)return
    !deallocate(work)

    ! Evaluate fit coefficients.
    c=matmul(Minv,c)

    ! Put parameters in output arrays.
    ip=0
    ! Loop over shared parameters.
    do ipoly=1,fit%npoly
      if(.not.fit%share(ipoly))cycle
      ip=ip+1
      a(ipoly,1:ndataset)=c(ip)
    enddo ! ipoly
    ! Loop over independent parameters.
    do ipoly=1,fit%npoly
      if(fit%share(ipoly))cycle
      do iset=1,ndataset
        ip=ip+1
        a(ipoly,iset)=c(ip)
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
        chi2=chi2+(y(ixy,iset)-e_fit)**2*w(ixy,iset)
      enddo ! ixy
    enddo ! iset

  END SUBROUTINE perform_multifit


  SUBROUTINE eval_multifit_monte_carlo(ndataset,dlist,drange,fit,mcparams,&
     &deval,ierr,fmean,ferr,chi2mean,chi2err,rmsymean,rmsyerr,amean,aerr,&
     &fmean_1s,ferr_1s,fmean_2s,ferr_2s,fmed,fskew,fkurt)
    !------------------------------------------------------!
    ! Perform Monte Carlo sampling of data space to obtain !
    ! fit values or derivatives at specified points with   !
    ! error bars.                                          !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_type),INTENT(in) :: fit
    TYPE(mc_params_type),INTENT(in) :: mcparams
    TYPE(eval_type),INTENT(in) :: deval
    INTEGER,INTENT(inout) :: ierr
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: &
       &fmean(deval%n,ndataset),ferr(deval%n,ndataset),&
       &chi2mean,chi2err,rmsymean,rmsyerr,&
       &amean(fit%npoly,ndataset),aerr(fit%npoly,ndataset),&
       &fmean_1s(deval%n,ndataset),ferr_1s(deval%n,ndataset),&
       &fmean_2s(deval%n,ndataset),ferr_2s(deval%n,ndataset),&
       &fmed(deval%n,ndataset),fskew(deval%n,ndataset),&
       &fkurt(deval%n,ndataset)
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy,xy_orig
    ! Pointers.
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    ! Polynomials resulting from operations on fitting polynomial.
    INTEGER op_npoly
    DOUBLE PRECISION op_pow(fit%npoly),op_a(fit%npoly)
    ! Random sampling arrays.
    DOUBLE PRECISION w_vector(mcparams%nsample)
    DOUBLE PRECISION,ALLOCATABLE :: f_array(:,:,:),a_array(:,:,:),&
       &chi2_array(:),rmsy_array(:)
    ! Distribution analysis.
    DOUBLE PRECISION var,skew,kurt,f2s_lo,f1s_lo,f1s_hi,f2s_hi
    ! Misc.
    LOGICAL need_f,need_a,need_chi2,need_rmsy
    INTEGER nsample,ipoly,ideriv,irandom,ix,iset
    DOUBLE PRECISION chi2,t1,a(fit%npoly,ndataset),f0,sum_weight

    ! Initialize.
    ierr=0

    ! Make copy of datasets.
    call clone_dlist(dlist,tmp_dlist)

    ! Figure out what we need and allocate arrays.
    need_f=present(fmean).or.present(ferr).or.present(fmean_1s).or.&
       &present(ferr_1s).or.present(fmean_2s).or.present(ferr_2s).or.&
       &present(fmed).or.present(fskew).or.present(fkurt)
    need_f=need_f.and.deval%n>0
    need_chi2=present(chi2mean).or.present(chi2err)
    need_rmsy=present(rmsymean).or.present(rmsyerr)
    need_a=present(amean).or.present(aerr)
    ! Determine what to evaluate (function per set, shared portion of
    ! function, "unshared" portion of function per set, sum of functions).
    ! Allocate f_array.
    if(need_f)then
      select case(trim(deval%what))
      case('shared','sum')
        allocate(f_array(mcparams%nsample,deval%n,1))
      case default
        allocate(f_array(mcparams%nsample,deval%n,ndataset))
      end select
      f_array=0.d0
    endif
    if(need_a)allocate(a_array(mcparams%nsample,fit%npoly,ndataset))
    if(need_chi2)allocate(chi2_array(mcparams%nsample))
    if(need_rmsy)allocate(rmsy_array(mcparams%nsample))

    ! Compute total dataset weight.
    sum_weight=0.d0
    do iset=1,ndataset
      sum_weight=sum_weight+tmp_dlist(iset)%dataset%weight
    enddo ! iset

    ! Initialize.
    f0=0.d0
    nsample=mcparams%nsample

    ! Loop over random points.
    do irandom=1,nsample
      do iset=1,ndataset
        dataset=>tmp_dlist(iset)%dataset
        xy=>dataset%xy
        xy_orig=>dlist(iset)%dataset%xy
        if(xy%have_dx)xy%x=xy_orig%x+gaussian_random_number(xy_orig%dx)
        if(xy%have_dy)xy%y=xy_orig%y+gaussian_random_number(xy_orig%dy)
        call refresh_dataset(dataset,drange)
      enddo ! iset
      call perform_multifit(ndataset,tmp_dlist,fit,chi2,a,ierr)
      if(ierr/=0)then
        call kill_dlist(tmp_dlist)
        return
      endif
      w_vector(irandom)=1.d0
      if(need_chi2)chi2_array(irandom)=chi2
      if(need_rmsy)rmsy_array(irandom)=sqrt(chi2/sum_weight)
      if(need_a)a_array(irandom,1:fit%npoly,1:ndataset)=&
         &a(1:fit%npoly,1:ndataset)
      ! Evaluate requested function.
      if(need_f)then
        select case(trim(deval%what))
        case('shared')
          ! Evaluate shared component of requested function.
          op_npoly=fit%npoly
          op_pow=fit%pow
          op_a=a(:,1)
          where(.not.fit%share)op_a=0.d0
          do ideriv=1,deval%nderiv
            call deriv_poly(op_npoly,op_pow,op_a)
          enddo ! ideriv
          if(deval%rel)f0=eval_poly(op_npoly,op_pow,op_a,-fit%X0)
          do ix=1,deval%n
            t1=eval_poly(op_npoly,op_pow,op_a,deval%x(ix)-fit%X0)
            f_array(irandom,ix,1)=t1-f0
          enddo ! ix
        case default
          do iset=1,ndataset
            op_npoly=fit%npoly
            op_pow=fit%pow
            op_a=a(:,iset)
            if(trim(deval%what)=='unshared')then
              where(fit%share)op_a=0.d0
            endif
            do ideriv=1,deval%nderiv
              call deriv_poly(op_npoly,op_pow,op_a)
            enddo ! ideriv
            if(deval%rel)f0=eval_poly(op_npoly,op_pow,op_a,-fit%X0)
            do ix=1,deval%n
              t1=eval_poly(op_npoly,op_pow,op_a,deval%x(ix)-fit%X0)
              if(trim(deval%what)=='sum')then
                f_array(irandom,ix,1)=f_array(irandom,ix,1)+&
                   &tmp_dlist(iset)%dataset%weight*(t1-f0)
              else
                f_array(irandom,ix,iset)=t1-f0
              endif
            enddo ! ix
          enddo ! iset
        end select ! deval%what
      endif ! need_f
    enddo ! irandom

    ! Return coefficients and statistical properties of requested function
    ! of fit.
    if(need_chi2)then
      call characterize_dist(nsample,chi2_array,w_vector,mean=t1,var=var)
      if(present(chi2mean))chi2mean=t1
      if(present(chi2err))chi2err=sqrt(var)
    endif
    if(need_rmsy)then
      call characterize_dist(nsample,rmsy_array,w_vector,mean=t1,var=var)
      if(present(rmsymean))rmsymean=t1
      if(present(rmsyerr))rmsyerr=sqrt(var)
    endif
    do iset=1,ndataset
      if(need_a)then
        do ipoly=1,fit%npoly
          call characterize_dist(nsample,a_array(:,ipoly,iset),w_vector,&
             &mean=t1,var=var)
          if(present(amean))amean(ipoly,iset)=t1
          if(present(aerr))aerr(ipoly,iset)=sqrt(var)
        enddo ! ipoly
      endif ! need_a
      if(need_f)then
        if((trim(deval%what)=='shared'.or.trim(deval%what)=='sum').and.&
           &iset>1)then
          if(present(fmean))fmean(:,iset)=fmean(:,1)
          if(present(ferr))ferr(:,iset)=ferr(:,1)
          if(present(fskew))fskew(:,iset)=fskew(:,1)
          if(present(fkurt))fkurt(:,iset)=fkurt(:,1)
          if(present(fmed))fmed(:,iset)=fmed(:,1)
          if(present(fmean_1s))fmean_1s(:,iset)=fmean_1s(:,1)
          if(present(ferr_1s))ferr_1s(:,iset)=ferr_1s(:,1)
          if(present(fmean_2s))fmean_2s(:,iset)=fmean_2s(:,1)
          if(present(ferr_2s))ferr_2s(:,iset)=ferr_2s(:,1)
        else ! per-set "what" or first set in global "what"
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
        endif ! global "what".and.iset==1 or not
      endif ! need_f
    enddo ! iset

    ! Clean up.
    call kill_dlist(tmp_dlist)

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


  SUBROUTINE read_file(fname,icol_x,icol_y,icol_dx,icol_dy,icol_w,&
     &nsearch,fsearch,search,ndiscr,fdiscr,ndataset,dlist,ierr)
    !-------------------------------------!
    ! Read in the data in the input file. !
    !-------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(in) :: icol_x,icol_y,icol_dx,icol_dy,icol_w,nsearch,&
       &fsearch(nsearch),ndiscr,fdiscr(ndiscr)
    CHARACTER(*),INTENT(in) :: search(nsearch)
    INTEGER,INTENT(inout) :: ndataset,ierr
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    TYPE(dataset_type),POINTER :: dataset
    CHARACTER(8192) line
    CHARACTER(8192),POINTER :: discr(:,:)
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
    nullify(discr)

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
        if(ndataset>0)then
          allocate(tmp_dlist(ndataset))
          tmp_dlist=dlist(1:ndataset)
          deallocate(dlist)
        endif
        allocate(dlist(ndataset+1))
        if(ndataset>0)then
          dlist(1:ndataset)=tmp_dlist
          deallocate(tmp_dlist)
        endif
        ndataset=ndataset+1
        allocate(dlist(ndataset)%dataset)
        dataset=>dlist(ndataset)%dataset
        iset=ndataset
        ! Store discriminators.
        call resize_pointer_char2(8192,(/ndiscr,ndataset/),discr)
        do idiscr=1,ndiscr
          discr(idiscr,iset)=trim(field(fdiscr(idiscr),line))
        enddo ! idiscr
      else
        dataset=>dlist(iset)%dataset
      endif
      ! Enlarge arrays.
      call increment_xy_size(icol_dx>0,icol_dy>0,icol_w>0,dataset%xy)
      i=dataset%xy%nxy
      ! Read data point from string.
      if(icol_x>0)then
        dataset%xy%x(i)=dble_field(icol_x,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of x in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
      else
        dataset%xy%x(i)=dble(i)
      endif
      if(icol_y>0)then
        dataset%xy%y(i)=dble_field(icol_y,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of y in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
      else
        dataset%xy%y(i)=dble(i)
      endif
      if(icol_dx>0)then
        dataset%xy%dx(i)=dble_field(icol_dx,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of dx in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
        if(dataset%xy%dx(i)<0.d0.and..not.are_equal(dataset%xy%dx(i),0.d0))then
          write(6,'(a)')'Found negative dx in "'//trim(fname)//'".'
          write(6,'()')
          ierr=-5
          exit
        endif
      endif ! have_dx
      if(icol_dy>0)then
        dataset%xy%dy(i)=dble_field(icol_dy,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of dy in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
        if(dataset%xy%dy(i)<0.d0.and..not.are_equal(dataset%xy%dy(i),0.d0))then
          write(6,'(a)')'Found negative dy in "'//trim(fname)//'".'
          write(6,'()')
          ierr=-6
          exit
        endif
      endif ! have_dy
      if(icol_w>0)then
        dataset%xy%w(i)=dble_field(icol_w,line,ierr)
        if(ierr/=0)then
          write(6,'(a)')'Failed to parse value of w in "'//trim(fname)//'".'
          write(6,'()')
          exit
        endif
        if(dataset%xy%w(i)<0.d0.or.are_equal(dataset%xy%w(i),0.d0))then
          write(6,'(a)')'Found non-positive w in "'//trim(fname)//'".'
          write(6,'()')
          ierr=-7
          exit
        endif
      endif ! have_w
    enddo ! i

    ! Close file.
    close(io)

    ! Clean up.
    if(associated(discr))deallocate(discr)

  END SUBROUTINE read_file


  ! POLYNOMIAL HANDLING UTILITIES.


  FUNCTION print_poly_sym(fit) RESULT(polystr)
    !---------------------------------------------------!
    ! Returns a string with the symbolic version of the !
    ! fitted polynomial.                                !
    !---------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),INTENT(in) :: fit
    CHARACTER(3+fit%npoly*54) :: polystr
    INTEGER j,ipow
    CHARACTER(40) pwstr
    CHARACTER(6) xstr
    if(abs(fit%X0)<tol_zero)then
      xstr='X'
    else
      xstr='(X-X0)'
    endif
    polystr='Y ='
    do j=1,fit%npoly
      ipow=nint(fit%pow(j))
      if(abs(dble(ipow)-fit%pow(j))<tol_zero)then
        if(ipow==0)then
          pwstr=''
        elseif(ipow==1)then
          pwstr='*'//trim(xstr)
        else
          pwstr='*'//trim(xstr)//'^'//trim(i2s(ipow))
        endif ! fit%pow=0
      else
        write(pwstr,*)fit%pow(j)
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


  FUNCTION print_poly_num(fit,a) RESULT(polystr)
    !----------------------------------------------------!
    ! Returns a string with the numerical version of the !
    ! fitted polynomial in a suitable format for pasting !
    ! into xmgrace.                                      !
    !----------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),INTENT(in) :: fit
    DOUBLE PRECISION,INTENT(in) :: a(fit%npoly)
    CHARACTER(2+fit%npoly*105) :: polystr
    INTEGER j,ipow
    CHARACTER(1) plusstr
    CHARACTER(32) coeffstr
    CHARACTER(72) pwstr
    CHARACTER(36) xstr
    if(abs(fit%X0)<tol_zero)then
      xstr='x'
    else
      write(xstr,*)fit%X0
      xstr='(x-'//trim(adjustl(xstr))//')'
    endif
    polystr='y='
    do j=1,fit%npoly
      ipow=nint(fit%pow(j))
      if(abs(dble(ipow)-fit%pow(j))<tol_zero)then
        if(ipow==0)then
          pwstr=''
        elseif(ipow==1)then
          pwstr='*'//trim(xstr)
        else
          pwstr='*'//trim(xstr)//'^'//trim(i2s(ipow))
        endif ! fit%pow=0
      else
        write(pwstr,*)fit%pow(j)
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


  SUBROUTINE increment_xy_size(have_dx,have_dy,have_w,xy)
    !--------------------------------------!
    ! Increment the size of pointers in an !
    ! xy_type-type pointer by 1.           !
    !--------------------------------------!
    IMPLICIT NONE
    LOGICAL,INTENT(in) :: have_dx,have_dy,have_w
    TYPE(xy_type),POINTER :: xy
    INTEGER nxy
    if(.not.associated(xy))then
      nxy=0
      allocate(xy)
    else
      nxy=xy%nxy
    endif
    nxy=nxy+1
    xy%nxy=nxy
    xy%have_dx=have_dx
    xy%have_dy=have_dy
    xy%have_w=have_w
    call resize_pointer_dble1((/nxy/),xy%x)
    call resize_pointer_dble1((/nxy/),xy%y)
    call resize_pointer_dble1((/nxy/),xy%dx)
    call resize_pointer_dble1((/nxy/),xy%dy)
    call resize_pointer_dble1((/nxy/),xy%w)
    ! Initialize x and y to index by default.
    xy%x(nxy)=dble(nxy)
    xy%y(nxy)=dble(nxy)
    xy%w(nxy)=1.d0
  END SUBROUTINE increment_xy_size


  SUBROUTINE kill_xy(xy)
    !----------------------------------!
    ! Destroy an xy_type-type pointer. !
    !----------------------------------!
    IMPLICIT NONE
    TYPE(xy_type),POINTER :: xy
    if(.not.associated(xy))return
    if(associated(xy%x))deallocate(xy%x)
    if(associated(xy%y))deallocate(xy%y)
    if(associated(xy%dx))deallocate(xy%dx)
    if(associated(xy%dy))deallocate(xy%dy)
    deallocate(xy)
    nullify(xy)
  END SUBROUTINE kill_xy


  SUBROUTINE kill_dataset(dataset)
    !--------------------------------------!
    ! Destroy a dataset_type-type pointer. !
    !--------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),POINTER :: dataset
    if(.not.associated(dataset))return
    call kill_xy(dataset%xy)
    call kill_xy(dataset%txy)
    call kill_xy(dataset%rtxy)
    deallocate(dataset)
    nullify(dataset)
  END SUBROUTINE kill_dataset


  SUBROUTINE kill_dlist(dlist)
    !-------------------------------------------!
    ! Destroy a dataset_list_type-type pointer. !
    !-------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_list_type),POINTER :: dlist(:)
    INTEGER ndataset,iset
    if(.not.associated(dlist))return
    ndataset=size(dlist)
    do iset=1,ndataset
      call kill_dataset(dlist(iset)%dataset)
    enddo ! iset
    deallocate(dlist)
    nullify(dlist)
  END SUBROUTINE kill_dlist


  SUBROUTINE clone_xy(xy1,xy2)
    !-----------------------------------------!
    ! Make an independent copy of xy1 as xy2. !
    !-----------------------------------------!
    IMPLICIT NONE
    TYPE(xy_type),POINTER :: xy1,xy2
    nullify(xy2)
    if(.not.associated(xy1))return
    allocate(xy2)
    xy2%nxy=xy1%nxy
    xy2%have_dx=xy1%have_dx
    xy2%have_dy=xy1%have_dy
    allocate(xy2%x(xy2%nxy),xy2%dx(xy2%nxy),xy2%y(xy2%nxy),xy2%dy(xy2%nxy),&
       &xy2%w(xy2%nxy))
    xy2%x=xy1%x
    xy2%dx=xy1%dx
    xy2%y=xy1%y
    xy2%dy=xy1%dy
    xy2%w=xy1%w
  END SUBROUTINE clone_xy


  SUBROUTINE clone_dataset(dataset1,dataset2)
    !---------------------------------------------------!
    ! Make an independent copy of dataset1 as dataset2. !
    !---------------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_type),POINTER :: dataset1,dataset2
    nullify(dataset2)
    if(.not.associated(dataset1))return
    allocate(dataset2)
    dataset2%itransfx=dataset1%itransfx
    dataset2%itransfy=dataset1%itransfy
    dataset2%weight=dataset1%weight
    dataset2%wexp=dataset1%wexp
    call clone_xy(dataset1%xy,dataset2%xy)
    call clone_xy(dataset1%txy,dataset2%txy)
    call clone_xy(dataset1%rtxy,dataset2%rtxy)
  END SUBROUTINE clone_dataset


  SUBROUTINE clone_dlist(dlist1,dlist2)
    !-----------------------------------------------!
    ! Make an independent copy of dlist1 as dlist2. !
    !-----------------------------------------------!
    IMPLICIT NONE
    TYPE(dataset_list_type),POINTER :: dlist1(:),dlist2(:)
    INTEGER iset,ndataset
    nullify(dlist2)
    if(.not.associated(dlist1))return
    ndataset=size(dlist1)
    allocate(dlist2(ndataset))
    do iset=1,ndataset
      call clone_dataset(dlist1(iset)%dataset,dlist2(iset)%dataset)
    enddo ! iset
  END SUBROUTINE clone_dlist


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
    !--------------------------------------------------------!
    ! Return the N-th field in string LINE, where the fields !
    ! are separated by one or more spaces.  An empty string  !
    ! is returned if N<1.                                    !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    CHARACTER(*),INTENT(in) :: line
    CHARACTER(len(line)) :: field
    CHARACTER(len(line)) remainder
    INTEGER i,k
    ! Initialize.
    field=''
    if(n<1)return
    remainder=adjustl(line)
    ! Loop over fields.
    i=0
    do
      i=i+1
      if(remainder(1:1)=='"')then
        ! Potentially start of a multi-word field.
        ! Locate end of field (quote followed by <space> or <EOL>).
        k=index(trim(remainder(2:)),'" ')+1
        if(k==1)then
          ! quote-space not found, so see if there is a quote at EOL.
          k=len_trim(remainder)
          if(remainder(k:k)/='"')k=1
        endif
        if(k>1)then
          ! Found end of field.
          if(i==n)then
            field=trim(remainder(2:k-1))
          else
            remainder=adjustl(remainder(k+1:))
          endif
          cycle
        endif
      endif
      ! Single-word field.
      ! Locate end of field.
      k=scan(trim(remainder),' ')
      if(k==0)then
        ! End is EOL.
        if(i==n)field=trim(remainder)
        return
      elseif(i==n)then
        field=trim(remainder(1:k-1))
        return
      else
        remainder=adjustl(remainder(k:))
      endif
    enddo
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
    int_field=parse_int(field(ifield,command),ierr)
  END FUNCTION int_field


  INTEGER FUNCTION parse_int(string,ierr)
    !--------------------------------------!
    ! Parse a string to obtain an integer. !
    !--------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: ierr
    read(string,*,iostat=ierr)parse_int
  END FUNCTION parse_int


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
    dble_field=parse_dble(field(ifield,command),ierr)
  END FUNCTION dble_field


  DOUBLE PRECISION FUNCTION parse_dble(string,ierr)
    !-----------------------------------------------------!
    ! Parse a string to obtain a double-precision number. !
    ! Supports fractions, e.g., string="1/3".             !
    !-----------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: ierr
    INTEGER ipos
    DOUBLE PRECISION t1,t2
    ipos=scan(string,'/')
    if(ipos==0)then
      read(string,*,iostat=ierr)parse_dble
    else
      ierr=-1
      if(ipos<2)return
      if(ipos>=len_trim(string))return
      read(string(1:ipos-1),*,iostat=ierr)t1
      if(ierr/=0)return
      read(string(ipos+1:),*,iostat=ierr)t2
      if(ierr/=0)return
      if(are_equal(t2,0.d0))then
        ierr=-1
        return
      endif
      parse_dble=t1/t2
    endif
  END FUNCTION parse_dble


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
    x=parse_dble(cx,ierr)
    if(ierr/=0)return
    y=parse_dble(cy,ierr)
    if(ierr/=0)return
    are_equal_string=are_equal(x,y,tol)
  END FUNCTION are_equal_string


  SUBROUTINE parse_range(string,drange)
    !------------------------------------------------!
    ! Parse a data range out of a string. The string !
    ! can take the form:                             !
    ! * <variable><operation><threshold>             !
    !   e.g., X>1, y<=3, X]4 (the latter means "the  !
    !   four largest values of X").                  !
    !------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    TYPE(range_type),INTENT(inout) :: drange
    CHARACTER(len_trim(string)) remainder
    INTEGER ierr

    ! Initialize.
    drange%var='X'
    drange%op=''
    drange%thres=0.d0
    drange%size=0
    drange%no_rhs=.false.
    remainder=trim(adjustl(string))

    ! Extract variable.
    select case(remainder(1:1))
    case('x','y','X','Y')
      drange%var=remainder(1:1)
    case default
      return
    end select
    remainder=remainder(2:)

    ! Check sort operation and set range parameters.
    select case(remainder(1:1))
    case('<','>')
      if(remainder(2:2)=='=')then
        drange%op=remainder(1:2)
        remainder=remainder(3:)
      else
        drange%op=remainder(1:1)
        remainder=remainder(2:)
      endif
      drange%thres=parse_dble(remainder,ierr)
      if(ierr/=0)drange%no_rhs=.true.
    case('[',']')
      drange%op=remainder(1:1)
      remainder=remainder(2:)
      drange%size=parse_int(remainder,ierr)
      if(ierr/=0)drange%no_rhs=.true.
    case default
      continue
    end select

  END SUBROUTINE parse_range


  SUBROUTINE parse_xeval(string,deval)
    !---------------------------------------------!
    ! Parse a list of x values from a string. The !
    ! string can take the form:                   !
    ! * <variable>=<comma-separated-list>, e.g.,  !
    !   X=1,2,3  gives  X={1,2,3}                 !
    ! * <variable>=<first>:<last>:<count>, e.g.,  !
    !   X=1:2:3  gives  X={1,1.5,2}               !
    !---------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    TYPE(eval_type),INTENT(inout) :: deval
    CHARACTER(len_trim(string)) remainder
    INTEGER ipos,it1,i,ierr
    DOUBLE PRECISION t1,t2

    ! Initialize.
    if(associated(deval%x))deallocate(deval%x)
    nullify(deval%x)
    deval%var='X'
    deval%n=0
    remainder=trim(adjustl(string))

    ! Extract variable.
    select case(remainder(1:1))
    case('X') ! FIXME - 'x' too?
      deval%var=remainder(1:1)
    case default
      return
    end select
    remainder=remainder(2:)

    ! Verify = sign.
    select case(remainder(1:1))
    case('=')
    case default
     return
    end select
    remainder=remainder(2:)

    ! Get whether this is a list or a range.
    ipos=scan(remainder,':')
    if(ipos>0)then
      t1=parse_dble(remainder(1:ipos-1),ierr)
      if(ierr/=0)return
      remainder=remainder(ipos+1:)
      ipos=scan(remainder,':')
      if(ipos<1)return
      t2=parse_dble(remainder(1:ipos-1),ierr)
      if(ierr/=0)return
      remainder=remainder(ipos+1:)
      it1=parse_int(remainder,ierr)
      if(ierr/=0)return
      if(it1<1)return
      if(it1==1.and..not.are_equal(t1,t2))return
      deval%n=it1
      allocate(deval%x(deval%n))
      if(deval%n==1)then
        deval%x(1)=t1
      else
        do i=1,deval%n
          deval%x(i)=t1+(t2-t1)*(dble(i-1)/dble(deval%n-1))
        enddo ! i
      endif
    else
      ! Assume comma-separated list.
      do while(len_trim(remainder)>0)
        ipos=scan(remainder,',')
        if(ipos<1)ipos=len_trim(remainder)+1
        t1=parse_dble(remainder(1:ipos-1),ierr)
        if(ierr/=0)then
          if(associated(deval%x))then
            deallocate(deval%x)
            nullify(deval%x)
          endif
          return
        endif
        deval%n=deval%n+1
        call resize_pointer_dble1((/deval%n/),deval%x)
        deval%x(deval%n)=t1
        remainder=remainder(ipos+1:)
      enddo
    endif

  END SUBROUTINE parse_xeval


  CHARACTER(7) FUNCTION type_string(have_dx,have_dy,have_w)
    !------------------------------------------------!
    ! Produce a dataset type string "x[dx]y[dy][w]". !
    !------------------------------------------------!
    IMPLICIT NONE
    LOGICAL,INTENT(in) :: have_dx,have_dy,have_w
    type_string='x'
    if(have_dx)type_string=trim(type_string)//'dx'
    type_string=trim(type_string)//'y'
    if(have_dy)type_string=trim(type_string)//'dy'
    if(have_w)type_string=trim(type_string)//'w'
  END FUNCTION type_string


  SUBROUTINE parse_type_string(string,ipos_x,ipos_dx,ipos_y,ipos_dy,ipos_w,&
     &ierr)
    !---------------------------------------------------------!
    ! Parse a dataset type string "x[dx]y[dy][w]" into flags. !
    !---------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: ipos_x,ipos_dx,ipos_y,ipos_dy,ipos_w,ierr
    ! Local variables.
    CHARACTER(len_trim(string)) remainder
    INTEGER ipos
    ierr=0
    ipos_x=0
    ipos_dx=0
    ipos_y=0
    ipos_dy=0
    ipos_w=0
    remainder=adjustl(trim(string))
    ipos=0
    do while(len_trim(remainder)>0)
      ipos=ipos+1
      if(remainder(1:1)=='x')then
        if(ipos_x>0)ierr=1
        ipos_x=ipos
        remainder=adjustl(remainder(2:))
      elseif(remainder(1:1)=='y')then
        if(ipos_y>0)ierr=1
        ipos_y=ipos
        remainder=adjustl(remainder(2:))
      elseif(remainder(1:2)=='dx')then
        if(ipos_dx>0)ierr=1
        ipos_dx=ipos
        remainder=adjustl(remainder(3:))
      elseif(remainder(1:2)=='dy')then
        if(ipos_dy>0)ierr=1
        ipos_dy=ipos
        remainder=adjustl(remainder(3:))
      elseif(remainder(1:1)=='w')then
        if(ipos_w>0)ierr=1
        ipos_w=ipos
        remainder=adjustl(remainder(2:))
      else
        ierr=1
      endif
      if(ierr/=0)return
    enddo
  END SUBROUTINE parse_type_string


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


  ! PRETTY-PRINT ROUTINE.


  SUBROUTINE pprint(text,indent1,indent)
    !-------------------------------------------------------------!
    ! Print TEXT to stdout, folding lines at column 79 and using  !
    ! an indentation of INDENT1 spaces on the first line and      !
    ! INDENT on the rest.  INDENT1 defaults to INDENT, and INDENT !
    ! defaults to zero.                                           !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: text
    INTEGER,INTENT(in),OPTIONAL :: indent1,indent
    INTEGER,PARAMETER :: line_width=79
    CHARACTER(line_width) line
    CHARACTER(len(text)) remainder,word
    INTEGER ind1,indn,ind,ipos
    if(len_trim(text)==0)then
      write(6,'()')
      return
    endif
    indn=0
    if(present(indent))indn=max(indent,0)
    ind1=indn
    if(present(indent1))ind1=max(indent1,0)
    ind=ind1
    line=''
    remainder=adjustl(text)
    do while(len_trim(remainder)>0)
      ipos=scan(trim(remainder),' ')
      if(ipos==0)then
        ipos=len_trim(remainder)+1
        word=trim(remainder)
        remainder=''
      else
        word=remainder(1:ipos-1)
        remainder=adjustl(remainder(ipos:))
      endif
      if(len_trim(line)==0)then
        ! Ensure overlong words get flushed without passing through buffer.
        do while(ind+len_trim(word)>line_width)
          write(6,'(a)')repeat(' ',ind)//word(1:line_width-ind)
          word=word(line_width-ind+1:)
          ind=indn
        enddo
        ! Start new line buffer.
        line=repeat(' ',ind)//trim(word)
        ind=indn
      else
        if(len_trim(line)+1+len_trim(word)>line_width)then
          ! Flush current line buffer.
          write(6,'(a)')trim(line)
          ! Ensure overlong words get flushed without passing through buffer.
          do while(ind+len_trim(word)>line_width)
            write(6,'(a)')repeat(' ',ind)//word(1:line_width-ind)
            word=word(line_width-ind+1:)
            ind=indn
          enddo
          ! Start new line buffer.
          line=repeat(' ',ind)//trim(word)
          ind=indn
        else
          ! Add line to buffer.
          line=trim(line)//' '//trim(word)
        endif
      endif
    enddo
    if(len_trim(line)>0)write(6,'(a)')trim(line)
  END SUBROUTINE pprint


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
