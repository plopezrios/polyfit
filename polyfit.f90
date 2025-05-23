PROGRAM polyfit
  !----------------------------------------------------!
  ! POLYFIT                                            !
  ! =======                                            !
  ! Toolbox for performing and manipulating polynomial !
  ! fits of data with statistical uncertainties        !
  !                                                    !
  ! PLR 09.2015                                        !
  !----------------------------------------------------!
  IMPLICIT NONE

  ! Global variables
  ! ================

  ! * Precision for real-to-integer conversions and exponent comparisons.
  DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-10

  ! * Scale transformations.
  INTEGER,PARAMETER :: NTRANSF=3
  INTEGER,PARAMETER :: ITRANSF_NONE=0,ITRANSF_REC=1,ITRANSF_LOG=2,ITRANSF_EXP=3
  CHARACTER(32),PARAMETER :: TRANSF_FORMX(0:NTRANSF)=&
     &(/'x     ','1/x   ','log(x)','exp(x)'/)
  CHARACTER(32),PARAMETER :: TRANSF_FORMY(0:NTRANSF)=&
     &(/'y     ','1/y   ','log(y)','exp(y)'/)
  LOGICAL,PARAMETER :: TRANSF_REQ_NONZERO(0:NTRANSF)=&
     &(/.false.,.true.,.true.,.false./)
  LOGICAL,PARAMETER :: TRANSF_REQ_POSITIVE(0:NTRANSF)=&
     &(/.false.,.false.,.true.,.false./)

  ! Derived types
  ! =============

  ! * xy[dx][dy][w] data.  NB, all pointers are allocated regardless of
  !   have_{dx,dy,w}, to simplify logic.
  TYPE xy_type
    INTEGER :: nxy=0
    LOGICAL :: have_dx=.false.,have_dy=.false.,have_w=.false.
    DOUBLE PRECISION,POINTER :: x(:)=>null(),y(:)=>null(),dx(:)=>null(),&
       &dy(:)=>null(),w(:)=>null()
  END TYPE xy_type

  ! * Dataset - comprising original data, dataset weight, weight exponent,
  !   x/y transformation info, transformed data, range-restricted
  !   transformed data, and fit index.
  !   FIXME - compute txy / rtxy on the fly and drop from derived type.
  TYPE dataset_type
    INTEGER :: itransfx=ITRANSF_NONE,itransfy=ITRANSF_NONE,ifit=1
    DOUBLE PRECISION :: wexp=1.d0
    DOUBLE PRECISION :: weight=1.d0
    TYPE(xy_type),POINTER :: xy=>null() ! original data
    TYPE(xy_type),POINTER :: txy=>null() ! transformed data
    TYPE(xy_type),POINTER :: rtxy=>null() ! range-restricted transformed data
    ! Quasi-random noise handling.
    LOGICAL :: apply_qrandom=.false.
    DOUBLE PRECISION :: qrandom_exp=0.d0,qrandom_centre=0.d0
  END TYPE dataset_type

  ! * Dataset list item - this enables defining "arrays" of datasets.
  TYPE dataset_list_type
    TYPE(dataset_type),POINTER :: dataset=>null()
  END TYPE dataset_list_type

  ! * Fit form.
  ! FIXME - allow setting fixed coefficients.  E.g., uncertainty
  ! as a function of MC sample size in log-log scale should have a -1/2
  ! slope; this could be done in linear scale by fitting with an (1/X)^1/2
  ! form, but Y is allowed to be negative which is unphysical, so log-log
  ! with fixed slope would be far better.
  TYPE fit_form_type
    ! Polynomial order, powers of X in polynomial, shared parameter indices.
    INTEGER :: npoly=0
    DOUBLE PRECISION,ALLOCATABLE :: pow(:)
    LOGICAL,ALLOCATABLE :: share(:)
  END TYPE fit_form_type

  ! * Fit form list item - this enables defining "arrays" of fit forms.
  TYPE fit_form_list_type
    TYPE(fit_form_type),POINTER :: fit=>null()
  END TYPE fit_form_list_type

  ! * Global range restriction info (applies to all datasets).
  !   FIXME - implement per-dataset ranges?
  TYPE range_type
    CHARACTER :: var='X'
    CHARACTER(2) :: op=''
    LOGICAL :: no_rhs=.false.
    CHARACTER(3) :: thres_op='non'
    DOUBLE PRECISION :: thres=0.d0
    INTEGER :: size=0
  END TYPE range_type

  ! * Information defining functions of fit to evaluate.
  TYPE eval_type
    ! Target function of fit (function/shared).
    CHARACTER(10) :: what='function'
    ! Relative to value at X0?
    LOGICAL :: rel=.false.
    DOUBLE PRECISION :: Xrel=0.d0
    ! Derivative (0 for function, 1 for first, etc.).
    INTEGER :: nderiv=0
    ! Plot range definition.
    INTEGER :: n=0 ! number of points
    CHARACTER :: var='X' ! variable defining grid
    DOUBLE PRECISION,POINTER :: x(:)=>null() ! grid point values
  END TYPE eval_type

  ! * Intersect ranges.
  TYPE intersect_range_type
    DOUBLE PRECISION :: x1=0.d0,x2=0.d0,xmid=0.d0
    LOGICAL :: have_x1=.false.,have_x2=.false.,have_xmid=.false.
  END TYPE intersect_range_type

  ! * Global parameters.
  TYPE global_params_type
    ! Number of Monte Carlo samples
    INTEGER :: nsample=5000
    ! x-offset parameter.
    DOUBLE PRECISION :: X0=0.d0
    CHARACTER(64) :: X0_string='0'
  END TYPE global_params_type

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
    ! Fit form list.
    INTEGER nfit
    TYPE(fit_form_list_type),POINTER :: flist(:),tmp_flist(:)
    ! All-set settings.
    TYPE(range_type) drange,tdrange
    ! Function evaluation.
    TYPE(eval_type) deval
    ! Monte Carlo parameters.
    TYPE(global_params_type) glob
    ! Defaults.
    TYPE(dataset_type) dataset_model,dataset_default
    TYPE(global_params_type) glob_default
    TYPE(intersect_range_type) :: intersect_range_default
    INTEGER itransfx_default,itransfy_default
    DOUBLE PRECISION wexp_default
    ! Input file information.
    INTEGER ncolumn,icol_x,icol_y,icol_dx,icol_dy,icol_w
    CHARACTER(256) fname,dname
    ! Input file search facility.
    LOGICAL neg
    INTEGER,PARAMETER :: search_size=8192
    INTEGER nsearch,ndiscr
    INTEGER,POINTER :: fsearch(:),fdiscr(:)
    CHARACTER(search_size),POINTER :: search(:)
    ! Input echo.
    LOGICAL input_echo
    ! Intersect command.
    LOGICAL use_mix,use_Xrel
    TYPE(intersect_range_type) :: intersect_range
    ! Misc variables.
    CHARACTER(7) settype
    CHARACTER(8192) command,token
    LOGICAL, ALLOCATABLE :: lmask(:)
    INTEGER i,j,ifield,ierr,ipos_x,ipos_y,ipos_dx,ipos_dy,ipos_w,ierr1,&
       &ierr2,nxy,itransf,iset,jset,npoly,ifit,jfit,max_npoly,ci_choice
    INTEGER, POINTER :: ilist(:)
    DOUBLE PRECISION t1,t2,wexp,find_target,Xrel
    DOUBLE PRECISION, POINTER :: rlist(:)

    ! Write header.
    write(6,'(a)')'===================================='
    write(6,'(a)')'POLYFIT - polynomial fitting toolbox'
    write(6,'(a)')'===================================='
    write(6,'()')
    write(6,'(a)')'Type "help" for a list of commands.'
    write(6,'()')

    ! Initialize.
    ndataset=0
    nullify(dlist)
    itransfx_default=ITRANSF_NONE
    itransfy_default=ITRANSF_NONE
    wexp_default=1.d0
    nfit=1
    allocate(flist(nfit))
    allocate(flist(1)%fit)
    flist(1)%fit%npoly=2
    allocate(flist(1)%fit%pow(flist(1)%fit%npoly),&
       &flist(1)%fit%share(flist(1)%fit%npoly))
    flist(1)%fit%pow(1:2)=(/0.d0,1.d0/)
    flist(1)%fit%share=.false.
    nullify(fsearch,fdiscr,search)
    input_echo=.false.
    call init_random(1)

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
      write(6,'(a)',advance='no')'POLYFIT> '
      read(5,'(a)',iostat=ierr)command
      if(ierr/=0)then
        write(6,'()')
        call quit()
      endif
      command=adjustl(command)
      if(input_echo)write(6,'(a)')trim(command)
      if(command(1:1)=='#')cycle
      write(6,'()')

      ! Execute command.
      select case(trim(field(1,command)))

      case('inspect')
        ! Report number of lines and columns in data file.
        fname=trim(field(2,command))
        call check_file(fname,nxy,ncolumn)
        select case(nxy)
        case(-1)
          call msg('Could not open file "'//trim(fname)//'".')
          cycle user_loop
        case(-2)
          call msg('Problem reading file "'//trim(fname)//'".')
          cycle user_loop
        case(-3)
          call msg('Column count problem in file "'//trim(fname)//'".')
          cycle user_loop
        case(0)
          call msg('File "'//trim(fname)//'" contains no useful data.')
          cycle user_loop
        end select
        call msg('File "'//trim(fname)//'" contains '//trim(i2s(nxy))//&
           &' lines with '//trim(i2s(ncolumn))//' columns.')

      case('load','wload')
        ! Load data from specified columns of data file.
        fname=trim(field(2,command))
        nxy=0
        ncolumn=0
        if(fname(1:2)/='<<')then
          call check_file(fname,nxy,ncolumn)
          select case(nxy)
          case(-1)
            call msg('Could not open file "'//trim(fname)//'".')
            cycle user_loop
          case(-2)
            call msg('Problem reading file "'//trim(fname)//'".')
            cycle user_loop
          case(-3)
            call msg('Column count problem in file "'//trim(fname)//'".')
            cycle user_loop
          case(0)
            call msg('File "'//trim(fname)//'" contains no useful data.')
            cycle user_loop
          end select
        endif

        ! Initialize search facility.
        nsearch=0
        ndiscr=0

        ! Get column selection.
        ! Initialize to default column indices.
        if(trim(field(1,command))=='wload')then
          settype='y'
        else
          select case(ncolumn)
          case(0)
            settype=''
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
              call msg('Syntax error in load command: "using" without &
                 &preceeding "type" not allowed.')
              cycle user_loop
            endif
            nullify(ilist)
            ifield=ifield+1
            call parse_ilist(field(ifield,command),ilist)
            if(.not.associated(ilist))then
              call msg('Problem parsing "using" index.')
              cycle user_loop
            endif
            if(size(ilist,1)/=1)then
              call msg('Expected one column index after "using".')
              deallocate(ilist)
              cycle user_loop
            endif
            icol_y=ilist(1)
            deallocate(ilist)
            ! Check column index.
            if((ncolumn>0.and.icol_y>ncolumn).or.icol_y<0)then
              call msg('Column index out of range.')
              cycle user_loop
            endif
          case('type')
            if(trim(field(1,command))=='wload')then
              call msg('Syntax error in wload command: "type" is not &
                 &an allowed subcommand.')
              cycle user_loop
            endif
            ifield=ifield+1
            call parse_type_string(field(ifield,command),ipos_x,ipos_dx,&
               &ipos_y,ipos_dy,ipos_w,ierr)
            if(ierr/=0)then
              call msg('Syntax error in load command: unrecognized &
                 &dataset type.')
              cycle user_loop
            endif
            if(trim(field(ifield+1,command))=='using')then
              ifield=ifield+2
              nullify(ilist)
              call parse_ilist(field(ifield,command),ilist,sortuniq=.false.)
              if(.not.associated(ilist))then
                call msg('Problem parsing "using" indices.')
                cycle user_loop
              endif
              if(count((/ipos_x,ipos_dx,ipos_y,ipos_dy,ipos_w/)/=0)/=&
                 &size(ilist,1))then
                call msg('Wrong number of indices provided for "using" &
                   &subcommmand.')
                cycle user_loop
              endif
              icol_x=0
              icol_dx=0
              icol_y=0
              icol_dy=0
              icol_w=0
              if(ipos_x>0)icol_x=ilist(ipos_x)
              if(ipos_dx>0)icol_dx=ilist(ipos_dx)
              if(ipos_y>0)icol_y=ilist(ipos_y)
              if(ipos_dy>0)icol_dy=ilist(ipos_dy)
              if(ipos_w>0)icol_w=ilist(ipos_w)
              deallocate(ilist)
            else
              icol_x=ipos_x
              icol_dx=ipos_dx
              icol_y=ipos_y
              icol_dy=ipos_dy
              icol_w=ipos_w
            endif
            ! Check column indices.
            if((ncolumn>0.and.(icol_x>ncolumn.or.icol_dx>ncolumn.or.&
               &               icol_y>ncolumn.or.icol_dy>ncolumn.or.&
               &               icol_w>ncolumn)).or.&
               &icol_x<0.or.icol_dx<0.or.icol_y<0.or.icol_dy<0.or.&
               &icol_w<0)then
              call msg('Column indices out of range.')
              cycle user_loop
            endif

          case('where')
            ! Add a search clause.
            ifield=ifield+1
            call parse_search_clause(field(ifield,command),i,token,neg)
            if(i<1)then
              call msg('Syntax error: could not parse argument to "where" &
                 &subcommand.')
              cycle user_loop
            endif
            if(i<1.or.(i>ncolumn.and.ncolumn>0))then
              call msg('Column index out of range in "where" subcommand.')
              cycle user_loop
            endif
            nsearch=nsearch+1
            call resize_pointer_int1((/nsearch/),fsearch)
            call resize_pointer_char1(search_size,(/nsearch/),search)
            fsearch(nsearch)=i
            if(neg)fsearch(nsearch)=-fsearch(nsearch)
            search(nsearch)=token

          case('by')
            ! Add a discriminator clause.
            if(trim(field(1,command))=='wload')then
              call msg('Syntax error: "by" is not an allowed subcommand of &
                 &"wload".')
              cycle user_loop
            endif
            ifield=ifield+1
            call parse_colnum(field(ifield,command),i)
            if(i<1)then
              call msg('Syntax error: could not parse argument to "by" &
                 &subcommand.')
              cycle user_loop
            endif
            if(i<1.or.(i>ncolumn.and.ncolumn>0))then
              call msg('Column index out of range in "by" subcommand.')
              cycle user_loop
            endif
            ndiscr=ndiscr+1
            call resize_pointer_int1((/ndiscr/),fdiscr)
            fdiscr(ndiscr)=i
            ifield=ifield+1
          case('')
            exit
          case default
            call msg('Syntax error in load command: unknown subcommand "'&
               &//trim(field(ifield,command))//'".')
            cycle user_loop
          end select
        enddo ! ifield

        ! Read data.
        if(icol_y<1)then
          call msg('Must provide a set type (e.g., "type xydy").')
          cycle user_loop
        endif
        if(.not.associated(search))allocate(fsearch(nsearch),search(nsearch))
        if(.not.associated(fdiscr))allocate(fdiscr(ndiscr))
        call read_file(fname,icol_x,icol_y,icol_dx,icol_dy,icol_w,nsearch,&
           &fsearch,search,ndiscr,fdiscr,dataset_model,file_ndataset,&
           &file_dlist,ierr)
        if(ierr/=0)cycle user_loop
        if(file_ndataset<1)then
          call msg('No data loaded.')
          cycle user_loop
        endif

        if(field(1,command)=='wload')then
          ! We are loading dataset weights, so lets see if we have the
          ! right number of them.
          if(file_ndataset/=1)then
            call msg('Problem with wload: loaded too many datasets.')
            call kill_dlist(file_dlist)
            cycle user_loop
          endif
          if(file_dlist(1)%dataset%xy%nxy/=ndataset)then
            call msg('Problem with wload: loaded '//&
               &trim(i2s(file_dlist(1)%dataset%xy%nxy))//' weights but &
               &expected '//trim(i2s(ndataset))//'.')
            call kill_dlist(file_dlist)
            cycle user_loop
          endif
          ! Check that dataset weights are > 0.
          if(any(file_dlist(1)%dataset%xy%y<=0.d0))then
            call msg('Problem with wload: found negative weights.')
            call kill_dlist(file_dlist)
            cycle user_loop
          endif
          ! Set weights.
          do iset=1,ndataset
            dlist(iset)%dataset%weight=file_dlist(1)%dataset%xy%y(iset)
          enddo ! iset
          ! Clean up and loop here to avoid giant if-block below.
          call kill_dlist(file_dlist)
          if(fname(1:2)/='<<')then
            call msg('Loaded data from "'//trim(fname)//&
               &'" as dataset weights.')
          else
            call msg('Loaded data from stdin as dataset weights.')
          endif
          cycle user_loop
        endif

        ! Check data are compatible with current transformations.
        do iset=1,file_ndataset
          dataset=>file_dlist(iset)%dataset
          xy=>dataset%xy
          dataset%itransfx=itransfx_default
          dataset%itransfy=itransfy_default
          dataset%wexp=wexp_default
          if(TRANSF_REQ_NONZERO(dataset%itransfy))then
            if(any(eq_dble(xy%x,0.d0)))then
              dataset%itransfx=ITRANSF_NONE
              call msg('Note: using linear X for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains x=0.')
            endif
          endif
          if(TRANSF_REQ_NONZERO(dataset%itransfy))then
            if(any(eq_dble(xy%y,0.d0)))then
              dataset%itransfy=ITRANSF_NONE
              call msg('Note: using linear Y for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains y=0.')
            endif
          endif
          if(TRANSF_REQ_POSITIVE(dataset%itransfx))then
            if(any(xy%x<0.d0))then
              dataset%itransfx=ITRANSF_NONE
              call msg('Note: using linear X for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains x<0.')
            endif
          endif
          if(TRANSF_REQ_POSITIVE(dataset%itransfy))then
            if(any(xy%y<0.d0))then
              dataset%itransfy=ITRANSF_NONE
              call msg('Note: using linear Y for set '//&
                 &trim(i2s(iset+ndataset))//' since it contains y<0.')
            endif
          endif
          if(dataset%xy%have_w)then
            if(any(lt_dble(xy%w,0.d0)))then
              if(neq_dble(dataset%wexp,0.d0))then
                dataset%wexp=0.d0
                call msg('Note: zeroing weight exponent for set '//&
                   &trim(i2s(iset+ndataset))//' since it contains w<0.')
              endif
            elseif(any(eq_dble(xy%w,0.d0)))then
              if(lt_dble(dataset%wexp,0.d0))then
                dataset%wexp=1.d0
                call msg('Note: setting weight exponent to 1 for set '//&
                   &trim(i2s(iset+ndataset))//' since it contains w=0.')
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
          if(fname(1:2)/='<<')then
            write(6,'(a)')'Loaded data from "'//trim(fname)//'" as dataset #'//&
               &trim(i2s(iset))//', type '//trim(type_string(xy%have_dx,&
               &xy%have_dy,xy%have_w))//', '//trim(i2s(xy%nxy))//' data.'
          else
            write(6,'(a)')'Loaded data from stdin as dataset #'//&
               &trim(i2s(iset))//', type '//trim(type_string(xy%have_dx,&
               &xy%have_dy,xy%have_w))//', '//trim(i2s(xy%nxy))//' data.'
          endif
        enddo ! iset
        call refresh_X0(ndataset,dlist,glob)
        write(6,'()')

      case('unload')
        if(ndataset<1)then
          call msg('No datasets loaded.')
          cycle user_loop
        endif
        allocate(lmask(ndataset))
        if(nfield(command)==1)then
          ! Unload all.
          lmask=.true.
        else
          ! Unload specified sets.
          lmask=.false.
          ifield=1
          do
            ifield=ifield+1
            if(ifield>nfield(command))exit
            i=parse_int(field(ifield,command),ierr)
            if(ierr/=0)then
              call msg('Invalid dataset index.')
              deallocate(lmask)
              cycle user_loop
            endif
            if(i<1.or.i>ndataset)then
              call msg('Dataset index out of range.')
              deallocate(lmask)
              cycle user_loop
            endif
            lmask(i)=.true.
          enddo
        endif
        ! Delete datasets backwards.
        if(.not.all(lmask))then
          allocate(tmp_dlist(count(.not.lmask)))
          tmp_dlist=pack(dlist,.not.lmask)
        endif
        do i=1,ndataset
          if(.not.lmask(i))cycle
          call kill_dataset(dlist(i)%dataset)
        enddo
        deallocate(dlist)
        nullify(dlist)
        if(.not.all(lmask))then
          allocate(dlist(count(.not.lmask)))
          dlist=tmp_dlist
          deallocate(tmp_dlist)
        endif
        ndataset=count(.not.lmask)
        call refresh_X0(ndataset,dlist,glob)
        call msg(trim(i2s(count(lmask)))//' datasets unloaded.')
        deallocate(lmask)

      case('fit')
        if(ndataset<1)then
          call msg('No datasets loaded.')
          cycle user_loop
        endif
        if(nfield(command)>1)then
          call msg('Syntax error: subcommand "'//trim(field(2,command))//&
             &'" not recognized.')
          cycle user_loop
        endif ! nfield>1
        call show_multipoly(ndataset,dlist,drange,nfit,flist,glob)

      case('plot')
        if(ndataset<1)then
          call msg('No datasets loaded.')
          cycle user_loop
        endif
        ! Initialize components.
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
          call msg('Syntax error: unknown function "'//&
             &trim(field(2,command))//'".')
          cycle user_loop
        end select
        ! Parse sub-commands.
        ci_choice=0
        fname='fit.plot'
        ifield=2
        do
          ifield=ifield+1
          if(ifield>nfield(command))exit
          select case(trim(field(ifield,command)))
          case('to')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "to" subcommand must be followed &
                 &by a filename.')
              cycle user_loop
            endif
            fname=field(ifield+1,command)
            ifield=ifield+1
          case('at')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "at" subcommand must be followed &
                 &by a data range.')
              cycle user_loop
            endif
            call parse_xeval(trim(field(ifield+1,command)),deval)
            if(.not.associated(deval%x))then
              call msg('Syntax error: could not parse range.')
              cycle user_loop
            endif
            ifield=ifield+1
          case('wrt')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "wrt" subcommand must be followed &
                 &by X value.')
              cycle user_loop
            endif
            t1=dble_field(ifield+1,command,ierr)
            if(ierr/=0)then
              call msg('Syntax error: could not parse "wrt" value.')
              cycle user_loop
            endif
            deval%rel=.true.
            deval%Xrel=t1
            ifield=ifield+1
          case('one-sigma')
            ci_choice=1
          case('two-sigma')
            ci_choice=2
          case default
            call msg('Syntax error: unknown subcommand "'//&
               &trim(field(ifield,command))//'".')
            cycle user_loop
          end select
        enddo
        call plot_multipoly(ndataset,dlist,drange,nfit,flist,deval,glob,&
           &trim(fname),ci_choice=ci_choice)

      case('probability')
        if(ndataset<1)then
          call msg('No datasets loaded.')
          cycle user_loop
        endif
        ! Initialize components.
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
          call msg('Syntax error: unknown function "'//&
             &trim(field(2,command))//'".')
          cycle user_loop
        end select
        ! Get condition.
        select case(field(3,command))
        case('positive','negative')
          token=field(3,command)
        case default
          call msg('Syntax error: unknown condition "'//&
             &trim(field(3,command))//'".')
          cycle user_loop
        end select
        ! Parse sub-commands.
        ifield=3
        do
          ifield=ifield+1
          if(ifield>nfield(command))exit
          select case(trim(field(ifield,command)))
          case('at')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "at" subcommand must be followed &
                 &by a data range.')
              cycle user_loop
            endif
            call parse_xeval(trim(field(ifield+1,command)),deval)
            if(.not.associated(deval%x))then
              call msg('Syntax error: could not parse range.')
              cycle user_loop
            endif
            ifield=ifield+1
          case('wrt')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "wrt" subcommand must be followed &
                 &by X value.')
              cycle user_loop
            endif
            t1=dble_field(ifield+1,command,ierr)
            if(ierr/=0)then
              call msg('Syntax error: could not parse "wrt" value.')
              cycle user_loop
            endif
            deval%rel=.true.
            deval%Xrel=t1
            ifield=ifield+1
          case default
            call msg('Syntax error: unknown subcommand "'//&
               &trim(field(ifield,command))//'".')
            cycle user_loop
          end select
        enddo
        call prob_multipoly(ndataset,dlist,drange,nfit,flist,deval,glob,&
           &token)

      case('assess')
        if(ndataset<1)then
          call msg('No datasets loaded.')
          cycle user_loop
        endif
        ! Check assessment target.
        select case(trim(field(2,command)))
        case('fit','range','range,fit','fit,range')
        case default
          call msg('Unknown variable to assess "'//&
             &trim(field(2,command))//'".')
          cycle user_loop
        end select
        ! Initialize.
        deval%n=0
        deval%rel=.false.
        tdrange%var='X'
        tdrange%op='<='
        tdrange%thres_op='non'
        tdrange%thres=0.d0
        tdrange%size=0
        tdrange%no_rhs=.true.
        nullify(ilist)
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
              call msg('Syntax error: unknown function "'//&
                 &trim(field(ifield+1,command))//'".')
              cycle user_loop
            end select
            ! Get where to evaluate it at.
            if(trim(field(ifield+2,command))/='at')then
              call msg('Syntax error: missing "at" subcommand.')
              cycle user_loop
            endif
            call parse_xeval(field(ifield+3,command),deval)
            if(.not.associated(deval%x))then
              call msg('Syntax error: problem parsing list of X values.')
              cycle user_loop
            endif
            ifield=ifield+3
          case('by')
            call parse_range(field(ifield+1,command),tdrange)
            if(tdrange%op=='')then
              call msg('Syntax error: problem parsing "by" argument.')
              cycle user_loop
            endif
            ifield=ifield+1
          case('for')
            call parse_ilist(field(ifield+1,command),ilist)
            if(.not.associated(ilist))then
              call msg('Syntax error in <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              if(ilist(i)<1.or.ilist(i)>ndataset)then
                call msg('Set index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
            ifield=ifield+1
          case default
            call msg('Syntax error: unknown subcommand "'//&
                 &trim(field(ifield,command))//'".')
            cycle user_loop
          end select
        enddo ! ifield
        ! Turn ilist into mask.
        allocate(lmask(ndataset))
        lmask=.false.
        if(.not.associated(ilist))then
          lmask=.true.
        else
          lmask(ilist(:))=.true.
          deallocate(ilist)
        endif
        ! Evaluate.
        select case(trim(field(2,command)))
        case('fit')
          call assess_fit(ndataset,dlist,drange,nfit,flist,glob,deval,lmask)
        case('range')
          call assess_range(ndataset,dlist,tdrange,nfit,flist,glob,deval,lmask)
        case('range,fit','fit,range')
          call assess_fit_range(ndataset,dlist,tdrange,nfit,flist,glob,&
             &deval,lmask)
        end select

      case('report')
        if(ndataset<1)then
          call msg('No datasets loaded.')
          cycle user_loop
        endif
        select case(trim(field(2,command)))
        case('range')
          call report_statistics(ndataset,dlist)
        case default
          call msg('Unknown report name "'//trim(field(2,command))//'".')
          cycle user_loop
        end select

      case('evaluate')
        if(ndataset<1)then
          call msg('No datasets loaded.')
          cycle user_loop
        endif
        ! Initialize components.
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
          call msg('Syntax error: unknown function "'//&
             &trim(field(2,command))//'".')
          cycle user_loop
        end select
        ! Get where to evaluate it at.
        ifield=2
        do
          ifield=ifield+1
          if(ifield>nfield(command))exit
          select case(trim(field(ifield,command)))
          case('at')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "at" subcommand must be followed &
                 &by X range.')
              cycle user_loop
            endif
            call parse_xeval(field(ifield+1,command),deval)
            if(.not.associated(deval%x))then
              call msg('Syntax error: could not parse range.')
              cycle user_loop
            endif
            ifield=ifield+1
          case('wrt')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "wrt" subcommand must be followed &
                 &by X value.')
              cycle user_loop
            endif
            t1=dble_field(ifield+1,command,ierr)
            if(ierr/=0)then
              call msg('Syntax error: could not parse "wrt" value.')
              cycle user_loop
            endif
            deval%rel=.true.
            deval%Xrel=t1
            ifield=ifield+1
          case default
            call msg('Syntax error: unknown subcommand "'//&
               &trim(field(ifield,command))//'".')
            cycle user_loop
          end select
        enddo
        ! Check that we have a range.
        if(.not.associated(deval%x))then
          call msg('Must provide "at" subcommand to define range.')
          cycle user_loop
        endif
        ! Perform evaluation.
        call evaluate_fit(ndataset,dlist,nfit,flist,glob,drange,deval)

      case('find')

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
          call msg('Syntax error: unknown function "'//&
             &trim(field(2,command))//'".')
          cycle user_loop
        end select

        ! Get target value.
        find_target=dble_field(3,command,ierr1)
        if(ierr1/=0)then
          call msg('Syntax error: could not parse target value.')
          cycle user_loop
        endif

        ! Parse options.
        intersect_range=intersect_range_default
        use_Xrel=.false.
        ifield=4
        do
          select case(field(ifield,command))
          case("between")
            if(intersect_range%have_x1.or.intersect_range%have_x2)then
              call msg('Left and/or right range limits specified more than &
                 &once.')
              cycle user_loop
            endif
            ifield=ifield+1
            if(trim(field(ifield,command))=='data')then
              t1=minval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t1=min(t1,minval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
              t2=maxval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t2=max(t2,maxval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
            else
              t1=dble_field(ifield,command,ierr1)
              ifield=ifield+1
              t2=dble_field(ifield,command,ierr2)
              if(ierr1/=0.or.ierr2/=0)then
                call msg('Syntax error: could not parse arguments of &
                   &"between" subcommand.')
                cycle user_loop
              endif
              if(t1>=t2)then
                call msg('Syntax error: intersection range has non-positive &
                   &length.')
                cycle user_loop
              endif
            endif
            intersect_range%x1=t1
            intersect_range%x2=t2
            intersect_range%have_x1=.true.
            intersect_range%have_x2=.true.
          case("rightof")
            if(intersect_range%have_x1)then
              call msg('Left range limit specified more than once.')
              cycle user_loop
            endif
            ifield=ifield+1
            if(trim(field(ifield,command))=='data')then
              t1=maxval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t1=max(t1,maxval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
            else
              t1=dble_field(ifield,command,ierr1)
              if(ierr1/=0)then
                call msg('Syntax error: could not parse argument of &
                   &"rightof" subcommand.')
                cycle user_loop
              endif
            endif
            if(intersect_range%have_x2)then
              if(t1>=intersect_range%x2)then
                call msg('Syntax error: intersection range has non-positive &
                   &length.')
                cycle user_loop
              endif
            endif
            intersect_range%x1=t1
            intersect_range%have_x1=.true.
          case("leftof")
            if(intersect_range%have_x2)then
              call msg('Right range limit specified more than once.')
              cycle user_loop
            endif
            ifield=ifield+1
            if(trim(field(ifield,command))=='data')then
              t2=minval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t2=min(t2,minval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
            else
              t2=dble_field(ifield,command,ierr2)
              if(ierr2/=0)then
                call msg('Syntax error: could not parse argument of &
                   &"leftof" subcommand.')
                cycle user_loop
              endif
            endif
            if(intersect_range%have_x1)then
              if(t2<=intersect_range%x1)then
                call msg('Syntax error: intersection range has non-positive &
                   &length.')
                cycle user_loop
              endif
            endif
            intersect_range%x2=t2
            intersect_range%have_x2=.true.
          case("near")
            if(intersect_range%have_xmid)then
              call msg('"near" reference point specified more than once.')
              cycle user_loop
            endif
            ifield=ifield+1
            t1=dble_field(ifield,command,ierr1)
            if(ierr1/=0)then
              call msg('Syntax error: could not parse argument of &
                 &"near" subcommand.')
              cycle user_loop
            endif
            intersect_range%xmid=t1
            intersect_range%have_xmid=.true.
          case('wrt')
            if(nfield(command)<ifield+1)then
              call msg('Syntax error: "wrt" subcommand must be followed &
                 &by X value.')
              cycle user_loop
            endif
            t1=dble_field(ifield+1,command,ierr)
            if(ierr/=0)then
              call msg('Syntax error: could not parse "wrt" value.')
              cycle user_loop
            endif
            use_Xrel=.true.
            Xrel=t1
            ifield=ifield+1
          case('')
            exit
          case default
            call msg('Syntax error: "'//trim(field(ifield,command))//&
               &'" subcommand missing.')
            cycle user_loop
          end select
          ifield=ifield+1
        enddo

        ! Perform intersection.
        call find_fit_value(ndataset,dlist,nfit,flist,glob,drange,&
           &intersect_range,deval,find_target,use_Xrel,Xrel)

      case('intersect')

        ! This is only useful when two or more datasets are loaded.
        if(ndataset<2)then
          call msg('Need at least two datasets to find intersections.')
          cycle user_loop
        endif

        ! Initialize.
        deval%rel=.false.
        deval%what='function'
        deval%nderiv=0
        dname=''
        fname=''

        ! Parse options.
        intersect_range=intersect_range_default
        use_mix=.false.
        ifield=2
        do
          select case(field(ifield,command))
          case("mix","mixpairs")
            use_mix=.true.
          case("between")
            if(intersect_range%have_x1.or.intersect_range%have_x2)then
              call msg('Left and/or right range limits specified more than &
                 &once.')
              cycle user_loop
            endif
            ifield=ifield+1
            if(trim(field(ifield,command))=='data')then
              t1=minval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t1=min(t1,minval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
              t2=maxval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t2=max(t2,maxval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
            else
              t1=dble_field(ifield,command,ierr1)
              ifield=ifield+1
              t2=dble_field(ifield,command,ierr2)
              if(ierr1/=0.or.ierr2/=0)then
                call msg('Syntax error: could not parse arguments of &
                   &"between" subcommand.')
                cycle user_loop
              endif
              if(t1>=t2)then
                call msg('Syntax error: intersection range has non-positive &
                   &length.')
                cycle user_loop
              endif
            endif
            intersect_range%x1=t1
            intersect_range%x2=t2
            intersect_range%have_x1=.true.
            intersect_range%have_x2=.true.
          case("rightof")
            if(intersect_range%have_x1)then
              call msg('Left range limit specified more than once.')
              cycle user_loop
            endif
            ifield=ifield+1
            if(trim(field(ifield,command))=='data')then
              t1=maxval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t1=max(t1,maxval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
            else
              t1=dble_field(ifield,command,ierr1)
              if(ierr1/=0)then
                call msg('Syntax error: could not parse argument of &
                   &"rightof" subcommand.')
                cycle user_loop
              endif
            endif
            if(intersect_range%have_x2)then
              if(t1>=intersect_range%x2)then
                call msg('Syntax error: intersection range has non-positive &
                   &length.')
                cycle user_loop
              endif
            endif
            intersect_range%x1=t1
            intersect_range%have_x1=.true.
          case("leftof")
            if(intersect_range%have_x2)then
              call msg('Right range limit specified more than once.')
              cycle user_loop
            endif
            ifield=ifield+1
            if(trim(field(ifield,command))=='data')then
              t2=minval(dlist(1)%dataset%rtxy%x)
              do iset=2,ndataset
                t2=min(t2,minval(dlist(iset)%dataset%rtxy%x))
              enddo ! iset
            else
              t2=dble_field(ifield,command,ierr2)
              if(ierr2/=0)then
                call msg('Syntax error: could not parse argument of &
                   &"leftof" subcommand.')
                cycle user_loop
              endif
            endif
            if(intersect_range%have_x1)then
              if(t2<=intersect_range%x1)then
                call msg('Syntax error: intersection range has non-positive &
                   &length.')
                cycle user_loop
              endif
            endif
            intersect_range%x2=t2
            intersect_range%have_x2=.true.
          case("near")
            if(intersect_range%have_xmid)then
              call msg('"near" reference point specified more than once.')
              cycle user_loop
            endif
            ifield=ifield+1
            t1=dble_field(ifield,command,ierr1)
            if(ierr1/=0)then
              call msg('Syntax error: could not parse argument of &
                 &"near" subcommand.')
              cycle user_loop
            endif
            intersect_range%xmid=t1
            intersect_range%have_xmid=.true.
          case('plot')
            if(len_trim(fname)>0)then
              call msg('Syntax error: "plot" subcommand found twice.')
              cycle user_loop
            endif
            fname='fitmix.plot'
            do
              select case(trim(field(ifield+1,command)))
              case('to')
                if(nfield(command)<ifield+2)then
                  call msg('Syntax error: "to" subcommand must be followed &
                     &by a filename.')
                  cycle user_loop
                endif
                fname=field(ifield+2,command)
                ifield=ifield+2
              case('at')
                if(nfield(command)<ifield+2)then
                  call msg('Syntax error: "at" subcommand must be followed &
                     &by a data range.')
                  cycle user_loop
                endif
                call parse_xeval(trim(field(ifield+2,command)),deval)
                if(.not.associated(deval%x))then
                  call msg('Syntax error: could not parse range.')
                  cycle user_loop
                endif
                ifield=ifield+2
              case default
                exit
              end select
            enddo
            if(.not.associated(deval%x))then
              call msg('Syntax error: no range given for "plot" subcommand.')
              cycle user_loop
            endif
          case('dump')
            if(len_trim(dname)>0)then
              call msg('Syntax error: "dump" subcommand found twice.')
              cycle user_loop
            endif
            dname='intersection.dump'
            do
              select case(trim(field(ifield+1,command)))
              case('to')
                if(nfield(command)<ifield+2)then
                  call msg('Syntax error: "to" subcommand must be followed &
                     &by a filename.')
                  cycle user_loop
                endif
                dname=field(ifield+2,command)
                ifield=ifield+2
              case default
                exit
              end select
            enddo
          case('')
            exit
          case default
            call msg('Syntax error: "'//trim(field(ifield,command))//&
               &'" subcommand missing.')
            cycle user_loop
          end select
          ifield=ifield+1
        enddo

        ! Check plot requirements.
        if(.not.use_mix.and.len_trim(fname)>0)then
          call msg('"plot" subcommand only available for "intersect mix".')
          cycle user_loop
        endif
        if(use_mix.and.len_trim(dname)>0)then
          call msg('"dump" subcommand not available for "intersect mix".')
          cycle user_loop
        endif

        ! Perform intersection.
        if(.not.use_mix)then
          call intersect_fit(ndataset,dlist,nfit,flist,glob,drange,&
             &intersect_range,dname)
        else
          call intersect_mix_fit(ndataset,dlist,nfit,flist,glob,drange,&
             &intersect_range,deval,fname)
        endif

      case('use')
        if(trim(field(2,command))/='fit'.or.trim(field(4,command))/='for'.or.&
           &nfield(command)/=5)then
          call msg('Syntax: use fit <i> for <set-list>')
          cycle user_loop
        endif
        ifit=int_field(3,command,ierr)
        if(ierr/=0)then
          call msg('Fit index should be an integer.')
          cycle user_loop
        endif
        if(ifit<1.or.ifit>nfit)then
          call msg('Fit index out of range.')
          cycle user_loop
        endif
        nullify(ilist)
        call parse_ilist(field(5,command),ilist)
        if(.not.associated(ilist))then
          if(ierr/=0)then
            call msg('Problem parsing <set-list>.')
            cycle user_loop
          endif
        endif
        do jset=1,size(ilist,1)
          iset=ilist(jset)
          if(iset<1.or.iset>ndataset)then
            call msg('Set index out of range.')
            deallocate(ilist)
            cycle user_loop
          endif
        enddo ! iset
        ! See if any fits become unused.
        allocate(lmask(nfit))
        lmask=.false.
        do iset=1,ndataset
          jfit=dlist(iset)%dataset%ifit
          if(any(ilist==iset))jfit=ifit
          lmask(jfit)=.true.
        enddo ! iset
        ! Report removals in multiple-fit cases.
        if(any(.not.lmask))then
          write(6,'(a)',advance='no')'Note: removing fits:'
          do jfit=1,nfit
            if(lmask(jfit))cycle
            write(6,'(a)',advance='no')' '//trim(i2s(jfit))
          enddo ! jfit
          write(6,'()')
        endif
        ! Apply 'use'.
        do iset=1,ndataset
          if(any(ilist==iset))dlist(iset)%dataset%ifit=ifit
        enddo ! iset
        ! Sort out removals.
        if(any(.not.lmask))then
          ! Reindex fits.
          do iset=1,ndataset
            jfit=dlist(iset)%dataset%ifit
            dlist(iset)%dataset%ifit=count(lmask(1:jfit))
          enddo ! iset
          ! Remove fit forms.
          allocate(tmp_flist(count(lmask)))
          do jfit=1,nfit
            if(.not.lmask(jfit))cycle
            tmp_flist(count(lmask(1:jfit)))%fit=>flist(jfit)%fit
          enddo
          deallocate(flist)
          nfit=count(lmask)
          allocate(flist(nfit))
          do jfit=1,nfit
            flist(jfit)%fit=>tmp_flist(jfit)%fit
          enddo
          deallocate(tmp_flist)
        endif ! any removals
        ! Report.
        write(6,'(a)',advance='no')'Using fit #'//trim(i2s(ifit))
        jfit=count(lmask(1:ifit))
        if(ifit/=jfit)write(6,'(a)',advance='no')' (reindexed as fit #'//&
           &trim(i2s(jfit))//')'
        write(6,'(a)',advance='no')' for sets:'
        do jset=1,size(ilist,1)
          iset=ilist(jset)
          write(6,'(a)',advance='no')' '//trim(i2s(iset))
        enddo ! iset
        write(6,'()')
        write(6,'()')
        ! Clean up.
        deallocate(ilist,lmask)

      case('set')

        ! Set variables.
        select case(trim(field(2,command)))

        case('X','Y')
          ! Set scale transformation.
          do itransf=0,NTRANSF
            select case(trim(field(2,command)))
            case('X')
              if(trim(field(3,command))==trim(TRANSF_FORMX(itransf)))exit
            case('Y')
              if(trim(field(3,command))==trim(TRANSF_FORMY(itransf)))exit
            end select
          enddo ! itransf
          if(itransf>NTRANSF)then
            call msg('Unknown value "'//trim(field(3,command))//'" for &
               &variable "'//trim(field(2,command))//'".')
            cycle user_loop
          endif
          if(nfield(command)>3)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(4,command))/='for')then
              call msg('Syntax error in set command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(5,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
          else ! 3 fields
            ! All-set setting.
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
          endif ! more than 3 fields or not
          ! Check transformation is applicable.
          do i=1,size(ilist,1)
            iset=ilist(i)
            if(iset<1.or.iset>ndataset)then
              call msg('Set index out of range.')
              deallocate(ilist)
              cycle user_loop
            endif
            xy=>dlist(iset)%dataset%xy
            if(TRANSF_REQ_NONZERO(itransf))then
              select case(trim(field(2,command)))
              case('X')
                if(any(eq_dble(xy%x,0.d0)))then
                  call msg('Cannot apply axis transformation: set #'//&
                     &trim(i2s(iset))//' contains x=0.')
                  deallocate(ilist)
                  cycle user_loop
                endif
              case('Y')
                if(any(eq_dble(xy%y,0.d0)))then
                  call msg('Cannot apply axis transformation: set #'//&
                     &trim(i2s(iset))//' contains y=0.')
                  deallocate(ilist)
                  cycle user_loop
                endif
              end select
            endif
            if(TRANSF_REQ_POSITIVE(itransf))then
              select case(trim(field(2,command)))
              case('X')
                if(any(xy%x<0.d0))then
                  call msg('Cannot apply axis transformation: set #'//&
                     &trim(i2s(iset))//' contains x<0.')
                  deallocate(ilist)
                  cycle user_loop
                endif
              case('Y')
                if(any(xy%y<0.d0))then
                  call msg('Cannot apply axis transformation: set #'//&
                     &trim(i2s(iset))//' contains y<0.')
                  deallocate(ilist)
                  cycle user_loop
                endif
              end select
            endif
          enddo ! i
          ! Update default itransf if no sets specified.
          if(nfield(command)>3)then
            select case(trim(field(2,command)))
            case('X')
              itransfx_default=itransf
              write(6,'(a)')'Set default X='//trim(TRANSF_FORMX(itransf))//&
                 &' for new sets.'
            case('Y')
              itransfy_default=itransf
              write(6,'(a)')'Set default Y='//trim(TRANSF_FORMY(itransf))//&
                 &' for new sets.'
            end select
          endif
          ! Store transformation.
          do i=1,size(ilist,1)
            iset=ilist(i)
            dataset=>dlist(iset)%dataset
            select case(trim(field(2,command)))
            case('X')
              dataset%itransfx=itransf
              write(6,'(a)')'Set X='//trim(TRANSF_FORMX(itransf))//&
                 &' for set #'//trim(i2s(iset))//'.'
            case('Y')
              dataset%itransfy=itransf
              write(6,'(a)')'Set Y='//trim(TRANSF_FORMY(itransf))//&
                 &' for set #'//trim(i2s(iset))//'.'
            end select
          enddo ! i
          write(6,'()')
          deallocate(ilist)
          ! Apply transformations.
          do iset=1,ndataset
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! iset
          ! Update X0.
          call refresh_X0(ndataset,dlist,glob)

        case('wexp')
          ! Set weight exponent.
          wexp=parse_dble(field(3,command),ierr)
          if(ierr/=0)then
            call msg('Problem parsing value of weight exponent.')
            cycle user_loop
          endif
          if(nfield(command)>3)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(4,command))/='for')then
              call msg('Syntax error in set command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)<5)then
              call msg('Syntax error in set command: "for" subcommand &
                 &requires an argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(5,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              iset=ilist(i)
              if(iset<1.or.iset>ndataset)then
                call msg('Set index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! 3 fields
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
          endif ! more than 3 fields or not
          ! Check exponent is applicable.
          do i=1,size(ilist,1)
            iset=ilist(i)
            if(lt_dble(wexp,0.d0))then
              xy=>dlist(iset)%dataset%xy
              if(any(eq_dble(xy%w,0.d0)))then
                call msg('Cannot apply weight exponent: set #'//&
                   &trim(i2s(iset))//' contains w=0.')
                deallocate(ilist)
                cycle user_loop
              endif
            endif
          enddo ! i
          ! Update default wexp if no sets specified.
          if(nfield(command)>3)then
            wexp_default=wexp
            write(6,'(a)')'Set default wexp for new sets.'
          endif
          ! Set exponent.
          do i=1,size(ilist,1)
            iset=ilist(i)
            dataset=>dlist(iset)%dataset
            dataset%wexp=wexp
            write(6,'(a)')'Set wexp for set #'//trim(i2s(iset))//'.'
          enddo ! i
          write(6,'()')

        case('fit')
          ! Set fit exponents.
          ! Figure out syntax used.
          nullify(ilist)
          select case(nfield(command))
          case(2)
            call msg('Syntax error: no value given for "fit".')
            cycle user_loop
          case(3)
            continue
          case(5)
            if(trim(field(4,command))/='for')then
              call msg('Syntax error: only "for <set-list>" allowed to &
                 &follow fit form.')
              cycle user_loop
            endif
            call parse_ilist(field(5,command),ilist)
            if(.not.associated(ilist))then
              call msg('Syntax error in <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              if(ilist(i)<1.or.ilist(i)>ndataset)then
                call msg('One or more dataset indices are out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          case default
            call msg('Syntax error: either three or five fields expected in &
               &"set fit".')
            cycle user_loop
          end select
          ! Parse exponents.
          select case(field(3,command))
          case('constant')
            allocate(rlist(1))
            rlist=(/(dble(i),i=0,0)/)
          case('linear')
            allocate(rlist(2))
            rlist=(/(dble(i),i=0,1)/)
          case('quadratic')
            allocate(rlist(3))
            rlist=(/(dble(i),i=0,2)/)
          case('cubic')
            allocate(rlist(4))
            rlist=(/(dble(i),i=0,3)/)
          case('quartic')
            allocate(rlist(5))
            rlist=(/(dble(i),i=0,4)/)
          case('quintic')
            allocate(rlist(6))
            rlist=(/(dble(i),i=0,5)/)
          case('sextic')
            allocate(rlist(7))
            rlist=(/(dble(i),i=0,6)/)
          case('septic')
            allocate(rlist(8))
            rlist=(/(dble(i),i=0,7)/)
          case default
            nullify(rlist)
            call parse_rlist(field(3,command),rlist)
            if(.not.associated(rlist))then
              if(associated(ilist))deallocate(ilist)
              call msg('Syntax error in <fit-form>.')
              cycle user_loop
            endif
          end select
          npoly=size(rlist,1)
          ! Decide what to do with new fit form.
          if(associated(ilist))then
            ! Define additional fit and point at it from provided datasets.
            ! See if any fits become unused.
            allocate(lmask(nfit+1))
            lmask=.false.
            do iset=1,ndataset
              if(any(ilist==iset))then
                ifit=nfit+1
              else
                ifit=dlist(iset)%dataset%ifit
              endif
              lmask(ifit)=.true.
            enddo ! iset
            ! Report removals in multiple-fit cases.
            if(any(.not.lmask))then
              write(6,'(a)',advance='no')'Note: removing fits:'
              do ifit=1,nfit
                if(lmask(ifit))cycle
                write(6,'(a)',advance='no')' '//trim(i2s(ifit))
              enddo ! ifit
              write(6,'()')
            endif
            ! Reindex fits, both pointing at new fit and removing unused fits.
            do iset=1,ndataset
              if(any(ilist==iset))then
                ifit=nfit+1
              else
                ifit=dlist(iset)%dataset%ifit
              endif
              dlist(iset)%dataset%ifit=count(lmask(1:ifit))
            enddo ! iset
            ! Remove unused fit forms and add new fit.
            allocate(tmp_flist(count(lmask)))
            do ifit=1,nfit
              if(.not.lmask(ifit))cycle
              tmp_flist(count(lmask(1:ifit)))%fit=>flist(ifit)%fit
            enddo
            deallocate(flist)
            nfit=count(lmask)
            allocate(tmp_flist(nfit)%fit)
            tmp_flist(nfit)%fit%npoly=npoly
            allocate(tmp_flist(nfit)%fit%pow(npoly),&
               &tmp_flist(nfit)%fit%share(npoly))
            tmp_flist(nfit)%fit%pow(1:npoly)=rlist(1:npoly)
            tmp_flist(nfit)%fit%share(1:npoly)=.false.
            allocate(flist(nfit))
            do ifit=1,nfit
              flist(ifit)%fit=>tmp_flist(ifit)%fit
            enddo
            deallocate(tmp_flist)
            ! Clean up.
            deallocate(lmask,rlist,ilist)
          else
            ! Replace all fits with this one and point at it form all datasets.
            ! Report removals in multiple-fit cases.
            if(nfit>1)then
              write(6,'(a)')'Note: replacing all '//trim(i2s(nfit))//' defined &
                 &fit forms with provided form.'
              if(ndataset>0)write(6,'(a)')'      All datasets use fit form #1.'
              write(6,'()')
            endif
            ! Replace all defined fit forms with this one.
            call kill_flist(flist)
            nfit=1
            allocate(flist(nfit))
            allocate(flist(nfit)%fit)
            flist(nfit)%fit%npoly=npoly
            allocate(flist(nfit)%fit%pow(npoly),flist(nfit)%fit%share(npoly))
            flist(nfit)%fit%pow(1:npoly)=rlist(1:npoly)
            flist(nfit)%fit%share(1:npoly)=.false.
            ! Make sure all datasets now point at this fit.
            do iset=1,ndataset
              dlist(iset)%dataset%ifit=1
            enddo ! iset
            ! Clean up.
            deallocate(rlist)
          endif
          ! Report newly added fit.
          write(6,'(a)')'Fit #'//trim(i2s(nfit))//' set to:'
          write(6,'(2x,a)')trim(print_poly_sym(flist(nfit)%fit,glob%X0))
          write(6,'(2x,a)')'Shared coefficients reset to: none'
          write(6,'()')

        case('range')
          ! Set fit range.
          ! Check sort variable.
          call parse_range(trim(field(3,command)),drange)
          if(drange%op=='')then
            call msg('Syntax error parsing range string.')
            cycle user_loop
          endif
          if(drange%no_rhs)then
            call msg('Syntax error parsing right-hand side of range.')
            drange%op=''
            drange%no_rhs=.false.
            cycle user_loop
          endif
          ! Apply transformations.
          do iset=1,ndataset
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! iset
          ! Update X0.
          call refresh_X0(ndataset,dlist,glob)
          ! Report.
          call msg('Range set.')

        case('shared')
          if(nfield(command)<3)then
            call msg('Syntax error: no value given for "shared".')
            cycle user_loop
          endif
          ! Get maximum fit size.
          max_npoly=0
          do ifit=1,nfit
            max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
          enddo ! ifit
          allocate(lmask(max_npoly))
          lmask=.false.
          ! See which coefficients need to be shared.
          if(trim(field(3,command))=='all')then
            lmask=.true.
          elseif(trim(field(3,command))=='none')then
            continue
          else
            nullify(ilist)
            call parse_ilist(field(3,command),ilist)
            if(.not.associated(ilist))then
              call msg('Syntax error: could not parse <coeff-list>.')
              deallocate(lmask)
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              if(ilist(i)<1.or.ilist(i)>max_npoly)then
                call msg('Coefficient index out of range.')
                deallocate(lmask,ilist)
                cycle user_loop
              endif
              lmask(ilist(i))=.true.
            enddo ! i
            deallocate(ilist)
          endif ! all/none/<list>
          ! See which fits this applies to.
          if(nfield(command)==3)then
            ! Not specified, so build corresponding ilist.
            allocate(ilist(nfit))
            ilist=(/(i,i=1,nfit)/)
          elseif(nfield(command)==5)then
            ! Specified, so parse ilist and check.
            if(trim(field(4,command))/='in')then
              call msg('Syntax error: only "in <fit-list>" can be provided &
                 &after <coeff-list>.')
              deallocate(lmask)
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(5,command),ilist)
            if(.not.associated(ilist))then
              call msg('Syntax error: could not parse <fit-list>.')
              deallocate(lmask)
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              if(ilist(i)<1.or.ilist(i)>nfit)then
                call msg('Fit index out of range.')
                deallocate(lmask,ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! nfield is neither 3 or 5
            call msg('Syntax error: either 3 or 5 fields expected.')
            deallocate(lmask)
            cycle user_loop
          endif
          ! Apply and report.
          do i=1,size(ilist,1)
            ifit=ilist(i)
            flist(ifit)%fit%share=lmask(1:flist(ifit)%fit%npoly)
            write(6,'(a)',advance='no')'Shared coefficients in fit #'//&
               &trim(i2s(ifit))//' set to:'
            if(.not.any(flist(ifit)%fit%share))then
              write(6,'(a)')' none'
            else
              do j=1,flist(ifit)%fit%npoly
                if(flist(ifit)%fit%share(j))write(6,'(a)',advance='no')' k'//&
                   &trim(i2s(j))
              enddo ! j
              write(6,'()')
            endif
          enddo ! i
          write(6,'()')
          ! Clean up.
          deallocate(ilist,lmask)

        case('centre')
          ! Check value.
          select case(field(3,command))
          case('left','right','max','min','centre','mean','median')
            continue
          case default
            t1=dble_field(3,command,ierr)
            if(ierr/=0)then
              call msg('Invalid value "'//trim(field(3,command))//'" for &
                 &variable "'//trim(field(2,command))//'".')
              cycle user_loop
            endif
          end select
          ! Set value.
          glob%X0_string=field(3,command)
          ! Update X0.
          call refresh_X0(ndataset,dlist,glob)
          call msg('centre set to '//trim(field(3,command))//'.')

        case('nsample')
          ! Check value.
          i=int_field(3,command,ierr)
          if(ierr/=0)then
            call msg('Invalid value "'//trim(field(3,command))//'" for &
               &variable nsample.')
            cycle user_loop
          endif
          if(i<10)then
            call msg('nsample value too small.')
            cycle user_loop
          endif
          glob%nsample=i
          call msg('nsample set to '//trim(i2s(i))//'.')

        case('random')
          ! Random number seed.
          select case(trim(field(3,command)))
          case('timer')
            i=0
          case default
            i=int_field(3,command,ierr)
            if(ierr/=0)then
              call msg('Invalid value "'//trim(field(3,command))//'" for &
                 &variable "random".  Must be "timer" or an integer.')
              cycle user_loop
            endif
          end select
          call init_random(i)
          call msg('Set random seed to '//trim(field(3,command))//'.')

        case('qrandom')
          ! Quasi-random noise handling.
          if(nfield(command)>2)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(3,command))/='for')then
              call msg('Syntax error in set command: unkwown subcommand &
                 &"'//trim(field(3,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)<4)then
              call msg('Syntax error in set command: "for" subcommand &
                 &requires an argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(4,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
          else ! 2 fields
            ! All-set setting.
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Set default.
            dataset_model%apply_qrandom=.true.
            write(6,'(a,es12.4,a)')'Enabled quasi-random error handling for &
               &new sets.'
          endif ! more than 2 fields or not
          do i=1,size(ilist,1)
            iset=ilist(i)
            dlist(iset)%dataset%apply_qrandom=.true.
            write(6,'(a)')'Enabled quasi-random error handling for set '//&
               &trim(i2s(iset))//'.'
          enddo ! i
          write(6,'()')
          deallocate(ilist)

        case('qrandom_exp')
          ! Quasi-random noise handling (model exponent).
          t1=dble_field(3,command,ierr)
          if(ierr/=0)then
            call msg('Invalid value "'//trim(field(3,command))//'" for &
               &variable qrandom_exp.')
            cycle user_loop
          endif
          if(nfield(command)>3)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(4,command))/='for')then
              call msg('Syntax error in set command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)<5)then
              call msg('Syntax error in set command: "for" subcommand &
                 &requires an argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(5,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
          else ! 3 fields
            ! All-set setting.
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Set default too.
            dataset_model%qrandom_exp=t1
            write(6,'(a,es12.4,a)')'Set qrandom model exponent to ',t1,&
               &' for new sets.'
          endif ! more than 3 fields or not
          do i=1,size(ilist,1)
            iset=ilist(i)
            dlist(iset)%dataset%qrandom_exp=t1
            write(6,'(a,es12.4,a)')'Set qrandom model exponent to ',t1,&
               &' for set '//trim(i2s(iset))//'.'
          enddo ! i
          write(6,'()')
          deallocate(ilist)

        case('qrandom_centre')
          ! Quasi-random noise handling (model centre).
          t1=dble_field(3,command,ierr)
          if(ierr/=0)then
            call msg('Invalid value "'//trim(field(3,command))//'" for &
               &variable qrandom_centre.')
            cycle user_loop
          endif
          if(nfield(command)>3)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(4,command))/='for')then
              call msg('Syntax error in set command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)<5)then
              call msg('Syntax error in set command: "for" subcommand &
                 &requires an argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(5,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
          else ! 3 fields
            ! All-set setting.
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Set default too.
            dataset_model%qrandom_centre=t1
            write(6,'(a,es12.4,a)')'Set qrandom model centre to ',t1,&
               &' for new sets.'
          endif ! more than 3 fields or not
          do i=1,size(ilist,1)
            iset=ilist(i)
            dlist(iset)%dataset%qrandom_centre=t1
            write(6,'(a,es12.4,a)')'Set qrandom model centre to ',t1,&
               &' for set '//trim(i2s(iset))//'.'
          enddo ! i
          write(6,'()')
          deallocate(ilist)

        case('echo')
          input_echo=.true.
          call msg('Enabled input echo.')

        case default
          call msg('Unknown variable "'//trim(field(2,command))//'".')
          cycle user_loop
        end select ! variable to set

      case('unset')
        ! Set variables.
        select case(trim(field(2,command)))
        case('X','Y')
          if(nfield(command)>2)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(3,command))/='for')then
              call msg('Syntax error in unset command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)/=4)then
              call msg('Syntax error in unset command: "for" subcommand &
                 &requires one argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(4,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              iset=ilist(i)
              if(iset<1.or.iset>ndataset)then
                call msg('Set index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! 3 fields
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Reset default while at it.
            select case(trim(field(2,command)))
            case('X')
              itransfx_default=ITRANSF_NONE
              call msg('Reset default X='//trim(TRANSF_FORMX(ITRANSF_NONE))//&
                 &' for new sets.')
            case('Y')
              itransfy_default=ITRANSF_NONE
              call msg('Reset default Y='//trim(TRANSF_FORMY(ITRANSF_NONE))//&
                 &' for new sets.')
            end select
          endif ! more than 2 fields or not
          ! Apply.
          do i=1,size(ilist,1)
            iset=ilist(i)
            select case(trim(field(2,command)))
            case('X')
              dlist(iset)%dataset%itransfx=ITRANSF_NONE
              call msg('Reset X='//trim(TRANSF_FORMX(ITRANSF_NONE))//&
                 &' for set #'//trim(i2s(iset))//'.')
            case('Y')
              dlist(iset)%dataset%itransfy=ITRANSF_NONE
              call msg('Reset Y='//trim(TRANSF_FORMY(ITRANSF_NONE))//&
                 &' for set #'//trim(i2s(iset))//'.')
            end select
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! i
          deallocate(ilist)
          call refresh_X0(ndataset,dlist,glob)
        case('wexp')
          if(nfield(command)>2)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(3,command))/='for')then
              call msg('Syntax error in unset command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)/=4)then
              call msg('Syntax error in unset command: "for" subcommand &
                 &requires one argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(4,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              iset=ilist(i)
              if(iset<1.or.iset>ndataset)then
                call msg('Set index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! 3 fields
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Reset default while at it.
            wexp_default=1.d0
          endif ! more than 2 fields or not
          ! Apply.
          do i=1,size(ilist,1)
            iset=ilist(i)
            dlist(iset)%dataset%wexp=1.d0
            call msg('Unset wexp for set #'//trim(i2s(iset))//'.')
          enddo ! ifield
          deallocate(ilist)
        case('fit')
          ! Report if we are 'undefining' existing fit forms.
          if(nfit>1)then
            write(6,'(a)')'Note: replacing all '//trim(i2s(nfit))//' defined &
               &fit forms with default.'
            if(ndataset>0)write(6,'(a)')'      All datasets use fit form #1.'
            write(6,'()')
          endif
          call kill_flist(flist)
          nfit=1
          allocate(flist(nfit))
          flist(nfit)%fit%npoly=2
          allocate(flist(nfit)%fit%pow(2),flist(nfit)%fit%share(2))
          flist(nfit)%fit%pow=(/0.d0,1.d0/)
          flist(nfit)%fit%share=.false.
          ! Make sure all datasets now point at this fit.
          do iset=1,ndataset
            dlist(iset)%dataset%ifit=1
          enddo ! iset
          ! Report.
          write(6,'(a)')'Fit #1 reset to:'
          write(6,'(2x,a)')trim(print_poly_sym(flist(nfit)%fit,glob%X0))
          write(6,'(2x,a)')'Shared coefficients reset to: none'
          write(6,'()')
        case('range')
          drange%var='X'
          drange%op=''
          drange%thres_op='non'
          drange%thres=0.d0
          drange%size=0
          drange%no_rhs=.false.
          do iset=1,ndataset
            call refresh_dataset(dlist(iset)%dataset,drange)
          enddo ! iset
          call refresh_X0(ndataset,dlist,glob)
          ! Report.
          call msg('Range unset.')
        case('shared')
          ! See which fits this applies to.
          if(nfield(command)==2)then
            ! Not specified, so build corresponding ilist.
            allocate(ilist(nfit))
            ilist=(/(i,i=1,nfit)/)
          elseif(nfield(command)==4)then
            ! Specified, so parse ilist and check.
            if(trim(field(3,command))/='in')then
              call msg('Syntax error: only "in <fit-list>" can be provided &
                 &after "unset shared".')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(4,command),ilist)
            if(.not.associated(ilist))then
              call msg('Syntax error: could not parse <fit-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              if(ilist(i)<1.or.ilist(i)>nfit)then
                call msg('Fit index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! nfield is neither 2 or 4
            call msg('Syntax error: either 2 or 4 fields expected.')
            cycle user_loop
          endif
          ! Apply and report.
          do i=1,size(ilist,1)
            ifit=ilist(i)
            flist(ifit)%fit%share=.false.
            write(6,'(a)')'Shared coefficients in fit #'//&
               &trim(i2s(ifit))//' set to: none'
          enddo ! i
          write(6,'()')
          ! Clean up.
          deallocate(ilist)
        case('centre')
          glob%X0_string=glob_default%X0_string
          call refresh_X0(ndataset,dlist,glob)
          call msg('centre reset to 0.')
        case('nsample')
          glob%nsample=glob_default%nsample
          call msg('nsample reset to '//&
             &trim(i2s(glob_default%nsample))//'.')
        case('random')
          call init_random(1)
          call msg('Reset random seed to 1.')
        case('qrandom')
          if(nfield(command)>2)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(3,command))/='for')then
              call msg('Syntax error in unset command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)/=4)then
              call msg('Syntax error in unset command: "for" subcommand &
                 &requires one argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(4,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              iset=ilist(i)
              if(iset<1.or.iset>ndataset)then
                call msg('Set index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! 3 fields
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Reset default while at it.
            dataset_model%apply_qrandom=dataset_default%apply_qrandom
            write(6,'(a)')'Disabled quasi-random error handling for new sets.'
          endif ! more than 2 fields or not
          ! Apply.
          do i=1,size(ilist,1)
            iset=ilist(i)
            dlist(iset)%dataset%apply_qrandom=.false.
            write(6,'(a)')'Disabled quasi-random error handling for set #'//&
               &trim(i2s(iset))//'.'
          enddo ! ifield
          write(6,'()')
          deallocate(ilist)
        case('qrandom_exp')
          if(nfield(command)>2)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(3,command))/='for')then
              call msg('Syntax error in unset command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)/=4)then
              call msg('Syntax error in unset command: "for" subcommand &
                 &requires one argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(4,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              iset=ilist(i)
              if(iset<1.or.iset>ndataset)then
                call msg('Set index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! 3 fields
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Reset default while at it.
            dataset_model%qrandom_exp=dataset_default%qrandom_exp
            write(6,'(a)')'Reset qrandom_exp for new sets.'
          endif ! more than 2 fields or not
          ! Apply.
          do i=1,size(ilist,1)
            iset=ilist(i)
            dlist(iset)%dataset%qrandom_exp=0.d0
            write(6,'(a)')'Reset qrandom_exp for set #'//trim(i2s(iset))//'.'
          enddo ! ifield
          deallocate(ilist)
        case('qrandom_centre')
          if(nfield(command)>2)then
            ! Set-by-set setting.  Check syntax.
            if(trim(field(3,command))/='for')then
              call msg('Syntax error in unset command: unkwown subcommand &
                 &"'//trim(field(4,command))//'".')
              cycle user_loop
            endif
            if(nfield(command)/=4)then
              call msg('Syntax error in unset command: "for" subcommand &
                 &requires one argument.')
              cycle user_loop
            endif
            nullify(ilist)
            call parse_ilist(field(4,command),ilist)
            if(.not.associated(ilist))then
              call msg('Could not parse <set-list>.')
              cycle user_loop
            endif
            do i=1,size(ilist,1)
              iset=ilist(i)
              if(iset<1.or.iset>ndataset)then
                call msg('Set index out of range.')
                deallocate(ilist)
                cycle user_loop
              endif
            enddo ! i
          else ! 3 fields
            allocate(ilist(ndataset))
            ilist=(/(i,i=1,ndataset)/)
            ! Reset default while at it.
            dataset_model%qrandom_centre=dataset_default%qrandom_centre
            write(6,'(a)')'Reset qrandom_centre for new sets.'
          endif ! more than 2 fields or not
          ! Apply.
          do i=1,size(ilist,1)
            iset=ilist(i)
            dlist(iset)%dataset%qrandom_centre=0.d0
            write(6,'(a)')'Reset qrandom_centre for set #'//trim(i2s(iset))//'.'
          enddo ! ifield
          deallocate(ilist)
        case('echo')
          input_echo=.false.
          call msg('Disabled input echo.')
        case default
          call msg('Unknown variable "'//trim(field(2,command))//'".')
          cycle user_loop
        end select

      case('status')
        ! Report status.
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
          select case(drange%thres_op)
          case('min','max')
            write(6,'(a,es11.4)')'  Data range: '//trim(drange%var)//' '//&
               &trim(drange%op)//' ',trim(drange%thres_op)//&
               &trim(i2s(drange%size))
          case default
            write(6,'(a,es11.4)')'  Data range: '//trim(drange%var)//' '//&
               &trim(drange%op)//' ',drange%thres
          end select
        end select
        write(6,'(a)')'  Number of Monte Carlo samples: '//&
           &trim(i2s(glob%nsample))
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
             &'X='//trim(TRANSF_FORMX(dataset%itransfx))//', '//&
             &'Y='//trim(TRANSF_FORMY(dataset%itransfy))
          write(6,'(a)',advance='no')', fit #'//trim(i2s(dataset%ifit))
          write(6,'(a)')'.'
          write(6,'(a,es10.4,a,es10.4)')'  wexp=',dataset%wexp,', sw=',&
             &dataset%weight
          if(dlist(iset)%dataset%apply_qrandom)then
            write(6,'(a,es12.4,a,es12.4)')'  Using quasi-random noise &
               &handling with exp=',dlist(iset)%dataset%qrandom_exp,&
               &', centre=',dlist(iset)%dataset%qrandom_centre
          else
            write(6,'(a)')'  Quasi-random noise handling disabled.'
          endif
        enddo ! i
        write(6,'()')
        do ifit=1,nfit
          write(6,'(a)')'Fit form #'//trim(i2s(ifit))//':'
          write(6,'(2x,a)')trim(print_poly_sym(flist(ifit)%fit,glob%X0))
          write(6,'(a)',advance='no')'  Shared coefficients:'
          if(.not.any(flist(ifit)%fit%share))then
            write(6,'(a)')' none'
          else
            do i=1,flist(ifit)%fit%npoly
              if(flist(ifit)%fit%share(i))write(6,'(a)',advance='no')' k'//&
                 &trim(i2s(i))
            enddo ! i
            write(6,'()')
          endif
        enddo ! ifit
        write(6,'("  X0 = ",a," =",es12.4)')trim(glob%X0_string),glob%X0
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
          call pprint('* load <file> [type <type> [using <column-list>]] &
             &[where $<column><eq><value>] [by $<column>]',0,2)
          call pprint('* wload <file> [using <column>] [where $<column><eq>&
             &<value>]',0,2)
          call pprint('* unload <set-index>',0,2)
          call pprint('* fit',0,2)
          call pprint('* plot <function> [at <xvalues>] [wrt <X>] &
             &[to <file-name>] [one-sigma|two-sigma]',0,2)
          call pprint('* assess <variables> [using <function> at X <X> &
             &[for <set-list>]]',0,2)
          call pprint('* report <report>',0,2)
          call pprint('* evaluate <function> at X <X>',0,2)
          call pprint('* find <function> <value> [between <X1> <X2>] &
             &[rightof <X1>] [leftof <X2>] [near <Xmid>]',0,2)
          call pprint('* intersect [mix] [between <X1> <X2>] [rightof <X1>] &
             &[leftof <X2>] [near <Xmid>] [plot [at <xvalues>] &
             &[to <file-name>]]',0,2)
          call pprint('* probability <function> [<condition>] [at <xvalues>] &
             &[wrt <offset>]',0,2)
          call pprint('* set <variable> <value> [for <set-list>]',0,2)
          call pprint('* use fit <fit-index> for <set-list>',0,2)
          call pprint('* unset <variable>',0,2)
          call pprint('* status',0,2)
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
          call pprint('Command: load <file> [type <type> &
             &[using <column-list>]] [where $<column><eq><value>] &
             &[by $<column>]',0,9)
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
          call pprint('Using "<<EOF" in place of <file> will cause standard &
             &input to be read until a line containing string "EOF" is found &
             &("EOF" can be replaced with another string; matching is case-&
             &sensitive).  This requires <type> to be explicitly specified in &
             &the "load" command.',2,2)
          call pprint('')
          call pprint('The "where" subcommand restricts file parsing to lines &
             &where column <column> is either equal (when <eq> is ==) or &
             &unequal (when <eq> is !=) to <value>.  <value> can be a &
             &string.  If multiple "where" subcommands are specified, only &
             &lines for which ALL specified columns take the required values &
             &are loadedi (i.e., logical AND operations combine multiple &
             &"where" clauses).',2,2)
          call pprint('')
          call pprint('The "by" subcommand allows loading multiple datasets &
             &at once, each corresponding to a different value of column &
             &<column>.  If multiple "by" commands are specified, two lines &
             &belong to different datasets if ANY of the specified columns &
             &differs in value.',2,2)
          call pprint('')
          call pprint('Example:',2,2)
          call pprint('')
          call pprint('load "../data.dat" type ywdy using 3,5,4 &
             &where $2!=bad where $7==1/4 by $1',4,6)
          call pprint('')
          call pprint('This will read file "../data.dat", loading y, dy, and &
             &w from colums 3, 4, and 5, setting x to the data-point index. &
             &Lines whose column 2 contains "bad" or whose column 7 does not &
             &contain the value 0.25 are skipped, and the loaded data are &
             &split into individual datasets for each distinct value of &
             &column 1.  Note that data-point indices for each dataset run &
             &independently, i.e., the first value of x in each dataset will &
             &be 1.',2,2)
          call pprint('')

        case('wload')
          call pprint('')
          call pprint('Command: wload <file> [using <column>] &
             &[where $<column><eq><value>]',0,9)
          call pprint('')
          call pprint('Loads global dataset weights from column <column> &
             &(column 1 by default) of <file>.  These weights are applied to &
             &all data in a dataset during fitting, and are used in the "sum" &
             &operation as linear coefficients.  E.g., "evaluate sumf at X=0" &
             &on two weighted datasets will compute w1*f1(0) + w2*f2(0).',2,2)
          call pprint('')
          call pprint('Using "<<EOF" in place of <file> will cause standard &
             &input to be read until a line containing string "EOF" is found &
             &("EOF" can be replaced with another string; matching is case-&
             &sensitive).',2,2)
          call pprint('')
          call pprint('The "where" subcommand restricts file parsing to lines &
             &where column <column> is either equal (when <eq> is ==) or &
             &unequal (when <eq> is !=) to <value>.  <value> can be a &
             &string.  If multiple "where" subcommands are specified, only &
             &lines for which ALL specified columns take the required values &
             &are loadedi (i.e., logical AND operations combine multiple &
             &"where" clauses).',2,2)
          call pprint('')

        case('unload')
          call pprint('')
          call pprint('Command: unload [<set1> [<set2> [...]]]',0,9)
          call pprint('')
          call pprint('Unload datasets from memory.  If no sets are &
             &specified, all datasets are unloaded.',2,2)
          call pprint('')

        case('fit')
          call pprint('')
          call pprint('Command: fit',0,9)
          call pprint('')
          call pprint('Perform fit of currently loaded datasets.',2,2)
          call pprint('')

        case('plot')
          call pprint('')
          call pprint('Command: plot <function> [at <xvalues>] [wrt <X>] &
             &[to <filename>] [one-sigma|two-sigma]',0,9)
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
          call pprint('The "wrt" option causes <function> to be evaluated &
             &relative to its value at <X>.  This is useful for assessing the &
             &statistics of the difference between values of the function at &
             &different X; note that it is the uncertainty on this difference &
             &that is written to the plot file, so, e.g., the plot value at X &
             &will always have zero uncertainty.  With the "wrt" option, if &
             &<function> is f or sharedf the original data are written with &
             &their y values are offset by the *mean* value of <function> at &
             &X.',2,2)
          call pprint('')
          call pprint('If any fit parameters are shared among datasets, the &
             &fit function is split into a shared part and a set-specific &
             &part -- data points are offset by the value of the set-specific &
             &part, and the shared part of the fit is plotted.',2,2)
          call pprint('')
          call pprint('By default, data are plotted as "x f df", where f is &
             &the mean value of the function and df its standard error.  If &
             &"one-sigma" or "two-sigma" are specified, data are instead &
             &written as "x f_lo f f_hi", where f_lo and f_hi are the lower &
             &and upper bounds of the requested (one-sigma or two-sigma) &
             &confidence interval on the value of f.',2,2)
          call pprint('')

        case('assess')
          call pprint('')
          call pprint('Command: assess <variables> [by <criterion>] [using &
             &<function> at <x-values> [for <set-list>]]',0,9)
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
          call pprint('* <minN',2,4)
          call pprint('* <=minN',2,4)
          call pprint('* >=maxN',2,4)
          call pprint('* >maxN',2,4)
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

        case('evaluate')
          call pprint('')
          call pprint('Command: evaluate <function> at <xvalues> [wrt <X>]',&
             &0,9)
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
          call pprint('The "wrt" option causes <function> to be evaluated &
             &relative to its value at <X>.  This is useful for assessing the &
             &statistics of the difference between values of the function at &
             &different X; note that it is the uncertainty on this difference &
             &that is reported, so, e.g., the value reported at X will always &
             &have zero uncertainty.',2,2)
          call pprint('')

        case('find')
          call pprint('')
          call pprint('Command: find <function> <value> [between <X1> <X2>] &
             &[rightof <X1>] [leftof <X2>] [near <Xmid>] [wrt <X0>]',0,9)
          call pprint('')
          call pprint('Report the location X at which <function> takes the &
             &value <value>.',2,2)
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
          call pprint('If <function> is a linear polynomial, limits X1 and/or &
             &X2 can be provided optionally, which will cause an error to be &
             &flagged if <function> is not <value> in the range [X1,X2].  &
             &Xmid is ignored in this case.',2,2)
          call pprint('')
          call pprint('If <function> is a quadratic polynomial, at least one &
             &of X1, X2, and Xmid must be provided, which selects which of &
             &the (potentially two) locations to choose.',2,2)
          call pprint('')
          call pprint('For all other fit forms, X1 and X2 must be provided, &
             &and Xmid is ignored.  <function> minus <value> is required to &
             &be of opposite signs at X1 and X2, and is assumed to change &
             &sign only once in the interval.  An error will be flagged if &
             &this is not consistently the case during random sampling.',2,2)
          call pprint('')
          call pprint('It is possible to specify "rightof data" or "leftof &
             &data", which respectively set X1 to the maximum value of X or &
             &X2 to the minimum value of X found in the datasets.  Also &
             &accepted is "between data", which sets X1 and X2 to the minimum &
             &and maximum value of X found in the datasets, respectively.',2,2)
          call pprint('')
          call pprint('The implementation tolerates up to 1% of the random &
             &resample to yield no (or out-of-range) locations.  The &
             &fraction of failures is reported as "missfrac".  The failures &
             &are simply ignored in computing the results.',2,2)
          call pprint('')
          call pprint('The "wrt" option causes the value of Y at X to be &
             &reported relative to its value at the specified position X0.',&
             &2,2)
          call pprint('')

        case('intersect')
          call pprint('')
          call pprint('Command: intersect [mix] [between <X1> <X2>] &
             &[rightof <X1>] [leftof <X2>] [near <Xmid>] [plot [at <xvalues>] &
             &[to <file-name>]] [dump [to <file-name>]]',0,9)
          call pprint('')
          call pprint('Evaluate the average location of the intersections &
             &between each pair of datasets.',2,2)
          call pprint('')
          call pprint('The intersect command finds the intersection between &
             &the fits to each pair of datasets and computes the average &
             &location of the intersection.  This operation requires two or &
             &more datasets to be loaded.',2,2)
          call pprint('')
          call pprint('For linear fits, limits X1 and/or X2 can be provided &
             &optionally, which will cause an error to be flagged if the &
             &intersection is not in the range [X1,X2].  Xmid is ignored &
             &in this case.',2,2)
          call pprint('')
          call pprint('For quadratic fits, at least one of X1, X2, and Xmid &
             &must be provided, which selects which of the (potentially two) &
             &intersections between each pair of curves to choose.',2,2)
          call pprint('')
          call pprint('For all other fit forms, X1 and X2 must be provided, &
             &and Xmid is ignored.  The difference between each pair of fits &
             &is required to be of opposite signs at X1 and X2, and is &
             &assumed to change sign only once in the interval.  An error &
             &will be flagged if this is not consistently the case during &
             &random sampling.',2,2)
          call pprint('')
          call pprint('It is possible to specify "rightof data" or "leftof &
             &data", which respectively set X1 to the maximum value of X or &
             &X2 to the minimum value of X found in the datasets.  Also &
             &accepted is "between data", which sets X1 and X2 to the minimum &
             &and maximum value of X found in the datasets, respectively.',2,2)
          call pprint('')
          call pprint('The implementation tolerates up to 1% of the random &
             &resample to yield no (or out-of-range) intersections.  The &
             &fraction of intersection failures is reported as "missfrac".  &
             &The failures are simply ignored in computing the results.',2,2)
          call pprint('')
          call pprint('Specifying "mix" triggers the use of an &
             &alternative approach in which two linear combinations of each &
             &pair of datasets are constructed, and the uncertainty in the &
             &intersection abcissa is minimized with respect to the two &
             &parameters BETA1 and BETA2 that define the linear &
             &combinations.',2,2)
          call pprint('')
          call pprint('The "mix" method requires datasets to have the same &
             &number of data and the same x values, and that these x values &
             &have no uncertainty.  This method is advantageous in the &
             &presence of quasirandom fluctuations which are correlated &
             &across different datasets; in this case it is advisable to set &
             &"qrandom" in addition to using the "mix" intersection method.',&
             &2,2)
          call pprint('')
          call pprint('The "plot" subcommand causes the "mix" datasets to &
             &be plotted.  See "help plot" for details on the "at" and "to" &
             &subcommands.',2,2)
          call pprint('')
          call pprint('The "dump" subcommand causes a scatter plot of the &
             &intersections found during random sampling to be dumped to &
             &<file-name> ("intersection.dump" by default), with plots for &
             &different dataset-pair intersections separated by two blank &
             &lines.',2,2)
          call pprint('')

        case('probability')
          call pprint('')
          call pprint('Command: probability <function> <condition> &
             &at <xvalues> [wrt <X>]',0,9)
          call pprint('')
          call pprint('Evaluate the probability that <function> satisfies &
             &<condition> (relative to <X>) at <xvalues>.',2,2)
          call pprint('')
          call pprint('<function> can be f, f'', or f'''' for the value, &
             &first, and second derivative, respectively, or sharedf, &
             &sharedf'', or sharedf'''' for the value, first, and second &
             &derivative of the shared part of the fit function (for &
             &multi-dataset fits containing shared parameters).',2,2)
          call pprint('')
          call pprint('<condition> can be either "positive" or "negative".',&
             &2,2)
          call pprint('')
          call pprint('<xvalues> is specified as &
             &"<variable>=<comma-separated-list>" (e.g., "X=0,1,2,3") or as &
             &"<variable>=<first>:<last>:<count>" (e.g., "X=0:3:4").',2,2)
          call pprint('')
          call pprint('The "wrt" option causes <function> to be evaluated &
             &relative to its value at <X>.  This is useful for assessing &
             &whether the function is ever above or below the value at <X>.',&
             &2,2)
          call pprint('')
          call pprint('If any fit parameters are shared among datasets, the &
             &fit function is split into a shared part and a set-specific &
             &part -- data points are offset by the value of the set-specific &
             &part, and the shared part of the fit is assessed.',2,2)
          call pprint('')

        case('use')
          call pprint('')
          call pprint('Command: use fit <fit-index> for <set-list>',0,9)
          call pprint('')
          call pprint('Use fit <fit-index> for sets listed in <set-list>.')
          call pprint('')

        case('set')
          if(nfield(command)==2)then
            call pprint('')
            call pprint('Command: set <variable> <value> &
               &[for <set-list>|in <fit-list>]',0,9)
            call pprint('')
            call pprint('Sets <variable> to <value>, either globally or for &
               &selected sets or fits (for per-set and per-fit variables, &
               &respectively).  The list of available variables is:',2,2)
            call pprint('')
            call pprint('* X',2,4)
            call pprint('* Y',2,4)
            call pprint('* wexp',2,4)
            call pprint('* fit',2,4)
            call pprint('* range',2,4)
            call pprint('* shared',2,4)
            call pprint('* centre',2,4)
            call pprint('* nsample',2,4)
            call pprint('* random',2,4)
            call pprint('* qrandom',2,4)
            call pprint('* qrandom_exp',2,4)
            call pprint('* qrandom_centre',2,4)
            call pprint('* echo',2,4)
            call pprint('')
            call pprint('Type "help set <variable>" for detailed &
               &information.',2,2)
            call pprint('')
          else
            select case(field(3,command))
            case('X','Y')
              call pprint('')
              call pprint('Variable: X, Y',0,10)
              call pprint('')
              call pprint('"X" and "Y" set scale transformations for the &
                 &independent and dependent variables, respectively.  &
                 &In POLYFIT notation, the original variables are called x &
                 &and y, and the transformed variables are called X and Y.',&
                 &2,2)
              call pprint('')
              call pprint('The allowed values for "X" are:',2,2)
              call pprint('* x',2,4)
              call pprint('* 1/x     [x==0 forbidden]',2,4)
              call pprint('* log(x)  [x<=0 forbidden]',2,4)
              call pprint('* exp(x)',2,4)
              call pprint('and similarly for "Y":',2,2)
              call pprint('* y',2,4)
              call pprint('* 1/y     [y==0 forbidden]',2,4)
              call pprint('* log(y)  [y<=0 forbidden]',2,4)
              call pprint('* exp(y)',2,4)
              call pprint('')
              call pprint('These variables can be set in a per-set manner or &
                 &globally.  Note that the global value applies to all loaded &
                 &datasets and becomes the default for new datasets.  POLYFIT &
                 &will refuse to apply a transformation to datasets &
                 &containing incompatible data, e.g., "set X log(x) &
                 &for 1" is not allowed if dataset #1 contains negative x &
                 &values.',2,2)
              call pprint('')
              call pprint('The default values of "X" and "Y" are "x" and "y", &
                 &respectively.',2,2)
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
                 &function.  The value can be specified as a comma-separated &
                 &list of tokens, each of which can either be a number or a &
                 &range, or as one of "constant", "linear", "quadratic", &
                 &..., "septic".  That is, "set fit 0,1,2:3" and "set fit &
                 &cubic" both generate the set of exponents {0,1,2,3}.  &
                 &(Exponents are not constrained to being integers; any real-&
                 &valued exponent can be used.)',2,2)
              call pprint('')
              call pprint('If specified globally (without <set-list>), any &
                 &previously defined fits associated with any loaded datasets &
                 &are replaced with the provided fit.  If "for <set-list>" is &
                 &specified, a new fit form is added and associated with the &
                 &provided datasets.  Note that POLYFIT automatically deletes &
                 &fit forms not associated with any datasets (except when no &
                 &datasets are loaded).',2,2)
              call pprint('')
              call pprint('See also the related "use fit" command which &
                 &modifies the dataset-fit association after having defined &
                 &the fit form..',2,2)
              call pprint('')
              call pprint('The default value of "fit" is "0,1", &
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
              call pprint('- <minN',4,6)
              call pprint('- <=minN',4,6)
              call pprint('- >=maxN',4,6)
              call pprint('- >maxN',4,6)
              call pprint('')
              call pprint('Note that in the above, T is a real-valued &
                 &threshold, and N is an integer.',2,4)
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
              call pprint('"shared" specifies which coefficients in each fit &
                 &are shared among datasets.  The value is a list of integer &
                 &coefficient indices (e.g. "set shared 2,3:4 in 2"); the &
                 &special values "none" and "all" are also allowed.  &
                 &Coefficients not flagged as shared take independent values &
                 &for each dataset that uses the same fit form.',2,2)
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
                 &dynamically changes when data are loaded/unloaded, when &
                 &x/y scales and ranges are redefined, etc.',2,2)
              call pprint('')
              call pprint('For technical reasons, "centre" is a global &
                 &variable (cannot be set on a per-set or per-fit basis).  &
                 &The default value of "centre" is "0".',2,2)
              call pprint('')
            case('nsample')
              call pprint('')
              call pprint('Variable: nsample',0,10)
              call pprint('')
              call pprint('"nsample" is a global integer variable which &
                 &determines the number of Monte Carlo samples to use in the &
                 &evaluation of uncertainties. Note that the uncertainty ddY &
                 &in an uncertainty dY is ddY = dY/sqrt(2*nsample).',2,2)
              call pprint('')
              call pprint('The default value of "nsample" is 5000, which &
                 &yields an uncertainty in the estimated uncertainty of 1% of &
                 &its value.',2,2)
              call pprint('')
            case('random')
              call pprint('')
              call pprint('Variable: random',0,10)
              call pprint('')
              call pprint('"random" sets the seed for the random number &
                 &generator.  Allowed values are "timer" (sets a seed based &
                 &on the current time and yields a non-reproducible random &
                 &number sequence) and any integer.  The default random seed &
                 &is 1.',2,2)
              call pprint('')
            case('qrandom')
              call pprint('')
              call pprint('Variable: qrandom',0,10)
              call pprint('')
              call pprint('"qrandom" is a Boolean variable which activates &
                 &handling of quasirandom noise in the data.  The reported &
                 &uncertainty of results (fit values, intersections, etc) &
                 &includes a contribution from the bias due to quasirandom &
                 &fluctuations of the data in addition to the &
                 &purely statistical uncertainty.',2,2)
              call pprint('')
              call pprint('Note that quasirandom fluctuations in X are not &
                 &accounted for by this facility.',2,2)
              call pprint('')
              call pprint('"qrandom" can be set on a per-set basis.',2,2)
              call pprint('')
              call pprint('See "help set qrandom_exp" for details of the &
                 &form used to model quasirandom noise.',2,2)
              call pprint('')
            case('qrandom_exp','qrandom_centre')
              call pprint('')
              call pprint('Variables: qrandom_exp / qrandom_centre',0,10)
              call pprint('')
              call pprint('Quasirandom noise in the data is modelled by &
                 &sigma_Y = <alpha> * |X-<qrandom_centre>|^<qrandom_exp>.  &
                 &Both qrandom_centre and qrandom_exponent default to zero. &
                 &Coefficient alpha is obtained by a preliminary fit to the &
                 &data.  Once alpha is found, sigma_Y is added in quadrature &
                 &to the original value of dY in the dataset, thus including &
                 &quasirandom noise as a contribution to the uncertainty.',2,2)
              call pprint('')
              call pprint('"qrandom_exp" and "qrandom_centre" can be set on a &
                 &per-set basis.',2,2)
              call pprint('')
            case('echo')
              call pprint('Variable: echo',0,10)
              call pprint('')
              call pprint('Setting "echo" causes input commands to be echoed &
                 &back to stdout, which is useful for scripts.')
              call pprint('')
            case default
              call pprint('No help for variable "'//trim(field(3,command))//&
                 &'".')
            end select
          endif

        case('unset')
          call pprint('')
          call pprint('Command: unset <variable>',0,9)
          call pprint('')
          call pprint('Sets <variable> to its default value.',2,2)
          call pprint('')
          call pprint('Type "help set <variable>" for detailed information on &
             &variables and their default values.',2,2)
          call pprint('')

        case('status')
          call pprint('')
          call pprint('Command: status',0,9)
          call pprint('')
          call pprint('Report currently loaded datasets and values of &
             &internal variables.',2,2)
          call pprint('')

        case default
          call pprint('No help for command "'//trim(field(2,command))//'".')
        end select

      case('quit','exit')
        call quit()

      case('')
        continue

      case default
        call msg('Command "'//trim(field(1,command))//'" not recognized.')
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
    LOGICAL mask(dataset%xy%nxy)

    ! Transform.
    xy=>dataset%xy
    call kill_xy(dataset%txy)
    call clone_xy(xy,txy)
    dataset%txy=>txy
    call scale_transform(xy%nxy,dataset%itransfx,xy%x,txy%x,xy%have_dx,xy%dx,&
       &txy%dx)
    call scale_transform(xy%nxy,dataset%itransfy,xy%y,txy%y,xy%have_dy,xy%dy,&
       &txy%dy)

    ! Get mask.
    call get_range_mask(drange,xy,txy,mask)

    ! Apply weight exponent.
    if(neq_dble(dataset%wexp,1.d0))then
      if(eq_dble(dataset%wexp,0.d0))then
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


  SUBROUTINE get_range_mask(drange,xy,txy,mask)
    !--------------------------------------!
    ! Get the range-restriction mask which !
    ! transforms TXY into RTXY.            !
    !--------------------------------------!
    IMPLICIT NONE
    TYPE(range_type), INTENT(in) :: drange
    TYPE(xy_type), POINTER :: xy,txy
    LOGICAL, INTENT(inout) :: mask(xy%nxy)
    INTEGER indx(xy%nxy),n,i1,i2
    DOUBLE PRECISION,POINTER :: sortvec(:)

    ! Initialize. default
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

    ! Get threshold.
    select case(trim(drange%thres_op))
    case('max','min')
      ! Sort data.
      call isort_dble(xy%nxy,sortvec,indx)
      ! Get range.
      n=drange%size
      select case(drange%op(1:1))
      case('<')
        i1=1
        select case(trim(drange%thres_op))
        case('min')
          i2=n-1
        case('max')
          i2=xy%nxy-n
        end select
        if(drange%op(2:2)=='=')i2=i2+1
      case('>')
        i2=xy%nxy
        select case(trim(drange%thres_op))
        case('min')
          i1=n+1
        case('max')
          i1=xy%nxy-n+2
        end select
        if(drange%op(2:2)=='=')i1=i1-1
      end select
      ! Get range within limits.
      i1=max(min(i1,xy%nxy),1)
      i2=max(min(i2,xy%nxy),1)
      ! Make mask.
      if(i1<=i2)mask(indx(i1:i2))=.true.
    case default
      ! Act on sort operation.
      select case(trim(drange%op))
      case('<')
        mask=lt_dble(sortvec,drange%thres)
      case('<=')
        mask=le_dble(sortvec,drange%thres)
      case('>')
        mask=gt_dble(sortvec,drange%thres)
      case('>=')
        mask=ge_dble(sortvec,drange%thres)
      end select
    end select

    ! Clean up (unnecessary).
    nullify(sortvec)

  END SUBROUTINE get_range_mask


  SUBROUTINE refresh_X0(ndataset,dlist,glob)
    !-----------------------------------------!
    ! Refresh value of X0 in following change !
    ! to X0_string, range or loaded data.     !
    !-----------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset
    TYPE(dataset_list_type),INTENT(in) :: dlist(ndataset)
    TYPE(global_params_type),INTENT(inout) :: glob
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    INTEGER iset,tot_nxy,ixy,ierr
    DOUBLE PRECISION t1,tx0
    DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:)

    select case(trim(glob%X0_string))

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
        dataset=>dlist(iset)%dataset
        xy=>dataset%rtxy
        if(xy%nxy==0)cycle
        x(tot_nxy+1:tot_nxy+xy%nxy)=xy%x(1:xy%nxy)
        y(tot_nxy+1:tot_nxy+xy%nxy)=xy%y(1:xy%nxy)
        tot_nxy=tot_nxy+xy%nxy
      enddo ! iset

      ! Locate X0.
      select case(trim(glob%X0_string))
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
      t1=parse_dble(glob%X0_string,ierr)
      tx0=t1
    end select

    ! Update X0.
    glob%X0=tx0

  END SUBROUTINE refresh_X0


  SUBROUTINE qrandom_apply(ndataset,dlist,drange,nfit,flist,glob_X0,report)
    !---------------------------------------------!
    ! Apply a quasirandom noise correction to the !
    ! standard error dy on all datasets.          !
    !---------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_list_type),INTENT(in) :: flist(nfit)
    DOUBLE PRECISION,INTENT(in) :: glob_X0
    LOGICAL,INTENT(in),OPTIONAL :: report
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    ! Local variables.
    LOGICAL any_apply_qrandom
    INTEGER max_npoly,ifit,iset,ixy,ipoly,ierr
    DOUBLE PRECISION sexp,sx0,alpha2,dy2,e_fit,chi2,t0,t1,t2,&
       &alpha(ndataset),alpha_prime(ndataset),dy2_save(ndataset)
    DOUBLE PRECISION,ALLOCATABLE :: a(:,:)

    ! See if we have anything to do here.
    any_apply_qrandom=.false.
    do iset=1,ndataset
      if(dlist(iset)%dataset%apply_qrandom)any_apply_qrandom=.true.
    enddo ! iset
    if(.not.any_apply_qrandom)return

    ! Allocate parameter vector.
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit
    allocate(a(max_npoly,ndataset))
    a=0.d0

    ! Initialize.
    alpha=0.d0

    ! Get un-resampled fit.
    call perform_multifit(ndataset,dlist,glob_X0,nfit,flist,chi2,a,ierr)
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      ifit=dataset%ifit
      if(dataset%rtxy%nxy<=flist(ifit)%fit%npoly)cycle
      sexp=dlist(iset)%dataset%qrandom_exp
      sx0=dlist(iset)%dataset%qrandom_centre
      if(.not.dataset%rtxy%have_dy)then
        ! Toggle have_dy on dataset.
        dataset%rtxy%have_dy=.true.
        dataset%rtxy%dy=0.d0
        dataset%txy%have_dy=.true.
        dataset%txy%dy=0.d0
        dataset%xy%have_dy=.true.
        dataset%xy%dy=0.d0
      endif
      ! Get alpha^2.
      alpha2=0.d0
      dy2=0.d0
      do ixy=1,dataset%rtxy%nxy
        e_fit=0.d0
        do ipoly=1,flist(ifit)%fit%npoly
          e_fit=e_fit+a(ipoly,iset)*&
             &(dataset%rtxy%x(ixy)-glob_X0)**flist(ifit)%fit%pow(ipoly)
        enddo ! ipoly
        t0=dataset%rtxy%y(ixy)-e_fit
        t1=t0**2
        t2=dataset%rtxy%dy(ixy)**2
        if(neq_dble(sexp,0.d0))then
          t1=t1/abs(dataset%rtxy%x(ixy)-sx0)**(2*sexp)
          t2=t2/abs(dataset%rtxy%x(ixy)-sx0)**(2*sexp)
        endif
        alpha2=alpha2+t1
        dy2=dy2+t2
      enddo ! ixy
      alpha2=alpha2/dble(dataset%rtxy%nxy-flist(ifit)%fit%npoly)
      ! alpha_prime is alpha pre dy correction.
      dy2_save(iset)=sqrt(dy2/dble(dataset%rtxy%nxy))
      alpha_prime(iset)=sqrt(alpha2)
      alpha2=alpha2-dy2/dble(dataset%rtxy%nxy)
      if(le_dble(alpha2,0.d0))alpha2=0.d0
      ! Store locally for reporting.
      alpha(iset)=sqrt(alpha2)
      ! Adjust stderrs.
      if(eq_dble(sexp,0.d0))then
        dataset%rtxy%dy=sqrt(dataset%rtxy%dy**2+alpha2)
      else
        dataset%rtxy%dy=sqrt(dataset%rtxy%dy**2+&
           &alpha2*abs(dataset%rtxy%x-sx0)**(2*sexp))
      endif
      ! Unrestrict range to get txy.
      call back_transform(drange,dataset)
    enddo ! iset

    ! Report.
    if(present(report))then
      if(report)then
        write(6,'(a)')'Quasi-random noise'
        write(6,'(a)')'------------------'
        write(6,'(2x,a5,1x,3(1x,a20))')'Set','alpha       ','stddev       ',&
           &'dy         '
        do iset=1,ndataset
          write(6,'(2x,i5,1x,3(1x,es20.12))')iset,alpha(iset),&
             &alpha_prime(iset),dy2_save(iset)
        enddo ! iset
        write(6,'()')
      endif
    endif

  END SUBROUTINE qrandom_apply


  SUBROUTINE show_multipoly(ndataset,dlist,drange,nfit,flist,glob)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(global_params_type),INTENT(in) :: glob
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    ! Local variables.
    TYPE(eval_type) deval ! (dummy arg.)
    LOGICAL fit_used(nfit)
    INTEGER tot_nparam,tot_nxy,iset,ifit,i,ierr,max_npoly
    DOUBLE PRECISION chi2,chi2err,rmsy,rmsyerr
    DOUBLE PRECISION,ALLOCATABLE :: a(:,:),da(:,:)

    ! Allocate parameter vector.
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit
    allocate(a(max_npoly,ndataset),da(max_npoly,ndataset))
    a=0.d0
    da=0.d0

    ! Initialize.
    tot_nparam=0
    tot_nxy=0
    fit_used=.false.
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      tot_nxy=tot_nxy+dataset%rtxy%nxy
      ifit=dataset%ifit
      tot_nparam=tot_nparam+count(.not.flist(ifit)%fit%share)
      if(.not.fit_used(ifit))tot_nparam=tot_nparam+count(flist(ifit)%fit%share)
      fit_used(ifit)=.true.
    enddo ! iset

    ! Perform fit.
    call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,glob,&
       &deval,ierr,chi2mean=chi2,chi2err=chi2err,amean=a,aerr=da,&
       &rmsymean=rmsy,rmsyerr=rmsyerr)
    if(ierr/=0)then
      call msg('Could not perform fit.')
      return
    endif

    ! Print table header.
    write(6,'(a)')'Fit parameters:'
    write(6,'()')
    write(6,'(3x)',advance='no')
    write(6,'(2x,a5,2x,a5,1x,2(1x,a20))')'Set','i','ki       ','dki       '
    write(6,'(3x)',advance='no')
    write(6,'(2x,'//trim(i2s(55))//'("-"))')
    ! Print table.
    do iset=1,ndataset
      ifit=dlist(iset)%dataset%ifit
      do i=1,flist(ifit)%fit%npoly
        write(6,'(a)',advance='no')'FIT'
        write(6,'(2x,i5,2x,i5,1x,2(1x,es20.12))')iset,i,a(i,iset),da(i,iset)
      enddo ! i
    enddo ! iset
    ! Print table footer.
    write(6,'(3x)',advance='no')
    write(6,'(2x,'//trim(i2s(55))//'("-"))')
    write(6,'()')

    ! Report chi-squared.
    write(6,'(a)')'Fit assessment:'
    write(6,'()')
    write(6,'(3x)',advance='no')
    write(6,'(2x,a12,1x,2(1x,a20))')'Measure  ','Value       ','Stderr      '
    write(6,'(3x)',advance='no')
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
    write(6,'(3x)',advance='no')
    write(6,'(2x,'//trim(i2s(55))//'("-"))')
    write(6,'()')

    ! Write out fit in XMGRACE format.
    ! NB, xmgrace has a string length limit of 256 characters, so this may not
    ! be very useful in some cases.
    write(6,'(a)')'Fit in XMGRACE format:'
    do iset=1,ndataset
      ifit=dlist(iset)%dataset%ifit
      write(6,'(a)')'  Set #'//trim(i2s(iset))//': '//&
         &trim(print_poly_num(flist(ifit)%fit,glob%X0,a(:,iset)))
    enddo ! iset
    write(6,'()')

  END SUBROUTINE show_multipoly


  SUBROUTINE evaluate_fit(ndataset,dlist,nfit,flist,glob,drange,deval)
    !-------------------------------------------------!
    ! Evaluate value or derivative of fit at provided !
    ! points.                                         !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_list_type),INTENT(in) :: flist(:)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(range_type),INTENT(in) :: drange
    TYPE(eval_type),INTENT(in) :: deval
    ! Local variables.
    INTEGER iset,ix,ierr
    DOUBLE PRECISION,ALLOCATABLE :: fmean(:,:),ferr(:,:)

    ! Evaluate.
    allocate(fmean(deval%n,ndataset),ferr(deval%n,ndataset))
    call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
       &glob,deval,ierr,fmean,ferr)
    if(ierr/=0)then
      call msg('Could not perform fit.')
      return
    endif

    ! Report.
    write(6,'(4x)',advance='no')
    write(6,'(2x,a4,1x,3(1x,a20))')'set','X       ','f       ','df       '
    write(6,'(4x)',advance='no')
    write(6,'(2x,68("-"))')
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
    write(6,'(2x,68("-"))')
    write(6,'()')
    deallocate(fmean,ferr)

  END SUBROUTINE evaluate_fit


  SUBROUTINE find_fit_value(ndataset,dlist,nfit,flist,glob,drange,&
    &intersect_range,deval,find_target,use_Xrel,Xrel)
    !----------------------------------------------!
    ! Find the location of the requested value and !
    ! report to stdout.                            !
    !----------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(range_type),INTENT(in) :: drange
    TYPE(intersect_range_type),INTENT(in) :: intersect_range
    TYPE(eval_type),INTENT(in) :: deval
    DOUBLE PRECISION,INTENT(in) :: find_target
    LOGICAL,INTENT(in) :: use_Xrel
    DOUBLE PRECISION,INTENT(in) :: Xrel
    ! Monte Carlo sample storage.
    DOUBLE PRECISION,ALLOCATABLE :: w_vector(:)
    DOUBLE PRECISION,ALLOCATABLE :: x0_array(:,:),y0_array(:,:)
    INTEGER,ALLOCATABLE :: err_array(:,:)
    ! Combined fit.
    INTEGER,ALLOCATABLE :: map1(:,:),map2(:,:)
    TYPE(fit_form_list_type),POINTER :: tmp_flist(:)
    ! Parameter vector.
    INTEGER max_npoly
    DOUBLE PRECISION,ALLOCATABLE :: a(:,:),a_diff(:),a_target(:)
    ! Pointers.
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy,xy_orig
    ! Local variables.
    TYPE(fit_form_type),POINTER :: fit_form_target,tmp_fit_form
    LOGICAL report_relative_y,is_consecutive,all_are_poly1,all_are_poly2
    INTEGER ifit,i,iset,irandom,nsample,ierr,ideriv
    DOUBLE PRECISION x0,y0,dx0,dy0,errfrac,chi2
    INTEGER op_npoly
    DOUBLE PRECISION,ALLOCATABLE :: op_a(:),op_pow(:)

    ! Build fit corresponding to target value.
    allocate(fit_form_target)
    allocate(fit_form_target%pow(1),fit_form_target%share(1),a_target(1))
    fit_form_target%npoly=1
    fit_form_target%pow(1)=0.d0
    fit_form_target%share(1)=.false.
    a_target(1)=find_target

    ! Allocate op_* vectors, tmp_fit_form, and parameter vectors.
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit
    allocate(op_a(max_npoly),op_pow(max_npoly))
    op_a=0.d0
    op_pow=0.d0
    allocate(tmp_fit_form)
    allocate(tmp_fit_form%pow(max_npoly),tmp_fit_form%share(max_npoly))
    tmp_fit_form%npoly=0
    tmp_fit_form%pow=0.d0
    tmp_fit_form%share=.false.
    allocate(a(max_npoly,ndataset),a_diff(1+max_npoly))
    a=0.d0
    a_diff=0.d0

    ! Range checks.
    all_are_poly1=.true.
    all_are_poly2=.true.
    do ifit=1,nfit
      op_npoly=flist(ifit)%fit%npoly
      op_pow(1:op_npoly)=flist(ifit)%fit%pow(1:op_npoly)
      do ideriv=1,deval%nderiv
        call deriv_poly(op_npoly,op_pow,op_a)
      enddo ! ideriv
      is_consecutive=all(eq_dble(op_pow(1:op_npoly),&
         &(/(dble(i-1),i=1,op_npoly)/)))
      all_are_poly1=all_are_poly1.and.is_consecutive.and.&
         &op_npoly==2.and.neq_dble(op_pow(op_npoly),0.d0)
      all_are_poly2=all_are_poly2.and.is_consecutive.and.&
         &op_npoly==3.and.neq_dble(op_pow(op_npoly),0.d0)
    enddo
    if(all_are_poly2)then
      if(.not.intersect_range%have_x1.and..not.intersect_range%have_x2.and.&
         &.not.intersect_range%have_xmid)then
        call msg('Need at least one reference point to intersect &
           &second-order polynomials.')
        call kill_fit_form(tmp_fit_form)
        call kill_fit_form(fit_form_target)
        return
      endif
    elseif(.not.all_are_poly1)then
      if(.not.intersect_range%have_x1.or..not.intersect_range%have_x2)then
        call msg('Need explicit range to intersect generic polynomials.')
        call kill_fit_form(tmp_fit_form)
        call kill_fit_form(fit_form_target)
        return
      endif
    endif

    ! Build combination fit forms.
    allocate(tmp_flist(nfit),map1(1+max_npoly,nfit),map2(1+max_npoly,nfit))
    map1=0
    map2=0
    do ifit=1,nfit
      op_npoly=flist(ifit)%fit%npoly
      op_pow(1:op_npoly)=flist(ifit)%fit%pow(1:op_npoly)
      do ideriv=1,deval%nderiv
        call deriv_poly(op_npoly,op_pow,op_a)
      enddo ! ideriv
      tmp_fit_form%npoly=op_npoly
      tmp_fit_form%pow(1:op_npoly)=op_pow(1:op_npoly)
      call combine_fit_form(tmp_fit_form,fit_form_target,&
         &tmp_flist(ifit)%fit,map1(1,ifit),map2(1,ifit))
    enddo ! ifit

    ! Make copy of datasets and apply qrandom.
    call clone_dlist(dlist,tmp_dlist)
    call qrandom_apply(ndataset,tmp_dlist,drange,nfit,flist,glob%X0)

    ! Allocate storage for location of intersection.
    allocate(w_vector(glob%nsample))
    allocate(x0_array(glob%nsample,ndataset),y0_array(glob%nsample,ndataset),&
       &err_array(glob%nsample,ndataset))
    x0_array=0.d0
    y0_array=0.d0
    err_array=0
    w_vector=1.d0

    ! Initialize.
    nsample=glob%nsample

    ! Loop over random points.
    irandom=0
    do irandom=1,nsample
      do iset=1,ndataset
        dataset=>tmp_dlist(iset)%dataset
        xy=>dataset%xy
        xy_orig=>dlist(iset)%dataset%xy
        ! Using xy%dx and xy%dy so we get qrandom accounted for.
        if(xy%have_dx)xy%x=xy_orig%x+gaussian_random_number(xy%dx)
        if(xy%have_dy)xy%y=xy_orig%y+gaussian_random_number(xy%dy)
        call refresh_dataset(dataset,drange)
      enddo ! iset
      call perform_multifit(ndataset,tmp_dlist,glob%X0,nfit,flist,chi2,a,ierr)
      if(ierr/=0)then
        call kill_fit_form(tmp_fit_form)
        call kill_fit_form(fit_form_target)
        call kill_dlist(tmp_dlist)
        call kill_flist(tmp_flist)
        call msg('Could not perform fit.')
        return
      endif
      ! Loop over datasets.
      do iset=1,ndataset
        ifit=tmp_dlist(iset)%dataset%ifit
        op_npoly=flist(ifit)%fit%npoly
        op_pow(1:op_npoly)=flist(ifit)%fit%pow(1:op_npoly)
        op_a(1:op_npoly)=a(1:op_npoly,iset)
        do ideriv=1,deval%nderiv
          call deriv_poly(op_npoly,op_pow,op_a)
        enddo ! ideriv
        tmp_fit_form%npoly=op_npoly
        tmp_fit_form%pow=op_pow(1:op_npoly)
        call combine_fit_coeffs(tmp_fit_form,fit_form_target,&
           &tmp_flist(ifit)%fit,map1(1,ifit),map2(1,ifit),op_a,&
           &a_target,1.d0,-1.d0,a_diff)
        ! Find intersection between sets ISET and JSET.
        call intersect_zero(tmp_flist(ifit)%fit,glob%X0,a_diff,&
           &intersect_range,x0,ierr)
        y0=eval_poly(flist(ifit)%fit%npoly,flist(ifit)%fit%pow,a(1,iset),&
           &x0-glob%X0)
        if (use_Xrel) y0=y0-eval_poly(flist(ifit)%fit%npoly,&
           &flist(ifit)%fit%pow,a(1,iset),Xrel-glob%X0)
        x0_array(irandom,iset)=x0
        y0_array(irandom,iset)=y0
        err_array(irandom,iset)=ierr
      enddo ! iset
    enddo ! irandom

    ! Compute and report location of value.
    write(6,'(a)')'Location of values:'
    write(6,'(4x)',advance='no')
    write(6,'(1x,a3,3(1x,a20))')'Set','X0          ','DX0         ',&
       &'Y0          ','DY0         ','missfrac      '
    do iset=1,ndataset
      call forgiving_analysis(nsample,x0_array(1,iset),y0_array(1,iset),&
         &err_array(1,iset),x0,dx0,y0,dy0,errfrac,ierr)
      if(ierr==0)then
        write(6,'(a4)',advance='no')'FIND'
        write(6,'(1x,i3,5(1x,es20.12))')iset,x0,dx0,y0,dy0,errfrac
      else
        write(6,'(a4)',advance='no')'FIND'
        write(6,'(1x,i3,5(1x,es20.12))')iset,0.d0,0.d0,0.d0,0.d0,errfrac
      endif
    enddo ! iset
    write(6,'()')

    ! Clean up.
    call kill_fit_form(tmp_fit_form)
    call kill_fit_form(fit_form_target)
    call kill_dlist(tmp_dlist)
    call kill_flist(tmp_flist)

  END SUBROUTINE find_fit_value


  SUBROUTINE intersect_fit(ndataset,dlist,nfit,flist,glob,drange,&
     &intersect_range,dname)
    !------------------------------------------------!
    ! Find the intersection of all dataset pairs and !
    ! report to stdout.                              !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(range_type),INTENT(in) :: drange
    TYPE(intersect_range_type),INTENT(in) :: intersect_range
    CHARACTER(*),INTENT(in) :: dname
    ! Monte Carlo sample storage.
    DOUBLE PRECISION,ALLOCATABLE :: w_vector(:)
    DOUBLE PRECISION,ALLOCATABLE :: x0_array(:,:),y0_array(:,:)
    INTEGER,ALLOCATABLE :: err_array(:,:)
    ! Combined fit.
    INTEGER cfmap(nfit,nfit)
    INTEGER,ALLOCATABLE :: map1(:,:),map2(:,:)
    TYPE(fit_form_list_type),POINTER :: tmp_flist(:)
    ! Parameter vector.
    INTEGER max_npoly
    DOUBLE PRECISION,ALLOCATABLE :: a(:,:),a_diff(:)
    ! Pointers.
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy,xy_orig
    ! Local variables.
    INTEGER,PARAMETER :: io=10
    LOGICAL is_consecutive,all_are_poly1,all_are_poly2,file_exists
    INTEGER ifit,jfit,icfit,ncfit,i,iset,jset,iyy,nyy,irandom,nsample,ierr
    INTEGER ngood
    DOUBLE PRECISION x0,y0,dx0,dy0,errfrac,chi2
    DOUBLE PRECISION x0l,x0r,x02sl,x02sr,x01sl,x01sr,x0med,y0l,y0r,y02sl,&
       &y02sr,y01sl,y01sr,y0med
    DOUBLE PRECISION, ALLOCATABLE :: xpack(:),ypack(:)

    ! Range checks.
    all_are_poly1=.true.
    all_are_poly2=.true.
    do ifit=1,nfit
      is_consecutive=all(eq_dble(flist(ifit)%fit%pow,&
         &(/(dble(i-1),i=1,flist(ifit)%fit%npoly)/)))
      all_are_poly1=all_are_poly1.and.is_consecutive.and.&
         &flist(ifit)%fit%npoly==2.and.&
         &neq_dble(flist(ifit)%fit%pow(flist(ifit)%fit%npoly),0.d0)
      all_are_poly2=all_are_poly2.and.is_consecutive.and.&
         &flist(ifit)%fit%npoly==3.and.&
         &neq_dble(flist(ifit)%fit%pow(flist(ifit)%fit%npoly),0.d0)
    enddo
    if(all_are_poly2)then
      if(.not.intersect_range%have_x1.and..not.intersect_range%have_x2.and.&
         &.not.intersect_range%have_xmid)then
        call msg('Need at least one reference point to intersect &
           &second-order polynomials.')
        return
      endif
    elseif(.not.all_are_poly1)then
      if(.not.intersect_range%have_x1.or..not.intersect_range%have_x2)then
        call msg('Need explicit range to intersect generic polynomials.')
        return
      endif
    endif

    ! Compute maximum polynomial size and allocate parameter vector.
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit
    allocate(a(max_npoly,ndataset),a_diff(2*max_npoly))
    a=0.d0
    a_diff=0.d0

    ! Build combination fits.
    ! NB, i,j and j,i stored separately for simplicity.
    ncfit=nfit*nfit
    allocate(tmp_flist(ncfit),map1(2*max_npoly,ncfit),map2(2*max_npoly,ncfit))
    map1=0
    map2=0
    icfit=0
    do ifit=1,nfit
      do jfit=1,nfit
        icfit=icfit+1
        cfmap(ifit,jfit)=icfit
        call combine_fit_form(flist(ifit)%fit,flist(jfit)%fit,&
           &tmp_flist(icfit)%fit,map1(1,icfit),map2(1,icfit))
      enddo ! jfit
    enddo ! ifit

    ! Make copy of datasets and apply qrandom.
    call clone_dlist(dlist,tmp_dlist)
    call qrandom_apply(ndataset,tmp_dlist,drange,nfit,flist,glob%X0)

    ! Allocate storage for location of intersection.
    nyy=(ndataset*(ndataset-1))/2
    allocate(w_vector(glob%nsample))
    allocate(x0_array(glob%nsample,nyy),y0_array(glob%nsample,nyy),&
       &err_array(glob%nsample,nyy))
    err_array=0
    w_vector=1.d0

    ! Initialize.
    nsample=glob%nsample

    ! Loop over random points.
    irandom=0
    do irandom=1,nsample
      do iset=1,ndataset
        dataset=>tmp_dlist(iset)%dataset
        xy=>dataset%xy
        xy_orig=>dlist(iset)%dataset%xy
        ! Using xy%dx and xy%dy so we get qrandom accounted for.
        if(xy%have_dx)xy%x=xy_orig%x+gaussian_random_number(xy%dx)
        if(xy%have_dy)xy%y=xy_orig%y+gaussian_random_number(xy%dy)
        call refresh_dataset(dataset,drange)
      enddo ! iset
      call perform_multifit(ndataset,tmp_dlist,glob%X0,nfit,flist,chi2,a,ierr)
      if(ierr/=0)then
        call kill_dlist(tmp_dlist)
        call kill_flist(tmp_flist)
        call msg('Could not perform fit.')
        return
      endif
      ! Loop over pairs of datasets.
      iyy=0
      do iset=1,ndataset
        ifit=tmp_dlist(iset)%dataset%ifit
        do jset=iset+1,ndataset
          iyy=iyy+1
          jfit=tmp_dlist(jset)%dataset%ifit
          icfit=cfmap(ifit,jfit)
          call combine_fit_coeffs(flist(ifit)%fit,flist(jfit)%fit,&
             &tmp_flist(icfit)%fit,map1(1,icfit),map2(1,icfit),a(1,iset),&
             &a(1,jset),1.d0,-1.d0,a_diff)
          ! Find intersection between sets ISET and JSET.
          call intersect_zero(tmp_flist(icfit)%fit,glob%X0,a_diff,&
             &intersect_range,x0,ierr)
          select case(ierr)
          case(0)
            y0=eval_poly(flist(ifit)%fit%npoly,flist(ifit)%fit%pow,a(1,iset),&
               &x0-glob%X0)
          case default
            x0=0.d0
            y0=0.d0
          end select
          x0_array(irandom,iyy)=x0
          y0_array(irandom,iyy)=y0
          err_array(irandom,iyy)=ierr
        enddo ! jset
      enddo ! iset
    enddo ! irandom

    ! Dump scatter plot if requested.
    if(len_trim(dname)/=0)then
      ! Delete existing dump file.
      inquire(file=trim(dname),exist=file_exists)
      if(file_exists)then
        open(io,file=trim(dname),status='old')
        close(io,status='delete')
      endif
      open(unit=io,file=dname,status='new',iostat=ierr)
      ! Loop over pairs of datasets.
      iyy=0
      do iset=1,ndataset
        do jset=iset+1,ndataset
          iyy=iyy+1
          ! Dump scatter plot.
          write(io,'("# Intersection of datasets ",i3," and ",i3)')iset,jset
          do irandom=1,nsample
            write(io,'(es20.12,1x,es20.12,1x,i1)')x0_array(irandom,iyy),&
               &y0_array(irandom,iyy),err_array(irandom,iyy)
          enddo ! irandom
          write(io,'()')
          write(io,'()')
        enddo ! jset
      enddo ! iset
      ! Close dump file.
      close(io)
    endif

    ! Compute and report intersection(s).
    write(6,'(a)')'Intersections:'
    write(6,'(4x)',advance='no')
    write(6,'(1x,a7,5(1x,a20))')'Sets ','X0          ','DX0         ',&
       &'Y0          ','DY0         ','missfrac      '
    iyy=0
    do iset=1,ndataset
      do jset=iset+1,ndataset
        iyy=iyy+1
        call forgiving_analysis(nsample,x0_array(1,iyy),y0_array(1,iyy),&
           &err_array(1,iyy),x0,dx0,y0,dy0,errfrac,ierr)
        if(ierr==0)then
          write(6,'(a4)',advance='no')'INTR'
          write(6,'(2(1x,i3),5(1x,es20.12))')iset,jset,&
             &x0,dx0,y0,dy0,errfrac
        else
          write(6,'(a4)',advance='no')'INTR'
          write(6,'(2(1x,i3),5(1x,es20.12))')iset,jset,&
             &0.d0,0.d0,0.d0,0.d0,errfrac
        endif
      enddo ! jset
    enddo ! iset
    write(6,'()')

    ! Report intersection quantiles:
    write(6,'(a)')'Intersection quantiles:'
    write(6,'(4x)',advance='no')
    write(6,'(1x,a7,1x,a16,2(1x,a20))')'Sets ','Quantile    ',&
       &'X0          ','Y0          '
    iyy=0
    do iset=1,ndataset
      do jset=iset+1,ndataset
        iyy=iyy+1
        call forgiving_analysis(nsample,x0_array(1,iyy),y0_array(1,iyy),&
           &err_array(1,iyy),x0,dx0,y0,dy0,errfrac,ierr)
        if(ierr==0)then
          ngood=count(err_array(:,iyy)==0)
          allocate(xpack(ngood),ypack(ngood))
          xpack=pack(x0_array(:,iyy),err_array(:,iyy)==0)
          ypack=pack(y0_array(:,iyy),err_array(:,iyy)==0)
          x0l=find_pth_smallest(nint(dble(ngood)*0.001d0),ngood,xpack)
          x0r=find_pth_smallest(nint(dble(ngood)*0.999d0),ngood,xpack)
          x02sl=find_pth_smallest(nint(dble(ngood)*0.022750131948d0),ngood,&
             &xpack)
          x02sr=find_pth_smallest(nint(dble(ngood)*0.977249868052d0),ngood,&
             &xpack)
          x01sl=find_pth_smallest(nint(dble(ngood)*0.158655254d0),ngood,&
             &xpack)
          x01sr=find_pth_smallest(nint(dble(ngood)*0.841344746d0),ngood,&
             &xpack)
          x0med=find_pth_smallest(nint(dble(ngood)*0.5d0),ngood,xpack)
          y0l=find_pth_smallest(nint(dble(ngood)*0.001d0),ngood,ypack)
          y0r=find_pth_smallest(nint(dble(ngood)*0.999d0),ngood,ypack)
          y02sl=find_pth_smallest(nint(dble(ngood)*0.022750131948d0),ngood,&
             &ypack)
          y02sr=find_pth_smallest(nint(dble(ngood)*0.977249868052d0),ngood,&
             &ypack)
          y01sl=find_pth_smallest(nint(dble(ngood)*0.158655254d0),ngood,ypack)
          y01sr=find_pth_smallest(nint(dble(ngood)*0.841344746d0),ngood,ypack)
          y0med=find_pth_smallest(nint(dble(ngood)*0.5d0),ngood,ypack)
          deallocate(xpack,ypack)
        else
          x0l=0.d0
          x0r=0.d0
          x02sl=0.d0
          x02sr=0.d0
          x01sl=0.d0
          x01sr=0.d0
          x0med=0.d0
          y0l=0.d0
          y0r=0.d0
          y02sl=0.d0
          y02sr=0.d0
          y01sl=0.d0
          y01sr=0.d0
          y0med=0.d0
        endif
        write(6,'(a4)',advance='no')'INTQ'
        write(6,'(2(1x,i3),1x,a16,2(1x,es20.12))')iset,jset,&
           &' 0.10%          ',x0l,y0l
        write(6,'(a4)',advance='no')'INTQ'
        write(6,'(2(1x,i3),1x,a16,4(1x,es20.12))')iset,jset,&
           &' 2.27% (2-sigma)',x02sl,y02sl
        write(6,'(a4)',advance='no')'INTQ'
        write(6,'(2(1x,i3),1x,a16,4(1x,es20.12))')iset,jset,&
           &'15.87% (1-sigma)',x01sl,y01sl
        write(6,'(a4)',advance='no')'INTQ'
        write(6,'(2(1x,i3),1x,a16,4(1x,es20.12))')iset,jset,&
           &'50.00% (median) ',x0med,y0med
        write(6,'(a4)',advance='no')'INTQ'
        write(6,'(2(1x,i3),1x,a16,4(1x,es20.12))')iset,jset,&
           &'84.13% (1-sigma)',x01sr,y01sr
        write(6,'(a4)',advance='no')'INTQ'
        write(6,'(2(1x,i3),1x,a16,4(1x,es20.12))')iset,jset,&
           &'97.72% (2-sigma)',x02sr,y02sr
        write(6,'(a4)',advance='no')'INTQ'
        write(6,'(2(1x,i3),1x,a16,2(1x,es20.12))')iset,jset,&
           &'99.9%           ',x0r,y0r
      enddo ! jset
    enddo ! iset
    write(6,'()')

    ! Clean up.
    call kill_dlist(tmp_dlist)
    call kill_flist(tmp_flist)

  END SUBROUTINE intersect_fit


  SUBROUTINE intersect_mix_fit(ndataset,dlist,nfit,flist,glob,drange,&
     &intersect_range,deval,fname)
    !-------------------------------------------------------------!
    ! Find the intersection of linear combinations of all dataset !
    ! pairs, minimizing the intersection abcissa uncertainty with !
    ! respect to the linear parameters, and report to stdout.     !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(range_type),INTENT(in) :: drange
    TYPE(intersect_range_type),INTENT(in) :: intersect_range
    TYPE(eval_type),INTENT(inout) :: deval
    CHARACTER(*),INTENT(in) :: fname
    ! Combined fit.
    INTEGER cfmap(nfit,nfit)
    INTEGER,ALLOCATABLE :: map_i(:,:),map_j(:,:)
    TYPE(fit_form_list_type),POINTER :: tmp_flist(:),flist_mix(:)
    ! Local variables.
    TYPE(dataset_list_type),POINTER :: dlist_mix(:)
    LOGICAL is_consecutive,all_are_poly1,all_are_poly2,any_have_dx,any_have_dy,&
       &file_exists
    INTEGER i,ifit,jfit,icfit,ncfit,nsample,ierr,ierr1,ierr2,max_nxy,&
       &iset,jset,max_npoly
    INTEGER,PARAMETER :: io=10
    DOUBLE PRECISION beta1,beta2,x0,dx0,y0,dy0,errfrac
    DOUBLE PRECISION,ALLOCATABLE :: ranx1(:),rany1(:),ranx2(:),rany2(:)

    ! Range checks.
    all_are_poly1=.true.
    all_are_poly2=.true.
    do ifit=1,nfit
      is_consecutive=all(eq_dble(flist(ifit)%fit%pow,&
         &(/(dble(i-1),i=1,flist(ifit)%fit%npoly)/)))
      all_are_poly1=all_are_poly1.and.is_consecutive.and.&
         &flist(ifit)%fit%npoly==2.and.&
         &neq_dble(flist(ifit)%fit%pow(flist(ifit)%fit%npoly),0.d0)
      all_are_poly2=all_are_poly2.and.is_consecutive.and.&
         &flist(ifit)%fit%npoly==3.and.&
         &neq_dble(flist(ifit)%fit%pow(flist(ifit)%fit%npoly),0.d0)
    enddo
    if(all_are_poly2)then
      if(.not.intersect_range%have_x1.and..not.intersect_range%have_x2.and.&
         &.not.intersect_range%have_xmid)then
        call msg('Need at least one reference point to intersect &
           &second-order polynomials.')
        return
      endif
    elseif(.not.all_are_poly1)then
      if(.not.intersect_range%have_x1.or..not.intersect_range%have_x2)then
        call msg('Need explicit range to intersect generic polynomials.')
        return
      endif
    endif

    ! Prepare storage for random numbers.
    nsample=glob%nsample
    max_nxy=0
    any_have_dx=.false.
    any_have_dy=.false.
    do iset=1,ndataset
      max_nxy=max(max_nxy,dlist(iset)%dataset%rtxy%nxy)
      any_have_dx=any_have_dx.or.dlist(iset)%dataset%rtxy%have_dx
      any_have_dy=any_have_dy.or.dlist(iset)%dataset%rtxy%have_dy
    enddo ! iset
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit
    allocate(ranx1(max_nxy*nsample),rany1(max_nxy*nsample),&
       &ranx2(max_nxy*nsample),rany2(max_nxy*nsample))

    ! Build combination fits.
    ! NB, i,j and j,i stored separately for simplicity.
    ncfit=nfit*nfit
    allocate(tmp_flist(ncfit),map_i(2*max_npoly,ncfit),&
       &map_j(2*max_npoly,ncfit))
    map_i=0
    map_j=0
    icfit=0
    do ifit=1,nfit
      do jfit=1,nfit
        icfit=icfit+1
        cfmap(ifit,jfit)=icfit
        call combine_fit_form(flist(ifit)%fit,flist(jfit)%fit,&
           &tmp_flist(icfit)%fit,map_i(1,icfit),map_j(1,icfit))
      enddo ! jfit
    enddo ! ifit

    ! Delete existing plot file.
    if(len_trim(fname)>0)then
      inquire(file=trim(fname),exist=file_exists)
      if(file_exists)then
        open(io,file=trim(fname),status='old')
        close(io,status='delete')
      endif
    endif

    ! Loop over pairs of datasets.
    write(6,'(a)')'Intersections:'
    write(6,'(4x)',advance='no')
    write(6,'(1x,a7,7(1x,a20))')'Sets ','X0          ','DX0         ',&
       &'Y0          ','DY0         ','missfrac      ','BETA1        ',&
       &'BETA2        '
    do iset=1,ndataset
      ifit=dlist(iset)%dataset%ifit
      do jset=iset+1,ndataset
        jfit=dlist(jset)%dataset%ifit
        icfit=cfmap(ifit,jfit)
        ! Generate random numbers.
        if(any_have_dx)then
          ranx1=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
          ranx2=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
        endif
        if(any_have_dy)then
          rany1=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
          rany2=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
        endif
        ! Find best intersection.
        beta1=0.d0
        beta2=0.d0
        call intersect_mix_search(tmp_flist(icfit)%fit,flist(ifit)%fit,&
           &flist(jfit)%fit,map_i(:,icfit),map_j(:,icfit),&
           &glob,drange,intersect_range,ranx1,rany1,beta1,&
           &dlist(iset)%dataset,dlist(jset)%dataset,ranx2,rany2,beta2,dy0,ierr)
        if(ierr/=0)then
          beta1=0.d0
          beta2=0.d0
        endif
        ! Generate another set of random numbers.
        if(any_have_dx)then
          ranx1=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
          ranx2=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
        endif
        if(any_have_dy)then
          rany1=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
          rany2=gaussian_random_number((/( 1.d0, i=1,max_nxy*nsample )/))
        endif
        ! Reevaluate at beta1,beta2.
        call intersect_mix_eval(tmp_flist(icfit)%fit,flist(ifit)%fit,&
           &flist(jfit)%fit,map_i(:,icfit),map_j(:,icfit),glob,drange,&
           &intersect_range,ranx1,rany1,beta1,dlist(iset)%dataset,&
           &dlist(jset)%dataset,ranx2,rany2,beta2,x0,dx0,y0,dy0,errfrac,ierr)
        ! Make plot if requested.
        if(len_trim(fname)>0)then
          allocate(dlist_mix(2))
          call construct_beta_dataset(drange,dlist(iset)%dataset,&
             &dlist(jset)%dataset,beta1,dlist_mix(1)%dataset,ierr1)
          call construct_beta_dataset(drange,dlist(jset)%dataset,&
             &dlist(iset)%dataset,beta2,dlist_mix(2)%dataset,ierr2)
          allocate(flist_mix(1))
          flist_mix(1)%fit=>tmp_flist(icfit)%fit
          call plot_multipoly(2,dlist_mix,drange,nfit,flist,deval,&
             &glob,trim(fname),append=.true.)
          call kill_dlist(dlist_mix)
          nullify(flist_mix(1)%fit)
          deallocate(flist_mix)
        endif
        ! Report.
        if(ierr==0)then
          write(6,'(a4)',advance='no')'INTR'
          write(6,'(2(1x,i3),7(1x,es20.12))')iset,jset,&
             &x0,dx0,y0,dy0,errfrac,beta1,beta2
        else
          write(6,'(a4)',advance='no')'INTR'
          write(6,'(2(1x,i3),7(1x,es20.12))')iset,jset,&
             &0.d0,0.d0,0.d0,0.d0,errfrac,beta1,beta2
        endif
      enddo ! jset
    enddo ! iset
    write(6,'()')

    ! Report having made plot.
    if(len_trim(fname)>0)call msg('Plotted "mixed" datasets to "'//&
       &trim(fname)//'".')

    ! Clean up.
    call kill_flist(tmp_flist)

  END SUBROUTINE intersect_mix_fit


  SUBROUTINE intersect_mix_search(cfit,fit_i,fit_j,map_i,map_j,glob,drange,&
     &intersect_range,ranx1,rany1,beta1,dataset_i,dataset_j,ranx2,rany2,&
     &beta2,dy,ierr)
    !---------------------------------------------------------!
    ! Find the values of beta1 and beta2 that yield the best- !
    ! resolved intersection between yi+beta1*(yj-yi) and      !
    ! yj+beta2*(yi-yj).                                       !
    !---------------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),POINTER :: cfit,fit_i,fit_j
    INTEGER,INTENT(in) :: map_i(fit_i%npoly),map_j(fit_j%npoly)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(range_type),INTENT(in) :: drange
    TYPE(intersect_range_type),INTENT(in) :: intersect_range
    TYPE(dataset_type),POINTER :: dataset_i,dataset_j
    DOUBLE PRECISION,INTENT(in) :: ranx1(*),rany1(*),ranx2(*),rany2(*)
    DOUBLE PRECISION,INTENT(inout) :: beta1,beta2,dy
    INTEGER,INTENT(inout) :: ierr
    ! Local variables.
    INTEGER ierr1,ierr2
    DOUBLE PRECISION dy1,dy2,dy_prev
    DOUBLE PRECISION, PARAMETER :: dy_rel_tol=1.d-5
    ! Initialize.
    beta1=0.d0
    beta2=0.d0
    ierr=1
    ! Search along (beta1,0) and (0,beta2)
    call intersect_mix_search_1d(cfit,fit_i,fit_j,map_i,map_j,glob,&
       &drange,intersect_range,ranx1,rany1,beta1,dataset_i,dataset_j,&
       &ranx2,rany2,0.d0,dy1,ierr1)
    call intersect_mix_search_1d(cfit,fit_j,fit_i,map_j,map_i,glob,&
       &drange,intersect_range,ranx2,rany2,beta2,dataset_j,dataset_i,&
       &ranx1,rany1,0.d0,dy2,ierr2)
    if(ierr1/=0.and.ierr2/=0)return
    ! Alternate search along each direction.
    dy_prev=dy1
    if(dy2<dy1.or.ierr1/=0)then
      call intersect_mix_search_1d(cfit,fit_i,fit_j,map_i,map_j,glob,&
         &drange,intersect_range,ranx1,rany1,beta1,dataset_i,dataset_j,&
         &ranx2,rany2,beta2,dy,ierr)
      if(ierr/=0)return
      if(abs(dy_prev-dy)<=dy_rel_tol*dy)return
      dy_prev=dy
    endif
    do
      call intersect_mix_search_1d(cfit,fit_j,fit_i,map_j,map_i,glob,&
         &drange,intersect_range,ranx2,rany2,beta2,dataset_j,dataset_i,&
         &ranx1,rany1,beta1,dy,ierr)
      if(ierr/=0)return
      if(abs(dy_prev-dy)<=dy_rel_tol*dy)return
      dy_prev=dy
      call intersect_mix_search_1d(cfit,fit_i,fit_j,map_i,map_j,glob,&
         &drange,intersect_range,ranx1,rany1,beta1,dataset_i,dataset_j,&
         &ranx2,rany2,beta2,dy,ierr)
      if(ierr/=0)return
      if(abs(dy_prev-dy)<=dy_rel_tol*dy)return
      dy_prev=dy
    enddo
  END SUBROUTINE intersect_mix_search


  SUBROUTINE intersect_mix_search_1d(cfit,fit_i,fit_j,map_i,map_j,glob,&
     &drange,intersect_range,ranx1,rany1,beta1,dataset_i,dataset_j,ranx2,rany2,&
     &beta2,dy,ierr)
    !-------------------------------------------------------!
    ! Find the value of beta1 that yields the best-resolved !
    ! intersection between yi+beta1*(yj-yi) and             !
    ! yj+beta2*(yi-yj), where beta2 is fixed.               !
    !-------------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),POINTER :: cfit,fit_i,fit_j
    INTEGER,INTENT(in) :: map_i(fit_i%npoly),map_j(fit_j%npoly)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(range_type),INTENT(in) :: drange
    TYPE(intersect_range_type),INTENT(in) :: intersect_range
    DOUBLE PRECISION,INTENT(in) :: ranx1(*),rany1(*),ranx2(*),rany2(*)
    TYPE(dataset_type),POINTER :: dataset_i,dataset_j
    DOUBLE PRECISION,INTENT(inout) :: beta1,dy
    DOUBLE PRECISION,INTENT(in) :: beta2
    INTEGER,INTENT(inout) :: ierr
    ! Local variables.
    INTEGER,PARAMETER :: init_ngrid=23 ! 4*integer+3
    DOUBLE PRECISION,PARAMETER :: grid_base=2.d0
    LOGICAL rejected,is_left_edge,is_right_edge
    INTEGER i,n,errvec(init_ngrid),imin
    DOUBLE PRECISION bu,bv,bw,fu,fv,fw,b,f,x0,dx0,y0,bvec(init_ngrid),&
       &fvec(init_ngrid),errfrac

    ! Construct grid for initial evaluation.
    n=0
    do i=(init_ngrid-1)/2,1,-1
      n=n+1
      bvec(n)=beta1-grid_base**dble(i-(init_ngrid+1)/4)
    enddo ! i
    n=n+1
    bvec(n)=beta1
    do i=1,(init_ngrid-1)/2
      n=n+1
      bvec(n)=beta1+grid_base**dble(i-(init_ngrid+1)/4)
    enddo ! i

    ! Perform an initial evaluation at several beta1.
    do i=1,init_ngrid
      b=bvec(i)
      call intersect_mix_eval(cfit,fit_i,fit_j,map_i,map_j,glob,drange,&
         &intersect_range,ranx1,rany1,b,dataset_i,dataset_j,ranx2,rany2,&
         &beta2,x0,dx0,y0,f,errfrac,ierr)
      fvec(i)=f
      errvec(i)=ierr
    enddo ! i

    ! Locate minimum within grid and set bracket middle point.
    imin=sum(minloc(fvec,errvec==0))
    if(imin==0)then
      ierr=1
      return
    endif
    bv=bvec(imin)
    fv=fvec(imin)

    ! Determine if minimum is the left edge, else set left backet point.
    is_left_edge=imin==1
    if(imin>1)is_left_edge=count(errvec(1:imin-1)==0)==0
    if(.not.is_left_edge)then
      do i=imin-1,1,-1
        if(errvec(i)/=0)cycle
        bu=bvec(i)
        fu=fvec(i)
        exit
      enddo ! i
    endif

    ! Determine if minimum is the right edge, else set left backet point.
    is_right_edge=imin==init_ngrid
    if(imin<init_ngrid)is_right_edge=count(errvec(imin+1:init_ngrid)==0)==0
    if(.not.is_right_edge)then
      do i=imin+1,init_ngrid
        if(errvec(i)/=0)cycle
        bw=bvec(i)
        fw=fvec(i)
        exit
      enddo ! i
    endif

    ! Extend range to left if minimum is the left edge.
    if(is_left_edge)then
      b=bv
      do
        b=b*grid_base
        call intersect_mix_eval(cfit,fit_i,fit_j,map_i,map_j,glob,drange,&
           &intersect_range,ranx1,rany1,b,dataset_i,dataset_j,ranx2,rany2,&
           &beta2,x0,dx0,y0,f,errfrac,ierr)
        if(ierr==0)then
          if(f<fv)then
            is_right_edge=.false.
            bw=bv
            fw=fv
            bv=b
            fv=f
          else
            bu=bv
            fu=fv
            exit
          endif
        endif
      enddo
    endif

    ! Extend range to right if minimum is the right edge.
    if(is_right_edge)then
      b=bv
      do
        b=b*grid_base
        call intersect_mix_eval(cfit,fit_i,fit_j,map_i,map_j,glob,drange,&
           &intersect_range,ranx1,rany1,b,dataset_i,dataset_j,ranx2,rany2,&
           &beta2,x0,dx0,y0,f,errfrac,ierr)
        if(ierr==0)then
          if(f<fv)then
            bu=bv
            fu=fv
            bv=b
            fv=f
          else
            bw=bv
            fw=fv
            exit
          endif
        endif
      enddo
    endif

    ! Now zoom in on minimum.
    do
      call parabolic_min(bu,bv,bw,fu,fv,fw,b,f,rejected)
      if (rejected) exit
      if (b<=bu.or.b>=bw) exit
      call intersect_mix_eval(cfit,fit_i,fit_j,map_i,map_j,glob,drange,&
         &intersect_range,ranx1,rany1,b,dataset_i,dataset_j,ranx2,rany2,&
         &beta2,x0,dx0,y0,f,errfrac,ierr)
      if(ierr/=0)return ! FIXME - ideally never happens, but could be handled
      if (f<fv) then
        if (b<bv) then
          bw = bv
          fw = fv
        elseif (b>bv) then
          bu = bv
          fu = fv
        else
          exit
        endif
        bv = b
        fv = f
      elseif (f>fv) then
        if (b<bv) then
          bu = b
          fu = f
        elseif (b>bv) then
          bw = b
          fw = f
        else
          exit
        endif
      else
        exit
      endif
      if (bw-bu<1.d-7) exit
    enddo

    ! Return minimum.
    beta1=bv
    dy=fv

  END SUBROUTINE intersect_mix_search_1d


  SUBROUTINE intersect_mix_eval(cfit,fit_i,fit_j,map_i,map_j,glob,drange,&
     &intersect_range,ranx1,rany1,beta1,dataset_i,dataset_j,ranx2,rany2,beta2,&
     &x0,dx0,y0,dy0,errfrac,ierr)
    !---------------------------------------------------!
    ! Evaluate the intersection of yi+beta1*(yj-yi) and !
    ! yj+beta2*(yi-yj).                                 !
    !---------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),POINTER :: cfit,fit_i,fit_j
    INTEGER,INTENT(in) :: map_i(fit_i%npoly),map_j(fit_j%npoly)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(range_type),INTENT(in) :: drange
    TYPE(intersect_range_type),INTENT(in) :: intersect_range
    DOUBLE PRECISION,INTENT(in) :: ranx1(*),rany1(*),ranx2(*),rany2(*)
    TYPE(dataset_type),POINTER :: dataset_i,dataset_j
    DOUBLE PRECISION,INTENT(in) :: beta1,beta2
    DOUBLE PRECISION,INTENT(inout) :: x0,dx0,y0,dy0,errfrac
    INTEGER,INTENT(inout) :: ierr
    ! Local variables.
    INTEGER nsample,nxy,irandom,ierr1,ierr2,errvec(glob%nsample),nfit
    DOUBLE PRECISION a_z(max(fit_i%npoly,fit_j%npoly),2),chi2,&
       &w_vector(glob%nsample),x0vec(glob%nsample),&
       &y0vec(glob%nsample),a_diff(cfit%npoly,2)
    TYPE(xy_type),POINTER :: xy,xy_orig
    TYPE(dataset_type),POINTER :: dataset
    TYPE(dataset_list_type),POINTER :: dlist(:),tmp_dlist(:)
    TYPE(fit_form_list_type),POINTER :: flist(:)

    ! Initialize.
    ierr=1
    nsample=glob%nsample
    nxy=dataset_i%xy%nxy
    w_vector=1.d0
    errfrac=1.d0

    ! Transform datasets.
    allocate(dlist(2))
    call construct_beta_dataset(drange,dataset_i,dataset_j,beta1,&
       &dlist(1)%dataset,ierr1)
    call construct_beta_dataset(drange,dataset_j,dataset_i,beta2,&
       &dlist(2)%dataset,ierr2)
    if(ierr1/=0.or.ierr2/=0)then
      call kill_dlist(dlist)
      return
    endif

    ! Generate fit list.
    nfit=2
    if(associated(fit_i,fit_j))nfit=1
    allocate(flist(nfit))
    flist(1)%fit=>fit_i
    flist(nfit)%fit=>fit_j
    dlist(1)%dataset%ifit=1
    dlist(2)%dataset%ifit=nfit

    ! Make another copy and apply qrandom.
    call clone_dlist(dlist,tmp_dlist)
    call qrandom_apply(2,tmp_dlist,drange,nfit,flist,glob%X0)

    ! Loop over random points.
    do irandom=1,nsample
      ! Modify the datasets at this (ibeta,jbeta).
      dataset=>tmp_dlist(1)%dataset
      xy=>dataset%xy
      xy_orig=>dlist(1)%dataset%xy
      if(xy%have_dx)xy%x=xy_orig%x+&
         &xy%dx*ranx1((irandom-1)*nxy+1:irandom*nxy)
      if(xy%have_dy)xy%y=xy_orig%y+&
         &xy%dy*rany1((irandom-1)*nxy+1:irandom*nxy)
      call refresh_dataset(dataset,drange)
      dataset=>tmp_dlist(2)%dataset
      xy=>dataset%xy
      xy_orig=>dlist(2)%dataset%xy
      if(xy%have_dx)xy%x=xy_orig%x+&
         &xy%dx*ranx2((irandom-1)*nxy+1:irandom*nxy)
      if(xy%have_dy)xy%y=xy_orig%y+&
         &xy%dy*rany2((irandom-1)*nxy+1:irandom*nxy)
      call refresh_dataset(dataset,drange)
      call perform_multifit(2,tmp_dlist,glob%X0,nfit,flist,chi2,a_z,ierr)
      if(ierr/=0)then
        call kill_dlist(tmp_dlist)
        call kill_dlist(dlist)
        deallocate(flist)
        return
      endif
      call combine_fit_coeffs(fit_i,fit_j,cfit,map_i,map_j,a_z(1,1),a_z(1,2),&
         &1.d0-beta1,beta1,a_diff(1,1))
      call combine_fit_coeffs(fit_j,fit_i,cfit,map_j,map_i,a_z(1,2),a_z(1,1),&
         &1.d0-beta2,beta2,a_diff(1,2))
      call intersect_zero(cfit,glob%X0,a_diff(:,2)-a_diff(:,1),&
         &intersect_range,x0,ierr)
      y0=eval_poly(fit_i%npoly,fit_i%pow,a_z(1,1),x0-glob%X0)
      x0vec(irandom)=x0
      y0vec(irandom)=y0
      errvec(irandom)=ierr
    enddo ! irandom

    ! Compute output.
    call forgiving_analysis(nsample,x0vec,y0vec,errvec,x0,dx0,y0,dy0,&
       &errfrac,ierr)

    ! Clean up.
    call kill_dlist(tmp_dlist)
    call kill_dlist(dlist)
    deallocate(flist)

  END SUBROUTINE intersect_mix_eval


  SUBROUTINE intersect_zero(fit,glob_X0,a,intersect_range,x0,ierr)
    !------------------------------------------------------------------!
    ! Find zero of fit with parameters A(:).  IERR return values:      !
    ! 0: intersection found                                            !
    ! 1: one linear/one or two quadratic intersections left of range   !
    !    or no quadratic intersection in leftof search                 !
    ! 2: one linear/one or two quadratic intersections right of range  !
    !    or no quadratic intersection in rightof search                !
    ! 3: two quadratic intersections, each either side of the range    !
    ! 4: bisection error / no quadratic intersections in fully bounded !
    !    or fully unbounded search                                     !
    !------------------------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),INTENT(in) :: fit
    DOUBLE PRECISION,INTENT(in) :: glob_X0,a(fit%npoly)
    TYPE(intersect_range_type),INTENT(in) :: intersect_range
    DOUBLE PRECISION,INTENT(inout) :: x0
    INTEGER,INTENT(inout) :: ierr
    LOGICAL is_consecutive,is_poly1,is_poly2,x0a_in_interval,x0b_in_interval,&
       &lplus,rplus
    INTEGER i
    DOUBLE PRECISION t1,x0a,x0b,xmid,xl,xr,yl,yr

    ! Initialize output variables.
    ierr=4
    x0=0.d0

    ! See if we are dealing with an especially simple case.
    is_consecutive=all(eq_dble(fit%pow,(/(dble(i-1),i=1,fit%npoly)/)))
    is_poly1=is_consecutive.and.fit%npoly==2.and.&
       &neq_dble(fit%pow(fit%npoly),0.d0)
    is_poly2=is_consecutive.and.fit%npoly==3.and.&
       &neq_dble(fit%pow(fit%npoly),0.d0)

    if(is_poly1)then

      ! Linear fit.
      x0=-a(1)/a(2)+glob_X0

      ! Flag error if outside range.
      if(intersect_range%have_x1.and.x0<intersect_range%x1)then
        ierr=1
        return
      elseif(intersect_range%have_x2.and.x0>intersect_range%x2)then
        ierr=2
        return
      endif

    elseif(is_poly2)then

      ! Quadratic fit.  Get discriminant.
      t1=a(2)**2-4.d0*a(1)*a(3)
      if(eq_dble(t1,0.d0))then
        ! Discriminant is zero: one intersection.
        x0=-0.5d0*a(2)/a(3)+glob_X0
        ! Flag error if outside range.
        if(intersect_range%have_x1.and.x0<intersect_range%x1)then
          ierr=1
          return
        elseif(intersect_range%have_x2.and.x0>intersect_range%x2)then
          ierr=2
          return
        endif
      elseif(t1<0.d0)then
        ! Discriminant is negative: no intersections.
        if(intersect_range%have_x1.and..not.intersect_range%have_x2)then
          ierr=2
          return
        elseif(intersect_range%have_x2.and..not.intersect_range%have_x1)then
          ierr=1
          return
        else
          ierr=4
          return
        endif
      else
        ! Discriminant is positive: two intersections.
        x0a=0.5d0*(-a(2)+sqrt(t1))/a(3)+glob_X0
        x0b=0.5d0*(-a(2)-sqrt(t1))/a(3)+glob_X0
        x0a_in_interval=.true.
        x0b_in_interval=.true.
        if(intersect_range%have_x1)then
          x0a_in_interval=x0a_in_interval.and.ge_dble(x0a,intersect_range%x1)
          x0b_in_interval=x0b_in_interval.and.ge_dble(x0b,intersect_range%x1)
        endif
        if(intersect_range%have_x2)then
          x0a_in_interval=x0a_in_interval.and.le_dble(x0a,intersect_range%x2)
          x0b_in_interval=x0b_in_interval.and.le_dble(x0b,intersect_range%x2)
        endif
        ! Decide which intersection to report.
        if(x0a_in_interval.and.x0b_in_interval)then
          xmid=0.d0
          if(intersect_range%have_xmid)then
            xmid=intersect_range%xmid
          elseif(intersect_range%have_x1.and.intersect_range%have_x2)then
            xmid=0.5d0*(intersect_range%x1+intersect_range%x2)
          elseif(intersect_range%have_x1)then
            xmid=intersect_range%x1
          elseif(intersect_range%have_x2)then
            xmid=intersect_range%x2
          endif
          ! Pick closest to xmid.
          if(abs(x0a-xmid)<abs(x0b-xmid))then
            x0=x0a
          else
            x0=x0b
          endif
        elseif(x0a_in_interval)then
          x0=x0a
        elseif(x0b_in_interval)then
          x0=x0b
        else
          if(intersect_range%have_x1.and.x0a<intersect_range%x1.and.&
             &x0b<intersect_range%x1)then
            ierr=1
            return
          elseif(intersect_range%have_x2.and.x0a>intersect_range%x2.and.&
             &x0b>intersect_range%x2)then
            ierr=2
            return
          else
            ierr=3
            return
          endif
        endif ! x0a and/or x0b in interval
      endif ! sign of discriminant

    else

      ! Not a "simple" fit form, so use bisection.  Initialize.
      xl=intersect_range%x1
      yl=eval_poly(fit%npoly,fit%pow,a,xl-glob_X0)
      xr=intersect_range%x2
      yr=eval_poly(fit%npoly,fit%pow,a,xr-glob_X0)
      if(abs(yl)<=0.d0.or.abs(yr)<=0.d0)return
      lplus=yl>0.d0
      rplus=yr>0.d0
      if(lplus.eqv.rplus)return

      ! Loop over bisection iterations.
      do
        x0=0.5d0*(xl+xr)
        if(eval_poly(fit%npoly,fit%pow,a,x0-glob_X0)>0.d0.eqv.lplus)then
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

    endif

    ! Flag no error.
    ierr=0

  END SUBROUTINE intersect_zero


  SUBROUTINE forgiving_analysis(n,xvec,yvec,errvec,x0,dx0,y0,dy0,errfrac,ierr)
    !-------------------------------------------------------------------!
    ! Given a set of N random-sampled intersections XVEC(1:N),YVEC(1:N) !
    ! and their corresponding error signals ERRVEC(1:N), compute the    !
    ! average intersection allowing "a few" errors to occur.            !
    !-------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    DOUBLE PRECISION,INTENT(in) :: xvec(n),yvec(n)
    INTEGER,INTENT(in) :: errvec(n)
    DOUBLE PRECISION,INTENT(inout) :: x0,dx0,y0,dy0,errfrac
    INTEGER,INTENT(inout) :: ierr
    DOUBLE PRECISION,PARAMETER :: max_errfrac=1.d0
    DOUBLE PRECISION tvec(n),w_vector(n),var
    INTEGER nerr,nn
    ! Initialize.
    ierr=1
    x0=0.d0
    dx0=0.d0
    y0=0.d0
    dy0=0.d0
    if(n<1)return
    ! Count number of errors.
    nerr=count(errvec/=0)
    errfrac=dble(nerr)/dble(n)
    if(ge_dble(errfrac,max_errfrac))return
    ! Analyse.
    ierr=0
    w_vector=1.d0
    nn=n-nerr
    tvec(1:nn)=pack(xvec,errvec==0)
    call characterize_dist(nn,tvec,w_vector,mean=x0,var=var)
    dx0=sqrt(var)
    tvec(1:nn)=pack(yvec,errvec==0)
    call characterize_dist(nn,tvec,w_vector,mean=y0,var=var)
    dy0=sqrt(var)
  END SUBROUTINE forgiving_analysis


  SUBROUTINE plot_multipoly(ndataset,dlist,drange,nfit,flist,deval,glob,&
     &fname,append,ci_choice)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(eval_type),INTENT(inout) :: deval
    TYPE(global_params_type),INTENT(in) :: glob
    CHARACTER(*),INTENT(in) :: fname
    LOGICAL,INTENT(in),OPTIONAL :: append
    INTEGER,INTENT(in),OPTIONAL :: ci_choice
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    INTEGER iset,ifit,i,ierr,op_npoly,max_npoly,plot_ci
    LOGICAL do_append
    DOUBLE PRECISION tx0,tx1,t1,f0
    DOUBLE PRECISION,ALLOCATABLE :: fmean(:,:),fmean_xs(:,:),ferr(:,:),a(:,:),&
       &op_a(:),op_pow(:)
    INTEGER,PARAMETER :: io=10

    ! Allocate op_* vectors.
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit
    allocate(op_a(max_npoly),op_pow(max_npoly))
    op_a=0.d0
    op_pow=0.d0

    ! Handle optional arguments.
    do_append=.false.
    if(present(append))do_append=append
    plot_ci=0
    if(present(ci_choice))plot_ci=max(0,min(2,ci_choice))

    ! Open plot file.
    if(.not.do_append)then
      open(unit=io,file=fname,status='replace',iostat=ierr)
    else
      open(unit=io,file=fname,position='append',iostat=ierr)
    endif
    if(ierr/=0)then
      call msg('Problem opening file.')
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
       &a(max_npoly,ndataset))
    if(plot_ci>0)allocate(fmean_xs(deval%n,ndataset))
    if(trim(deval%what)/='shared')then
      select case(plot_ci)
      case(0)
        call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
           &glob,deval,ierr,fmean=fmean,ferr=ferr,amean=a,silent=.true.)
      case(1)
        call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
           &glob,deval,ierr,fmean=fmean,fmean_1s=fmean_xs,ferr_1s=ferr,&
           &amean=a,silent=.true.)
      case(2)
        call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
           &glob,deval,ierr,fmean=fmean,fmean_2s=fmean_xs,ferr_2s=ferr,&
           &amean=a,silent=.true.)
      end select
      if(ierr/=0)then
        call msg('Could not perform fit.')
        close(io,status='delete')
        return
      endif
      if(deval%nderiv==0)then
        f0=0.d0
        ! Plot data.
        do iset=1,ndataset
          ifit=dlist(iset)%dataset%ifit
          if(deval%rel)f0=eval_poly(flist(ifit)%fit%npoly,flist(ifit)%fit%pow,&
             &a(1,iset),deval%Xrel-glob%X0)
          xy=>dlist(iset)%dataset%txy
          if(xy%have_dx.and.xy%have_dy)then
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i)-f0,xy%dx(i),xy%dy(i)
            enddo ! i
          elseif(xy%have_dx)then
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i)-f0,xy%dx(i)
            enddo ! i
          elseif(xy%have_dy)then
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i)-f0,xy%dy(i)
            enddo ! i
          else
            do i=1,xy%nxy
              write(io,*)xy%x(i),xy%y(i)-f0
            enddo ! i
          endif
          write(io,'()')
          write(io,'()')
        enddo ! iset
      endif
      ! Plot fit functions.
      do iset=1,ndataset
        do i=1,deval%n
          select case(plot_ci)
          case(0)
            write(io,*)deval%x(i),fmean(i,iset),ferr(i,iset)
          case(1)
            write(io,*)deval%x(i),fmean_xs(i,iset)-ferr(i,iset),&
               &fmean(i,iset),fmean_xs(i,iset)+ferr(i,iset)
          case(2)
            write(io,*)deval%x(i),fmean_xs(i,iset)-2.d0*ferr(i,iset),&
               &fmean(i,iset),fmean_xs(i,iset)+2.d0*ferr(i,iset)
          end select
        enddo ! i
        write(io,'()')
        write(io,'()')
      enddo ! iset
    else
      select case(plot_ci)
      case(0)
        call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
           &glob,deval,ierr,fmean=fmean,ferr=ferr,amean=a,silent=.true.)
      case(1)
        call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
           &glob,deval,ierr,fmean=fmean,fmean_1s=fmean_xs,ferr_1s=ferr,&
           &amean=a,silent=.true.)
      case(2)
        call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
           &glob,deval,ierr,fmean=fmean,fmean_2s=fmean_xs,ferr_2s=ferr,&
           &amean=a,silent=.true.)
      end select
      if(ierr/=0)then
        call msg('Could not perform fit.')
        close(io,status='delete')
        return
      endif
      ! Plot data minus independent bit.
      if(deval%nderiv==0)then
        do iset=1,ndataset
          ifit=dlist(iset)%dataset%ifit
          xy=>dlist(iset)%dataset%txy
          op_npoly=flist(ifit)%fit%npoly
          op_pow(1:op_npoly)=flist(ifit)%fit%pow(1:op_npoly)
          op_a(1:op_npoly)=a(1:op_npoly,iset)
          where(flist(ifit)%fit%share(1:op_npoly))op_a(1:op_npoly)=0.d0
          if(xy%have_dx.and.xy%have_dy)then
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i)),xy%dx(i),&
                 &xy%dy(i)
            enddo ! i
          elseif(xy%have_dx)then
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i)),xy%dx(i)
            enddo ! i
          elseif(xy%have_dy)then
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i)),xy%dy(i)
            enddo ! i
          else
            do i=1,xy%nxy
              write(io,*)xy%x(i),&
                 &xy%y(i)-eval_poly(op_npoly,op_pow,op_a,xy%x(i))
            enddo ! i
          endif
          write(io,'()')
          write(io,'()')
        enddo ! iset
      endif
      ! Plot shared part of fit functions.
      do i=1,deval%n
        write(io,*)deval%x(i),fmean(i,1),ferr(i,1)
        select case(plot_ci)
        case(0)
          write(io,*)deval%x(i),fmean(i,1),ferr(i,1)
        case(1)
          write(io,*)deval%x(i),fmean_xs(i,1)-ferr(i,1),&
             &fmean(i,1),fmean_xs(i,1)+ferr(i,1)
        case(2)
          write(io,*)deval%x(i),fmean_xs(i,1)-2.d0*ferr(i,1),&
             &fmean(i,1),fmean_xs(i,1)+2.d0*ferr(i,1)
        end select
      enddo ! i
      write(io,'()')
      write(io,'()')
    endif
    nullify(xy)

    ! Close file.
    close(io)

    ! Report.
    if(.not.do_append)call msg('Plot saved to "'//trim(fname)//'".')

  END SUBROUTINE plot_multipoly


  SUBROUTINE prob_multipoly(ndataset,dlist,drange,nfit,flist,deval,glob,&
     &fcondition)
    !--------------------------------!
    ! Perform chosen fit and report. !
    !--------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(eval_type),INTENT(inout) :: deval
    TYPE(global_params_type),INTENT(in) :: glob
    CHARACTER(*),INTENT(in) :: fcondition
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    ! Local variables.
    INTEGER iset,i,ierr
    DOUBLE PRECISION fcount(ndataset),tx0,tx1,t1

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
    call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
       &glob,deval,ierr,fcount=fcount,fcondition=fcondition,silent=.true.)
    if(ierr/=0)then
      call msg('Could not perform fit.')
      return
    endif

    ! Report.
    write(6,'(4x)',advance='no')
    write(6,'(2x,a4,2x,a20)')'set','Probability'
    write(6,'(4x)',advance='no')
    write(6,'(2x,26("-"))')
    do iset=1,ndataset
      write(6,'(a4)',advance='no')'PROB'
      write(6,'(2x,i4,2x,es20.12)')iset,fcount(iset)
    enddo ! iset
    write(6,'(4x)',advance='no')
    write(6,'(2x,26("-"))')
    write(6,'()')

  END SUBROUTINE prob_multipoly


  SUBROUTINE assess_fit(ndataset,dlist,drange,nfit,flist,glob,deval,&
     &setmask)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! expansion order.                               !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(eval_type),INTENT(in) :: deval
    LOGICAL,INTENT(in) :: setmask(ndataset)
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    ! Local variables.
    TYPE(fit_form_list_type),POINTER :: tflist(:)
    LOGICAL fit_used(nfit)
    INTEGER max_npoly,tot_nxy,tot_nparam,iset,ifit,npoly,ix,prev_npoly,ierr
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION,ALLOCATABLE :: chi2_all(:),chi2err_all(:),&
       &fmean_all(:,:,:),ferr_all(:,:,:)

    ! Compute array sizes.
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit
    fit_used=.false.
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      ifit=dataset%ifit
      fit_used(ifit)=.true.
    enddo ! iset

    ! Allocate arrays.
    allocate(chi2_all(max_npoly),chi2err_all(max_npoly),&
       &fmean_all(deval%n,ndataset,max_npoly),&
       &ferr_all(deval%n,ndataset,max_npoly))
    call clone_flist(flist,tflist)

    ! Loop over fits.
    do ifit=1,nfit
      if(.not.fit_used(ifit))cycle
      write(6,'(a)')'Fit #'//trim(i2s(ifit))//':'

      ! Print table header.
      write(6,'(2x,a5,1x,2(1x,a12))',advance='no')'Order','chi^2/Ndf',&
         &'dchi^2/Ndf'
      do iset=1,ndataset
        do ix=1,deval%n
          if(setmask(iset))write(6,'(2(1x,a12))',advance='no')&
             &'f'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  ',&
             &'df'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  '
        enddo ! ix
      enddo ! iset
      write(6,'()')
      write(6,'(2x,'//trim(i2s(32+26*deval%n*count(setmask)))//'("-"))')

      ! Initialize test results to "invalid".
      chi2_all=-1.d0
      chi2err_all=-1.d0
      fmean_all=0.d0
      ferr_all=-1.d0

      ! Loop over expansion orders.
      do npoly=1,flist(ifit)%fit%npoly
        ! Set expansion order.
        tflist(ifit)%fit%npoly=npoly
        ! Adjust counters.
        tot_nparam=count(tflist(ifit)%fit%share(1:npoly))
        tot_nxy=0
        do iset=1,ndataset
          dataset=>dlist(iset)%dataset
          if(ifit/=dataset%ifit)cycle
          tot_nxy=tot_nxy+dataset%rtxy%nxy
          tot_nparam=tot_nparam+count(.not.tflist(ifit)%fit%share(1:npoly))
        enddo ! iset
        if(tot_nparam<tot_nxy)then
          ! Perform fit and report.
          call eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,tflist,&
             &glob,deval,ierr,fmean,ferr,chi2,chi2err)
          if(ierr==0)then
            chi2_all(npoly)=chi2/dble(tot_nxy-tot_nparam)
            chi2err_all(npoly)=chi2err/dble(tot_nxy-tot_nparam)
            fmean_all(:,:,npoly)=fmean
            ferr_all(:,:,npoly)=ferr
            write(6,'(2x,i5,1x,2(1x,es12.4))',advance='no')npoly-1,&
               &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
            do ix=1,deval%n
              do iset=1,ndataset
                if(setmask(iset))write(6,'(2(1x,es12.4))',advance='no')&
                   &fmean(ix,iset),ferr(ix,iset)
              enddo ! iset
            enddo ! ix
            write(6,'()')
          endif
        endif
      enddo ! npoly

      ! Print table footer.
      write(6,'(2x,'//trim(i2s(32+26*deval%n*count(setmask)))//'("-"))')
      write(6,'()')

      ! Report best choice of parameters.
      ! Locate minimum chi2/Ndf.
      min_chi2=-1.d0
      do npoly=1,flist(ifit)%fit%npoly
        if(chi2_all(npoly)<0.d0)cycle
        t1=chi2_all(npoly)+2.d0*chi2err_all(npoly)
        if(min_chi2<0.d0.or.t1<min_chi2)min_chi2=t1
      enddo ! npoly
      ! Converge function value with respect to npoly.
      prev_npoly=0
      do npoly=1,flist(ifit)%fit%npoly
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
      if(npoly<=flist(ifit)%fit%npoly)then
        call msg('Suggested fit: '//trim(i2s(prev_npoly-1)))
      else
        call msg('Could not find optimal fit by criteria.')
      endif

    enddo ! ifit

    ! Clean up.
    call kill_flist(tflist)

  END SUBROUTINE assess_fit


  SUBROUTINE assess_range(ndataset,dlist,drange,nfit,flist,glob,deval,&
     &setmask)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! number of data points.                         !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(range_type),INTENT(inout) :: drange
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(eval_type),INTENT(in) :: deval
    LOGICAL,INTENT(in) :: setmask(ndataset)
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Pointer-resizing pointers.
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    ! Local variables.
    LOGICAL fit_used(nfit)
    INTEGER ifit,ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,ix,prev_igrid,ierr
    INTEGER,ALLOCATABLE :: indx(:),tngrid(:),tot_nxy_grid(:)
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION,ALLOCATABLE :: chi2_all(:),chi2err_all(:),&
       &fmean_all(:,:,:),ferr_all(:,:,:)
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:)

    ! Compute array sizes.
    tot_nparam=0
    tot_nxy=0
    fit_used=.false.
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      tot_nxy=tot_nxy+dataset%txy%nxy
      ifit=dataset%ifit
      tot_nparam=tot_nparam+count(.not.flist(ifit)%fit%share)
      if(.not.fit_used(ifit))tot_nparam=tot_nparam+count(flist(ifit)%fit%share)
      fit_used(ifit)=.true.
    enddo ! iset

    ! Compute grid.
    allocate(tot_nxy_grid(tot_nxy))
    select case(drange%thres_op)
    case('min','max')
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
    case default
      allocate(txall(tot_nxy),indx(tot_nxy),txgrid(tot_nxy),tngrid(tot_nxy))
      tngrid=0
      tot_nxy=0
      do iset=1,ndataset
        xy=>dlist(iset)%dataset%txy
        txall(tot_nxy+1:tot_nxy+xy%nxy)=xy%x(1:xy%nxy)
        tot_nxy=tot_nxy+xy%nxy
      enddo ! iset
      call isort_dble(tot_nxy,txall,indx)
      if(drange%op(1:1)=='>')then
        do ixy=1,tot_nxy
          if(ixy>=tot_nxy-ixy+1)exit
          call iswap1(indx(ixy),indx(tot_nxy-ixy+1))
        enddo ! ixy
      endif
      ngrid=1
      txgrid(1)=txall(indx(1))
      do ixy=2,tot_nxy
        if(eq_dble(txall(indx(ixy)),txgrid(ngrid)))cycle
        ngrid=ngrid+1
        txgrid(ngrid)=txall(indx(ixy))
      enddo ! ixy
      deallocate(txall,indx)
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
    write(6,'(2x,a5,1x,2(1x,a12))',advance='no')'Ndata','chi^2/Ndf',&
       &'dchi^2/Ndf'
    do ix=1,deval%n
      do iset=1,ndataset
        if(setmask(iset))write(6,'(2(1x,a12))',advance='no')&
           &'f'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  ',&
           &'df'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  '
      enddo ! iset
    enddo ! ix
    write(6,'()')
    write(6,'(2x,'//trim(i2s(32+26*deval%n*count(setmask)))//'("-"))')

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
        call eval_multifit_monte_carlo(ndataset,tmp_dlist,drange,nfit,flist,&
           &glob,deval,ierr,fmean,ferr,chi2,chi2err)
        if(ierr/=0)cycle
        chi2_all(igrid)=chi2/dble(tot_nxy-tot_nparam)
        chi2err_all(igrid)=chi2err/dble(tot_nxy-tot_nparam)
        fmean_all(:,:,igrid)=fmean
        ferr_all(:,:,igrid)=ferr
        write(6,'(2x,i5,1x,2(1x,es12.4))',advance='no')tot_nxy,&
           &chi2/dble(tot_nxy-tot_nparam),chi2err/dble(tot_nxy-tot_nparam)
        do ix=1,deval%n
          do iset=1,ndataset
            if(setmask(iset))write(6,'(2(1x,es12.4))',advance='no')&
               &fmean(ix,iset),ferr(ix,iset)
          enddo ! iset
        enddo ! ix
        write(6,'()')
      endif
    enddo ! igrid

    ! Print table footer.
    write(6,'(2x,'//trim(i2s(32+26*deval%n*count(setmask)))//'("-"))')
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
      call msg('Suggested grid point: '//trim(i2s(tot_nxy_grid(prev_igrid))))
    else
      call msg('Could not find optimal range by criteria.')
    endif

    ! Clean up.
    call kill_dlist(tmp_dlist)
    deallocate(chi2_all,chi2err_all,fmean_all,ferr_all)

  END SUBROUTINE assess_range


  SUBROUTINE assess_fit_range(ndataset,dlist,drange,nfit,flist,glob,deval,&
     &setmask)
    !------------------------------------------------!
    ! Test convergence of chi^2/Ndf as a function of !
    ! expansion order and number of data points.     !
    !------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(inout) :: drange
    TYPE(fit_form_list_type),POINTER :: flist(:)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(eval_type),INTENT(in) :: deval
    LOGICAL,INTENT(in) :: setmask(ndataset)
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    TYPE(fit_form_list_type),POINTER :: tflist(:)
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    LOGICAL fit_used(nfit)
    INTEGER ifit,ixy,igrid,ngrid,tot_nxy,tot_nparam,iset,npoly,ix,prev_igrid,&
       &prev_npoly,prev_prev_npoly,ierr
    INTEGER,ALLOCATABLE :: indx(:),tngrid(:),tot_nxy_grid(:)
    DOUBLE PRECISION chi2,chi2err,fmean(deval%n,ndataset),&
       &ferr(deval%n,ndataset),t1,t2,min_chi2
    DOUBLE PRECISION,ALLOCATABLE :: chi2_all(:,:),chi2err_all(:,:),&
       &fmean_all(:,:,:,:),ferr_all(:,:,:,:)
    DOUBLE PRECISION,ALLOCATABLE :: txall(:),txgrid(:)

    ! Compute array sizes.
    tot_nparam=0
    tot_nxy=0
    fit_used=.false.
    do iset=1,ndataset
      dataset=>dlist(iset)%dataset
      tot_nxy=tot_nxy+dataset%txy%nxy
      ifit=dataset%ifit
      tot_nparam=tot_nparam+count(.not.flist(ifit)%fit%share)
      if(.not.fit_used(ifit))tot_nparam=tot_nparam+count(flist(ifit)%fit%share)
      fit_used(ifit)=.true.
    enddo ! iset

    ! Compute grid.
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
      call isort_dble(tot_nxy,txall,indx)
      if(drange%op(1:1)=='>')then
        do ixy=1,tot_nxy
          if(ixy>=tot_nxy-ixy+1)exit
          call iswap1(indx(ixy),indx(tot_nxy-ixy+1))
        enddo ! ixy
      endif
      ngrid=1
      txgrid(1)=txall(indx(1))
      do ixy=2,tot_nxy
        if(eq_dble(txall(indx(ixy)),txgrid(ngrid)))cycle
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

    ! Allocate work arrays.
    call clone_flist(flist,tflist)

    ! Loop over fits.
    do ifit=1,nfit
      if(.not.fit_used(ifit))cycle
      write(6,'(a)')'Fit #'//trim(i2s(ifit))//':'

      ! Allocate arrays to store test results and initialize to "invalid".
      allocate(chi2_all(flist(ifit)%fit%npoly,tot_nxy),&
         &chi2err_all(flist(ifit)%fit%npoly,tot_nxy),&
         &fmean_all(deval%n,ndataset,flist(ifit)%fit%npoly,tot_nxy),&
         &ferr_all(deval%n,ndataset,flist(ifit)%fit%npoly,tot_nxy))
      chi2_all=-1.d0
      chi2err_all=-1.d0
      fmean_all=0.d0
      ferr_all=-1.d0

      ! Make copy of datasets.
      call clone_dlist(dlist,tmp_dlist)

      ! Print table header.
      write(6,'(2x,a5,2x,a5,1x,2(1x,a12))',advance='no')'Ndata','Order',&
         &'chi^2/Ndf','dchi^2/Ndf'
      do ix=1,deval%n
        do iset=1,ndataset
          if(setmask(iset))write(6,'(2(1x,a12))',advance='no')&
             &'f'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  ',&
             &'df'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  '
        enddo ! iset
      enddo ! ix
      write(6,'()')
      write(6,'(2x,'//trim(i2s(39+26*deval%n*count(setmask)))//'("-"))')

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
        do npoly=1,flist(ifit)%fit%npoly
          ! Set expansion order.
          tflist(ifit)%fit%npoly=npoly
          ! Adjust counters.
          tot_nparam=count(tflist(ifit)%fit%share(1:npoly))
          do iset=1,ndataset
            dataset=>dlist(iset)%dataset
            if(ifit/=dataset%ifit)cycle
            tot_nparam=tot_nparam+count(.not.tflist(ifit)%fit%share(1:npoly))
          enddo ! iset
          ! Adjust counters.
          if(tot_nparam<tot_nxy)then
            ! Perform fit and report.
            call eval_multifit_monte_carlo(ndataset,tmp_dlist,drange,&
               &nfit,tflist,glob,deval,ierr,fmean,ferr,chi2,chi2err)
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
                if(setmask(iset))write(6,'(2(1x,es12.4))',advance='no')&
                   &fmean(ix,iset),ferr(ix,iset)
              enddo ! iset
            enddo ! ix
            write(6,'()')
          endif
        enddo ! npoly
      enddo ! igrid

      ! Print table footer.
      write(6,'(2x,'//trim(i2s(39+26*deval%n*count(setmask)))//'("-"))')
      write(6,'()')

      ! Print second table header.
      write(6,'()')
      write(6,'(2x,a5,2x,a5,1x,2(1x,a12))',advance='no')'Ndata','Order',&
         &'chi^2/Ndf','dchi^2/Ndf'
      do ix=1,deval%n
        do iset=1,ndataset
          if(setmask(iset))write(6,'(2(1x,a12))',advance='no')&
             &'f'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  ',&
             &'df'//trim(i2s(iset))//'(X'//trim(i2s(ix))//')  '
        enddo ! iset
      enddo ! ix
      write(6,'()')
      write(6,'(2x,'//trim(i2s(39+26*deval%n*count(setmask)))//'("-"))')

      ! Report best choice of parameters.
      ! Locate minimum chi2/Ndf.
      min_chi2=-1.d0
      do igrid=1,ngrid
        do npoly=1,flist(ifit)%fit%npoly
          if(chi2_all(npoly,igrid)<0.d0)cycle
          t1=chi2_all(npoly,igrid)+2.d0*chi2err_all(npoly,igrid)
          if(min_chi2<0.d0.or.t1<min_chi2)min_chi2=t1
        enddo ! npoly
      enddo ! igrid
      ! Converge function value with respect to igrid.
      prev_igrid=0
      do igrid=ngrid,1,-1
        prev_npoly=0
        do npoly=1,flist(ifit)%fit%npoly
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
        if(npoly>flist(ifit)%fit%npoly)cycle
        ! Write table entry.
        write(6,'(2x,i5,2x,i5,1x,2(1x,es12.4))',advance='no')&
           &tot_nxy_grid(igrid),prev_npoly-1,&
           &chi2_all(prev_npoly,igrid),chi2err_all(prev_npoly,igrid)
        do ix=1,deval%n
          do iset=1,ndataset
            if(setmask(iset))write(6,'(2(1x,es12.4))',advance='no')&
               &fmean_all(ix,iset,prev_npoly,igrid),&
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
      write(6,'(2x,'//trim(i2s(39+26*deval%n*count(setmask)))//'("-"))')
      write(6,'()')

      ! Report suggestion.
      if(igrid>=1)then
        write(6,'(a)')'Suggested fit: '//trim(i2s(prev_prev_npoly-1))
        write(6,'(a)')'Suggested grid point: '//&
           &trim(i2s(tot_nxy_grid(prev_igrid)))
        write(6,'()')
      else
        call msg('Could not find optimal range by criteria.')
      endif

      ! Clean up.
      call kill_dlist(tmp_dlist)
      deallocate(chi2_all,chi2err_all,fmean_all,ferr_all)

    enddo ! ifit

    ! Clean up.
    call kill_flist(tflist)

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
      tx=1.d0/abs(x)
      if(have_dx.and.present(dx).and.present(dtx))dtx=dx*tx**2
    case(ITRANSF_LOG)
      tx=log(abs(x))
      if(have_dx.and.present(dx).and.present(dtx))dtx=dx/abs(x)
    case(ITRANSF_EXP)
      tx=exp(x)
      if(have_dx.and.present(dx).and.present(dtx))dtx=dx*tx
    end select
  END SUBROUTINE scale_transform


  SUBROUTINE scale_untransform(nxy,itransfx,tx,x,have_dx,dtx,dx)
    !---------------------------------------!
    ! Apply a scale transformation x <- tx. !
    !---------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: nxy,itransfx
    DOUBLE PRECISION,INTENT(in) :: tx(nxy)
    DOUBLE PRECISION,INTENT(inout) :: x(nxy)
    LOGICAL,INTENT(in) :: have_dx
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: dtx(nxy)
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: dx(nxy)
    if(present(dx))dx=0.d0
    select case(itransfx)
    case(ITRANSF_NONE)
      x=tx
      if(have_dx.and.present(dtx).and.present(dx))dx=dtx
    case(ITRANSF_REC)
      x=1.d0/abs(tx)
      if(have_dx.and.present(dtx).and.present(dx))dx=dtx*x**2
    case(ITRANSF_LOG)
      x=exp(tx)
      if(have_dx.and.present(dtx).and.present(dx))dx=dtx*x
    case(ITRANSF_EXP)
      x=log(abs(tx))
      if(have_dx.and.present(dtx).and.present(dx))dx=dtx/abs(tx)
    end select
  END SUBROUTINE scale_untransform


  SUBROUTINE back_transform(drange,dataset)
    !-----------------------------------------------------!
    ! Copy/transform data in dataset%rtxy to %txy and %xy !
    ! (i.e., in reverse).                                 !
    !-----------------------------------------------------!
    IMPLICIT NONE
    TYPE(range_type),INTENT(in) :: drange
    TYPE(dataset_type),POINTER :: dataset
    LOGICAL,ALLOCATABLE :: mask(:)
    INTEGER ixy,jxy
    ! Unrestrict range to get txy.
    allocate(mask(dataset%xy%nxy))
    call get_range_mask(drange,dataset%xy,dataset%txy,mask)
    jxy=0
    do ixy=1,dataset%txy%nxy
      if(.not.mask(ixy))cycle
      jxy=jxy+1
      dataset%txy%y(ixy)=dataset%rtxy%y(jxy)
      dataset%txy%dy(ixy)=dataset%rtxy%dy(jxy)
    enddo ! ixy
    deallocate(mask)
    ! Now transform back to xy.
    call scale_untransform(dataset%txy%nxy,dataset%itransfy,dataset%txy%y,&
       &dataset%xy%y,.true.,dataset%txy%dy,dataset%xy%dy)
  END SUBROUTINE back_transform


  ! COMPUTATION ROUTINES


  SUBROUTINE perform_multifit(ndataset,dlist,glob_X0,nfit,flist,chi2,a,ierr)
    !----------------------------------------------------!
    ! Perform least-squares fit of sets of (weighted) xy !
    ! data to the polynomial of exponents pow(1:npoly),  !
    ! with equal/independent coefficients for each set   !
    ! depending on the value of share(1:npoly).          !
    ! Wrapper around dependent/independent versions.     !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),INTENT(in) :: dlist(ndataset)
    DOUBLE PRECISION,INTENT(in) :: glob_X0
    TYPE(fit_form_list_type),INTENT(in) :: flist(nfit)
    DOUBLE PRECISION,INTENT(inout) :: chi2,a(:,:)
    INTEGER,INTENT(inout) :: ierr
    DOUBLE PRECISION tchi2,ta(size(a,1),size(a,2))
    INTEGER ifit,iset,nuse
    chi2=0.d0
    a=0.d0
    do ifit=1,nfit
      nuse=0
      do iset=1,ndataset
        if(dlist(iset)%dataset%ifit==ifit)nuse=nuse+1
      enddo ! iset
      if(nuse<1)cycle
      if(count(flist(ifit)%fit%share)>0)then
        call perform_multifit_dep(ndataset,dlist,glob_X0,ifit,flist(ifit)%fit,&
           &tchi2,ta,ierr)
      else
        call perform_multifit_indep(ndataset,dlist,glob_X0,ifit,&
           &flist(ifit)%fit,tchi2,ta,ierr)
      endif
      chi2=chi2+tchi2
      do iset=1,ndataset
        if(dlist(iset)%dataset%ifit==ifit)then
          a(:,iset)=ta(:,iset)
        endif
      enddo ! iset
    enddo ! ifit
  END SUBROUTINE perform_multifit


  SUBROUTINE perform_multifit_indep(ndataset,dlist,glob_X0,ifit,fit,chi2,a,ierr)
    !----------------------------------------------------!
    ! Perform least-squares fit of sets of (weighted) xy !
    ! data to the polynomial of exponents pow(1:npoly),  !
    ! with equal/independent coefficients for each set   !
    ! depending on the value of share(1:npoly).          !
    ! This version performs each fit independently.      !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,ifit
    TYPE(dataset_list_type),INTENT(in) :: dlist(ndataset)
    DOUBLE PRECISION,INTENT(in) :: glob_X0
    TYPE(fit_form_type),INTENT(in) :: fit
    DOUBLE PRECISION,INTENT(inout) :: chi2,a(:,:)
    INTEGER,INTENT(inout) :: ierr
    ! Quick-access pointers.
    TYPE(xy_type),POINTER :: xy
    ! Local variables.
    DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:),w(:)
    INTEGER npoly,nxy,iset,ixy,ipoly,jpoly,lwork,i,j
    DOUBLE PRECISION e_fit,set_weight
    DOUBLE PRECISION,ALLOCATABLE :: M(:,:),Minv(:,:),c(:),work(:)
    INTEGER,ALLOCATABLE :: ipiv(:)

    ! Initialize.
    chi2=0.d0
    a=0.d0
    ierr=0

    ! Perpare storage for fit.
    npoly=fit%npoly
    allocate(M(npoly,npoly),Minv(npoly,npoly),ipiv(npoly),c(npoly))

    ! Loop over datasets.
    do iset=1,ndataset
      if(dlist(iset)%dataset%ifit/=ifit)cycle
      xy=>dlist(iset)%dataset%rtxy
      nxy=xy%nxy
      if(nxy<npoly)cycle
      set_weight=dlist(iset)%dataset%weight
      allocate(x(nxy),y(nxy),w(nxy))
      x=xy%x-glob_X0
      y=xy%y
      w=xy%w
      ! Construct c vector and M matrix.
      do ipoly=1,npoly
        c(ipoly)=sum(w(:)*y(:)*x(:)**fit%pow(ipoly))
        do jpoly=1,npoly
          M(jpoly,ipoly)=sum(w(:)*x(:)**(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! jpoly
      enddo ! ipoly
      ! Invert M.
      Minv=M
      allocate(work(1))
      lwork=-1
      call dsytrf('L',npoly,Minv,npoly,ipiv,work,lwork,ierr)
      if(ierr/=0)return
      lwork=nint(work(1))
      deallocate(work)
      allocate(work(lwork),stat=ierr)
      if(ierr/=0)return
      call dsytrf('L',npoly,Minv,npoly,ipiv,work,lwork,ierr)
      if(ierr/=0)return
      deallocate(work)
      allocate(work(npoly),stat=ierr)
      if(ierr/=0)return
      call dsytri('L',npoly,Minv,npoly,ipiv,work,ierr)
      if(ierr/=0)return
      deallocate(work)
      ! Complete Minv.
      do i=1,npoly
        do j=i+1,npoly
          Minv(i,j)=Minv(j,i)
        enddo ! j
      enddo ! i
      ! Evaluate fit coefficients.
      a(1:npoly,iset)=matmul(Minv,c)
      ! Add to chi^2 value.
      do ixy=1,nxy
        e_fit=0.d0
        do ipoly=1,npoly
          e_fit=e_fit+a(ipoly,iset)*x(ixy)**fit%pow(ipoly)
        enddo ! ipoly
        chi2=chi2+(y(ixy)-e_fit)**2*w(ixy)*set_weight
      enddo ! ixy
      ! Clean up for next set.
      deallocate(x,y,w)
    enddo ! iset

  END SUBROUTINE perform_multifit_indep


  SUBROUTINE perform_multifit_dep(ndataset,dlist,glob_X0,ifit,fit,chi2,a,ierr)
    !----------------------------------------------------!
    ! Perform least-squares fit of sets of (weighted) xy !
    ! data to the polynomial of exponents pow(1:npoly),  !
    ! with equal/independent coefficients for each set   !
    ! depending on the value of share(1:npoly).          !
    ! This version performs all fits simultaenously.     !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,ifit
    TYPE(dataset_list_type),INTENT(in) :: dlist(ndataset)
    DOUBLE PRECISION,INTENT(in) :: glob_X0
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
    chi2=0.d0
    a=0.d0
    ierr=0

    ! Extract fit properties.
    tot_nparam=count(fit%share)+ndataset*count(.not.fit%share)
    max_nxy=0
    tot_nxy=0
    do iset=1,ndataset
      if(dlist(iset)%dataset%ifit/=ifit)cycle
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
      if(dlist(iset)%dataset%ifit/=ifit)cycle
      set_weight=dlist(iset)%dataset%weight
      xy=>dlist(iset)%dataset%rtxy
      if(xy%nxy==0)cycle
      x(1:xy%nxy,iset)=xy%x(1:xy%nxy)-glob_X0
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
        if(dlist(iset)%dataset%ifit/=ifit)cycle
        c(ieq)=c(ieq)+sum(w(:,iset)*y(:,iset)*x(:,iset)**fit%pow(ipoly))
      enddo ! iset
      ! Coefficients of shared parameters.
      jp=0
      do jpoly=1,fit%npoly
        if(.not.fit%share(jpoly))cycle
        jp=jp+1
        M(jp,ieq)=0.d0
        do iset=1,ndataset
          if(dlist(iset)%dataset%ifit/=ifit)cycle
          M(jp,ieq)=M(jp,ieq)+sum(w(:,iset)*x(:,iset)**&
             &(fit%pow(jpoly)+fit%pow(ipoly)))
        enddo ! iset
      enddo ! jpoly
      ! Coefficients of independent parameters.
      do jpoly=1,fit%npoly
        if(fit%share(jpoly))cycle
        do iset=1,ndataset
          if(dlist(iset)%dataset%ifit/=ifit)cycle
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
        if(dlist(iset)%dataset%ifit/=ifit)cycle
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
            if(dlist(jset)%dataset%ifit/=ifit)cycle
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

    ! NB, for some reason the following does not work in the presence of
    ! shared parameters.  I would have expected this to work more reliably
    ! than the above, if anything..
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
        if(dlist(iset)%dataset%ifit/=ifit)cycle
        ip=ip+1
        a(ipoly,iset)=c(ip)
      enddo ! iset
    enddo ! ipoly

    ! Return chi^2 value.
    chi2=0.d0
    do iset=1,ndataset
      if(dlist(iset)%dataset%ifit/=ifit)cycle
      do ixy=1,max_nxy
        e_fit=0.d0
        do ipoly=1,fit%npoly
          e_fit=e_fit+a(ipoly,iset)*x(ixy,iset)**fit%pow(ipoly)
        enddo ! ipoly
        chi2=chi2+(y(ixy,iset)-e_fit)**2*w(ixy,iset)
      enddo ! ixy
    enddo ! iset

  END SUBROUTINE perform_multifit_dep


  SUBROUTINE eval_multifit_monte_carlo(ndataset,dlist,drange,nfit,flist,&
     &glob,deval,ierr,fmean,ferr,chi2mean,chi2err,rmsymean,rmsyerr,&
     &amean,aerr,fmean_1s,ferr_1s,fmean_2s,ferr_2s,fmed,fskew,fkurt,fcount,&
     &fcondition,silent)
    !------------------------------------------------------!
    ! Perform Monte Carlo sampling of data space to obtain !
    ! fit values or derivatives at specified points with   !
    ! error bars.                                          !
    !------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ndataset,nfit
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(range_type),INTENT(in) :: drange
    TYPE(fit_form_list_type),INTENT(in) :: flist(nfit)
    TYPE(global_params_type),INTENT(in) :: glob
    TYPE(eval_type),INTENT(in) :: deval
    INTEGER,INTENT(inout) :: ierr
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: &
       &fmean(deval%n,ndataset),ferr(deval%n,ndataset),&
       &chi2mean,chi2err,rmsymean,rmsyerr,&
       &fmean_1s(deval%n,ndataset),ferr_1s(deval%n,ndataset),&
       &fmean_2s(deval%n,ndataset),ferr_2s(deval%n,ndataset),&
       &fmed(deval%n,ndataset),fskew(deval%n,ndataset),&
       &fkurt(deval%n,ndataset),&
       &fcount(ndataset)
    ! Arguments with size (max_npoly,ndataset) - left implicit.
    DOUBLE PRECISION,INTENT(inout),OPTIONAL :: amean(:,:),aerr(:,:)
    CHARACTER(*),INTENT(in),OPTIONAL :: fcondition
    LOGICAL,INTENT(in),OPTIONAL :: silent
    ! Quick-access pointers.
    TYPE(dataset_type),POINTER :: dataset
    TYPE(xy_type),POINTER :: xy,xy_orig
    ! Pointers.
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    ! Polynomials resulting from operations on fitting polynomial.
    INTEGER op_npoly
    DOUBLE PRECISION,ALLOCATABLE :: op_pow(:),op_a(:)
    ! Random sampling arrays.
    DOUBLE PRECISION w_vector(glob%nsample)
    DOUBLE PRECISION,ALLOCATABLE :: f_array(:,:,:),a_array(:,:,:),&
       &chi2_array(:),rmsy_array(:)
    ! Distribution analysis.
    DOUBLE PRECISION var,skew,kurt,f2s_lo,f1s_lo,f1s_hi,f2s_hi
    ! Misc.
    LOGICAL need_f,need_a,need_chi2,need_rmsy,be_silent
    INTEGER ifit,max_npoly,nsample,ipoly,ideriv,irandom,ix,iset
    DOUBLE PRECISION chi2,t1,f0,sum_weight
    DOUBLE PRECISION,ALLOCATABLE :: a(:,:)

    ! Initialize.
    ierr=0
    be_silent=.false.
    if(present(silent))be_silent=silent

    ! Make copy of datasets.
    call clone_dlist(dlist,tmp_dlist)

    ! Apply qrandom.
    call qrandom_apply(ndataset,tmp_dlist,drange,nfit,flist,glob%X0,&
       &report=.not.be_silent)

    ! Figure out what we need.
    need_f=present(fmean).or.present(ferr).or.present(fmean_1s).or.&
       &present(ferr_1s).or.present(fmean_2s).or.present(ferr_2s).or.&
       &present(fmed).or.present(fskew).or.present(fkurt).or.&
       &(present(fcount).and.present(fcondition))
    need_f=need_f.and.deval%n>0
    need_chi2=present(chi2mean).or.present(chi2err)
    need_rmsy=present(rmsymean).or.present(rmsyerr)
    need_a=present(amean).or.present(aerr)
    max_npoly=0
    do ifit=1,nfit
      max_npoly=max(max_npoly,flist(ifit)%fit%npoly)
    enddo ! ifit

    ! Determine what to evaluate (function per set, shared portion of
    ! function, "unshared" portion of function per set, sum of functions) and
    ! allocate f_array accordingly.
    if(need_f)then
      select case(trim(deval%what))
      case('shared')
        allocate(f_array(glob%nsample,deval%n,nfit))
      case('sum')
        allocate(f_array(glob%nsample,deval%n,1))
      case default
        allocate(f_array(glob%nsample,deval%n,ndataset))
      end select
      f_array=0.d0
    endif
    ! Allocate things of size max_npoly.
    if(need_f)allocate(op_pow(max_npoly),op_a(max_npoly))
    allocate(a(max_npoly,ndataset))
    a=0.d0
    ! Allocate other things.
    if(need_a)allocate(a_array(glob%nsample,max_npoly,ndataset))
    if(need_chi2)allocate(chi2_array(glob%nsample))
    if(need_rmsy)allocate(rmsy_array(glob%nsample))

    ! Compute total dataset weight.
    sum_weight=0.d0
    do iset=1,ndataset
      sum_weight=sum_weight+tmp_dlist(iset)%dataset%weight
    enddo ! iset

    ! Initialize.
    f0=0.d0
    nsample=glob%nsample

    ! Loop over random points.
    do irandom=1,nsample
      do iset=1,ndataset
        dataset=>tmp_dlist(iset)%dataset
        xy=>dataset%xy
        xy_orig=>dlist(iset)%dataset%xy
        ! Using xy%dx and xy%dy so we get qrandom accounted for.
        if(xy%have_dx)xy%x=xy_orig%x+gaussian_random_number(xy%dx)
        if(xy%have_dy)xy%y=xy_orig%y+gaussian_random_number(xy%dy)
        call refresh_dataset(dataset,drange)
      enddo ! iset
      call perform_multifit(ndataset,tmp_dlist,glob%X0,nfit,flist,chi2,a,ierr)
      if(ierr/=0)then
        call kill_dlist(tmp_dlist)
        return
      endif
      w_vector(irandom)=1.d0
      if(need_chi2)chi2_array(irandom)=chi2
      if(need_rmsy)rmsy_array(irandom)=sqrt(chi2/sum_weight)
      if(need_a)a_array(irandom,1:max_npoly,1:ndataset)=&
         &a(1:max_npoly,1:ndataset)
      ! Evaluate requested function.
      if(need_f)then
        select case(trim(deval%what))
        case('shared')
          ! Evaluate shared component of requested function for each fit form.
          do ifit=1,nfit
            op_npoly=flist(ifit)%fit%npoly
            op_pow(1:op_npoly)=flist(ifit)%fit%pow(1:op_npoly)
            op_a(1:op_npoly)=a(1:op_npoly,1)
            where(.not.flist(ifit)%fit%share)op_a(1:op_npoly)=0.d0
            do ideriv=1,deval%nderiv
              call deriv_poly(op_npoly,op_pow,op_a)
            enddo ! ideriv
            if(deval%rel)f0=eval_poly(op_npoly,op_pow,op_a,deval%Xrel-glob%X0)
            do ix=1,deval%n
              t1=eval_poly(op_npoly,op_pow,op_a,deval%x(ix)-glob%X0)
              f_array(irandom,ix,ifit)=t1-f0
            enddo ! ix
          enddo ! ifit
        case default
          do iset=1,ndataset
            ifit=dlist(iset)%dataset%ifit
            op_npoly=flist(ifit)%fit%npoly
            op_pow(1:op_npoly)=flist(ifit)%fit%pow(1:op_npoly)
            op_a(1:op_npoly)=a(1:op_npoly,iset)
            if(trim(deval%what)=='unshared')then
              where(flist(ifit)%fit%share)op_a=0.d0
            endif
            do ideriv=1,deval%nderiv
              call deriv_poly(op_npoly,op_pow,op_a)
            enddo ! ideriv
            if(deval%rel)f0=eval_poly(op_npoly,op_pow,op_a,deval%Xrel-glob%X0)
            do ix=1,deval%n
              t1=eval_poly(op_npoly,op_pow,op_a,deval%x(ix)-glob%X0)
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
      ifit=dlist(iset)%dataset%ifit
      if(need_a)then
        if(present(amean))amean(:,iset)=0.d0
        if(present(aerr))aerr(:,iset)=0.d0
        do ipoly=1,flist(ifit)%fit%npoly
          call characterize_dist(nsample,a_array(:,ipoly,iset),w_vector,&
             &mean=t1,var=var)
          if(present(amean))amean(ipoly,iset)=t1
          if(present(aerr))aerr(ipoly,iset)=sqrt(var)
        enddo ! ipoly
      endif ! need_a
      if(need_f)then
        if(trim(deval%what)=='shared'.and.iset>nfit)then
          ! FIXME - why bother with this?
          if(present(fmean))fmean(:,iset)=fmean(:,1)
          if(present(ferr))ferr(:,iset)=ferr(:,1)
          if(present(fskew))fskew(:,iset)=fskew(:,1)
          if(present(fkurt))fkurt(:,iset)=fkurt(:,1)
          if(present(fmed))fmed(:,iset)=fmed(:,1)
          if(present(fmean_1s))fmean_1s(:,iset)=fmean_1s(:,1)
          if(present(ferr_1s))ferr_1s(:,iset)=ferr_1s(:,1)
          if(present(fmean_2s))fmean_2s(:,iset)=fmean_2s(:,1)
          if(present(ferr_2s))ferr_2s(:,iset)=ferr_2s(:,1)
        elseif(trim(deval%what)=='sum'.and.iset>1)then
          ! FIXME - why bother with this?
          if(present(fmean))fmean(:,iset)=fmean(:,1)
          if(present(ferr))ferr(:,iset)=ferr(:,1)
          if(present(fskew))fskew(:,iset)=fskew(:,1)
          if(present(fkurt))fkurt(:,iset)=fkurt(:,1)
          if(present(fmed))fmed(:,iset)=fmed(:,1)
          if(present(fmean_1s))fmean_1s(:,iset)=fmean_1s(:,1)
          if(present(ferr_1s))ferr_1s(:,iset)=ferr_1s(:,1)
          if(present(fmean_2s))fmean_2s(:,iset)=fmean_2s(:,1)
          if(present(ferr_2s))ferr_2s(:,iset)=ferr_2s(:,1)
        else ! statistics not yet computed
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
        endif ! "iset" index out of or in range for set/fit/global functions
        if(present(fcount).and.present(fcondition))then
          fcount(iset)=0.d0
          do irandom=1,nsample
            if(fcondition=='negative')then
              if(any(f_array(irandom,1:deval%n,iset)<0.d0))&
                 &fcount(iset)=fcount(iset)+1.d0
            else
              if(any(f_array(irandom,1:deval%n,iset)>0.d0))&
                 &fcount(iset)=fcount(iset)+1.d0
            endif
          enddo ! irandom
          fcount(iset)=fcount(iset)/dble(nsample)
        endif
      endif ! need_f
    enddo ! iset

    ! Clean up.
    call kill_dlist(tmp_dlist)

  END SUBROUTINE eval_multifit_monte_carlo


  SUBROUTINE construct_beta_dataset(drange,dataset_i,dataset_j,beta,dataset,&
     &ierr)
    !------------------------------------------!
    ! Construct linear combination of datasets !
    ! Y = yi + beta*(yj-yi).                   !
    !------------------------------------------!
    IMPLICIT NONE
    TYPE(range_type),INTENT(in) :: drange
    TYPE(dataset_type),POINTER :: dataset_i,dataset_j,dataset
    DOUBLE PRECISION,INTENT(in) :: beta
    INTEGER,INTENT(inout) :: ierr
    INTEGER ixy
    ierr=1
    call clone_dataset(dataset_i,dataset)
    if(dataset_i%rtxy%nxy/=dataset_j%rtxy%nxy)return
    if(any(neq_dble(dataset_i%rtxy%x,dataset_j%rtxy%x)))return
    ierr=0
    do ixy=1,dataset%rtxy%nxy
      dataset%rtxy%y(ixy)=dataset_i%rtxy%y(ixy)+&
        &beta*(dataset_j%rtxy%y(ixy)-dataset_i%rtxy%y(ixy))
      dataset%rtxy%dy(ixy)=sqrt(&
        &((1.d0-beta)*dataset_i%rtxy%dy(ixy))**2+&
        &(beta*dataset_j%rtxy%dy(ixy))**2)
    enddo ! ixy
    call back_transform(drange,dataset)
  END SUBROUTINE construct_beta_dataset


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
     &nsearch,fsearch,search,ndiscr,fdiscr,dataset_model,ndataset,dlist,ierr)
    !-------------------------------------!
    ! Read in the data in the input file. !
    !-------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(in) :: icol_x,icol_y,icol_dx,icol_dy,icol_w,nsearch,&
       &fsearch(nsearch),ndiscr,fdiscr(ndiscr)
    CHARACTER(*),INTENT(in) :: search(nsearch)
    INTEGER,INTENT(inout) :: ndataset,ierr
    TYPE(dataset_type),INTENT(in) :: dataset_model
    TYPE(dataset_list_type),POINTER :: dlist(:)
    TYPE(dataset_list_type),POINTER :: tmp_dlist(:)
    TYPE(dataset_type),POINTER :: dataset
    CHARACTER(8192) line,label,sname
    CHARACTER(8192),POINTER :: discr(:,:)
    INTEGER i,ipos,isearch,iset,idiscr,io
    LOGICAL neg,match
    ! Constants.
    INTEGER, PARAMETER :: io_stdin=5,io_file=10

    label=''
    if(fname(1:2)=='<<')then
      io=io_stdin
      sname='stdin'
      label=adjustl(fname(3:))
      if(len_trim(label)==0)then
        call msg('End-of-input label must not be an empty string.')
        return
      endif
    else
      io=io_file
      sname='"'//trim(fname)//'"'
      ! Open file.
      open(unit=io,file=trim(fname),status='old',iostat=ierr)
      if(ierr/=0)then
        call msg('Problem opening '//trim(sname)//'.')
        return
      endif
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
        call msg('Problem getting line from '//trim(sname)//'.')
        exit
      endif
      line=adjustl(line)
      ! Skip comments.
      ipos=scan(line,'#!')
      if(ipos==1)cycle
      if(ipos>1)line=line(1:ipos-1)
      ! Skip empty lines.
      if(len_trim(line)==0)cycle
      ! See if this is the end-of-stream label.
      if(len_trim(label)>0)then
        if(trim(line)==trim(label))exit
      endif
      ! Verify that this line contains all search strings.
      do isearch=1,nsearch
        neg=fsearch(isearch)<0
        match=eq_dble_string(trim(field(abs(fsearch(isearch)),line)),&
          &trim(search(isearch)))
        if(neg.eqv.match)exit
      enddo ! isearch
      if(isearch<=nsearch)cycle
      ! Decide which dataset this goes in.
      do iset=1,ndataset
        do idiscr=1,ndiscr
          if(.not.eq_dble_string(trim(field(fdiscr(idiscr),line)),&
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
        dataset%apply_qrandom=dataset_model%apply_qrandom
        dataset%qrandom_exp=dataset_model%qrandom_exp
        dataset%qrandom_centre=dataset_model%qrandom_centre
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
          call msg('Failed to parse value of x in '//trim(sname)//'.')
          exit
        endif
      else
        dataset%xy%x(i)=dble(i)
      endif
      if(icol_y>0)then
        dataset%xy%y(i)=dble_field(icol_y,line,ierr)
        if(ierr/=0)then
          call msg('Failed to parse value of y in '//trim(sname)//'.')
          exit
        endif
      else
        dataset%xy%y(i)=dble(i)
      endif
      if(icol_dx>0)then
        dataset%xy%dx(i)=dble_field(icol_dx,line,ierr)
        if(ierr/=0)then
          call msg('Failed to parse value of dx in '//trim(sname)//'.')
          exit
        endif
        if(lt_dble(dataset%xy%dx(i),0.d0))then
          call msg('Found negative dx in '//trim(sname)//'.')
          ierr=-5
          exit
        endif
      endif ! have_dx
      if(icol_dy>0)then
        dataset%xy%dy(i)=dble_field(icol_dy,line,ierr)
        if(ierr/=0)then
          call msg('Failed to parse value of dy in '//trim(sname)//'.')
          exit
        endif
        if(lt_dble(dataset%xy%dy(i),0.d0))then
          call msg('Found negative dy in '//trim(sname)//'.')
          ierr=-6
          exit
        endif
      endif ! have_dy
      if(icol_w>0)then
        dataset%xy%w(i)=dble_field(icol_w,line,ierr)
        if(ierr/=0)then
          call msg('Failed to parse value of w in '//trim(sname)//'.')
          exit
        endif
        if(le_dble(dataset%xy%w(i),0.d0))then
          call msg('Found non-positive w in '//trim(sname)//'.')
          ierr=-7
          exit
        endif
      endif ! have_w
    enddo ! i

    ! Close file.
    if(io==io_file)close(io)

    ! Clean up.
    if(associated(discr))deallocate(discr)

  END SUBROUTINE read_file


  ! POLYNOMIAL HANDLING UTILITIES.


  FUNCTION print_poly_sym(fit,glob_X0) RESULT(polystr)
    !---------------------------------------------------!
    ! Returns a string with the symbolic version of the !
    ! fitted polynomial.                                !
    !---------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),INTENT(in) :: fit
    DOUBLE PRECISION,INTENT(in) :: glob_X0
    CHARACTER(3+fit%npoly*54) :: polystr
    INTEGER j,ipow
    CHARACTER(40) pwstr
    CHARACTER(6) xstr
    if(abs(glob_X0)<tol_zero)then
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


  FUNCTION print_poly_num(fit,glob_X0,a) RESULT(polystr)
    !----------------------------------------------------!
    ! Returns a string with the numerical version of the !
    ! fitted polynomial in a suitable format for pasting !
    ! into xmgrace.                                      !
    !----------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),INTENT(in) :: fit
    DOUBLE PRECISION,INTENT(in) :: glob_X0,a(fit%npoly)
    CHARACTER(2+fit%npoly*105) :: polystr
    INTEGER j,ipow
    CHARACTER(1) plusstr
    CHARACTER(32) coeffstr
    CHARACTER(72) pwstr
    CHARACTER(36) xstr
    if(abs(glob_X0)<tol_zero)then
      xstr='x'
    else
      write(xstr,*)glob_X0
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
      elseif(neq_dble(x_target,0.d0))then
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
    if(associated(xy%w))deallocate(xy%w)
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
    dataset2%ifit=dataset1%ifit
    dataset2%weight=dataset1%weight
    dataset2%wexp=dataset1%wexp
    dataset2%apply_qrandom=dataset1%apply_qrandom
    dataset2%qrandom_exp=dataset1%qrandom_exp
    dataset2%qrandom_centre=dataset1%qrandom_centre
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


  SUBROUTINE clone_fit_form(fit1,fit2)
    !-------------------------------------------!
    ! Make an independent copy of fit1 as fit2. !
    !-------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),POINTER :: fit1,fit2
    nullify(fit2)
    if(.not.associated(fit1))return
    allocate(fit2)
    fit2%npoly=fit1%npoly
    allocate(fit2%pow(fit2%npoly),fit2%share(fit2%npoly))
    fit2%pow=fit1%pow(1:fit2%npoly)
    fit2%share=fit1%share(1:fit2%npoly)
  END SUBROUTINE clone_fit_form


  SUBROUTINE combine_fit_form(fit1,fit2,fit,map1,map2)
    !-------------------------------------------!
    ! Make an independent fit_form object which !
    ! combines fit1 and fit2.                   !
    !-------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),POINTER :: fit1,fit2,fit
    INTEGER,INTENT(inout) :: map1(fit1%npoly),map2(fit2%npoly)
    DOUBLE PRECISION pow(fit1%npoly+fit2%npoly)
    INTEGER ipoly,jpoly,npoly,indx(fit1%npoly+fit2%npoly)

    ! Initialize.
    nullify(fit)
    if(.not.associated(fit1))return
    if(.not.associated(fit2))return

    ! For same fit use clone routine so we get "shared" right.
    if(associated(fit1,fit2))then
      call clone_fit_form(fit1,fit)
      map1(1:fit1%npoly)=(/(ipoly,ipoly=1,fit1%npoly)/)
      map2(1:fit2%npoly)=(/(ipoly,ipoly=1,fit2%npoly)/)
      return
    endif

    ! Copy exponent lists, then eliminate repeated exponents.
    npoly=fit1%npoly+fit2%npoly
    pow(1:npoly)=(/fit1%pow(1:fit1%npoly),fit2%pow(1:fit2%npoly)/)
    call isort_dble(npoly,pow,indx)
    pow(1:npoly)=pow(indx(1:npoly))
    ipoly=0
    do
      if(ipoly==npoly)exit
      ipoly=ipoly+1
      jpoly=ipoly
      do
        if(jpoly==npoly)exit
        jpoly=jpoly+1
        if(eq_dble(pow(ipoly),pow(jpoly)))then
          pow(jpoly:npoly-1)=pow(jpoly+1:npoly)
          npoly=npoly-1
          jpoly=jpoly-1
        endif
      enddo ! jpoly
    enddo ! ipoly

    ! Build index maps.
    do ipoly=1,fit1%npoly
      do jpoly=1,npoly
        if(eq_dble(fit1%pow(ipoly),pow(jpoly)))exit
      enddo ! jpoly
      map1(ipoly)=jpoly
    enddo ! ipoly
    do ipoly=1,fit2%npoly
      do jpoly=1,npoly
        if(eq_dble(fit2%pow(ipoly),pow(jpoly)))exit
      enddo ! jpoly
      map2(ipoly)=jpoly
    enddo ! ipoly

    ! Make combined fit.
    allocate(fit)
    fit%npoly=npoly
    allocate(fit%pow(fit%npoly),fit%share(npoly))
    fit%pow=pow(1:fit%npoly)
    fit%share=.false.

  END SUBROUTINE combine_fit_form


  SUBROUTINE combine_fit_coeffs(fit1,fit2,cfit,map1,map2,a1,a2,c1,c2,a)
    !---------------------------------------------------------------!
    ! Given fit forms fit1, fit2 and their combination cfit, return !
    ! the coefficients of the linear combination c1*fit1 + c2*fit2. !
    !---------------------------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),POINTER :: fit1,fit2,cfit
    INTEGER,INTENT(in) :: map1(fit1%npoly),map2(fit2%npoly)
    DOUBLE PRECISION,INTENT(in) :: a1(fit1%npoly),a2(fit2%npoly),c1,c2
    DOUBLE PRECISION,INTENT(inout) :: a(cfit%npoly)
    a=0.d0
    a(map1(1:fit1%npoly))=a(map1(1:fit1%npoly))+c1*a1(1:fit1%npoly)
    a(map2(1:fit2%npoly))=a(map2(1:fit2%npoly))+c2*a2(1:fit2%npoly)
  END SUBROUTINE combine_fit_coeffs


  SUBROUTINE clone_flist(flist1,flist2)
    !-----------------------------------------------!
    ! Make an independent copy of flist1 as flist2. !
    !-----------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_list_type),POINTER :: flist1(:),flist2(:)
    INTEGER ifit,nfit
    nullify(flist2)
    if(.not.associated(flist1))return
    nfit=size(flist1)
    allocate(flist2(nfit))
    do ifit=1,nfit
      call clone_fit_form(flist1(ifit)%fit,flist2(ifit)%fit)
    enddo ! ifit
  END SUBROUTINE clone_flist


  SUBROUTINE kill_fit_form(fit)
    !---------------------------------------!
    ! Destroy a fit_form_type-type pointer. !
    !---------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_type),POINTER :: fit
    if(.not.associated(fit))return
    if(allocated(fit%pow))deallocate(fit%pow)
    if(allocated(fit%share))deallocate(fit%share)
    deallocate(fit)
    nullify(fit)
  END SUBROUTINE kill_fit_form


  SUBROUTINE kill_flist(flist)
    !--------------------------------------------!
    ! Destroy a fit_form_list_type-type pointer. !
    !--------------------------------------------!
    IMPLICIT NONE
    TYPE(fit_form_list_type),POINTER :: flist(:)
    INTEGER nfit,ifit
    if(.not.associated(flist))return
    nfit=size(flist)
    do ifit=1,nfit
      call kill_fit_form(flist(ifit)%fit)
    enddo ! ifit
    deallocate(flist)
    nullify(flist)
  END SUBROUTINE kill_flist


  ! GENERIC NUMERICAL UTILITIES.


  SUBROUTINE parabolic_min(x1,x2,x3,y1,y2,y3,x0,y0,rejected)
    !-----------------------------------------------------------------!
    ! Fit three points to a parabola and return (x,y) of the min/max. !
    !-----------------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x1,x2,x3,y1,y2,y3
    DOUBLE PRECISION,INTENT(out) :: x0,y0
    LOGICAL,INTENT(out) :: rejected
    DOUBLE PRECISION a,b,c,numa,numb,numc,den,x1_sq,x2_sq,x3_sq,&
       &x21,x32,x31,z1,z2,z3,invden
    ! Initialize.
    x0=x2
    y0=y2
    rejected=.false.
    ! Check that x and y values are distinguishable and in correct order.
    if(x1>=x2.or.x2>=x3.or.(y2>=y1.eqv.y3>=y2))then
      rejected=.true.
      return
    endif
    ! Compute squares.
    x1_sq=x1*x1
    x2_sq=x2*x2
    x3_sq=x3*x3
    ! Renormalize for better numerics.
    x31=x3-x1
    x21=(x2-x1)/x31
    x32=(x3-x2)/x31
    z1=y1*x32
    z2=y2
    z3=y3*x21
    ! Solve linear system.
    den=-x1_sq*x32+x2_sq-x3_sq*x21
    numa=-z1+z2-z3
    numb=z1*(x2+x3)-z2*(x1+x3)+z3*(x1+x2)
    numc=-z1*x2*x3+z2*x3*x1-z3*x1*x2
    ! Find x0 and y0.
    invden=1.d0/den
    a=numa*invden
    b=numb*invden
    c=numc*invden
    x0=-0.5d0*numb/numa
    y0=(a*x0+b)*x0+c
  END SUBROUTINE parabolic_min


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
      if(K2<=0.d0)return
    endif
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


  SUBROUTINE isort_dble(n,x,indx)
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
  END SUBROUTINE isort_dble


  SUBROUTINE isort_int(n,x,indx)
    !----------------------------------------------------------!
    ! Perform insertion sort on an integer vector X(1:N) so    !
    ! that X(I)<X(J) if I<J.  This sorting algorithm typically !
    ! costs ~ N^2, but is stable (preserves the order of       !
    ! entries with same X) and becomes order N when X(:) is    !
    ! nearly sorted.                                           !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    INTEGER,INTENT(in) :: x(n)
    INTEGER,INTENT(inout) :: indx(n)
    INTEGER i,j
    INTEGER xi
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
  END SUBROUTINE isort_int


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


  LOGICAL ELEMENTAL FUNCTION eq_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X == Y within a reasonable tolerance.             !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    DOUBLE PRECISION abs_x,abs_y,big,small
    ! Parameters.
    DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-12
    DOUBLE PRECISION,PARAMETER :: tol_rel=1.d-9
    if(present(tol))then
      eq_dble=abs(x-y)<=tol
    else
      abs_x=abs(x)
      abs_y=abs(y)
      if(abs_x<=tol_zero.and.abs_y<=tol_zero)then
        eq_dble=.true.
      elseif(x>0.d0.eqv.y>0.d0)then
        big=max(abs_x,abs_y)
        small=min(abs_x,abs_y)
        eq_dble=big-small<=big*tol_rel
      else
        eq_dble=.false.
      endif
    endif
  END FUNCTION eq_dble


  LOGICAL ELEMENTAL FUNCTION neq_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X /= Y within a reasonable tolerance.             !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    neq_dble=.not.eq_dble(x,y,tol)
  END FUNCTION neq_dble


  LOGICAL ELEMENTAL FUNCTION lt_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X < Y within a reasonable tolerance.              !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    lt_dble=x<y.and.neq_dble(x,y,tol)
  END FUNCTION lt_dble


  LOGICAL ELEMENTAL FUNCTION le_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X <= Y within a reasonable tolerance.             !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    le_dble=x<y.or.eq_dble(x,y,tol)
  END FUNCTION le_dble


  LOGICAL ELEMENTAL FUNCTION gt_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X > Y within a reasonable tolerance.              !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    gt_dble=x>y.and.neq_dble(x,y,tol)
  END FUNCTION gt_dble


  LOGICAL ELEMENTAL FUNCTION ge_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X >= Y within a reasonable tolerance.             !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    ge_dble=x>y.or.eq_dble(x,y,tol)
  END FUNCTION ge_dble


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
      if(eq_dble(t2,0.d0))then
        ierr=-1
        return
      endif
      parse_dble=t1/t2
    endif
  END FUNCTION parse_dble


  LOGICAL FUNCTION eq_dble_string(cx,cy,tol)
    !-------------------------------------------------------!
    ! Check if two strings are equal, or if their numerical !
    ! values are equal within a reasonable tolerance.       !
    !-------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: cx,cy
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    DOUBLE PRECISION x,y
    INTEGER ierr
    eq_dble_string=trim(cx)==trim(cy)
    if(eq_dble_string)return
    x=parse_dble(cx,ierr)
    if(ierr/=0)return
    y=parse_dble(cy,ierr)
    if(ierr/=0)return
    eq_dble_string=eq_dble(x,y,tol)
  END FUNCTION eq_dble_string


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
    drange%thres_op='non'
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

    ! Get operator.
    select case(remainder(1:1))
    case('<','>')
      if(remainder(2:2)=='=')then
        drange%op=remainder(1:2)
        remainder=remainder(3:)
      else
        drange%op=remainder(1:1)
        remainder=remainder(2:)
      endif
    case default
      return
    end select

    ! Get sort operation (if any) and threshold point.
    select case(remainder(1:3))
    case('max','min')
      drange%thres_op=remainder(1:3)
      remainder=remainder(4:)
      drange%size=parse_int(remainder,ierr)
    case default
      drange%thres=parse_dble(remainder,ierr)
    end select
    drange%no_rhs=ierr/=0

  END SUBROUTINE parse_range


  SUBROUTINE parse_xeval(string,deval)
    !--------------------------------------------------------!
    ! Parse a list of x values from a string of the form     !
    ! <variable>=<list>.  See parse_rlist for <list> format. !
    !--------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    TYPE(eval_type),INTENT(inout) :: deval
    CHARACTER(len_trim(string)) remainder

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

    ! Parse range.
    call parse_rlist(remainder,deval%x)
    if(.not.associated(deval%x))return
    deval%n=size(deval%x,1)

  END SUBROUTINE parse_xeval


  SUBROUTINE parse_ilist(string,xlist,sortuniq)
    !-------------------------------------------------------------!
    ! Parse a list of integers expressed as a comma-separated     !
    ! list of tokens, in which each token can be a single number  !
    ! or a uniform grid expressed as "X1:X2:N" or "X1:X2" [N is   !
    ! then assumed to be abs(X2-X1)].  Returns a null pointer if  !
    ! there are parsing errors.  If sortuniq=.true. (default),    !
    ! xlist is sorted and deduplicated.                           !
    ! E.g., "1,4,2,3", "1:4", and "2:4:3,1" all give {1,2,3,4}.   !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER, POINTER :: xlist(:)
    LOGICAL, INTENT(in), OPTIONAL :: sortuniq
    CHARACTER(len_trim(string)) remainder, token
    LOGICAL is_valid, do_sortuniq
    LOGICAL, ALLOCATABLE :: lmask(:)
    INTEGER ipos, ipos1, ipos2, it1, nx, i, ierr
    INTEGER t1, t2
    INTEGER, POINTER :: xx(:)
    INTEGER, ALLOCATABLE :: indx(:)

    ! Initialize.
    nullify(xlist)
    remainder=string
    nullify(xx)
    nx=0
    do_sortuniq=.true.
    if(present(sortuniq))do_sortuniq=sortuniq

    ! Loop over tokens in comma-separated list.
    do while(len_trim(remainder)>0)

      ! Mark as not valid until we are done, so we can just "exit" on errors.
      is_valid=.false.

      ! Get token before next comma.
      ipos=scan(remainder,',')
      if(ipos<1)then
        token=remainder
        remainder=''
      else
        token=remainder(:ipos-1)
        if(ipos==len_trim(remainder))then
          exit ! empty last item
        else
          remainder=remainder(ipos+1:)
        endif
      endif

      ! Determine type of token.
      ipos1=scan(token,':')
      ipos2=scan(token,':',back=.true.)
      if(ipos1==0)then
        ! Single number.
        t1=parse_int(token,ierr)
        if(ierr/=0)exit
        t2=t1
        it1=1
      else
        ! Automatic list.
        t1=parse_int(token(1:ipos1-1),ierr)
        if(ierr/=0)exit
        if(ipos1==ipos2)then
          t2=parse_int(token(ipos1+1:),ierr)
          if(ierr/=0)exit
          it1=abs(t2-t1)+1
        else
          t2=parse_int(token(ipos1+1:ipos2-1),ierr)
          if(ierr/=0)exit
          it1=parse_int(token(ipos2+1:),ierr)
          if(ierr/=0)exit
          if(it1<1)exit
          if(it1>1)then
            if(mod(abs(t2-t1),it1-1)/=0)exit
          endif
        endif
      endif

      ! Parsing done, so token is valid.
      is_valid=.true.

      ! Add list to vector.
      call resize_pointer_int1((/nx+it1/),xx)
      xx(nx+1)=t1
      if(it1>1)then
        t2=(t2-t1)/(it1-1)
        do i=2,it1
          t1=t1+t2
          xx(nx+i)=t1
        enddo ! i
      endif
      nx=nx+it1

    enddo ! while remainder/=''

    ! Decide whether to return a vector or not.
    if(.not.is_valid)then
      if(associated(xx))deallocate(xx)
      return
    endif
    xlist=>xx
    if(nx==0)return
    if(.not.do_sortuniq)return

    ! Sort vector.
    allocate(indx(nx))
    call isort_int(nx,xlist,indx)
    xlist=xlist(indx)
    deallocate(indx)

    ! Remove duplicates.
    allocate(lmask(nx))
    lmask(1)=.true.
    lmask(2:nx)=xlist(2:nx)/=xlist(1:nx-1)
    nx=count(lmask)
    xlist(1:nx)=pack(xlist,lmask)
    deallocate(lmask)
    call resize_pointer_int1((/nx/),xlist)

  END SUBROUTINE parse_ilist


  SUBROUTINE parse_rlist(string,xlist,sortuniq)
    !-------------------------------------------------------------!
    ! Parse a list of real numbers expressed as a comma-separated !
    ! list of tokens, in which each token can be a single number  !
    ! or a uniform grid expressed as "X1:X2:N" or "X1:X2" [N is   !
    ! then assumed to be abs(X2-X1)].  Returns a null pointer if  !
    ! there are parsing errors.  If sortuniq=.true. (default),    !
    ! xlist is sorted and deduplicated.                           !
    ! E.g., "1,2,3,4", "1:4", and "1,2:4:3" all give {1,2,3,4}.   !
    !-------------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    DOUBLE PRECISION, POINTER :: xlist(:)
    LOGICAL, INTENT(in), OPTIONAL :: sortuniq
    CHARACTER(len_trim(string)) remainder, token
    LOGICAL is_valid, do_sortuniq
    LOGICAL, ALLOCATABLE :: lmask(:)
    INTEGER ipos, ipos1, ipos2, it1, nx, i, ierr
    INTEGER, ALLOCATABLE :: indx(:)
    DOUBLE PRECISION t1, t2
    DOUBLE PRECISION, POINTER :: xx(:)

    ! Initialize.
    nullify(xlist)
    remainder=string
    nullify(xx)
    nx=0
    do_sortuniq=.true.
    if(present(sortuniq))do_sortuniq=sortuniq

    ! Loop over tokens in comma-separated list.
    do while(len_trim(remainder)>0)

      ! Mark as not valid until we are done, so we can just "exit" on errors.
      is_valid=.false.

      ! Get token before next comma.
      ipos=scan(remainder,',')
      if(ipos<1)then
        token=remainder
        remainder=''
      else
        token=remainder(:ipos-1)
        if(ipos==len_trim(remainder))then
          exit ! empty last item
        else
          remainder=remainder(ipos+1:)
        endif
      endif

      ! Determine type of token.
      ipos1=scan(token,':')
      ipos2=scan(token,':',back=.true.)
      if(ipos1==0)then
        ! Single number.
        t1=parse_dble(token,ierr)
        if(ierr/=0)exit
        t2=t1
        it1=1
      else
        ! Automatic list.
        t1=parse_dble(token(1:ipos1-1),ierr)
        if(ierr/=0)exit
        if(ipos1==ipos2)then
          t2=parse_dble(token(ipos1+1:),ierr)
          if(ierr/=0)exit
          if(.not.eq_dble(abs(t2-t1),dble(nint(abs(t2-t1)))))exit
          it1=nint(abs(t2-t1))+1
        else
          t2=parse_dble(token(ipos1+1:ipos2-1),ierr)
          if(ierr/=0)exit
          it1=parse_int(token(ipos2+1:),ierr)
          if(ierr/=0)exit
          if(it1<1.or.(it1==1.and..not.eq_dble(t1,t2)))exit
        endif
      endif

      ! Parsing done, so token is valid.
      is_valid=.true.

      ! Add list to vector.
      call resize_pointer_dble1((/nx+it1/),xx)
      xx(nx+1)=t1
      if(it1>1)then
        t2=(t2-t1)/dble(it1-1)
        do i=2,it1
          t1=t1+t2
          xx(nx+i)=t1
        enddo ! i
      endif
      nx=nx+it1

    enddo ! while remainder/=''

    ! Decide whether to return a vector or not.
    if(.not.is_valid)then
      if(associated(xx))deallocate(xx)
      return
    endif
    xlist=>xx
    if(nx==0)return
    if(.not.do_sortuniq)return

    ! Sort vector.
    allocate(indx(nx))
    call isort_dble(nx,xlist,indx)
    xlist=xlist(indx)
    deallocate(indx)

    ! Remove duplicates.
    allocate(lmask(nx))
    lmask(1)=.true.
    lmask(2:nx)=.not.eq_dble(xlist(2:nx),xlist(1:nx-1))
    nx=count(lmask)
    xlist(1:nx)=pack(xlist,lmask)
    deallocate(lmask)
    call resize_pointer_dble1((/nx/),xlist)

  END SUBROUTINE parse_rlist


  SUBROUTINE parse_search_clause(string,icol,search,neg)
    !---------------------------------------------------------!
    ! Given string "$n==search" or "$n!=search" return icol=n !
    ! (where n is a positive integer) and search (a string),  !
    ! or icol=0 if there is an error.                         !
    !---------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: icol
    CHARACTER(*),INTENT(inout) :: search
    LOGICAL,INTENT(inout) :: neg
    INTEGER ipos

    ! Initialize.
    icol=0
    search=''
    neg=.false.

    ! See if this is a negative clause.
    ipos=scan(string,'!')
    if(ipos<1)ipos=scan(string,'/')
    if(ipos>0)then
      neg=.true.
    else
      ! Find '=' and parse column number to its left.
      ipos=scan(string,'=')
      if(ipos<1)return
    end if
    call parse_colnum(string(1:ipos-1),icol)
    if(icol<1)return

    ! Find search string.
    if(string(ipos+1:ipos+1)/='=')then
      icol=0
      return
    endif
    search=string(ipos+2:)

  END SUBROUTINE parse_search_clause


  SUBROUTINE parse_colnum(string,icol)
    !--------------------------------------------------------!
    ! Given string "$n" return icol=n (where n is a positive !
    ! integer), or icol=0 if there is an error.              !
    !--------------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: icol
    INTEGER ierr
    icol=0
    if(string(1:1)/='$')return
    icol=parse_int(string(2:),ierr)
    if(ierr/=0)then
      icol=0
      return
    endif
    if(icol<1)icol=0
  END SUBROUTINE parse_colnum


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


  SUBROUTINE msg(messg)
    !----------------------------------------!
    ! Write a single-line message to stdout, !
    ! followed by a blank line.              !
    !----------------------------------------!
    IMPlICIT NONE
    CHARACTER(*), INTENT(in), OPTIONAL :: messg
    write(6,'(a)')messg
    write(6,'()')
  END SUBROUTINE msg


  ! RUN CONTROL UTILITIES.


  SUBROUTINE init_random(irandom)
    !-----------------------------------------------!
    ! Initialize built-in random seed so that       !
    ! different processess use different sequences. !
    !-----------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: irandom
    INTEGER i,j,k,l,n,iskip
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER, PARAMETER :: NSKIP=100
    DOUBLE PRECISION t1
    call random_seed(size=n)
    allocate(seed(n))
    do i=1,n
      if(irandom<1)then
        call system_clock(j,k,l)
        seed(i)=0*n+j
      else
        seed(i)=((irandom-1)*1+0)*n+i
      endif
    enddo ! i
    call random_seed(put=seed)
    ! Skip a few random numbers (ifort has trouble with first).
    do iskip=1,NSKIP
      call random_number(t1)
    enddo ! iskip
    deallocate(seed)
  END SUBROUTINE init_random


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
