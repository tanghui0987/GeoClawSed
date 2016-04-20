! ============================================================================
!  File:        sediment_module.f90
!  Program:     sediment_module
!  Author:      Hui Tang tanghui@vt.edu
!  This module is developed based on previous work from Kyle Mandli and 
!  David George. This module is designed to read sediment data and setup sediment
!  condition to Geoclaw
! ============================================================================
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================
! ============================================================================
! Read sediment files as specified in sediment.data
! including two files:
! sediment thickness file and grainsize distribution file
! Each sediment thickness file has a type stored in sedtype(i).
!   sedtype = 1:  standard GIS format: 3 columns: lon,lat,thickness(m)
!   sedtype = 2:  Header as in DEM file, thickness(m) one value per line
!   sedtype = 3:  Header as in DEM file, thickness(m) one row per line
! For other formats modify readtopo routine.
! Grainsize distribution file is similar:
!   sedtype = 1:  standard GIS format: 2+gmax columns: lon,lat,P1...Pgmax
!   sedtype = 2:  Header as in DEM file, percentages for one location per line
!   sedtype = 3:  Header as in DEM file, percentages for one location per line
! advancing northwest to northeast then from north to south. Values should
! be uniformly spaced.
!
! Sediment data will update during computatio
! =============================================================================
    module sediment_module

        use Set_Precision, only: Prec
        use set_variable
        use amr_module, only: tstart_thisrun,node,ndihi,ndilo,ndjhi,ndjlo
        use topo_module, only: rectintegral, intersection

        implicit none
        save
        integer, parameter :: SED_PARM_UNIT = 78

        ! Work array for Sediment for all t
        real(kind=8), allocatable :: sedtwork(:),sedpwork(:,:)

        ! Sediment file data
        integer :: test_sediment !don't need in this case actually, only for test
        character(len=150), allocatable :: sedtfname(:),sedpfname(:)
        integer :: mtsedfiles,mtsedsize,mpsedsize,mpsedfiles
        real(kind=Prec), allocatable :: xlowsed(:), ylowsed(:), tlowsed(:)
        real(kind=Prec), allocatable :: xhised(:), yhised(:), thised(:)
        real(kind=Prec), allocatable :: dxsed(:), dysed(:)
        real(kind=Prec), allocatable :: sedtime(:)
        integer, allocatable ::  mxsed(:), mysed(:)
        real(kind=Prec), allocatable :: totalthick_temp(:),pbbed_temp(:,:)
        integer, allocatable :: i0sed(:), msed(:), msedorder(:)
        integer, allocatable :: minlevelsed(:), maxlevelsed(:), itsedtype(:),ipsedtype(:)
        integer, allocatable :: sedID(:),sed0save(:)
        logical :: sed_finalized ! don't need here 

        ! sediment transport support
        integer :: isedtran, aux_finalized

        !Analytic sediment add later

        ! Initial sediment
        ! Work array for initial topography (only arrays where topo evolves)
        real(kind=Prec), allocatable :: sed0twork(:),sed0pwork(:,:)
        integer, allocatable :: i0sed0(:),sed0ID(:)
        integer :: mtsed0size,mtsed0files,mpsed0size,mpsed0files


        ! ========================================================================
        !  Constants and parameter
        ! ========================================================================
        Real(kind=Prec) ::      rhos,rho,por,cmax,facDc,hcr,thick
        integer :: nd_var
        Integer :: gmax,lmax,mgc
        Real(kind=Prec), dimension(:),allocatable :: D
        Real(kind=Prec) ::      g,vis,Te,Trep,eps,k0,k1,m0,tsfac,Tsmin,smax,cf
        Real(kind=Prec) ::      nuh,nuhfac,gammaWs,a1,hswitch,wetslp,dryslp,toler
        Real(kind=Prec) ::      facsl
        Integer :: sourcesink,struct,morfac,sws,avalanching,thetanum
        Real(kind=Prec) ::      morstart
        Real(kind=Prec) ::       vareps,split,merge,beta
        CHARACTER(120)  ::      limit_method,trim,method
        Logical         ::      aval
        integer :: nyc = 50,nxc = 50,ngc

        contains

            subroutine read_sed_settings(file_name1,file_name2)

                use geoclaw_module
                use amr_module, only: xlower,xupper,ylower,yupper
                use Set_Precision, only: Prec


                implicit none

                !Input arguments
                character(len=*), intent(in), optional :: file_name1,file_name2

                !locals
                integer, parameter :: iunit = 7,SED_PARM_UNIT=78
                integer :: i,j,k,ised,finer_than,rank,nl
                real(kind=Prec) :: area_i,area_j,x_junk,y_junk
                real(kind=Prec) :: area, area_domain

                open(unit=SED_PARM_UNIT,file='fort.geo',status="unknown",action="write")
                ! Open and begin parameter file output
                write(SED_PARM_UNIT,*) ' '
                write(SED_PARM_UNIT,*) '--------------------------------------------'
                write(SED_PARM_UNIT,*) 'SETSED:'
                write(SED_PARM_UNIT,*) '---------'

                if (present(file_name1)) then
                    call opendatafile(iunit, file_name1)
                else
                    call opendatafile(iunit, 'thick.data')
                endif

                ! Read in sediment specification type
                read(iunit,"(i1)") test_sediment
                ! Primary sediment type, read in sediment files specified
                if (test_sediment == 0) then
                    read(iunit,*) mtsedfiles
                    if (mtsedfiles == 0) then
                        write(SED_PARM_UNIT,*) '   mtsedfiles = 0'
                        write(SED_PARM_UNIT,*) '   No sediment thickness files specified, '
                        write(SED_PARM_UNIT,*) '          will set totalthick(x,y) = 0 in setaux'
                        return
                    endif
                    write(SED_PARM_UNIT,*) '   msedfiles = ',mtsedfiles

                    ! Read and allocate data parameters for each file, as thickness file and grainsize file is same, we use one to put, two file share the same data may not good
                    allocate(mxsed(mtsedfiles),mysed(mtsedfiles))
                    allocate(xlowsed(mtsedfiles),ylowsed(mtsedfiles))
                    allocate(tlowsed(mtsedfiles),xhised(mtsedfiles),yhised(mtsedfiles))
                    allocate(thised(mtsedfiles),dxsed(mtsedfiles),dysed(mtsedfiles))
                    allocate(sedtfname(mtsedfiles),itsedtype(mtsedfiles))
                    allocate(minlevelsed(mtsedfiles),maxlevelsed(mtsedfiles))
                    allocate(i0sed(mtsedfiles),msed(mtsedfiles),msedorder(mtsedfiles))
                    allocate(sedID(mtsedfiles),sedtime(mtsedfiles),sed0save(mtsedfiles))
                    allocate(i0sed0(mtsedfiles),sed0ID(mtsedfiles))



                    do i=1,mtsedfiles
                        read(iunit,*) sedtfname(i)
                        read(iunit,*) itsedtype(i),minlevelsed(i), maxlevelsed(i), &
                            tlowsed(i),thised(i)

                        write(SED_PARM_UNIT,*) '   '
                        write(SED_PARM_UNIT,*) '   ',sedtfname(i)
                        write(SED_PARM_UNIT,*) '  itsedtype = ', itsedtype(i)
                        write(SED_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                            minlevelsed(i), maxlevelsed(i)
                        write(SED_PARM_UNIT,*) '  tlow, thi = ', tlowsed(i),thised(i)
                        call read_sed_header(sedtfname(i),itsedtype(i),mxsed(i), &
                                mysed(i),xlowsed(i),ylowsed(i),xhised(i),yhised(i), &
                                dxsed(i),dysed(i))
                        sedID(i) = i
                        msed(i) = mxsed(i)*mysed(i)
                    enddo
                    ! Indexing into work array
                    i0sed(1)=1
                    if (mtsedfiles > 1) then
                        do i=2,mtsedfiles
                            i0sed(i)=i0sed(i-1) + msed(i-1)
                        enddo
                    endif

                    ! Read sediment information and allocate space for each file
                    mtsedsize = sum(msed)
                    allocate(sedtwork(mtsedsize))
                    do i=1,mtsedfiles
                        sedID(i) = i
                        sedtime(i) = -huge(1.0)
                        allocate(totalthick_temp(1:mxsed(i)*mysed(i)),pbbed_temp(1:mxsed(i)*mysed(i),gmax))
                        call read_tsed_file(mxsed(i),mysed(i),itsedtype(i),sedtfname(i), &
                            sedtwork(i0sed(i):i0sed(i)+msed(i)-1))
                        allocate(totalthick(1-mgc:mxsed(i)+mgc,1-mgc:mysed(i)+mgc))
                        !do j = 1, mysed(i)
                            !do k = 1, mxsed(i)
                                !totalthick(i,k) = totalthick_temp(i+k*mysed(i))
                                !if (totalthick(i,k)>thick) then
                                    !totalnum(i,k) = nint(totalthick(i,k)/thick)
                                    !if (totalthick(i,k)>totalnum(i,k)*thick) then
                                        !totalnum(i,k)=totalnum(i,k)+1
                                        !dzbed(i,k,1:int(totalnum(i,k)-1))= thick
                                        !dzbed(i,k,totalnum(i,k)) =totalthick(i,k)-(totalnum(i,k)-1)*thick
                                    !else
                                        !dzbed(i,k,1:totalnum(i,k)-1)= thick
                                        !dzbed(i,k,totalnum(i,k)) =totalthick(i,k)-(totalnum(i,k)-1)*thick
                                    !endif
                                !else
                                    !totalnum(i,k) = 1
                                    !dzbed(i,k,1) = totalthick(i,k)
                                !endif
                                !if(dzbed(i,j,totalnum(i,j))<toler*thick) then
                                    !totalnum(i,j) = totalnum(i,j) - 1
                                    !dzbed(i,j,totalnum(i,j)) = 0.0
                                !endif
                            !enddo
                        !enddo
                    enddo


                    ! Sediment order, which determines theorder to process sediment data
                    !
                    ! The finest one will be given priority in any region
                    ! msedorder(rank) = i means that i'th sediment file has rank,
                    ! where the file with rank=1 is the finest and considered first.
                    ! we check for both sediment thickness profile and grain size distribution file
                    do i=1,mtsedfiles
                        finer_than = 0
                        do j=1,mtsedfiles
                            if (j /= i) then
                                area_i=dxsed(i)*dysed(i)
                                area_j=dxsed(j)*dysed(j)
                                if (area_i < area_j) finer_than = finer_than + 1
                                ! if two files have the same resolution, order is
                                ! arbitrarily chosen
                                if ((area_i == area_j).and.(j < i)) then
                                    finer_than = finer_than + 1
                                endif
                            endif
                        enddo
                        ! ifinerthan tells how many other files, file i is finer than
                        rank = mtsedfiles - finer_than
                        msedorder(rank) = i
                    enddo

                    write(SED_PARM_UNIT,*) ' '
                    write(SED_PARM_UNIT,*) '  Ranking of sediment files', &
                                            '  finest to coarsest: ', &
                                (msedorder(rank),rank=1,mtsedfiles)
                    write(SED_PARM_UNIT,*) ' '
                    i0sed0(1) = 1
                    mtsed0size = dot_product(msed,sed0save)
                    allocate(sed0twork(mtsed0size))
                    do i = 2,mtsedfiles
                        i0sed0(i)= i0sed0(i-1) + msed(i-1)*sed0save(i-1)
                    enddo

                    do i = 1,mtsedfiles
                        if (sed0save(i)>0) then
                            sed0twork(i0sed0(i):i0sed0(i)+msed(i)-1) = &
                            sedtwork(i0sed(i):i0sed(i)+msed(i)-1)
                        endif
                    enddo
                endif
                if (present(file_name2)) then
                    call opendatafile(iunit, file_name2)
                else
                    call opendatafile(iunit, 'perc.data')
                endif
                ! Read in sediment specification type
                read(iunit,"(i1)") test_sediment
                ! Primary sediment type, read in sediment files specified
                if (test_sediment == 0) then
                    read(iunit,*) mpsedfiles

                    if (mpsedfiles == 0) then
                        write(SED_PARM_UNIT,*) '   mpsedfiles = 0'
                        write(SED_PARM_UNIT,*) '   No sediment grainsize distribution files specified, '
                        write(SED_PARM_UNIT,*) '          will set pebbed(x,y,n) = 1/gmax in setaux'
                        return
                    endif
                    if (mtsedfiles .ne. mtsedfiles) then
                        write(SED_PARM_UNIT,*) 'Warning:'
                        write(SED_PARM_UNIT,*) '   The number of thickness and grainsize distribution file is not same! '
                        write(SED_PARM_UNIT,*) '   Check your data! '
                        return
                    endif

                    write(SED_PARM_UNIT,*) '   psedfiles = ',mpsedfiles

                    allocate(sedpfname(mpsedfiles),ipsedtype(mpsedfiles)) !add to file

                    do i=1,mpsedfiles
                        read(iunit,*) sedpfname(i)
                        read(iunit,*) ipsedtype(i),minlevelsed(i), maxlevelsed(i), &
                            tlowsed(i),thised(i)
                        write(SED_PARM_UNIT,*) '   '
                        write(SED_PARM_UNIT,*) '   ',sedpfname(i)
                        write(SED_PARM_UNIT,*) '  ipsedtype = ', ipsedtype(i)
                        write(SED_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                            minlevelsed(i), maxlevelsed(i)
                        write(SED_PARM_UNIT,*) '  tlow, thi = ', tlowsed(i),thised(i)
                        call read_per_header(sedpfname(i),ipsedtype(i),mxsed(i), &
                            mysed(i),xlowsed(i),ylowsed(i),xhised(i),yhised(i), &
                            dxsed(i),dysed(i),gmax)
                    enddo
                    mpsedsize = sum(msed)
                    allocate(sedpwork(mpsedsize,gmax))


                    do i=1,mpsedfiles
                        sedID(i) = i
                        sedtime(i) = -huge(1.0)
                        allocate(pbbed(1-mgc:mxsed(i),1-mgc:mysed(i),lmax,gmax))
                        call read_psed_file(mxsed(i),mysed(i),ipsedtype(i),sedpfname(i), &
                            pbbed_temp(i0sed(i):i0sed(i)+msed(i)-1,:))
                        do j = 1, mysed(i)
                            do k = 1, mxsed(i)
                                do nl = 1, lmax
                                    pbbed(k,j,nl,:) = pbbed_temp(i+k*mysed(i)-1,:) !set uniform grain size distribution
                                enddo
                            enddo
                        enddo
                        pbbed(0,:,:,:) = pbbed(1,:,:,:)
                        pbbed(-1,:,:,:) = pbbed(0,:,:,:)
                        pbbed(:,0,:,:) = pbbed(:,1,:,:)
                        pbbed(:,-1,:,:) = pbbed(:,0,:,:)
                        pbbed(mxsed(i),:,:,:) = pbbed(mxsed(i)-1,:,:,:)
                        pbbed(mxsed(i)+1,:,:,:) = pbbed(mxsed(i),:,:,:)
                        pbbed(:,mysed(i),:,:) = pbbed(:,mysed(i)-1,:,:)
                        pbbed(:,mysed(i)+1,:,:) = pbbed(:,mysed(i),:,:)
                    enddo
                    do i=1,mpsedfiles
                        finer_than = 0
                        do j=1,mpsedfiles
                            if (j /= i) then
                                area_i=dxsed(i)*dysed(i)
                                area_j=dxsed(j)*dysed(j)
                                if (area_i < area_j) finer_than = finer_than + 1
                                        ! if two files have the same resolution, order is
                                        ! arbitrarily chosen
                                if ((area_i == area_j).and.(j < i)) then
                                        finer_than = finer_than + 1
                                endif
                            endif
                        enddo
                        ! ifinerthan tells how many other files, file i is finer than
                        rank = mpsedfiles - finer_than
                        msedorder(rank) = i
                    enddo

                    write(SED_PARM_UNIT,*) ' '
                    write(SED_PARM_UNIT,*) '  Ranking of sediment grain size distribution files', &
                        '  finest to coarsest: ', &
                                (msedorder(rank),rank=1,mpsedfiles)
                    write(SED_PARM_UNIT,*) ' '
                    i0sed0(1) = 1
                    mpsed0size = dot_product(msed,sed0save)
                    allocate(sed0pwork(mpsed0size,gmax))
                    do i = 2,mpsedfiles
                        i0sed0(i)= i0sed0(i-1) + msed(i-1)*sed0save(i-1)
                    enddo

                    do i = 1,mpsedfiles
                        if (sed0save(i)>0) then
                            sed0pwork(i0sed0(i):i0sed0(i)+msed(i)-1,:) = &
                            sedpwork(i0sed(i):i0sed(i)+msed(i)-1,:)
                        endif
                    enddo


                    ! Check that topo arrays cover full domain: need to check the area but not now!
                    call sedarea(xlower,xupper,ylower,yupper,1,area,mtsedfiles,mpsedfiles)
                    area_domain = (yupper-ylower)*(xupper-xlower)
                    if (abs(area - area_domain) > 1e-12*area_domain) then
                        write(6,*) '**** sediment arrays do not cover domain'
                        write(6,*) '**** area of overlap = ', area
                        write(6,*) '**** area of domain  = ', area_domain
                        stop
                    endif
                end if
            end subroutine read_sed_settings
            ! ========================================================================
            !  read_tsed_file(mx,my,sed_type,fname,totalthick)
            !
            !  Read sediment thicness file.
            ! ========================================================================

            subroutine read_tsed_file(mx,my,sed_type,fname,totalthick_temp)

                use geoclaw_module
                use Set_Precision, only: Prec

                implicit none

                ! Arguments
                integer, intent(in) :: mx,my,sed_type
                character(len=150), intent(in) :: fname
                real(kind=Prec), intent(inout) :: totalthick_temp(1:mx*my)

                ! Locals
                integer, parameter :: iunit = 19, miss_unit = 17
                real(kind=Prec), parameter :: sed_missing = -150.d0
                logical, parameter :: maketype2 = .false.
                integer :: i,j,num_points,missing,status,sed_start
                real(kind=Prec) :: no_data_value,x,y,z,th_temp

                print *, ' '
                print *, 'Reading sediment thickness file  ', fname

                open(unit=iunit, file=fname, status='unknown',form='formatted')
                select case(abs(sed_type))
                ! ASCII file with x,y,z values on each line.
                ! (progressing from upper left corner across rows, then down)
                ! Assumes a uniform rectangular grid of data values.
                    case(1)
                        i = 0
                        status = 0
                        do i=1,mx*my
                            read(iunit,fmt=*,iostat=status) x,y,th_temp
                            if ((i > mx * my) .and. (status == 0)) then
                                print *,'*** Error: i > mx*my = ',mx*my
                                print *,'*** i, mx, my: ',i,mx,my
                                print *,'*** status = ',status
                                stop
                            endif

                            if (status /= 0) then
                                print *,"Error reading sediment thickness file, reached EOF."
                                print *,"  File = ",fname
                                stop
                            else
                                totalthick_temp(i) = th_temp
                            endif
                        enddo

                    ! ================================================================
                    ! ASCII file with header followed by thickness data
                    ! (progressing from upper left corner across rows, then down)
                    ! one value per line if sed_type=2 or
                    ! mx values per line if sed_type=3
                    ! ================================================================
                    case(2:3)
                    ! Read header
                        do i=1,5
                            read(iunit,*)
                        enddo
                        read(iunit,*) no_data_value

                        ! Read in data
                        missing = 0
                        select case(abs(sed_type))
                            case(2)
                                do i=1,mx*my
                                    read(iunit,*) totalthick_temp(i)
                                    if (totalthick_temp(i) == no_data_value) then
                                        missing = missing + 1
                                        totalthick_temp(i) = sed_missing
                                    endif
                                enddo
                            case(3)
                                do j=1,my
                                    read(iunit,*) (totalthick_temp((j-1)*mx + i),i=1,mx)
                                    do i=1,mx
                                        if (totalthick_temp((j-1)*mx + i) == no_data_value) then
                                            missing = missing + 1
                                            totalthick_temp((j-1)*mx + i) = sed_missing
                                        endif
                                    enddo
                                enddo
                        end select

                        ! Write a warning if we found and missing values
                        if (missing > 0)  then
                            print *, '   WARNING...some missing data values this file'
                            print *, '       ',missing,' missing data values'
                            print *, '              (see fort.missing)'
                            print *, '   These values have arbitrarily been set to ',&
                                    sed_missing
                        endif
                    end select

                    close(unit=iunit)

            end subroutine read_tsed_file


            ! ========================================================================
            !  read_psed_file(mx,my,sed_type,fname,totalthick)
            !
            !  Read sediment grain size distribution file.
            ! ========================================================================

            subroutine read_psed_file(mx,my,sed_type,fname,pbbed_temp)

                use geoclaw_module
                use Set_Precision, only: Prec

                implicit none

                ! Arguments
                integer, intent(in) :: mx,my,sed_type
                character(len=150), intent(in) :: fname
                real(kind=Prec), intent(inout) :: pbbed_temp(1:mx*my,gmax)

                ! Locals
                integer, parameter :: iunit = 19, miss_unit = 17,punit = 79
                real(kind=Prec), parameter :: sed_missing = -150.d0
                logical, parameter :: maketype2 = .false.
                integer :: i,j,k,num_points,missing,status,sed_start
                real(kind=Prec) :: no_data_value,x,y,z
                real(kind=Prec) :: p_temp(gmax)

                print *, ' '
                print *, 'Reading sediment grain size distribution file  ', fname

                open(unit=iunit, file=fname, status='unknown',form='formatted')

                select case(abs(sed_type))
                    ! ASCII file with x,y,z values on each line.
                    ! (progressing from upper left corner across rows, then down)
                    ! Assumes a uniform rectangular grid of data values.
                    case(1)
                        i = 0
                        status = 0
                        do i=1,mx*my
                            read(iunit,fmt=*,iostat=status) x,y
                            read(punit,fmt=*,iostat=status) p_temp !may have problem
                            if ((i > mx * my) .and. (status == 0)) then
                                print *,'*** Error: i > mx*my = ',mx*my
                                print *,'*** i, mx, my: ',i,mx,my
                                print *,'*** status = ',status
                                stop
                            endif

                            if (status /= 0) then
                                print *,"Error reading sediment thickness file, reached EOF."
                                print *,"  File = ",fname
                                stop
                            else
                                pbbed_temp(i,:) = p_temp
                            endif
                        enddo

                        ! ================================================================
                        ! ASCII file with header followed by thickness data
                        ! (progressing from upper left corner across rows, then down)
                        ! one value per line if topo_type=2 or
                        ! mx values per line if topo_type=3
                        ! ================================================================
                        case(2:3)
                            ! Read header
                            do i=1,6
                                read(iunit,*)
                            enddo
                            read(iunit,*) no_data_value

                            ! Read in data
                            missing = 0
                            select case(abs(sed_type))
                                case(2)
                                    do i=1,mx*my
                                        read(iunit,*) pbbed_temp(i,:)!1),pbbed_temp(i,2)
                                        do j = 1, gmax
                                            if (pbbed_temp(i,j) == no_data_value) then
                                                missing = missing + 1
                                                pbbed_temp(i,j) = sed_missing
                                            endif
                                        enddo
                                    enddo
                                case(3)
                                    do j=1,my
                                        read(punit,*) (pbbed_temp((j-1)*mx + i,:),i=1,mx)
                                        do i=1,mx
                                            do k = 1, gmax
                                                if (pbbed_temp((j-1)*mx + i,k) == no_data_value) then
                                                    missing = missing + 1
                                                    pbbed_temp((j-1)*mx + i,k) = sed_missing
                                                endif
                                            enddo
                                        enddo
                                    enddo
                            end select

                            ! Write a warning if we found and missing values
                        if (missing > 0)  then
                            print *, '   WARNING...some missing data values this file'
                            print *, '       ',missing,' missing data values'
                            print *, '              (see fort.missing)'
                            print *, '   These values have arbitrarily been set to ',&
                                            sed_missing
                        endif
                    end select
                    close(unit=iunit)

            end subroutine read_psed_file


            ! ========================================================================
            ! subroutine read_sed_header(fname,sed_type,mx,my,xll,yll,xhi,yhi,dx,dy)
            ! ========================================================================
            !  Read sediment file header to determine space needed in allocatable array
            !
            !  :Input:
            !   - fname - (char) Name of file
            !   - sed_type - (int) Type of sediment file (-3 < topo_type < 3)
            !
            !  :Output:
            !   - mx,my - (int) Number of grid points
            !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
            !   - dx,dy - (float) Spatial resolution of grid
            ! ========================================================================
            subroutine read_sed_header(fname,sed_type,mx,my,xll,yll,xhi,yhi,dx,dy)

                use geoclaw_module

                implicit none

                ! Input and Output
                character(len=150), intent(in) :: fname
                integer, intent(in) :: sed_type
                integer, intent(out) :: mx,my
                real(kind=8), intent(out) :: xll,yll,xhi,yhi,dx,dy

                ! Local
                integer, parameter :: iunit = 19
                integer :: sed_size, status
                real(kind=8) :: x,y,th,no_data_value
                logical :: found_file

                inquire(file=fname,exist=found_file)
                if (.not. found_file) then
                    print *, 'Missing sediment file:'
                    print *, '   ', fname
                    stop
                endif

                open(unit=iunit, file=fname, status='unknown',form='formatted')


                select case(abs(sed_type))
                    ! ASCII file with 3 columns
                    ! determine data size
                    case(1)
                    ! Initial size variables
                        sed_size = 0
                        mx = 0

                        ! Read in first values, determines xlow and yhi
                        read(iunit,*) xll,yhi
                        sed_size = sed_size + 1
                        mx = mx + 1

                        ! Go through first row figuring out mx, continue to count
                        y = yhi
                        do while (yhi == y)
                            read(iunit,*) x,y,th
                            sed_size = sed_size + 1
                            mx = mx + 1
                        enddo
                        mx = mx - 1
                        ! Continue to count the rest of the lines
                        do
                            read(iunit,fmt=*,iostat=status) x,y,th
                            if (status /= 0) exit
                            sed_size = sed_size + 1
                        enddo
                        if (status > 0) then
                            print *,"ERROR:  Error reading header of sediment file ",fname
                            stop
                        endif

                        ! Calculate remaining values
                        my = sed_size / mx
                        xhi = x
                        yll = y
                        dx = (xhi-xll) / (mx-1)
                        dy = (yhi-yll) / (my-1)

                        ! ASCII file with header followed by z data
                    case(2:3)
                        read(iunit,*) mx
                        read(iunit,*) my
                        read(iunit,*) xll
                        read(iunit,*) yll
                        read(iunit,*) dx
                        read(iunit,*) no_data_value
                        dy = dx
                        xhi = xll + (mx-1)*dx
                        yhi = yll + (my-1)*dy

                    case default
                        print *, 'ERROR:  Unrecognized sed_type'
                        print *, '    sed_type = ',sed_type
                        print *, '  for sediment file:'
                        print *, '   ', fname
                        stop
                end select

                close(iunit)

                write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
                write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
                write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

            end subroutine read_sed_header

    ! ========================================================================
    ! subroutine read_sed_header(fname,sed_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read sediment file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - sed_type - (int) Type of sediment file (-3 < topo_type < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================

            subroutine read_per_header(fname,per_type,mx,my,xll,yll,xhi,yhi,dx,dy,mclasses)

                use geoclaw_module

                implicit none

                ! Input and Output
                character(len=150), intent(in) :: fname
                integer, intent(in) :: per_type
                integer, intent(out) :: mx,my,mclasses
                real(kind=8), intent(out) :: xll,yll,xhi,yhi,dx,dy

                ! Local
                integer, parameter :: iunit = 19
                integer :: sed_size, status
                real(kind=8) :: x,y,th,no_data_value
                logical :: found_file

                inquire(file=fname,exist=found_file)
                if (.not. found_file) then
                    print *, 'Missing sediment file:'
                    print *, '   ', fname
                    stop
                endif

                open(unit=iunit, file=fname, status='unknown',form='formatted')


                select case(abs(per_type))
                    ! ASCII file with 3 columns
                    ! determine data size
                    case(1)
                        ! Initial size variables
                        sed_size = 0
                        mx = 0

                        ! Read in first values, determines xlow and yhi
                        read(iunit,*) xll,yhi
                        sed_size = sed_size + 1
                        mx = mx + 1

                        ! Go through first row figuring out mx, continue to count
                        y = yhi
                        do while (yhi == y)
                            read(iunit,*) x,y
                            sed_size = sed_size + 1
                            mx = mx + 1
                        enddo
                        mx = mx - 1
                        ! Continue to count the rest of the lines
                        do
                            read(iunit,fmt=*,iostat=status) x,y
                            if (status /= 0) exit
                                sed_size = sed_size + 1
                        enddo
                        if (status > 0) then
                            print *,"ERROR:  Error reading header of sediment file ",fname
                            stop
                        endif

                        ! Calculate remaining values
                        my = sed_size / mx
                        xhi = x
                        yll = y
                        dx = (xhi-xll) / (mx-1)
                        dy = (yhi-yll) / (my-1)

                        ! ASCII file with header followed by z data
                    case(2:3)
                        read(iunit,*) mx
                        read(iunit,*) my
                        read(iunit,*) mclasses
                        read(iunit,*) xll
                        read(iunit,*) yll
                        read(iunit,*) dx
                        read(iunit,*) no_data_value
                        dy = dx
                        xhi = xll + (mx-1)*dx
                        yhi = yll + (my-1)*dy

                    case default
                        print *, 'ERROR:  Unrecognized sed_type'
                        print *, '    sed_type = ',per_type
                        print *, '  for sediment file:'
                        print *, '   ', fname
                        stop
                end select

                close(iunit)

                write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
                write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
                write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

            end subroutine read_per_header

            ! ========================================================================
            !  set_geo(fname)
            !  Reads in user parameters from the given file name if provided
            subroutine set_sediment(file_name)

                use geoclaw_module, only : pi
                !use set_variable

                implicit none

                ! Input
                character(len=*), intent(in), optional :: file_name
                integer, parameter :: unit = 19

                !nxc      = 50!node(ndihi,2) - node(ndilo,2) + 1
                !nyc      = 50!node(ndjhi,2) - node(ndjlo,2) + 1
                !print *, node
                open(unit=SED_PARM_UNIT,file='fort.geo',status="unknown",action="write")
                call test
                write(SED_PARM_UNIT,*) ' '
                write(SED_PARM_UNIT,*) '--------------------------------------------'
                write(SED_PARM_UNIT,*) 'Sediment Parameters:'
                write(SED_PARM_UNIT,*) '-------------------'

                ! Read user parameters from setsed.data or sediment.data, see detail in these documents
                if (present(file_name)) then
                    call opendatafile(unit, file_name)
                else
                    call opendatafile(unit, 'sediment.data')
                endif
                ! Sediment parameters
                read(unit,"(1d16.8)") rhos
                read(unit,"(1d16.8)") rho
                read(unit,"(1d16.8)") por
                read(unit,"(1d16.8)") cmax
                read(unit,"(1d16.8)") facDc
                read(unit,"(1d16.8)") thick
                read(unit,*) nd_var
                read(unit,*) gmax
                read(unit,*) lmax
                read(unit,*) mgc
                read(unit,*) hcr
                if (gmax>0) then
                    allocate(D(gmax))
                    read(unit,*) D(:)
                endif
                !physical parameter
                read(unit,*) g
                read(unit,*) vis
                read(unit,*) Te
                read(unit,*) Trep
                read(unit,*) eps
                read(unit,*) k0
                read(unit,*) m0
                read(unit,*) tsfac
                read(unit,*) Tsmin
                read(unit,*) smax
                read(unit,*) cf
                read(unit,*) nuh
                read(unit,*) nuhfac
                read(unit,*) hswitch
                read(unit,*) wetslp
                read(unit,*) dryslp
                read(unit,*) toler
                read(unit,*) facsl
                !processes control paramter 
                read(unit,*) sourcesink
                read(unit,*) aval
                read(unit,*) avalanching
                read(unit,*) struct
                read(unit,*) morfac
                read(unit,*) thetanum
                read(unit,*) sws
                read(unit,*) morstart
                read(unit,*) split
                read(unit,*) merge
                !algorithm parameter
                read(unit,*) vareps
                read(unit,*) k1
                !method
                read(unit,*) limit_method
                read(unit,*) trim
                read(unit,*) method
                close(unit)
                if (gmax<1) then
                    print *, 'ERROR in setsed:  must set up grain size for sediment'
                    stop
                endif
                write(SED_PARM_UNIT,*) '   Sediment density:',rhos
                write(SED_PARM_UNIT,*) '   Water density:',rho
                write(SED_PARM_UNIT,*) '   Porosity:',por
                write(SED_PARM_UNIT,*) '   Maximum allowed sediment concentration:',cmax
                write(SED_PARM_UNIT,*) '   Control sediment diffusion coefficient:',facDc
                write(SED_PARM_UNIT,*) '   Sediment thickness for each layer:', thick
                write(SED_PARM_UNIT,*) '   Initial active sediment layer:', nd_var
                write(SED_PARM_UNIT,*) '   Number of Grain size classes:',gmax
                write(SED_PARM_UNIT,*) '   Number of sediment layers:',lmax
                write(SED_PARM_UNIT,*) '   Number of ghost cell:',mgc
                write(SED_PARM_UNIT,*) '   Water depth consider sediment transport:',hcr
                if (gmax>0) then
                    write(SED_PARM_UNIT,*) '   Sediment grain size classes:', D(:)
                endif
                write(SED_PARM_UNIT,*) '   gravity:',g
                write(SED_PARM_UNIT,*) '   kinematic viscosity:',vis
                write(SED_PARM_UNIT,*) '   Water temperature:',Te
                write(SED_PARM_UNIT,*) '   Representative wave period:',Trep
                write(SED_PARM_UNIT,*) '   Threshold water depth above which cells are considered wet:', eps
                write(SED_PARM_UNIT,*) '   von kaman coefficient:',k0
                write(SED_PARM_UNIT,*) '   mining coeffcient:',m0
                write(SED_PARM_UNIT,*) '   Coefficient determining Ts = tsfac * h/ws in sediment source term:',tsfac
                write(SED_PARM_UNIT,*) '   Minimum adaptation time scale in advection diffusion equation sediment:', Tsmin
                write(SED_PARM_UNIT,*) '   maximum Shields parameter for ceq Diane Foster:',smax
                write(SED_PARM_UNIT,*) '   Friction coefficient flow:',cf
                write(SED_PARM_UNIT,*) '   horizontal background viscosity:',nuh
                write(SED_PARM_UNIT,*) '   Viscosity switch for roller induced turbulent horizontal viscosity:',nuhfac
                write(SED_PARM_UNIT,*) '   Water depth at which is switched from wetslp to dryslp',hswitch
                write(SED_PARM_UNIT,*) '   Critical avalanching slope under water:',wetslp
                write(SED_PARM_UNIT,*) '   Critical avalanching slope above water:',dryslp
                write(SED_PARM_UNIT,*) '   toler for sediment flux limitor:',toler
                write(SED_PARM_UNIT,*) '   Slope affect factor:',facsl
                write(SED_PARM_UNIT,*) '    use source-sink terms to calculate bed level change:',sourcesink
                write(SED_PARM_UNIT,*) '    have avalanched or not:',aval
                write(SED_PARM_UNIT,*) '    Include avalanching:',avalanching
                write(SED_PARM_UNIT,*) '    Switch for hard structures:',struct
                write(SED_PARM_UNIT,*) '    morphological acceleration factor:',morfac
                write(SED_PARM_UNIT,*) '    Coefficient determining scheme:',thetanum
                write(SED_PARM_UNIT,*) '    short wave & roller stirring and undertow:',sws
                write(SED_PARM_UNIT,*) '    Start time:',morstart
                write(SED_PARM_UNIT,*) '    Split threshold for variable sediment layer:',split
                write(SED_PARM_UNIT,*) '    Mergethreshold for variable sediment layer ',merge
                write(SED_PARM_UNIT,*) '    Coefficient determining order of accuracy:',vareps
                write(SED_PARM_UNIT,*) '    Coefficient determining whether fully upwind',k1
                write(SED_PARM_UNIT,*) '    method to caculate sediment flux limiter',limit_method
                write(SED_PARM_UNIT,*) '    method to caculate equibrium sediment concentration',trim
                write(SED_PARM_UNIT,*) '    method to caculate sediment flux',method
            end subroutine set_sediment

            recursive subroutine sedarea(x1,x2,y1,y2,m,area,mtsedfiles,mpsedfiles)

                use Set_Precision, only: Prec

            ! Compute the area of overlap of sediment file  with the rectangle (x1,x2) x (y1,y2)
            ! using sediment arrays indexed mtsedorder(mtopofiles) through mtsedorder(m)
            ! (coarse to fine).
            ! The main call to this subroutine has corners of a physical domain for
            ! the rectangle and m = 1 in order to compute the area of overlap of
            ! domain by all topo arrays.  Used to check inputs and insure topo
            ! covers domain.

            ! similar to topoarea

                implicit none

                ! arguments
                real (kind=Prec), intent(in) :: x1,x2,y1,y2
                integer, intent(in) :: m,mtsedfiles,mpsedfiles
                real (kind=Prec), intent(out) :: area

                ! local
                real(kind=Prec) :: xmlo,xmhi,ymlo,ymhi,x1m,x2m, &
                    y1m,y2m, area1,area2,area_m
                integer :: mfid, indicator, i0

                mfid = msedorder(m)
                i0=i0sed(mfid)
                if ((mtsedfiles /= mpsedfiles)) then
                    write(SED_PARM_UNIT,*) 'Warning:'
                    write(SED_PARM_UNIT,*) '   The number of thickness and grainsize distribution file is not same! '
                    write(SED_PARM_UNIT,*) '   Check your data! '
                endif
                if (m == mtsedfiles) then
                    ! innermost step of recursion reaches this point.
                    ! only using coarsest topo grid -- compute directly...
                    call intersection(indicator,area,xmlo,xmhi, &
                            ymlo,ymhi, x1,x2,y1,y2, &
                            xlowsed(mfid),xhised(mfid),ylowsed(mfid),yhised(mfid))

                else
                    ! recursive call to compute area using one fewer topo grids:
                    call sedarea(x1,x2,y1,y2,m+1,area1,mtsedfiles,mpsedfiles)

                    ! region of intersection of cell with new topo grid:
                    call intersection(indicator,area_m,x1m,x2m, &
                        y1m,y2m, x1,x2,y1,y2, &
                        xlowsed(mfid),xhised(mfid),ylowsed(mfid),yhised(mfid))


                    if (area_m > 0) then

                        ! correction to subtract out from previous set of topo grids:
                        call sedarea(x1m,x2m,y1m,y2m,m+1,area2,mtsedfiles,mpsedfiles)

                        ! adjust integral due to corrections for new topo grid:
                        area = area1 - area2 + area_m
                    else
                        area = area1
                    endif
                endif
            end subroutine sedarea

    end module sediment_module



























