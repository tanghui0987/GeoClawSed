!
subroutine setauxnew(meqn,mbc,mx,my,xlow,ylow,dx,dy,dt,maux,aux,time,q)
!     ============================================
!
!     # set auxiliary arrays
!
!     aux(1,i,j) = Z(x,y) topography (negative below sea level for topoymetry)
!
!     If coordinate_system=2 then lat-lon coordinates on the sphere and
!        aux(2,i,j) = area ratio (capacity function -- set mcapa = 2)
!        aux(3,i,j) = length ratio for edge
!
!

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level
    use amr_module, only: mcapa, xupper, yupper, xlower, ylower, t0


    use topo_module

    !use sediment_module, only: pbbed
    use bedupdate_module, only: bed_update
    use sediment_module, only : eps
    use set_variable, only : zb,ero
    implicit none

    ! Arguments
    integer, intent(in) :: mbc,mx,my,maux,meqn
    real(kind=8), intent(in) :: xlow,ylow,dx,dy,dt,time
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    character(len=90) :: filename = "q0.txt"
    
    ! Locals
    integer :: ii,jj,m, iint,jint,i,j
    real(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
    integer :: skipcount,iaux,ilo,jlo
    real(kind=8) :: u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)
    ! Get primative variable
    do i = 1-mbc, mx+mbc
        do j = 1-mbc, my+mbc
            if (q(1,i,j) .lt. eps) then
                h(i,j) = 0.0
                u(i,j) = 0.0
                v(i,j) = 0.0
            else
                h(i,j) = q(1,i,j)
                u(i,j) = q(2,i,j)/h(i,j)
                v(i,j) = q(3,i,j)/h(i,j)
            end if
        enddo
    enddo

    ! Lat-Long coordinate system in use, check input variables
    if (coordinate_system == 2) then
        if (mcapa /= 2 .or. maux < 3) then
            print *,'ERROR in setaux:  for coordinate_system=2'
            print *,'     need mcapa = 2 and maux >= 3'
            print *,'     have mcapa = ',mcapa,'  maux = ',maux
            stop
        endif
    endif

    ! Set default values for aux variables
    !aux(1,:,:) = 0.d0 ! Bathymetry
    aux(2,:,:) = 1.d0 ! Grid cell area
    aux(3,:,:) = 1.d0 ! Length ratio for edge

    ! Set analytical bathymetry here if requested
    if (test_topography > 0) then
        forall (ii=1-mbc:mx+mbc,jj=1-mbc:my+mbc)
            aux(1,ii,jj) = test_topo(xlow + (ii - 0.5d0) * dx,       &
                                     ylow + (jj - 0.5d0) * dy)
        end forall
    endif

! test:  compute integer indices based off same corner of domain 
!        to reduce round off discrepancies
    ilo = floor((xlow - xlower + .05d0*dx)/dx)
    jlo = floor((ylow - ylower + .05d0*dy)/dy)

    ! Set bathymetry
    skipcount = 0

    do jj=1-mbc,my+mbc

        ym = ylower + (jlo+jj-1.d0) * dy
        yp = ylower + (jlo+jj) * dy
        y = 0.5d0*(ym+yp)

        do ii=1-mbc,mx+mbc

            xm = xlower + (ilo+ii-1.d0) * dx
            xp = xlower + (ilo+ii) * dx
            x = 0.5d0*(xm+xp)


            !write(*,444)ii,jj,aux(1,ii,jj)
444         format("in setaux ",2i4,e12.5)

            ! Set lat-long cell info
            if (coordinate_system == 2) then
                aux(2,ii,jj) = deg2rad * earth_radius**2 * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                aux(3,ii,jj) = ym * deg2rad
            endif

            ! skip setting aux(1,ii,jj) in ghost cell if outside physical domain
            ! since topo files may not cover ghost cell, and values
            ! should be extrapolated, which is done in next set of loops.
            if ((y>yupper) .or. (y<ylower) .or. &
                (x>xupper) .or. (x<xlower)) cycle

!           ### parameter NEEDS_TO_BE_SET initialized in amr_module.f90
!           ### saves time by otherwise copying instead of reinitializing
            !if (aux(1,ii,jj) .ne. NEEDS_TO_BE_SET) then
            !   skipcount = skipcount + 1
            !   cycle  ! new system copies bathy where possible
            !endif


            ! Use input topography files if available
            if (mtopofiles > 0 .and. test_topography == 0 .and. time == t0) then

                topo_integral = 0.d0

                call cellgridintegrate(topo_integral,xm,x,xp,ym,y,yp, &
                    xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo, &
                    mxtopo,mytopo,mtopo,i0topo,mtopoorder, &
                    mtopofiles,mtoposize,topowork)

                    aux(1,ii,jj) = topo_integral / (dx * dy * aux(2,ii,jj))
            endif

        enddo
    enddo
    if (mtopofiles > 0 .and. test_topography == 0 .and. time /= t0 ) then

        !call bed_update(mbc,mx,my,time,dt,dx,dy,u,v,h,aux,maux)
        !ero = aux(1,:,:)
        aux(1,:,:)= zb(:,:)
        do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
                if (ero(i,j)- zb(i,j)>0.0) then
                    ero(i,j) = zb(i,j)
                endif
            end do
        end do
        open (unit=1,status='replace',position='Append',file="ero.txt",action="write")
        write (1,*),ero(1,:)
        close(1)
        open (unit=1,status='replace',position='Append',file="ero1.txt",action="write")
        write (1,*),ero(2,:)
        close(1)
        open (unit=1,status='replace',position='Append',file="ero2.txt",action="write")
        write (1,*),ero(3,:)
        close(1)
        open (unit=1,status='replace',position='Append',file="ero3.txt",action="write")
        write (1,*),ero(4,:)
        close(1)
        write (filename, '( "h",I2,".txt" )' ) int(time)
        close(1)
        open (unit=1,status='replace',position='Append',file=filename,action="write")
        write (1,*),q(1,2,:)
        close(1)
        write (filename, '( "u",I2,".txt" )' ) int(time)
        close(1)
        open (unit=1,status='replace',position='Append',file=filename,action="write")
        write (1,*),u(2,:)
        close(1)
    endif
    !write(*,*)" skipcount = ",skipcount
    !print *, q(1,15,10)

    ! Copy topo to ghost cells if outside physical domain

    do jj=1-mbc,my+mbc
        y = ylower + (jlo+jj-.5d0) * dy
        if ((y < ylower) .or. (y>yupper)) then
            do ii=1-mbc,mx+mbc
                x = xlower + (ilo+ii-.5d0) * dx
                iint = ii + max(0, ceiling((xlower-x)/dx)) &
                          - max(0, ceiling((x-xupper)/dx))
                jint = jj + max(0, ceiling((ylower-y)/dy)) &
                          - max(0, ceiling((y-yupper)/dy))
                aux(1,ii,jj) = aux(1,iint,jint)
            enddo
        endif
    enddo


    do ii=1-mbc,mx+mbc
        x =  xlower + (ilo+ii-.5d0) * dx
        if ((x < xlower) .or. (x > xupper)) then
            do jj=1-mbc,my+mbc
                y = ylower + (jlo+jj-.5d0) * dy
                iint = ii + max(0, ceiling((xlower-x)/dx)) &
                          - max(0, ceiling((x-xupper)/dx))
                jint = jj + max(0, ceiling((ylower-y)/dy)) &
                          - max(0, ceiling((y-yupper)/dy))
                aux(1,ii,jj) = aux(1,iint,jint)
            enddo
        endif
    enddo



    ! Output for debugging to fort.23
    if (.false.) then
        print *,'Writing out aux arrays'
        print *,' '
        write(23,230)  mbc,mx,my,dx,dy,xlow,ylow
 230    format('==> mbc, mx, my:  ',3i5,'  dx, dy:',2f10.6, &
                '  xlow,ylow:', 2f10.6)
        do jj=1-mbc,my+mbc
            do ii=1-mbc,mx+mbc
                x = xlow + (ii-0.5d0)*dx
                y = ylow + (jj-0.5d0)*dy
                if ((x>223) .and. (x<232) .and. (y<37)) &
                write(23,231) ii,jj,x,y,(aux(m,ii,jj),m=1,maux)
 231            format(2i4,2f10.3,3e20.10)
            enddo
        enddo
    endif

end subroutine setauxnew
