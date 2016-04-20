subroutine initial_sed(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)


    use sediment_module, only: lmax,gmax,thick
    use set_variable

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer :: i,j
    integer, parameter :: READ_UNIT = 101
    !real(kind=8), ::
    if (.not. allocated(Susg)) then
        allocate(Susg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),subg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(Svsg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),svbg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(cu(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cub(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(cv(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cvb(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(ccg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),ccbg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(Tsg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(ceqbg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),ceqsg(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(ws(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(z0(1-mbc:mx+mbc,1-mbc:my+mbc))
        allocate(zb(1-mbc:mx+mbc,1-mbc:my+mbc),ero(1-mbc:mx+mbc,1-mbc:my+mbc))
        allocate(ub_cr(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),us_cr1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),&
        us_cr2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(Sus(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Svs(1-mbc:mx+mbc,1-mbc:my+mbc,gmax), &
        Sub(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Svb(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),&
            percd(1-mbc:mx+mbc,1-mbc:my+mbc,gmax))
        allocate(dzbed(1-mbc:mx+mbc,1-mbc:my+mbc,lmax))
        allocate(totalnum(1-mbc:mx+mbc,1-mbc:my+mbc)) !totalthick(1-mbc:mx+mbc,1-mbc:my+mbc),
        allocate(z0bed(1-mbc:mx+mbc,1-mbc:my+mbc),sedero(1-mbc:mx+mbc,1-mbc:my+mbc))
        allocate(profile(24606,2))
    end if
    ub_cr = 0.0
    us_cr1 = 0.0
    us_cr2 = 0.0
    ws = 0.0
    cu = 0.0
    cub = 0.0
    cv = 0.0
    cvb = 0.0
    ccg = 0.0
    ccbg = 0.0
    Susg = 0.0
    Subg = 0.0
    Svbg = 0.0
    Svsg = 0.0
    ceqbg = 0.0
    ceqsg = 0.0
    Tsg = 0.0
    z0 = 0.0
    Sus = 0.0
    Sub = 0.0
    Svs = 0.0
    Svb = 0.0
    percd(:,:,1) = 1.0
    percd(:,:,2) = 0.0
    totalnum = 20
    do i = 1, totalnum
        dzbed(:,:,i) = thick
    end do
    totalthick = 1.0
    zb = aux(1,:,:)
    ero = aux(1,:,:)
    z0bed = zb-totalthick
    sedero = 0.01
    open(unit=READ_UNIT,file='Boudary_depth_data.txt',status="old",action="read")
    do i = 1, 24606
        read(READ_UNIT,*) profile(i,1),profile(i,2)
    enddo
    close(READ_UNIT)
end subroutine initial_sed