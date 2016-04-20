    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This Part is used to calculate bed roughness                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine bedrough(mbc,mx,my,u,v,h,z0)

                use sediment_module, only: hcr,D,g,m0,rhos,rho,gmax
                use Set_Precision, only: Prec

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc),&
                        v(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local

                integer :: i,j
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,hloc,Dmm,zon,a2,ustarc,ustarcrit, &
                                taub,taucrit,Tstar,delb,zos
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,ub_cr,us_cr1,us_cr2
                Real(kind=Prec) :: a1,gammaWs

                !output
                Real(kind=Prec),intent(out),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: z0

                call crtical_velocity1(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

                gammaWs = 0.056
                vmag2 = u**2.0+v**2.0
                hloc = max(h,hcr)
                Dmm = meansize(D,pbbed(:,:,1,:),mx,my,mbc,gmax)
                zon = Dmm/30.0
                a1 = 0.68
                a2 = 0.0204*(log(Dmm*100))**2.0+0.0220*log(Dmm*100)+0.0709
                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        ustarc(i,j) = sqrt(g*m0**2.0*vmag2(i,j)/((rhos-rho)*Dmm(i,j)*hloc(i,j)**(1.0/3.0)))
                        ustarcrit(i,j) = sqrt(g*m0**2.0*ub_cr(i,j,int(gmax/2.0+1.0))**2.0/((rhos-rho)*Dmm(i,j)&
                            *hloc(i,j)**(1.0/3.0)))
                    enddo
                enddo
                taub = rho*ustarc**2.0
                taucrit=rho*ustarcrit**2.0
                do i = 1-mbc,mx+mbc
                    do j =1-mbc, my+mbc
                        Tstar(i,j)=taub(i,j)/taucrit(i,j)
                        delb(i,j)=Dmm(i,j)*a1*Tstar(i,j)/(1.0+a2(i,j)*Tstar(i,j))
                    enddo
                enddo
                zos = gammaWs*delb
                z0 = zon+zos
            end subroutine bedrough
            !*******************************************************************
            !
            !
            !This function is used to caculate mean grain size
            !
            !
            !********************************************************************
            function meansize(D,fr,mx,my,mbc,gmax)



                implicit none

                integer :: mx,my,mbc,gmax
                Real(kind=Prec), Dimension(gmax) :: D
                Real(kind=Prec), Dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: fr
                Real(kind=Prec), Dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: meansize
                integer :: i,j,k
                meansize = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        do k = 1, gmax
                            meansize(i,j) = meansize(i,j)+D(k)*fr(i,j,k)
                        enddo
                    enddo
                enddo
            end function meansize
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate critical velocity for each grain size classes!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine crtical_velocity1(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)!,u,v)

                use sediment_module, only: rhos,rho,gmax,g,D,hcr,k0,m0,pbbed
                use Set_Precision, only: Prec

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc) &
                        ,v(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local
                integer :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: hloc,z0,Dmm
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,ub_crt, us_cr1t,us_cr2t
                !Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) ::u_sh
                Real(kind=Prec) :: delta
                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ub_cr, us_cr1, us_cr2

                frc = pbbed(:,:,1,:)
                call settling_velocity(mbc,mx,my,ws)
                !call bedrough(mbc,mx,my,u,v,h,z0)
                Dmm = meansize(D,pbbed(:,:,1,:),mx,my,mbc,gmax)

                z0 = Dmm/30.0
                delta  = (rhos-rho)/rho
                hloc = max(h,hcr)
                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        do k = 1, gmax
                            if (D(k) <= 0.0005) then
                                ub_crt(i,j,k) = 0.19*D(k)**0.1*dlog10(4*hloc(i,j)/D(k))
                            elseif (D(k) > 0.0005 .AND. D(k) <= 0.002) then
                                ub_crt(i,j,k) = 8.5*D(k)**0.6*dlog10(4*hloc(i,j)/D(k))
                            else
                                print *, "The Grain size is out of range"
                            end if
                            !u_sh(i,j,k) = ub_cr*k0*(hloc(i,j)-z0(i,j))/(z0(i,j)-hloc(i,j)+hloc(i,j)*log(hloc(i,j)/z0(i,j)
                            us_cr1t(i,j,k) = ws(i,j,k)*((z0(i,j)-hloc(i,j)+hloc(i,j)*log(hloc(i,j)/z0(i,j)))/ &
                                ((hloc(i,j)-z0(i,j))*2.5))!*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(2.5*k0*sqrt(g)*m0)
                            us_cr2t(i,j,k) = ws(i,j,k)*((z0(i,j)-hloc(i,j)+hloc(i,j)*log(hloc(i,j)/z0(i,j)))/ &
                                ((hloc(i,j)-z0(i,j))*1.2))!*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(1.2*k0*sqrt(g)*m0)
                            ub_cr(i,j,k) = min(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                            us_cr2(i,j,k) = max(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                            us_cr1(i,j,k) =ub_crt(i,j,k)+us_cr1t(i,j,k)+us_cr2t(i,j,k)-us_cr2(i,j,k)-ub_cr(i,j,k)

                        enddo
                    enddo
                enddo
            end subroutine crtical_velocity1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate critical velocity for each grain size classes!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine crtical_velocity2(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

                use sediment_module, only: rhos,rho,gmax,g,D,hcr,Trep,k0,m0,pbbed
                use Set_Precision, only: Prec

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc) &
                    ,v(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local
                integer :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: hloc,z0,Dmm
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,ub_crt, us_cr1t,us_cr2t
                Real(kind=Prec) :: delta
                !output
                Real(kind=Prec),intent(inout) :: ub_cr(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
                    , us_cr1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
                    , us_cr2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

                frc = pbbed(:,:,1,:)
                Dmm = meansize(D,pbbed(:,:,1,:),mx,my,mbc,gmax)
                call settling_velocity(mbc,mx,my,ws)
                z0 = Dmm/30.0
                delta  = (rhos-rho)/rho

                hloc = max(h,hcr)

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        do k = 1, gmax

                            if (D(k) <= 0.0005) then
                                ub_crt(i,j,k) = 0.24*(delta*g)**(2.0/3.0)*(D(k)*Trep)**(1.0/3.0)
                            else
                                ub_crt(i,j,k) = 0.95*(delta*g)**(0.57)*D(k)**0.43*Trep**0.14
                            end if
                            us_cr1t(i,j,k) = ws(i,j,k)*((z0(i,j)-hloc(i,j)+hloc(i,j)*log(hloc(i,j)/z0(i,j)))/ &
                                ((hloc(i,j)-z0(i,j))*2.5))
                            us_cr2t(i,j,k) = ws(i,j,k)*((z0(i,j)-hloc(i,j)+hloc(i,j)*log(hloc(i,j)/z0(i,j)))/ &
                                ((hloc(i,j)-z0(i,j))*1.2))
                            ub_cr(i,j,k) = min(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                            us_cr2(i,j,k) = max(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                            us_cr1(i,j,k) =ub_crt(i,j,k)+us_cr1t(i,j,k)+us_cr2t(i,j,k)-us_cr2(i,j,k)-ub_cr(i,j,k)
                            !uscr1(i,j,k) = ws(i,j,k)*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(2.5*k0*sqrt(g)*m0)
                            !uscr2(i,j,k) = ws(i,j,k)*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(1.2*k0*sqrt(g)*m0)
                        end do
                    end do
                end do
            end subroutine crtical_velocity2

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate settling velocity for each grain size classes!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine settling_velocity(mbc,mx,my,ws)

                use sediment_module, only: rhos,rho,gmax,g,D
                use Set_Precision, only: Prec

                !use params

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my

                !local
                integer :: i,j,k
                Real(kind=Prec) :: delta,Te,vis,Sster,c1,c2,wster
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: C
                Real(kind=Prec),dimension(gmax) :: w,R,alpha
                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ws

                delta  = (rhos-rho)/rho

                C = ccbg+ccg !how to pass to here

                do k = 1, gmax
                    Te    = 20.0
                    vis  = 4.0/(20.0+Te)*1e-5 ! Van rijn, 1993
                    Sster = D(k)/(4*vis)*sqrt((rhos/rho-1)*g*D(k))
                    c1    = 1.06*tanh(0.064*Sster*exp(-7.5/Sster**2.0))
                    c2    = 0.220*tanh(2.34*Sster**(-1.180)*exp(-0.0064*Sster**2.0))
                    wster = c1+c2*Sster
                    w(k) = wster*sqrt((rhos/rho-1.0)*g*D(k))
                    R(k) = w(k)*D(k)/vis
                    alpha(k) = 2.35*(2.0+0.175*R(k)**(3.0/4.0))/(1.0+0.175*R(k)**(3.0/4.0))
                    do i =1-mbc, mx+mbc
                        do j = 1-mbc, my+mbc
                            ws(i,j,k) = (1-C(i,j,k))**alpha(k)*w(k)
                        enddo
                    enddo
                enddo

            end subroutine settling_velocity

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Soulsby-VanRijn Method                                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine sb_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)

                use sediment_module, only: eps,rhos,rho,gmax,g,D,hcr,tsfac,Tsmin,cmax,sws
                use Set_Precision, only: Prec
                use Set_variable, only:pbbed

                implicit none

                ! Arguments
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local
                integer :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   wet,vmg,urms,urms2,hloc,Cd,Asb,term1,term2
                Real(kind=Prec),dimension(gmax) :: dster
                Real(kind=Prec) :: delta,Ass,perc
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Ts,ceq,ceqs,ceqb

                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Tsg,ceqbg,ceqsg


                wet = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        if (h(i,j)>eps) then
                            wet(i,j) = 1.0
                        endif
                    enddo
                enddo

                delta = rhos - rho

                do k=1,gmax
                    dster(k)=(delta*g/1e-12)**(1.0/3.0)*D(k)
                enddo

                urms = 0.0

                vmg  = dsqrt(u**2+v**2)

                urms2  = urms**2.0

                hloc = max(h,hcr)

                ! calculate threshold velocity Ucr for bedload

                call settling_velocity(mbc,mx,my,ws)

                call crtical_velocity1(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

                ! bed roughness

                call bedrough(mbc,mx,my,u,v,h,z0)



                do k = 1, gmax

                    Ts(:,:,k) = tsfac*hloc/ws(:,:,k)
                    Tsg(:,:,k) = max(Ts(:,:,k),Tsmin)

                    ! drag coefficient
                    do i = 1-mbc, mx+mbc
                        do j = 1-mbc, my+mbc
                            Cd(i,j)=(0.40/(log(max(hloc(i,j),10.0*z0(i,j))/z0(i,j))-1.0))**2.0
                        enddo
                    enddo

                    ! transport parameters
                    Asb=0.005*hloc*(D(k)/hloc/(delta*g*D(k)))**1.20         ! bed load coefficent
                    Ass=0.0120*D(k)*dster(k)**(-0.60)/(delta*g*D(k))**1.20  ! suspended load coeffient

                    term1 = (vmg**2.0+0.018/Cd*sws*urms2)
                    term1 = sqrt(term1)
                    term2 = 0.0*term1

                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            if(term1(i,j)>Ub_cr(i,j,k) .and. hloc(i,j)>eps) then
                                term2(i,j)=(term1(i,j)-Ub_cr(i,j,k))**2.40
                            endif
                            ceq(i,j,k) = (Asb(i,j)+Ass)*term2(i,j)/hloc(i,j)
                        enddo
                    enddo
                    do j = 1-mbc,my+mbc
                        do i = 1-mbc,mx+mbc
                            if(term1(i,j)<Us_cr1(i,j,k)) then
                                ceqb(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = 0.0
                            elseif(term1(i,j)>Us_cr2(i,j,k)) then
                                ceqb(i,j,k) = 0.0
                                ceqs(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                            else
                                perc = term1(i,j)/Us_cr2(i,j,k)
                                ceqb(i,j,k) = (1-perc)*min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = perc*min(ceq(i,j,k),cmax/gmax/2.0)
                            end if
                            ceqbg(i,j,k) = ceqb(i,j,k)*wet(i,j)*pbbed(i,j,1,k)
                            ceqsg(i,j,k) = ceqs(i,j,k)*wet(i,j)*pbbed(i,j,1,k)
                        enddo
                    enddo
                enddo

            end subroutine sb_vr

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Van Thiel-Van Rijn Method                                                          !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine vt_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)

                use sediment_module, only: eps,rhos,rho,gmax,g,D,hcr,tsfac,Tsmin,cmax,sws
                use Set_Precision, only: Prec
                use Set_variable, only:pbbed

                implicit none

                ! Arguments
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local

                integer :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   wet,vmg,urms,urms2,hloc,Cd,Asb,term1,term2,term3
                Real(kind=Prec),dimension(gmax) :: dster
                Real(kind=Prec) :: delta,Ass,perc
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Ts,ceq,ceqs,ceqb
                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Tsg,ceqbg,ceqsg



                wet = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        if (h(i,j)>eps) then
                            wet(i,j) = 1.0
                        endif
                    enddo
                enddo

                delta = rhos - rho

                do k=1,gmax
                    dster(k)=(delta*g/1e-12)**(1.0/3.0)*D(k)
                enddo

                urms = 0.0

                vmg  = dsqrt(u**2+v**2)

                urms2  = urms**2.0

                hloc = max(h,hcr)

                ! calculate threshold velocity Ucr for bedload

                call settling_velocity(mbc,mx,my,ws)

                call crtical_velocity2(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

                do k = 1,gmax

                    ! transport parameters

                    Asb=0.015*hloc*(D(k)/hloc)**1.20/(delta*g*D(k))**0.75        !bed load coefficent
                    Ass=0.012*D(k)*dster(k)**(-0.60)/(delta*g*D(k))**1.20        !suspended load coeffient

                    term1=vmg**2+0.640*sws*urms2

                    term1=sqrt(term1)

                    !initialize ceqs
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            if(term1(i,j)>Ub_cr(i,j,k) .and. h(i,j)>eps) then
                                term2(i,j)=(term1(i,j)-Ub_cr(i,j,k))**1.50
                                term3(i,j)=(term1(i,j)-Ub_cr(i,j,k))**2.40
                            endif
                        enddo
                    enddo
                    ceq(:,:,k) = (Asb*term2+Ass*term3)/hloc
                    do j = 1-mbc,my+mbc
                        do i = 1-mbc,mx+mbc
                            if(term1(i,j)<Us_cr1(i,j,k)) then
                                ceqb(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = 0.0
                            elseif(term1(i,j)>Us_cr2(i,j,k)) then
                                ceqb(i,j,k) = 0.0
                                ceqs(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                            else
   
                                ceqb(i,j,k) = Asb(i,j)*term2(i,j)/hloc(i,j)!(1-perc)*min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = Ass*term3(i,j)/hloc(i,j)!perc*min(ceq(i,j,k),cmax/gmax/2.0)
                            endif
                            ceqbg(i,j,k) = ceqb(i,j,k)*wet(i,j)*pbbed(i,j,1,k)
                            ceqsg(i,j,k) = ceqs(i,j,k)*wet(i,j)*pbbed(i,j,1,k)
                        enddo
                    enddo
                enddo
            end subroutine vt_vr
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate Source term for finite volume method                !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine transus(mbc,mx,my,dx,dy,time,u,v,h,dt,cut,cubt,cvt,cvbt,ccgt,ccbgt,Susgt,Subgt,Svbgt,Svsgt)

                use flux, only: Flux_vector
                use sediment_module, only: trim,gmax,morfac,por,D,thetanum,cmax,lmax,eps,facDC,nuh,nuhfac,rho !&
                    !,facsl
                use Set_Precision, only: Prec
                use Set_variable,only: zb

                implicit none

                ! Arguments
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: dx,dy,dt,time
                real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local

                integer :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,dcsdy,dcbdy,dcsdx,dcbdx,hold,DR,Dc,wet,dzbdx,&
                    pbbedu,pbbedv,dzbdy!hold?
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,fac,ero1,ero2,depo_ex1,depo_ex2, &
                                cc,ccb
                Real(kind=Prec) :: exp_ero,facsl

                !out
                Real(kind=Prec),intent(inout) :: cut(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cubt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
                Real(kind=Prec),intent(inout) :: cvt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cvbt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
                Real(kind=Prec),intent(inout) :: ccgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),ccbgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
                Real(kind=Prec),intent(inout) :: Susgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Subgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
                Real(kind=Prec),intent(inout) :: Svbgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Svsgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)




                vmag2     = u**2+v**2
                cut = 0.0
                cubt = 0.0
                cvt = 0.0
                cvbt = 0.0
                dcsdx = 0.0
                dcbdx = 0.0
                dcsdy = 0.0
                dcbdy = 0.0
                DR = 0.0 !turn off roller model
                wet = 0.0
                exp_ero = 0.0
                facsl = 1.6
                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        if (h(i,j)>eps) then
                            wet(i,j) = 1.0
                        endif
                    enddo
                enddo
                do i = 1-mbc, mx+mbc-1
                    do j = 1-mbc, my+mbc
                        dzbdx = (zb(i+1,j)-zb(i,j))/dx
                    enddo
                enddo
                dzbdx(mx+mbc,:) = dzbdx(mx+mbc-1,:)
                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc-1
                        dzbdy = (zb(i,j+1)-zb(i,j))/dy
                    enddo
                enddo
                dzbdy(:,my+mbc) = dzbdy(:,my+mbc-1)
                ! calculate equibrium sediment concentration
                if (trim=='soulsby_vanrijn') then           ! Soulsby van Rijn
                    call sb_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)
                elseif (trim=='vanthiel_vanrijn') then       ! Van Thiel de Vries & Reniers 2008
                    call vt_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)
                end if
                ! compute reduction factor for sediment sources due to presence of hard layers

                frc = pbbed(:,:,1,:)

                do k = 1,gmax
                    do j= 1-mbc,my+mbc
                        do i= 1-mbc,mx+mbc
                            exp_ero = morfac*dt/(1.0-por)*h(i,j)*(ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k) &
                                    + ceqbg(i,j,k)*frc(i,j,k)/dt) !ceqsg, ceqbg from vt_vr or sb_vr
                            fac(i,j,k) =min(1.0,dzbed(i,j,1)*frc(i,j,k)/max(tiny(0.0),exp_ero))! totalthick or dzbed TODO
                        enddo
                    enddo
                enddo

                ! compute diffusion coefficient
                Dc = facDc*(nuh+nuhfac*h*(DR/rho)**(1.0/3.0))
                cc = ccgt
                ccb = ccbgt
                Sus = 0.0
                Sub = 0.0
                do k = 1,gmax

                    if (D(k)>0.002) then
                        print *, "WARNING: Grain size is larger than 2 mm"
                    endif
                    ! x-direction
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc-1
                            if(u(i,j)>0.0) then
                                cut(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(min(i+1,mx+mbc-1),j,k)
                                cubt(i,j,k)=thetanum*pbbed(i,j,1,k)*ceqbg(i,j,k)+(1.0-thetanum)&
                                            *pbbed(min(i+1,mx+mbc-1),j,1,k)*ceqbg(min(i+1,mx+mbc-1),j,k)
                                !cubt(i,j,k)=thetanum*ccb(i,j,k)+(1.0-thetanum)*ccb(min(i+1,mx+mbc),j,k)
                            elseif(u(i,j)<0.0) then
                                cut(i,j,k)=thetanum*cc(i+1,j,k)+(1.0-thetanum)*cc(max(i,2-mbc),j,k)
                                cubt(i,j,k)=thetanum*pbbed(i+1,j,1,k)*ceqbg(i+1,j,k)+(1.0-thetanum)&
                                        *pbbed(max(i,2-mbc),j,1,k)*ceqbg(max(i,2-mbc),j,k)
                                !cubt(i,j,k)=thetanum*ccb(i+1,j,k)+(1.0-thetanum)*ccb(max(i,2-mbc),j,k)
                            else
                                cut(i,j,k)=0.50*(cc(i,j,k)+cc(i+1,j,k))
                                cubt(i,j,k)=0.50*(pbbed(i,j,1,k)*ceqbg(i,j,k)+pbbed(i+1,j,1,k)*ceqbg(i+1,j,k))
                                !cubt(i,j,k)=0.50*(ccb(i,j,k)+ccb(i+1,j,k))
                            endif
                            dcsdx(i,j)=(sum(cc(i+1,j,:))-sum(cc(i,j,:)))/dx
                            dcbdx(i,j)=(sum(ccb(i+1,j,:))-sum(ccb(i,j,:)))/dx
                        enddo
                    enddo
                    cut(mx+mbc,:,:) = cc(mx+mbc,:,:)
                    cubt(mx+mbc,:,:) = ccb(mx+mbc,:,:)

                    !Sus = 0.0
                    !Sub = 0.0
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                        ! suspended load
                            Sus(i,j,k)=(cut(i,j,k)*u(i,j)*h(i,j)-Dc(i,j)*h(i,j)*dcsdx(i,j) &
                                -facsl*cut(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdx(i,j))*wet(i,j)   !No bed slope term in suspended transport?
                        ! bed load
                            Sub(i,j,k)=(cubt(i,j,k)*u(i,j)*h(i,j)-facsl*cubt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdx(i,j))*wet(i,j)
                        enddo
                    enddo
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc-1
                            if(Sub(i,j,k)>0.0) then
                                pbbedu(i,j) = pbbed(i,j,1,k)
                            elseif(Sub(i,j,k)<0.0) then
                                pbbedu(i,j)= pbbed(i+1,j,1,k)
                            else
                                pbbedu(i,j)=0.50*(pbbed(i,j,1,k)+pbbed(i+1,j,1,k))
                            endif
                        enddo
                    enddo
                    pbbedu(mx+mbc,:) = pbbedu(mx+mbc-1,:)
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            Sub(i,j,k) = pbbedu(i,j)*Sub(i,j,k)
                            Sus(i,j,k) = pbbedu(i,j)*Sus(i,j,k)
                        enddo
                    enddo
                    !y-direction
                    if (my>0) then
                        do j=1-mbc,my+mbc-1
                            do i=1-mbc,mx+mbc
                                if(v(i,j)>0) then
                                    cvt(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(i,min(j+1,my+mbc),k)
                                    cvbt(i,j,k)=thetanum*pbbed(i,j,1,k)*ceqbg(i,j,k)+(1.0-thetanum)&
                                        *pbbed(i,min(j+1,my+mbc),1,k)*ceqbg(i,min(j+1,my+mbc),k)
                                    !cvbt(i,j,k)=thetanum*ccb(i,j,k)+(1.0-thetanum)*ccb(i,min(j+1,my+mbc),k)
                                elseif(v(i,j)<0) then
                                    cvt(i,j,k)=thetanum*cc(i,j+1,k)+(1.0-thetanum)*cc(i,max(j,2-mbc),k)
                                    cvbt(i,j,k)=thetanum*pbbed(i,j+1,1,k)*ceqbg(i,j+1,k)+(1.0-thetanum)&
                                        *pbbed(i,max(j,2),1,k)*ceqbg(i,max(j,2),k)
                                    !cvbt(i,j,k)=thetanum*ccb(i,j+1,k)+(1.0-thetanum)*ccb(i,max(j,2-mbc),k)
                                else
                                    cvt(i,j,k)=0.50*(cc(i,j,k)+cc(i,j+1,k)) !Jaap: cc instead of cv
                                    cvbt(i,j,k)=0.50*(pbbed(i,j,1,k)*ceqbg(i,j,k)+pbbed(i,j+1,1,k)*ceqbg(i,j+1,k))
                                    !cvbt(i,j,k)=0.50*(ccb(i,j,k)+ccb(i,j+1,k))
                                end if
                                dcsdy(i,j)=(sum(cc(i,j+1,:))-sum(cc(i,j,:)))/dy !Jaap
                                dcbdy(i,j)=(sum(ccb(i,j+1,:))-sum(ccb(i,j,:)))/dy
                            end do
                        end do
                        cvt(:,my+mbc,:) = cc(:,my+mbc,:)
                        cvbt(:,my+mbc,:) = ccb(:,my+mbc,:)
                    else
                        cvt = cc
                        cvbt = ceqbg
                    endif ! my>0

                    !Svs = 0.0
                    !Svb = 0.0
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            ! suspended load
                            Svs(i,j,k)=(cvt(i,j,k)*v(i,j)*h(i,j)-Dc(i,j)*h(i,j)*dcsdy(i,j) &
                                    -facsl*cvt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdy(i,j))*wet(i,j)   !No bed slope term in suspended transport?
                            !bed load
                            Svb(i,j,k)=(cvbt(i,j,k)*v(i,j)*h(i,j)-facsl*cvbt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdy(i,j))*wet(i,j)
                        enddo
                    enddo
                    do j=1-mbc,my+mbc-1
                        do i=1-mbc,mx+mbc
                            if(Svb(i,j,k)>0.0) then
                                pbbedv(i,j) = pbbed(i,j,1,k)
                            elseif(Svb(i,j,k)<0.0) then
                                pbbedv(i,j)= pbbed(i,j+1,1,k)
                            else
                                pbbedv(i,j)=0.50*(pbbed(i,j,1,k)+pbbed(i,j+1,1,k))
                            endif
                        enddo
                    enddo
                    pbbedv(:,my+mbc) = pbbedv(:,my+mbc-1)
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            Svs(i,j,k) = pbbedv(i,j)*Svs(i,j,k)
                            Svb(i,j,k) = pbbedv(i,j)*Svb(i,j,k)

                        enddo
                    enddo
                    !call Flux_vector(mbc,mx,my,u,v,h,cc,ccb)
                    if (my>0) then
                        do j=2-mbc,my+mbc-1
                            do i=2-mbc,mx+mbc-1
                                !suspended sediment transport
                                ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                                cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (h(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dx-Sus(i-1,j,k)*dx+&
                                        Svs(i,j,k)*dy-Svs(i,j-1,k)*dy)/(dx*dy)-ero1(i,j,k)))
                                cc(i,j,k)=max(cc(i,j,k),0.00)
                                cc(i,j,k)=min(cc(i,j,k),cmax/2.0*h(i,j))
                                depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)
                                !bed sediment tranpsort
                                ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                                ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (h(i,j)*cc(i,j,k)/dt -((Sub(i,j,k)*dx-Sub(i-1,j,k)*dx+&
                                        Svb(i,j,k)*dy-Svb(i,j-1,k)*dy)/(dx*dy)-ero2(i,j,k)))
                                ccb(i,j,k)=max(ccb(i,j,k),0.00)
                                ccb(i,j,k)=min(ccb(i,j,k),cmax/2.0*h(i,j))
                                depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                            enddo
                        enddo
                    else
                        j=1
                        do i=2-mbc,mx+mbc-1
                            !suspended sediment transport
                            ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                            cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (h(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dx-Sus(i-1,j,k)*dx)/(dx*dy)-&
                                        ero1(i,j,k)))
                            cc(i,j,k)=max(cc(i,j,k),0.00)
                            cc(i,j,k)=min(cc(i,j,k),cmax/2.0*h(i,j))
                            depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)
                            !bed sediment tranpsort
                            ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                            ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (h(i,j)*ccb(i,j,k)/dt -((Sub(i,j,k)*dx-Sub(i-1,j,k)*dx)/(dx*dy)-&
                                        ero2(i,j,k)))
                            ccb(i,j,k)=max(ccb(i,j,k),0.00)
                            ccb(i,j,k)=min(ccb(i,j,k),cmax/2.0*h(i,j))
                            depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                        enddo
                    endif
                    do j= 1-mbc,my+mbc
                        do i= 1-mbc,mx+mbc
                            if (h(i,j) .lt. eps) then
                                cc(i,j,k) = 0.0
                                ccb(i,j,k) = 0.0
                            else
                                cc(i,j,k) = cc(i,j,k)/h(i,j)
                                ccb(i,j,k) = ccb(i,j,k)/h(i,j)
                            end if
                        end do
                    end do
                    ccbgt = ccb
                    ccgt = cc
                end do
                Svsgt = Svs
                Susgt = Sus
                Svbgt = Svb
                Subgt = Sub
            end subroutine transus
        end module transport_module
