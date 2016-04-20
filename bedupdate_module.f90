        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Copyright (C) 2015 Virginia Tech, Sediment transport Processes Group    !
        ! Hui Tang, Wei Cheng and Robert Weiss                                    !
        !                                                                         !
        ! tanghui@vt.edu                                                          !
        !                                                                         !
        ! This library is free software; you can redistribute it and/or           !
        ! modify it under the terms of the GNU Lesser General Public              !
        ! License as published by the Free Software Foundation; either            !
        ! version 2.1 of the License, or (at your option) any later version.      !
        !                                                                         !
        ! This library is distributed in the hope that it will be useful,         !
        ! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
        ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
        ! Lesser General Public License for more details.                         !
        !                                                                         !
        ! You should have received a copy of the GNU Lesser General Public        !
        ! License along with this library; if not, write to the Free Software     !
        ! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
        ! USA                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        module bedupdate_module

            use Set_precision, only: Prec
            use transport_module, only: vt_vr,sb_vr,transus
            use update
            use Set_variable
            use sediment_module, only : morstart,morfac,struct,gmax,sourcesink,por, & !dzbdt add to sediment_module
                thick,toler,nd_var,trim,eps
            !use Flux, only: Flux_vector


            IMPLICIT NONE

            contains

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !This part is used for containting the bed updating part                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine bed_update(mbc,mx,my,t,dt,dx,dy,u,v,h,aux,naux)

                use Set_variable

                IMPLICIT NONE

                ! argument
                integer, intent(in) :: mbc,mx,my,naux
                real(kind=Prec), intent(in) :: t,dt,dx,dy
                real(kind=Prec), intent(in) :: u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc)
                real(kind=8), intent(inout) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)
                real(kind=Prec), intent(inout) :: aux(naux,1-mbc:mx+mbc,1-mbc:my+mbc)

                !local
                integer :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: dzbdt,wet,sederold
                integer,dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: indSus,indSub,indSvs,indSvb
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Sout,fre,dzg,edg
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ftem
                Real(kind=Prec) :: Savailable

                wet = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        if (h(i,j)>eps) then
                            wet(i,j) = 1.0
                        endif
                    enddo
                enddo
                zb = aux(1,:,:)

                call transus(mbc,mx,my,dx,dy,t,u,v,h,dt,cu,cub,cv,cvb,ccg,ccbg,Susg,Subg,Svbg,Svsg)


                dzbdt  = 0.0
                !sedero = 0.0

                if (t>=morstart .and. morfac > .9990) then


                    if (struct == 1) then
                        indSus = 0
                        indSub = 0
                        indSvs = 0
                        indSvb = 0
                        Sout   = 0.0
                        do k = 1,gmax
                            do j=2,my+mbc
                                do i=2,mx+mbc
                                    ! fluxes at i,j
                                    if (Subg(i,j,k) > 0.0) then      ! bed load x-direction
                                        indSub(i,j,k) = 1
                                        Sout(i,j,k) = Sout(i,j,k) + Subg(i,j,k)*dx
                                    endif
                                    ! fluxes at i-1,j
                                    if (Subg(i-1,j,k) < 0.0 ) then   ! bed load x-direction
                                        Sout(i,j,k) = Sout(i,j,k) - Subg(i-1,j,k)*dx
                                    endif
                                    if (sourcesink==0) then
                                    ! fluxes at i,j
                                        if (Susg(i,j,k) > 0.0 ) then     ! suspended load x-direction
                                            indSus(i,j,k) = 1
                                            Sout(i,j,k) = Sout(i,j,k) + Susg(i,j,k)*dx
                                        endif
                                        ! fluxes at i-1,j
                                        if (Susg(i-1,j,k) < 0.0 ) then   ! suspended load x-direction
                                            Sout(i,j,k) = Sout(i,j,k) - Susg(i-1,j,k)*dx
                                        endif
                                    endif !sourcesink==0
                                enddo! imax
                            enddo! jmax
                            if (my>0) then
                                do j=2,my+mbc
                                    do i=2,mx+mbc
                                        if (Svbg(i,j,k) > 0.0 ) then     ! bed load y-direction
                                            indSvb(i,j,k) = 1
                                            Sout(i,j,k) = Sout(i,j,k) + Svbg(i,j,k)*dy
                                        endif
                                        ! fluxes at i,j-1
                                        if (Svbg(i,j-1,k) < 0.0 ) then   ! bed load y-direction
                                            Sout(i,j,k) = Sout(i,j,k) - Svbg(i,j-1,k)*dy
                                        endif
                                        if (sourcesink==0) then
                                            if (Svsg(i,j,k) > 0.0 ) then     ! suspended load y-direction
                                                indSvs(i,j,k) = 1
                                                Sout(i,j,k) = Sout(i,j,k) + Svsg(i,j,k)*dy
                                            endif
                                            ! fluxes at i,j-1
                                            if (Svsg(i,j-1,k) < 0.0 ) then   ! suspended load y-direction
                                                Sout(i,j,k) = Sout(i,j,k) - Svsg(i,j-1,k)*dy
                                            endif
                                        endif ! sourcesink = 0
                                    enddo !imax
                                enddo !jmax
                            endif !jmax>0

                            ! reduce the sediment flux/ transport due to hard structure layers
                            do j=2,my+mbc
                                do i=2,mx+mbc
                                    ! reduction factor for cell outgoing sediment transports
                                    Savailable = dzbed(i,j,1)*pbbed(i,j,1,k)/morfac/dt*(1.0-por)*dx*dy
                                    fre(i,j,k)  = min(1.0,Savailable/max(Sout(i,j,k),tiny(0.0)) )
                                    ! fix sediment transports for the presence of a hard layer; remind indSus etc are 1 in cases of cell outgoing transports
                                        ! updated S         all outgoing transports                  cell incoming transports
                                    if (fre(i,j,k) < 1.0)then
                                        Subg(i,j,k)   = fre(i,j,k)*indSub(i,j,k)*Subg(i,j,k) &
                                                        + (1-indSub(i,j,k))*Subg(i,j,k)
                                        Subg(i-1,j,k) = fre(i-1,j,k)*(1-indSub(i-1,j,k))*Subg(i-1,j,k) &
                                                        + indSub(i-1,j,k)*Subg(i-1,j,k)
                                        if (my>0) then
                                            Svbg(i,j,k)   = fre(i,j,k)*indSvb(i,j,k)*Svbg(i,j,k) &
                                                        + (1-indSvb(i,j,k))*Svbg(i,j,k)
                                            Svbg(i,j-1,k) = fre(i,j-1,k)*(1-indSvb(i,j-1,k))*Svbg(i,j-1,k) &
                                                        + indSvb(i,j-1,k)*Svbg(i,j-1,k)
                                        endif
                                        if (sourcesink==0) then
                                            Susg(i,j,k)   = fre(i,j,k)*indSus(i,j,k)*Susg(i,j,k) &
                                                        + (1-indSus(i,j,k))*Susg(i,j,k)
                                            Susg(i-1,j,k) = fre(i-1,j,k)*(1-indSus(i-1,j,k))*Susg(i-1,j,k) &
                                                        + indSus(i-1,j,k)*Susg(i-1,j,k)
                                            if (my>0) then
                                                Svsg(i,j,k)   = fre(i,j,k)*indSvs(i,j,k)*Svsg(i,j,k) &
                                                        + (1-indSvs(i,j,k))*Svsg(i,j,k)
                                                Svsg(i,j-1,k) = fre(i,j,k)*(1-indSvs(i,j-1,k))*Svsg(i,j-1,k) &
                                                        + indSvs(i,j-1,k)*Svsg(i,j-1,k)
                                            endif !jmax > 0
                                        endif ! sourcesink = 0
                                    endif !fac<1.0
                                enddo ! imax
                            enddo !jmax
                        enddo !gmax
                    endif !struct == 1
                    if (trim=='soulsby_vanrijn') then           ! Soulsby van Rijn
                        call sb_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)
                    elseif (trim=='vanthiel_vanrijn') then       ! Van Thiel de Vries & Reniers 2008
                        call vt_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)
                    !elseif (trim=='flux') then
                    !    call Flux_vector(mbc,mx,my,u,v,h,cc,ccb)
                    end if
                    if (my>0) then
                        do j=2,my
                            do i=2,mx
                            ! bed level changes per fraction in this morphological time step in meters sand including pores
                            ! positive in case of erosion

                                if (sourcesink==0) then
                                    do k = 1, gmax
                                        dzg(i,j,k)=wet(i,j)*morfac*dt/(1.0-por)*( &
                                            ! dz from sus transport gradients
                                            (Susg(i,j,k)*dy-Susg(i-1,j,k)*dy +&
                                            Svsg(i,j,k)*dx-Svsg(i,j-1,k)*dx +&
                                            ! dz from bed load transport gradients
                                            Subg(i,j,k)*dy-Subg(i-1,j,k)*dy+&
                                            Svbg(i,j,k)*dx-Svbg(i,j-1,k)*dx))/(dx*dy)! +&
                                    end do
                                        !source term
                                elseif (sourcesink==1) then
                                    dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                        ! dz from sus transport gradients
                                        Susg(i,j,:)-Susg(i-1,j,:) +&
                                        Svsg(i,j,:)-Svsg(i,j-1,:) +&
                                        !source term
                                        (cu(i,j,:)+cv(i,j,:)-ceqsg(i,j,:))*h(i,j)/Tsg(i,j,:))
                                endif
                                if (gmax==1) then ! Simple bed update in case one fraction
                                    zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                                    dzbdt(i,j) = dzbdt(i,j)-sum(dzg(i,j,:))
                                    sedero(i,j) = sedero(i,j)-sum(dzg(i,j,:))
                                    totalthick(i,j) = max(0.0,totalthick(i,j)-sum(dzg(i,j,:)))
                                else
                                    edg(i,j,:) = dzg(i,j,:)*(1.0-por)/dt
                                    if (totalthick(i,j)>thick) then
                                        totalnum(i,j) = nint(totalthick(i,j)/thick)
                                        if (totalthick(i,j)>totalnum(i,j)*thick) then
                                            totalnum(i,j)=totalnum(i,j)+1
                                        endif
                                    else
                                        totalnum(i,j) = 1
                                    endif
                                    if(dzbed(i,j,totalnum(i,j))<toler*thick) then
                                        totalnum(i,j) = totalnum(i,j) - 1
                                    endif
                                    nd_var = totalnum(i,j)
                                    zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                                    sederold(i,j)=sedero(i,j)
                                    sedero(i,j) = sedero(i,j)+sum(dzg(i,j,:))
                                    if (abs(sedero(i,j)) > 1e-10) then
                                        do k=1, gmax
                                            pbbed(i,j,1,k) = abs((pbbed(i,j,1,k)*sederold(i,j)+dzg(i,j,k))/sedero(i,j))
                                        end do
                                    end if
                                    do k=1,gmax
                                        pbbed(i,j,1,k) = pbbed(i,j,1,k)/sum(pbbed(i,j,1,:))
                                    end do
                                    !call update_fractions(i,j,dzbed(i,j,:),pbbed(i,j,:,:),edg(i,j,:),sum(dzg(i,j,:)),dt)
                                endif
                            enddo ! imax+1
                            zb(2-mbc,j) = zb(1-mbc,j)
                        enddo ! jmax+1
                        zb(:,2-mbc) = zb(:,1-mbc)
                        i =2
                    else
                        j=1
                        do i=1-mbc,mx+mbc
                            ! bed level changes per fraction in this morphological time step in meters sand including pores
                            ! positive in case of erosion
                            if (sourcesink==0) then
                                dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                    ! dz from sus transport gradients
                                    Susg(i,j,:)-Susg(i-1,j,:) +&
                                    ! dz from bed load transport gradients
                                    Subg(i,j,:)-Subg(i-1,j,:) +&
                                    !source term
                                    (cu(i,j,:)+cub(i,j,:)-ceqsg(i,j,:)-ceqbg(i,j,:))*h(i,j)/Tsg(i,j,:))
                            elseif (sourcesink==1) then
                                dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                    ! dz from sus transport gradients
                                    Susg(i,j,:)*dx-Susg(i-1,j,:)*dx +&
                                    !source term
                                    (cu(i,j,:)-ceqsg(i,j,:))*h(i,j)/Tsg(i,j,:))
                            endif
                            if (gmax==1) then ! Simple bed update in case one fraction
                                zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                                dzbdt(i,j) = dzbdt(i,j)-sum(dzg(i,j,:))
                                sedero(i,j) = sedero(i,j)-sum(dzg(i,j,:))
                                totalthick(i,j) = max(0.0,totalthick(i,j)-sum(dzg(i,j,:)))
                            else ! multiple fractions...
                                edg(i,j,:) = dzg(i,j,:)*(1.0-por)/dt
                                if (totalthick(i,j)>thick) then
                                    totalnum(i,j) = nint(totalthick(i,j)/thick)
                                    if (totalthick(i,j)>totalnum(i,j)*thick) then
                                        totalnum(i,j)=totalnum(i,j)+1
                                    endif
                                else
                                    totalnum(i,j) = 1
                                endif
                                if(dzbed(i,j,totalnum(i,j))<toler*thick) then
                                    totalnum(i,j) = totalnum(i,j) - 1
                                endif
                                nd_var = totalnum(i,j)
                                !call update_fractions(i,j,dzbed(i,j,:),pbbed(i,j,:,:),edg(i,j,:),sum(dzg(i,j,:)),dt)
                            endif
                        enddo ! imax
                    endif !jmax = 1
                endif !t
                i = 1

                open (unit=1,status='replace',position='Append',file="zb.txt",action="write")
                write (1,*),zb(1,:)
                close(1)
                open (unit=1,status='replace',position='Append',file="zb1.txt",action="write")
                write (1,*),zb(2,:)
                close(1)
                open (unit=1,status='replace',position='Append',file="zb2.txt",action="write")
                write (1,*),zb(3,:)
                close(1)
                open (unit=1,status='replace',position='Append',file="zb3.txt",action="write")
                write (1,*),zb(4,:)
                open (unit=1,status='replace',position='Append',file="pbbed1.txt",action="write")
                write (1,*),pbbed(1,:,1,1)
                close(1)
                open (unit=1,status='replace',position='Append',file="pbbed2.txt",action="write")
                write (1,*),pbbed(2,:,1,1)
                close(1)
                open (unit=1,status='replace',position='Append',file="pbbed3.txt",action="write")
                write (1,*),pbbed(3,:,1,1)
                close(1)
                open (unit=1,status='replace',position='Append',file="pbbed4.txt",action="write")
                write (1,*),pbbed(4,:,1,1)
                if (t<1.0) then
                    open (unit=1,status='replace',position='Append',file="zb4.txt",action="write")
                    write (1,*),zb(1,:)
                    close(1)
                end if


            end subroutine bed_update

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This part is designed for avanlanching scheme
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine avalanch(mbc,mx,my,dx,dy,dt,naux,aux,h,t)

                use sediment_module, only : avalanching,morfac,gmax,hswitch,eps,wetslp,dryslp,&
                        por,nd_var,aval,trim!dzbed,sedcal

                IMPLICIT NONE
                ! argument
                integer, intent(in) :: mbc,mx,my,naux
                real(kind=Prec), intent(in) :: t,dt,dx,dy
                real(kind=Prec), intent(inout) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)
                real(kind=Prec), intent(inout) :: aux(naux,1-mbc:mx+mbc,1-mbc:my+mbc)
                !local
                Integer         :: ie,id,jdz,je,jd,i,j,k,ndz,ii,indx
                Real(kind=Prec) :: sdz,dzb,one=1.00,dzmax,dzleft,dzavt,dzt
                Real(kind=Prec),dimension(gmax) :: edg1,edg2
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: dzbdx,dzbdy,zb,dh,sedero,dzbdt,wet
                !
                zb = aux(1,:,:)

                wet = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        if (h(i,j)>eps) then
                            wet(i,j) = 1.0
                        endif
                    enddo
                enddo

                if (avalanching==1) then
                    do ii=1,morfac
                        aval=.false.
                        dzbdx=0.0
                        dzbdy=0.0
                        do j=1-mbc,my+mbc
                            do i=1-mbc,mx+mbc-1
                                dzbdx(i,j)=(zb(i+1,j)-zb(i,j))/dx
                            enddo
                        enddo
                        do j=1-mbc,my+mbc-1
                            do i=1-mbc,mx+mbc
                                dzbdy(i,j)=(zb(i,j+1)-zb(i,j))/dy
                            enddo
                        enddo

                        indx = mx+1
                        do i=1-mbc,mx+mbc-1
                            do j=1-mbc,my+mbc
                                !decide Maximum bedlevel change due to avalanching
                                if(max(h(i,j),h(i+1,j))>hswitch+eps) then ! Jaap instead of hh
                                    dzmax=wetslp
                                    if (wet(i,j)==1) then ! tricks: seaward of indx (transition from sand to structure) wetslope is set to 0.03; TODO
                                        dzmax = max(wetslp,abs(dzbdx(i,j))*0.990)
                                    endif
                                else
                                    dzmax=dryslp
                                endif
                                if(abs(dzbdx(i,j))>dzmax .and. h(i+nint(max(0.0,sign(one,dzbdx(i,j)))),j)>eps) then
                                    aval=.true.
                                    ! the total mass need to be moved
                                    dzb=sign(one,dzbdx(i,j))*(abs(dzbdx(i,j))-dzmax)*dx
                                    if (dzb >= 0.0) then
                                        ie = i+1                                        ! index erosion point
                                        id = i                                          ! index deposition point
                                        dzb=min(dzb,dzmax*dt/dx)                        ! make sure dzb is not in conflict with maximum erosion rate dzmax
                                        dzb=min(dzb,totalthick(i+1,j))                 ! make sure dzb is not larger than sediment layer thickness
                                    else
                                        ie = i                                          ! index erosion point
                                        id = i+1                                        ! index deposition point
                                        dzb=max(dzb,-dzmax*dt/dx)
                                        dzb=max(dzb,-totalthick(i,j))
                                    endif
                                    if (gmax == 1) then ! Simple bed update in case one fraction
                                        dzleft = abs(dzb)
                                        zb(id,j) = zb(id,j)+dzleft
                                        zb(ie,j) = zb(ie,j)-dzleft
                                        dzbdt(id,j) = dzbdt(id,j)+dzleft
                                        dzbdt(ie,j) = dzbdt(ie,j)-dzleft
                                        sedero(id,j) = sedero(id,j)+dzleft
                                        sedero(ie,j) = sedero(ie,j)-dzleft
                                        totalthick(id,j) = max(0.0,totalthick(id,j)+dzleft)
                                        totalthick(ie,j) = max(0.0,totalthick(ie,j)-dzleft)
                                        h(id,j)  = h(id,j)+dzleft
                                        h(ie,j)  = h(ie,j)-dzleft
                                        !dh(id,j)= dzav(id,j)+dzleft
                                        !dh(ie,j)= dzav(ie,j)-dzleft
                                    else ! multiple fractions...

                                        ! now fix fractions....
                                        !dz = dzbed(ie,j,:)
                                        !pb = pbbed(ie,j,:,:)
                                        ! figure out how many sediment layers (ndz) are eroded in point iii
                                        sdz = 0.0
                                        ndz = 0
                                        do while (sdz<abs(dzb))
                                            ndz = ndz+1
                                            sdz = sdz+dzbed(ie,j,ndz)
                                        enddo
                                        !print *, ndz
                                        ! now update bed and fractions by stepping through each layer seperately
                                        dzleft = abs(dzb)
                                        dzavt  = 0.0
                                        do jdz=1,ndz
                                            dzt = min(dzbed(ie,j,jdz),dzleft)
                                            dzleft = dzleft-dzt
                                            ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion
                                            do k=1,gmax
                                                edg2(k) =  dzt*pbbed(ie,j,jdz,k)*(1.0-por)/dt                ! erosion
                                                edg1(k) = -edg2(k)                                   ! deposition
                                            enddo
                                            dzavt = dzavt + sum(edg2)*dt/(1.0-por)
                                            nd_var = totalnum(ie,j)

                                            call update_fractions(ie,j,dzbed(ie,j,:),pbbed(ie,j,:,:),edg2,dzavt,dt)! update bed in eroding point
                                            nd_var = totalnum(id,j)
                                            call update_fractions(id,j,dzbed(id,j,:),pbbed(id,j,:,:),edg1,-dzavt,dt) ! update bed in deposition point
                                            nd_var = maxval(totalnum)
                                        enddo
                                        ! update water levels and dzav
                                        h(ie,j)  = h(ie,j)-dzavt
                                        dh(ie,j)= dh(ie,j)-dzavt
                                        h(id,j)  = h(id,j)+dzavt
                                        dh(id,j)= dh(id,j)+dzavt
                                    end if ! yes/no multiple fractions
                                end if !dzmax
                            end do !jmax
                        end do !imax
                        !JJ: update y slopes after avalanching in X-direction seems more appropriat

                        do j=1-mbc,my+mbc-1 !
                            do i=1-mbc,mx+mbc
                                if(max(h(i,j),h(i,j+1))>hswitch+eps) then
                                    dzmax=wetslp
                                else
                                    dzmax=dryslp
                                end if

                                if(abs(dzbdy(i,j))>dzmax .and. h(i,j+nint(max(0.0,sign(one,dzbdy(i,j)))))>eps) then

                                    aval=.true.
                                    dzb=sign(one,dzbdy(i,j))*(abs(dzbdy(i,j))-dzmax)*dy
                                    if (dzb >= 0.0) then
                                            je = j+1                                        ! index erosion point
                                            jd = j                                          ! index deposition point
                                            dzb=min(dzb,dzmax*dt/dy)
                                            dzb=min(dzb,totalthick(i,j+1))
                                    else
                                            je = j                                          ! index erosion point
                                            jd = j+1                                        ! index deposition point
                                            dzb=max(dzb,-dzmax*dt/dy)
                                            dzb=max(dzb,-totalthick(i,j))
                                    endif
                                    if (gmax == 1) then ! Simple bed update in case one fraction
                                        dzleft = abs(dzb)
                                        zb(i,jd) = zb(i,jd)+dzleft
                                        zb(i,je) = zb(i,je)-dzleft
                                        dzbdt(i,jd) = dzbdt(i,jd)+dzleft
                                        dzbdt(i,je) = dzbdt(i,je)-dzleft
                                        sedero(i,jd) = sedero(i,jd)+dzleft
                                        sedero(i,je) = sedero(i,je)-dzleft
                                        totalthick(i,jd) = max(0.0,totalthick(i,jd)+dzleft)
                                        totalthick(i,je) = max(0.0,totalthick(i,je)-dzleft)
                                        h(i,jd)  = h(i,jd)+dzleft
                                        h(i,je)  = h(i,je)-dzleft
                                        dh(i,jd)= dh(i,jd)+dzleft
                                        dh(i,je)= dh(i,je)-dzleft
                                    else ! multiple fractions...\

                                        ! figure out how many depth layers (ndz) are affected
                                        sdz = 0
                                        ndz = 0

                                        do while (sdz<abs(dzb))
                                            ndz = ndz+1
                                            sdz = sdz+dzbed(i,je,ndz)
                                        enddo
                                        ! now update bed and fractions by stepping through each layer seperately
                                        dzleft = abs(dzb)
                                        dzavt  = 0.0
                                        do jdz=1,ndz
                                            dzt = min(dzbed(i,je,jdz),dzleft)
                                            dzleft = dzleft-dzt;
                                            !print *, dzbed(i,je,jdz)
                                            ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion
                                            do k=1,gmax
                                                edg2(k) = dzt*pbbed(i,je,jdz,k)*(1.0-por)/dt        ! erosion
                                                edg1(k) = -edg2(k)                            ! deposition
                                            enddo
                                            dzavt = dzavt + sum(edg2)*dt/(1.0-por)
                                            nd_var = totalnum(i,je)
                                            call update_fractions(i,je,dzbed(i,je,:),pbbed(i,je,:,:),edg2,dzavt,dt)           ! upwind point
                                            nd_var = totalnum(i,jd)
                                            call update_fractions(i,jd,dzbed(i,jd,:),pbbed(i,jd,:,:),edg1,-dzavt,dt)    ! downwind point
                                            nd_var = maxval(totalnum)
                                        enddo
                                        ! update water levels and dzav
                                        h(i,je)  = h(i,je)-dzavt
                                        dh(i,je)= dh(i,je)-dzavt
                                        h(i,jd)  = h(i,jd)+dzavt
                                        dh(i,jd)= dh(i,jd)+dzavt
                                    endif !yes/no multiple fractions
                                end if !dzmax
                            end do !imax
                        end do !jmax
                        if (.not.aval) exit
                    end do !morfac
                end if ! avalanching
            end subroutine avalanch


        end module bedupdate_module









