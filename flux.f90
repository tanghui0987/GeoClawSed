!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Part is used to calculate sediment Flux for finite volume method              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module Flux

        !use params
        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used to caculate flux term                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine Flux_vector(mbc,mx,my,u,v,h,cc,ccb)

            use sediment_module,only: gmax,vareps,k1,method,toler,vareps,g
            use Set_Precision,only: Prec

            implicit none

            ! Arguments
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)
            real(kind=Prec), intent(in) ::cc(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),ccb(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            !local
            integer :: i,j,k,jg
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,3) :: V_new,Vx_l,Vx_r,Vy_l,Vy_r
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax,2) ::C_new_x,C_new_y,Cx_l,Cx_r,Cy_l,Cy_r
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3) :: psix1,psix2,psiy1,psiy2
            !flux limitor

            V_new(1-mbc:mx+mbc,1-mbc:my+mbc,1) = u
            V_new(1-mbc:mx+mbc,1-mbc:my+mbc,2) = v
            V_new(1-mbc:mx+mbc,1-mbc:my+mbc,3) = h
            C_new_x(1-mbc:mx+mbc,1-mbc:my+mbc,:,1)= cc
            C_new_x(1-mbc:mx+mbc,1-mbc:my+mbc,:,2)= ccb
            C_new_y(1-mbc:mx+mbc,1-mbc:my+mbc,:,1)= cc
            C_new_y(1-mbc:mx+mbc,1-mbc:my+mbc,:,2)= ccb
            V_new(1-mbc,:,:) = V_new(2-mbc,:,:)
            V_new(mx+mbc,:,:) = V_new(mx+mbc-1,:,:)
            V_new(:,1-mbc,:) = V_new(:,2-mbc,:)
            V_new(:,my+mbc,:) = V_new(:,my+mbc-1,:)
            C_new_x(1-mbc,:,:,:) = C_new_x(2-mbc,:,:,:)
            C_new_x(mx+mbc,:,:,:) = C_new_x(mx+mbc-1,:,:,:)
            C_new_x(:,1-mbc,:,:) = C_new_x(:,2-mbc,:,:)
            C_new_x(:,my+mbc,:,:) =C_new_x(:,my+mbc-1,:,:)
            C_new_y(1-mbc,:,:,:) = C_new_y(2-mbc,:,:,:)
            C_new_y(mx+mbc,:,:,:) = C_new_y(mx+mbc-1,:,:,:)
            C_new_y(:,1-mbc,:,:) = C_new_y(:,2-mbc,:,:)
            C_new_y(:,my+mbc,:,:) =C_new_y(:,my+mbc-1,:,:)

            !flux limitor
            call flux_limitor(mbc,mx,my,V_new,C_new_x,C_new_y,psix1,psix2,psiy1,psiy2)

            do i=1-mbc, mx+mbc
                do j=1-mbc, my+mbc
                    do jg = 1, 3
                        Vx_l(i,j,jg) = V_new(i-1,j,jg) + 0.25*vareps* &
                                        ((1-k1)*psix1(i-1,j,jg)*(V_new(i,j,jg)-V_new(i-1,j,jg))+ &
                                        (1+k1)*psix2(i,j,jg)*(V_new(i+1,j,jg)-V_new(i,j,jg))      )
                        Vx_r(i,j,jg) = V_new(i,j,jg)  + 0.25*vareps* &
                                        ((1-k1)*psix2(i+1,j,jg)*(V_new(i+2,j,jg)-V_new(i+1,j,jg)) + &
                                        (1+k1)*psix1(i,j,jg)*(V_new(i+1,j,jg)-V_new(i,j,jg))      )
                        Vy_l(i,j,jg) = V_new(i,j-1,jg) + 0.25*vareps* &
                                        ((1-k1)*psiy1(i,j-1,jg)*(V_new(i,j,jg)-V_new(i,j-1,jg))+ &
                                        (1+k1)*psiy2(i,j,jg)*(V_new(i+1,j,jg)-V_new(i,j,jg))      )
                        Vy_r(i,j,jg) = V_new(i,j,jg)  + 0.25*vareps* &
                                        ((1-k1)*psiy2(i,j+1,jg)*(V_new(i,j+2,jg)-V_new(i,j+1,jg)) + &
                                        (1+k1)*psiy1(i,j,jg)*(V_new(i,j+1,jg)-V_new(i,j,jg))      )
                    enddo
                    do jg=1,gmax
                            Cx_l(i,j,jg,:) = C_new_x(i-1,j,jg,:) + 0.25*vareps* &
                                        ((1-k1)*psix1(i-1,j,jg+3)*(C_new_x(i,j,jg,:)-C_new_x(i-1,j,jg,:))+ &
                                        (1+k1)*psix2(i,j,jg+3)*(C_new_x(i+1,j,jg,:)-C_new_x(i,j,jg,:))      )
                            Cx_r(i,j,jg,:) = C_new_x(i,j,jg,:)  + 0.25*vareps* &
                                        ((1-k1)*psix2(i+1,j,jg+3)*(C_new_x(i+2,j,jg,:)-C_new_x(i+1,j,jg,:)) + &
                                        (1+k1)*psix1(i,j,jg+3)*(C_new_x(i+1,j,jg,:)-C_new_x(i,j,jg,:))      )
                            Cy_l(i,j,jg,:) = C_new_y(i,j-1,jg,:) + 0.25*vareps* &
                                        ((1-k1)*psiy1(i,j-1,jg+3)*(C_new_y(i,j,jg,:)-C_new_y(i,j-1,jg,:))+ &
                                        (1+k1)*psiy2(i,j,jg+3)*(C_new_y(i+1,j,jg,:)-C_new_y(i,j,jg,:))      )
                            Cy_r(i,j,jg,:) = C_new_y(i,j,jg,:)  + 0.25*vareps* &
                                        ((1-k1)*psiy2(i,j+1,jg+3)*(C_new_y(i,j+2,jg,:)-C_new_y(i,j+1,jg,:)) + &
                                        (1+k1)*psiy1(i,j,jg+3)*(C_new_y(i,j+1,jg,:)-C_new_y(i,j,jg,:))      )
                    enddo
                enddo
            enddo
            do i=1-mbc, mx+mbc
                do j=1-mbc, my+mbc
                    if (method == 'SVL') then

                        Sus(i,j,:)  = flux_calc_SVL(Vx_l(i,j,1),Vx_l(i,j,2),Vx_l(i,j,3),Cx_l(i,j,:,:),Vx_r(i,j,1), &
                                        Vx_r(i,j,2),Vx_r(i,j,3),Cx_r(i,j,:,:),gmax)
                        Svs(i,j,:)  = flux_calc_SVL(Vy_l(i,j,1),Vy_l(i,j,2),Vy_l(i,j,3),Cy_l(i,j,:,:),Vy_r(i,j,1), &
                                        Vy_r(i,j,2),Vy_r(i,j,3),Cy_r(i,j,:,:),gmax)
                        Sub(i,j,:)  = flux_calc_SVL(Vx_l(i,j,1),Vx_l(i,j,2),Vx_l(i,j,3),Cx_l(i,j,:,:),Vx_r(i,j,1), &
                                        Vx_r(i,j,2),Vx_r(i,j,3),Cx_r(i,j,:,:),gmax)
                        Svb(i,j,:)  = flux_calc_SVL(Vy_l(i,j,1),Vy_l(i,j,2),Vy_l(i,j,3),Cy_l(i,j,:,:),Vy_r(i,j,1), &
                                        Vy_r(i,j,2),Vy_r(i,j,3),Cy_r(i,j,:,:),gmax)
                    else
                        print *, 'WARNING: DECIDE TURN ON OR OFF FLUX LIMITER'
                    endif
                enddo
            enddo

        end subroutine Flux_vector
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used to Standard Van Leer method
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function flux_calc_SVL(u_l,v_l,h_l,c_l,u_r,v_r,h_r,c_r,gmax)

            use sediment_module,only: g,eps

            implicit none
            integer :: gmax
            Real(kind=Prec)  :: u_l,v_l,h_l,u_r,v_r,h_r
            Real(kind=Prec)  :: speed_l,speed_r,fr_l,fr_r,Fr1,Fr2,Beta_r,Beta_l,Alpha1,Alpha2,C1,C2,one
            Real(kind=Prec), Dimension(gmax) :: flux_calc_SVL
            Real(kind=Prec), Dimension(gmax) :: c_l, c_r
            integer :: jg
            one = 1.0
            speed_l = dsqrt(dabs(g*h_l))
            speed_r = dsqrt(dabs(g*h_r))
            if (h_l .lt. eps) then
                fr_l = 0.0
            else
                fr_l     = u_l/speed_l
            end if
            if (h_r .lt. eps) then
                fr_l = 0.0
            else
                fr_r     = u_r/speed_r
            end if
            Fr1      = 0.25*(fr_l+1.0)**2.0
            Fr2      = -0.25*(fr_r-1.0)**2.0
            Beta_l  = -max(0,1-INT(dabs(fr_l)))
            Beta_r  = -max(0,1-INT(dabs(fr_r)))
            Alpha1  = 0.5*(1.0+SIGN(one,fr_l))
            Alpha2  = 0.5*(1.0-SIGN(one,fr_r))
            C1      = Alpha1*(1.0+Beta_l)*fr_l-Beta_l*Fr1
            C2      = Alpha2*(1.0+Beta_r)*fr_r-Beta_r*Fr2
            do jg = 1, gmax
                flux_calc_SVL(jg) = c_l(jg)*speed_l*C1+c_r(jg)*speed_r*C2
            end do
        end function flux_calc_SVL
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used to put flux limiter
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine flux_limitor(mbc,mx,my,V_new,C_new_x,C_new_y,psix1,psix2,psiy1,psiy2)

            use Set_Precision, only: Prec
            use sediment_module, only: toler,limit_method,gmax,vareps
            implicit none

            ! Arguments
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: V_new(1-mbc:mx+mbc,1-mbc:my+mbc,3),C_new_x(1-mbc:mx+mbc,1-mbc:my+mbc,gmax,2),&
                                        C_new_y(1-mbc:mx+mbc,1-mbc:my+mbc,gmax,2)
            !output
            real(kind=Prec), intent(inout) :: psix1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3),psix2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3), &
                         psiy1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3),psiy2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3)

            !local

            integer :: i,j,k,ii
            Real(kind=Prec),Dimension(1-mbc:mx+mbc-1,1-mbc:my+mbc-1,gmax+3) :: Rx1,Rx2 !
            Real(kind=Prec),Dimension(1-mbc:mx+mbc-1,1-mbc:my+mbc-1,gmax+3) :: Ry1,Ry2  !
            Real(kind=Prec),Dimension(1-mbc:mx+mbc-1,1-mbc:my+mbc-1,gmax+3) :: denx, deny
            Real(kind=Prec),Dimension(1-mbc:mx+mbc,1-mbc:my+mbc-1,gmax) :: C_new_x_t,C_new_y_t


            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do ii = 1, gmax
                        C_new_x_t(i,j,ii) = sum(C_new_x(i,j,ii,:))
                        C_new_y_t(i,j,ii) = sum(C_new_y(i,j,ii,:))
                    enddo
                enddo
            enddo

            do i = 1-mbc, mx+mbc-1
                do j = 1-mbc, my+mbc-1
                    do ii = 1, 3
                        denx(i,j,ii) = V_new(i+1,j,ii)-V_new(i,j,ii)
                        denx(i,j,ii) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))
                        deny(i,j,ii) = V_new(i,j+1,ii)-V_new(i,j,ii)
                        deny(i,j,ii) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))
                        Rx1(i,j,ii) = (V_new(i+2,j,ii)-V_new(i+1,j,ii))/denx(i,j,ii)
                        Rx2(i,j,ii) = (V_new(i,j,ii)-V_new(i-1,j,ii))/denx(i,j,ii)
                        Ry1(i,j,ii) = (V_new(i,j+2,ii)-V_new(i,j+1,ii))/deny(i,j,ii)
                        Ry2(i,j,ii) = (V_new(i,j,ii)-V_new(i,j-1,ii))/deny(i,j,ii)
                    enddo
                    do ii = 1, gmax
                        denx(i,j,ii+3) = C_new_x_t(i+1,j,ii)-C_new_x_t(i,j,ii)
                        denx(i,j,ii+3) = sign(max(abs(denx(i,j,ii)),toler),denx(i,j,ii))
                        deny(i,j,ii+3) = C_new_y_t(i,j+1,ii)-C_new_y_t(i,j,ii)
                        deny(i,j,ii+3) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))
                        Rx1(i,j,ii+3) = (C_new_x_t(i+2,j,ii)-C_new_x_t(i+1,j,ii))/denx(i,j,ii)
                        Rx2(i,j,ii+3) = (C_new_x_t(i,j,ii)-C_new_x_t(i-1,j,ii))/denx(i,j,ii)
                        Ry1(i,j,ii+3) = (C_new_y_t(i,j+2,ii)-C_new_y_t(i,j+1,ii))/deny(i,j,ii)
                        Ry2(i,j,ii+3) = (C_new_y_t(i,j,ii)-C_new_y_t(i,j-1,ii))/deny(i,j,ii)
                    enddo
                    if (limit_method == 'Vanleer') then
                        psix1(i,j,:) = (Rx1(i,j,:)+abs(Rx1(i,j,:)))/(1.0+Rx1(i,j,:))
                        psix2(i,j,:) = (Rx2(i,j,:)+abs(Rx2(i,j,:)))/(1.0+Rx2(i,j,:))
                        psiy1(i,j,:) = (Ry1(i,j,:)+abs(Ry1(i,j,:)))/(1.0+Ry1(i,j,:))
                        psiy2(i,j,:) = (Ry2(i,j,:)+abs(Ry2(i,j,:)))/(1.0+Ry2(i,j,:))
                    elseif (limit_method == 'VanAlbada') then
                        psix1(i,j,:) = (Rx1(i,j,:)+(Rx1(i,j,:))**2.0)/(1.0+(Rx1(i,j,:))**2.0)
                        psix2(i,j,:) = (Rx2(i,j,:)+(Rx2(i,j,:))**2.0)/(1.0+(Rx2(i,j,:))**2.0)
                        psiy1(i,j,:) = (Ry1(i,j,:)+(Ry1(i,j,:))**2.0)/(1.0+(Ry1(i,j,:))**2.0)
                        psiy2(i,j,:) = (Ry2(i,j,:)+(Ry2(i,j,:))**2.0)/(1.0+(Ry2(i,j,:))**2.0)
                    elseif (limit_method == 'minmod') then
                        do ii = 1, gmax+3
                            if (Rx1(i,j,ii)>0) then
                                psix1(i,j,ii) = min(1.0, Rx1(i,j,ii))
                            else
                                psix1(i,j,ii) = 0.0
                            endif
                            if (Rx2(i,j,ii)>0) then
                                psix2(i,j,ii) = min(1.0, Rx2(i,j,ii))
                            else
                                psix2(i,j,ii) = 0.0
                            endif
                            if (Ry1(i,j,ii)>0) then
                                psiy1(i,j,ii) = min(1.0, Ry1(i,j,ii))
                            else
                                psiy1(i,j,ii) = 0.0
                            endif
                            if (Ry2(i,j,ii)>0) then
                                psiy2(i,j,ii) = min(1.0, Ry2(i,j,ii))
                            else
                                psiy2(i,j,ii) = 0.0
                            endif
                        enddo
                    endif
                enddo
            enddo
            psix1(2-mbc,:,:)=psix1(3-mbc,:,:)
            psix1(1-mbc,:,:)=psix1(2-mbc,:,:)
            psix1(mx,:,:)=psix1(mx-1,:,:)
            psix1(mx+mbc-1,:,:)=psix1(mx,:,:)
            psix1(mx+mbc,:,:)=psix1(mx+mbc-1,:,:)
            psix2(2-mbc,:,:)=psix2(3-mbc,:,:)
            psix2(1-mbc,:,:)=psix2(2-mbc,:,:)
            psix2(mx,:,:)=psix2(mx-1,:,:)
            psix2(mx+mbc-1,:,:)=psix2(mx,:,:)
            psix2(mx+mbc,:,:)=psix1(mx+mbc-1,:,:)
            psix1(:,2-mbc,:)=psix1(:,3-mbc,:)
            psix1(:,1-mbc,:)=psix1(:,2-mbc,:)
            psix1(:,my,:)=psix1(:,my-1,:)
            psix1(:,my+mbc-1,:)=psix1(:,my,:)
            psix1(:,my+mbc,:)=psix1(:,my+mbc-1,:)
            psix2(:,2-mbc,:)=psix2(:,3-mbc,:)
            psix2(:,1-mbc,:)=psix2(:,2-mbc,:)
            psix2(:,my,:)=psix2(:,my-1,:)
            psix2(:,my+mbc-1,:)=psix2(:,my,:)
            psix2(:,my+mbc,:)=psix2(:,my+mbc,:)
            psiy1(2-mbc,:,:)=psiy1(3-mbc,:,:)
            psiy1(1-mbc,:,:)=psiy1(2-mbc,:,:)
            psiy1(mx,:,:)=psiy1(mx-1,:,:)
            psiy1(mx+mbc-1,:,:)=psiy1(mx,:,:)
            psiy1(mx+mbc,:,:)=psiy1(mx+mbc-1,:,:)
            psiy2(2-mbc,:,:)=psiy2(3-mbc,:,:)
            psiy2(1-mbc,:,:)=psiy2(2-mbc,:,:)
            psiy2(mx,:,:)=psiy2(mx-1,:,:)
            psiy2(mx+mbc-1,:,:)=psiy2(mx,:,:)
            psiy2(mx+mbc,:,:)=psiy2(mx+mbc-1,:,:)
            psiy1(:,2-mbc,:)=psiy1(:,3-mbc,:)
            psiy1(:,1-mbc,:)=psiy1(:,2-mbc,:)
            psiy1(:,my,:)=psiy1(:,my-1,:)
            psiy1(:,my+mbc-1,:)=psiy1(:,my,:)
            psiy1(:,my+mbc,:)=psiy1(:,my+mbc-1,:)
            psiy2(:,2-mbc,:)=psiy2(:,3-mbc,:)
            psiy2(:,1-mbc,:)=psiy2(:,2-mbc,:)
            psiy2(:,my,:)=psiy2(:,my-1,:)
            psiy2(:,my+mbc-1,:)=psiy2(:,my,:)
            psiy2(:,my+mbc,:)=psiy2(:,my+mbc,:)
            if (vareps ==0)  then
                psix1 = 0.0
                psix2 = 0.0
                psiy1 = 0.0
                psiy2 = 0.0
            end if
        end subroutine flux_limitor
    end module flux
