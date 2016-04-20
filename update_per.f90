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
    module update

        use Set_precision, only: Prec
        use set_variable
        !use sed

        IMPLICIT NONE

        contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This Part is used for updating the grain size distribution for each cell
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine update_fractions(is,js,dzs,pbs,edgs,dzbs,dt)

            use sediment_module, only: lmax,gmax,por,nd_var,split,merge,toler,thick
            use Set_precision, only: Prec

            IMPLICIT NONE

            !Argument
            integer, intent(in) :: is, js
            real(Kind= Prec), intent(in)                              :: dzbs,dt
            real(Kind= Prec), dimension(gmax),intent(in)          :: edgs
            real(Kind= Prec) , dimension(lmax),intent(inout)    :: dzs
            real(Kind= Prec) , dimension(lmax,gmax),intent(inout) :: pbs
            !local
            integer                                         :: jg,jd,ki,mi,nl
            real(Kind= Prec)                                :: ED,zbold,dzbt,dzb_loc,fc,zero
            real(Kind= Prec) , dimension(:),allocatable,save                :: Ap,b
            real(Kind= Prec) , dimension(:,:),allocatable,save              :: Sm,As

            if (.not. allocated(Ap)) then
                allocate(Ap(lmax))
                allocate(b(lmax))
                allocate(Sm(lmax,gmax))
                allocate(As(lmax,3))
            endif
            dzb_loc = dzbs ! sediment thickness change
            !initialize Sm
            ED = sum(edgs) ! total sediment without pore

            mi = 0

            zero = 0.0
            !print *, dzb_loc
            !print *, is,js
            do while (abs(dzb_loc) .gt. 0.0)
                !print *, dzb_loc
                dzbt     = min(dzb_loc,dzs(1))                 ! make sure erosion (dzg is positive) is limited to thickness of variable layer
                dzbt     = max(dzbt,-dzs(1))                 ! make sure deposition (dzg is negative) is limited to thickness of first layer below variable layer
                fc = dzbt/dzbs                                    ! factor over mass change in cells to limit erosion and deposition
                dzb_loc = dzb_loc-dzbt                                 ! update dzb
                if (dzb_loc .LT. 0.0) then
                    dzb_loc = 0.0
                end if
                do ki=1,gmax

                    As=0.0
                    Ap=0.0
                    b=0.0

                    Sm(:,ki)=pbs(:,ki)*dzs*(1.0-por)                        ! Sm is the sediment mass per fraction per layer

                    select case(nd_var)
                    case(1)
                    !in this case: As=0
                    case(2)
                        As (1,1:3)= (/zero           ,  min(ED,zero) ,  max(ED,zero) /)
                        As (2,1:3)= (/-min(ED,zero)  , -max(ED,zero) ,  zero        /)
                    case(3:1000)
                        As (1,1:3)= (/zero           ,  min(ED,zero) ,  max(ED,zero) /)
                        As (2:nd_var-1,1)=-min(ED,zero)
                        As (2:nd_var-1,2)=-abs(ED)
                        As (2:nd_var-1,3)=max(ED,zero)
                        As (nd_var,1:3)= (/-min(ED,zero) , -max(ED,zero) ,   zero    /)
                    end select
                    !!!determine RHS

                    ! Ap makes sure that layer nd_var varies in thickness in time
                    ! Ap = 0 with single fraction (in case(1))
                    Ap(1) = sum(As(1,2:3)*pbs(1:2,ki)) !decide the distirbution used by depostion or erosion
                    do jd = 2,nd_var
                        Ap(jd) = sum(As (jd,:)*pbs(jd-1:jd+1,ki))
                    enddo
                    Ap(nd_var+1) = sum(As(nd_var+1,1:2)*pbs(nd_var:nd_var+1,ki))
                    ! b represents the actual erosion and deposition in the top layer.
                    ! However, the thickness of the top layer remains constant and instead layer nd_var will breath in thickness
                    b(1) = -edgs(ki)
                    !!!update Sm
                    ! Sm is the sediment mass per fraction per layer
                    Sm(1:nd_var+1,ki) = Sm(1:nd_var+1,ki) + dt*fc*(Ap(1:nd_var+1)+b(1:nd_var+1))
                enddo !fractions
                ! From Sm we can compute the new fraction ratios per layer and the layer thickness...

                do jd=1,nd_var+1
                    if (sum(Sm(jd,:))>0) then
                        pbs(jd,:) = Sm(jd,:)/sum(Sm(jd,:))
                    else
                        pbs(jd,:) = pbs(jd,:)
                    endif
                    dzs(jd) = sum(Sm(jd,:))/(1.0-por)
                enddo
                !!! modify grid
                !merge two upper layers in case of erosion
                if (dzs(nd_var) .lt. merge*dzs(nd_var+1)) then
                    forall (ki=1:gmax)
                        pbs(nd_var,ki) = (dzs(nd_var)*pbs(nd_var,ki) + dzs(nd_var+1)* &
                        pbs(nd_var+1,ki))/(dzs(nd_var)+dzs(nd_var+1))
                        pbs(nd_var+1:lmax-1,ki) = pbs(nd_var+2:lmax,ki)
                        pbs(lmax,ki) = pbs(lmax,ki)
                    endforall
                    z0bed(is,js) = z0bed(is,js)-dzs(nd_var+1)
                    dzs(nd_var) = dzs(nd_var+1)+dzs(nd_var)
                endif
                !split upper layer in case of sedimentation
                if (dzs(nd_var)>split*dzs(nd_var+1)) then
                    mi = mi + 1
                    pbs(nd_var+mi,:) = pbs(nd_var,:)
                    z0bed(is,js) = z0bed(is,js)+dzs(nd_var+mi)
                    dzs(nd_var) = dzs(nd_var)-dzs(nd_var+mi)
                endif
            enddo ! nt_sub
            pbs = max(0.0,min(pbs,1.0))
            zbold = zb(is,js)
            zb(is,js) = z0bed(is,js)+sum(dzs)
            sedero(is,js) = sedero(is,js)+(zb(is,js)-zbold)
            totalthick(is,js) = max(0.0,totalthick(is,js)+(zb(is,js)-zbold))
            pbbed(is,js,:,:)=pbs
            dzbed(is,js,:)= 0.0
            if (totalthick(is,js)>thick) then
                totalnum(is,js) = nint(totalthick(is,js)/thick)
                do nl= 1, totalnum(is,js)
                    dzbed(is,js,nl) = thick
                enddo
                if (totalthick(is,js)>totalnum(is,js)*thick) then
                    dzbed(is,js,totalnum(is,js)+1)=totalthick(is,js)-totalnum(is,js)*thick
                    totalnum(is,js)=totalnum(is,js)+1
                else
                    dzbed(is,js,totalnum(is,js))=totalthick(is,js)-(totalnum(is,js)-1)*thick
                endif
            else
                dzbed(is,js,1) =totalthick(is,js)
            endif
            if(dzbed(is,js,totalnum(is,js))<toler*thick) then
                dzbed(is,js,totalnum(is,js)) = 0.0
                totalnum(is,js) = totalnum(is,js) - 1
            endif
            do nl= 1, lmax
                if (dzbed(is,js,nl)<toler*thick) then
                    pbbed(is,is,nl,:) = 0.0
                endif
            enddo
        end subroutine update_fractions
    end module update










