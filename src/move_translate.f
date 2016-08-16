c***********************************************************************
      subroutine move_translate(init,x,y,xv,yv,xold,yold,xole,yole,ug,vg,jd,kd)
c

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd), xv(jmax,kmax), yv(jmax,kmax)
      real xold(jmax,kmax),yold(jmax,kmax),xole(jmax,kmax),yole(jmax,kmax)

      ! local variables

      real xnew,ynew,dx,xo,yo
      integer j,k
c***********************************************************************


        dx=rotf*dt

      if (init.eq.1) dx = 0
      print*,istep,dx

c.. rotate background grid in global sense about origin


!...  !update cell center co-ordinates
      !-------------------------------
      do k=1,kd
        do j=1,jd
          xnew = x(j,k) + dx 
          ynew = y(j,k)

          !grid velocities cell centers
          if(ntac.eq.1.or.init.eq.1) then
            ug(j,k)=(xnew-x(j,k))/dt
            vg(j,k)=(ynew-y(j,k))/dt
          else
            xo=x(j,k)-dx
            yo=y(j,k)
            ug(j,k)=(1.5*xnew-2*x(j,k)+0.5*xo)/dt
            vg(j,k)=(1.5*ynew-2*y(j,k)+0.5*yo)/dt
          endif

          x(j,k)=xnew
          y(j,k)=ynew
          
        end do !jd
      end do !kd
      !-------------------------------
      

!...  !update cell vertex co-ordinates
      !-------------------------------
      do k = 1,kmax
      do j = 1,jmax

c..store the old mesh for time metrics

        xole(j,k) = xold(j,k)
        yole(j,k) = yold(j,k)

        xold(j,k) = xv(j,k)
        yold(j,k) = yv(j,k)

c..update mesh

        xv(j,k) = xold(j,k) + dx 
        yv(j,k) = yold(j,k) 

      enddo
      enddo
!      if(istep.eq.50) then
!       write(3000) 2
!       write(3000) jmax,kmax,jd,kd
!       write(3000) ((xv(j,k),j=1,jmax),k=1,kmax),
!     >             ((yv(j,k),j=1,jmax),k=1,kmax)
!       write(3000) ((x(j,k),j=1,jd),k=1,kd),
!     >             ((y(j,k),j=1,jd),k=1,kd)
!      end if


      return
      end subroutine move_translate

