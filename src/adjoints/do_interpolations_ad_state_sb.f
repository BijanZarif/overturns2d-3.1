c***********************************************************************
      subroutine do_interpolations_ad_state(qg,sg,jmx,kmx,ibcg,imeshg,idonorg,
     c        fracg,nfringeg,ndonorg,iisptrg,iieptrg,idsize,qsize,ssize,nmesh)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer idsize,qsize,ssize,nmesh
      integer jmx(nmesh),kmx(nmesh)
      integer nfringeg(nmesh,nspec),ndonorg(nmesh,nspec)
      integer iisptrg(nmesh,nspec),iieptrg(nmesh,nspec)
      real qg(qsize),sg(ssize) !sg :global adjoint vector/state
      integer imeshg(idsize,2,nmesh,nspec),idonorg(idsize,2,nmesh,nspec)
      integer ibcg(idsize,nmesh,nspec)
      real fracg(idsize,2,nmesh,nspec)

c..local variables

      integer bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im
      integer qptr,sptr,nfringe,ndonor,iisptr,iieptr,nsp
      real djm1,dj0,djp1,dkm1,dk0,dkp1
      real w1,w2,w3,w4,w5,w6,w7,w8,w9,onefourth
      real,allocatable :: sbc(:,:)
      real,allocatable :: q(:,:,:),s(:,:,:)

      onefourth = 1./4

!!!$OMP PARALLEL IF (NSPEC > 1)
!!
!!!$OMP DO
!!!$OMP& PRIVATE(bcdim,ii,jj,kk,iim,jjm,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im)
!!!$OMP& PRIVATE(qptr,sptr,nfringe,ndonor,iisptr,iieptr)
!!!$OMP& PRIVATE(djm1,dj0,djp1,dkm1,dk0,dkp1)
!!!$OMP& PRIVATE(w1,w2,w3,w4,w5,w6,w7,w8,w9)
!!!$OMP& PRIVATE(q,s,sbc)

      spectralloop:do nsp = 1,nspec

        bcdim=iieptrg(nmesh,nsp)
        allocate(sbc(bcdim,nv))

c...LOOP THROUGH ALL THE MESHES AND COLLECT GLOBAL SBC (or QBC)

        qptr = 1
        sptr = 1
        do im = 1,nmesh

          ndonor = ndonorg(im,nsp)
          iisptr = iisptrg(im,nsp)
          iieptr = iieptrg(im,nsp)
          jmax = jmx(im); kmax = kmx(im)
          
          allocate(q(jmax,kmax,nq))
          allocate(s(jmax,kmax,nv))
           
c.....  assign local q from global values

          do n = 1,nq
            do k=1,kmax
              do j=1,jmax
                q(j,k,n) = qg(qptr-1 + jmax*kmax*nspec*(n-1) + 
     &                     jmax*kmax*(nsp-1) + jmax*(k-1) + j)
              enddo
            enddo
          enddo
c.....  assign local s from global values

          do n = 1,nv
            do k=1,kmax
              do j=1,jmax
                s(j,k,n) = sg(sptr-1 + jmax*kmax*nspec*(n-1) + 
     &                     jmax*kmax*(nsp-1) + jmax*(k-1) + j)
              enddo
            enddo
          enddo

          qptr = qptr + jmax*kmax*nspec*nq
          sptr = sptr + jmax*kmax*nspec*nv
   
          do id=1,ndonor

            ii=idonorg(id,1,im,nsp)
            jj=idonorg(id,2,im,nsp)

            iim=ii-1
            jjm=jj-1

            iip=ii+1
            jjp=jj+1

            djm1 = 1.+fracg(id,1,im,nsp)
            dkm1 = 1.+fracg(id,2,im,nsp)
            dj0  = fracg(id,1,im,nsp)
            dk0  = fracg(id,2,im,nsp)
            djp1 = 1.-fracg(id,1,im,nsp)
            dkp1 = 1.-fracg(id,2,im,nsp)

            w1   =    (dj0  * djp1 * dk0  * dkp1)*onefourth
            w2   = -2*(djm1 * djp1 * dk0  * dkp1)*onefourth
            w3   = -  (djm1 * dj0  * dk0  * dkp1)*onefourth
            w4   = -2*(dj0  * djp1 * dkm1 * dkp1)*onefourth
            w5   =  4*(djm1 * djp1 * dkm1 * dkp1)*onefourth
            w6   =  2*(djm1 * dj0  * dkm1 * dkp1)*onefourth
            w7   = -  (dj0  * djp1 * dkm1 * dk0 )*onefourth
            w8   =  2*(djm1 * djp1 * dkm1 * dk0 )*onefourth
            w9   =    (djm1 * dj0  * dkm1 * dk0 )*onefourth

c.....collect in global pointer sbc from pointer iisptr->iieptr

            do n=1,nmv
               sbc(iisptr-1+id,n)= 
     &                     w1*s(iim,jjm,n)*q(iim,jjm,nq) 
     &                   + w2*s(ii ,jjm,n)*q(ii ,jjm,nq)
     &                   + w3*s(iip,jjm,n)*q(iip,jjm,nq) 
     &                   + w4*s(iim,jj ,n)*q(iim,jj ,nq)
     &                   + w5*s(ii ,jj ,n)*q(ii ,jj ,nq)
     &                   + w6*s(iip,jj ,n)*q(iip,jj ,nq)
     &                   + w7*s(iim,jjp,n)*q(iim,jjp,nq)
     &                   + w8*s(ii ,jjp,n)*q(ii ,jjp,nq)
     &                   + w9*s(iip,jjp,n)*q(iip,jjp,nq)
            enddo

            do n=nmv+1,nv
               sbc(iisptr-1+id,n)= 
     &                     w1*s(iim,jjm,n)
     &                   + w2*s(ii ,jjm,n)
     &                   + w3*s(iip,jjm,n)
     &                   + w4*s(iim,jj ,n)
     &                   + w5*s(ii ,jj ,n)
     &                   + w6*s(iip,jj ,n)
     &                   + w7*s(iim,jjp,n)
     &                   + w8*s(ii ,jjp,n)
     &                   + w9*s(iip,jjp,n)
            enddo

!amm            do n=1,nv
!amm               sbc(iisptr-1+id,n)= 
!amm     &                     w1*q(iim,jjm,n)
!amm     &                   + w2*q(ii ,jjm,n)
!amm     &                   + w3*q(iip,jjm,n)
!amm     &                   + w4*q(iim,jj ,n)
!amm     &                   + w5*q(ii ,jj ,n)
!amm     &                   + w6*q(iip,jj ,n)
!amm     &                   + w7*q(iim,jjp,n)
!amm     &                   + w8*q(ii ,jjp,n)
!amm     &                   + w9*q(iip,jjp,n)
!amm            enddo

          enddo

          deallocate(q)
          deallocate(s)
        enddo

c...LOOP THROUGH ALL MESHES AND UPDATE VALUES FROM QBC ARRAY

        qptr = 1
        sptr = 1
        do im = 1,nmesh

          nfringe = nfringeg(im,nsp)
          jmax = jmx(im); kmax = kmx(im)

          allocate(q(jmax,kmax,nq)) !state has nq>nv
          allocate(s(jmax,kmax,nv))

c.....re-assign local q from global values
          do n=1,nq
            do k=1,kmax
              do j=1,jmax
                q(j,k,n) = qg(qptr-1 + jmax*kmax*nspec*(n-1) + 
     &                      jmax*kmax*(nsp-1) + jmax*(k-1) + j)
              enddo
            enddo
          enddo

c.....re-assign local s from global values
          do n=1,nv
            do k=1,kmax
              do j=1,jmax
                s(j,k,n) = sg(sptr-1 + jmax*kmax*nspec*(n-1) + 
     &                      jmax*kmax*(nsp-1) + jmax*(k-1) + j)
              enddo
            enddo
          enddo
c.....over write fringe points solution w/ donor global sbc solution

          do id=1,nfringe
            j = ibcg(id,im,nsp)
            ii = imeshg(id,1,im,nsp)
            jj = imeshg(id,2,im,nsp)

            do n = 1,nmv
               s(ii,jj,n) = sbc(j,n)/q(ii,jj,nq)
            enddo

            do n = nmv+1,nv
               s(ii,jj,n) = sbc(j,n)
            enddo

!amm            do n = 1,nv
!amm               q(ii,jj,n) = qbc(j,n)
!amm            enddo

          enddo
       
c.....reassign sbc to global s (containing all Ng meshes)
          do n=1,nv
            do k=1,kmax
              do j=1,jmax
                sg(sptr-1 + jmax*kmax*nspec*(n-1) + jmax*kmax*(nsp-1) + 
     &               jmax*(k-1) + j) = s(j,k,n)
              enddo
            enddo
          enddo

          qptr = qptr + jmax*kmax*nspec*nq
          sptr = sptr + jmax*kmax*nspec*nv

          deallocate(q)
          deallocate(s)
        enddo 

        deallocate(sbc)

      enddo spectralloop
!!!$OMP END DO

!!!$OMP END PARALLEL

      return
      end subroutine do_interpolations_ad_state
