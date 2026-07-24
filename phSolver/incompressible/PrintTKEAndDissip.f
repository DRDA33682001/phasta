      subroutine PrintTKEAndDissip (y, x, shp, shgl, istp, t)

      use pointer_data          ! brings in mien(iblk)%p (connectivity)

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension y(nshg,ndof),               x(numnp,nsd),
     &          shp(MAXTOP,maxsh,MAXQPT),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT)

      integer istp
      real*8  t

      real*8, allocatable :: vol_e(:), Ek_e(:), Enstr_e(:),
     &                       Viscdissip_e(:), divu_e(:),
     &                       graddivdissip_e(:)
      real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)

      real*8  vol, Ek, Enstr, Viscdissip, divu, graddivdissip
      real*8  volTotal, EkTotal, EnstrTotal, ViscdissipTotal, dissip
      real*8  divuTotal, SGSdissipTotal, fdotudissipTotal, graddivdissipTotal
      real*8  nuu

      vol        = zero
      Ek         = zero
      Enstr      = zero
      Viscdissip = zero
      divu       = zero
      graddivdissip = zero
      nuu        = datmat(1,2,1)   ! molecular viscosity (rho=1 => = nu)

c
c.... loop over the element blocks (same pattern as ElmGMR)
c
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         lcsyst = lcblk(3,iblk)
         nenl   = lcblk(5,iblk)        ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         npro   = lcblk(1,iblk+1) - iel
         ngauss = nint(lcsyst)

         allocate ( vol_e(npro)        )
         allocate ( Ek_e(npro)         )
         allocate ( Enstr_e(npro)      )
         allocate ( Viscdissip_e(npro) )
         allocate ( divu_e(npro)       )
         allocate ( graddivdissip_e(npro) )
         
         ! Allocate contiguous arrays for shape functions
         allocate ( tmpshp(nshl,ngauss) )
         allocate ( tmpshgl(nsd,nshl,ngauss) )

         ! Perform copy-in
         tmpshp(1:nshl,1:ngauss) = shp(lcsyst,1:nshl,1:ngauss)
         tmpshgl(:,1:nshl,1:ngauss) = shgl(lcsyst,:,1:nshl,1:ngauss)

         call ComputeElmTKEAndDissip ( y, x,
     &              tmpshp, tmpshgl,
     &              mien(iblk)%p,
     &              vol_e, Ek_e, Enstr_e, Viscdissip_e, divu_e,
     &              graddivdissip_e )

         do i = 1, npro
            vol        = vol        + vol_e(i)
            Ek         = Ek         + Ek_e(i)
            Enstr      = Enstr      + Enstr_e(i)
            Viscdissip = Viscdissip + Viscdissip_e(i)
            divu       = divu       + divu_e(i)
            graddivdissip = graddivdissip + graddivdissip_e(i)
         enddo

         deallocate ( tmpshp       )
         deallocate ( tmpshgl      )
         deallocate ( vol_e        )
         deallocate ( Ek_e         )
         deallocate ( Enstr_e      )
         deallocate ( Viscdissip_e )
         deallocate ( divu_e       )
         deallocate ( graddivdissip_e )
      enddo

c
c.... reduce across MPI ranks
c
      volTotal        = vol
      EkTotal         = Ek
      EnstrTotal      = Enstr
      ViscdissipTotal = Viscdissip
      divuTotal       = divu
      graddivdissipTotal = graddivdissip

      if (numpe > 1) then
         call drvAllreducesclr ( vol,        volTotal        )
         call drvAllreducesclr ( Ek,         EkTotal         )
         call drvAllreducesclr ( Enstr,      EnstrTotal      )
         call drvAllreducesclr ( Viscdissip, ViscdissipTotal )
         call drvAllreducesclr ( divu,       divuTotal       )
         call drvAllreducesclr ( graddivdissip, graddivdissipTotal )
      endif

c
c.... volume-average and form dissipation
c
      EkTotal         = EkTotal         / volTotal
      EnstrTotal      = EnstrTotal      / volTotal
      dissip          = 2.0d0 * nuu * EnstrTotal
      ViscdissipTotal = ViscdissipTotal / volTotal
      divuTotal       = divuTotal       / volTotal
      graddivdissipTotal = graddivdissipTotal / volTotal
      
c.... Variables unavailable in this minimal build are set to 0.0 
c.... to match the 8-column format safely.
      SGSdissipTotal     = zero
      
c.... Calculate fdotudissipTotal based on linear forcing A_force * (2*EkTotal)
      if (EnableForcing .eq. 1) then
         if (current_TKE .gt. 1.0d-8) then
            fdotudissipTotal = (PowerInput / current_TKE) * EkTotal
         else
            fdotudissipTotal = zero
         endif
      else
         fdotudissipTotal = zero
      endif
c.... Update the global TKE for the next time step's forcing calculation
      current_TKE = EkTotal

c
c.... master rank appends one line
c
      if (myrank .eq. master) then
         open (unit=195, file="TKE-dissip.dat", access="append",
     &         status="unknown")
         write(195,*) t, EkTotal, dissip, ViscdissipTotal, 
     &                SGSdissipTotal, fdotudissipTotal, 
     &                graddivdissipTotal, divuTotal
         close(195)
      endif

      return
      end



c-----------------------------------------------------------------------
c
c  ComputeElmTKEAndDissip : per-element-block work.
c  Gathers the local solution + coordinates, loops the quadrature points,
c  and accumulates the volume, kinetic energy, enstrophy, and strain-based
c  viscous dissipation for every element in the block.
c
c-----------------------------------------------------------------------
      subroutine ComputeElmTKEAndDissip ( y, x, shp, shgl, ien,
     &                        vol_e, Ek_e, Enstr_e, Viscdissip_e, divu_e,
     &                        graddivdissip_e )

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension y(nshg,ndof),               x(numnp,nsd),
     &          shp(nshl,ngauss),           shgl(nsd,nshl,ngauss),
     &          ien(npro,nshl)

      dimension yl(npro,nshl,ndof),         xl(npro,nenl,nsd),
     &          sgn(npro,nshl),             shape(npro,nshl),
     &          shdrv(npro,nsd,nshl)

      dimension g1yi(npro,nflow),           g2yi(npro,nflow),
     &          g3yi(npro,nflow),           shg(npro,nshl,nsd),
     &          dxidx(npro,nsd,nsd),        WdetJ(npro)

      dimension u1(npro),  u2(npro),  u3(npro),  Ek_tmp(npro),
     &          vort(npro,nsd),             Enstr_tmp(npro),
     &          Sij(npro,nsd,nsd),          Snorm(npro),
     &          divu_tmp(npro)

      dimension vol_e(npro), Ek_e(npro), Enstr_e(npro),
     &          Viscdissip_e(npro), divu_e(npro),
     &          graddivdissip_e(npro)

      dimension gijd(npro,6), tauM(npro), tauC(npro), fact(npro)
      real*8    dts, rho, rnu, fff, ff
      real*8    nuu

c
c.... mode-sign matrix for the hierarchic basis (as in AsIq)
c
      do i = 1, nshl
         where ( ien(:,i) < 0 )
            sgn(:,i) = -one
         elsewhere
            sgn(:,i) =  one
         endwhere
      enddo

c
c.... gather local solution and coordinates
c
      call localy (y, yl, ien, ndof, 'gather  ')
      call localx (x, xl, ien, nsd,  'gather  ')

      nuu          = datmat(1,2,1)
      rho          = datmat(1,1,1)
      rnu          = nuu / rho
      dts          = Dtgl * dtsfct
      ff           = taucfct / dtsfct
      
      if(lcsyst.eq.1) then
         fff = 36.0d0
      elseif(lcsyst.eq.2) then
         fff = 60.0d0
      else
         fff = 128.0d0
      endif

      vol_e        = zero
      Ek_e         = zero
      Enstr_e      = zero
      Viscdissip_e = zero
      divu_e       = zero
      graddivdissip_e = zero

c
c.... loop over the integration points (intp is the common counter that
c     getshp / e3qvar read)
c
      do intp = 1, ngauss
         if (Qwt(lcsyst,intp) .eq. zero) cycle      ! precaution

c....    shape functions / local gradients at this quad point
         call getshp ( shp, shgl, sgn, shape, shdrv )

c....    velocity gradients (g*yi), global shape grads, weighted Jacobian
         call e3qvar ( yl, shdrv, xl, g1yi, g2yi, g3yi,
     &                 shg, dxidx, WdetJ )

c....    volume
         vol_e(:) = vol_e(:) + WdetJ(:)

c....    velocity at the integration point  (y(:,2:4) = u1,u2,u3)
         u1 = zero
         u2 = zero
         u3 = zero
         do n = 1, nshl
            u1 = u1 + shape(:,n) * yl(:,n,2)
            u2 = u2 + shape(:,n) * yl(:,n,3)
            u3 = u3 + shape(:,n) * yl(:,n,4)
         enddo
         Ek_tmp = u1*u1 + u2*u2 + u3*u3
         Ek_e   = Ek_e + 0.5d0 * Ek_tmp * WdetJ

c....    vorticity  w = curl(u)  from the velocity gradients
c        g1yi(:,k)=d(y_k)/dx1, g2yi=d/dx2, g3yi=d/dx3 ; k=2,3,4 -> u1,u2,u3
         vort(:,1) = g2yi(:,4) - g3yi(:,3)   ! du3/dx2 - du2/dx3
         vort(:,2) = g3yi(:,2) - g1yi(:,4)   ! du1/dx3 - du3/dx1
         vort(:,3) = g1yi(:,3) - g2yi(:,2)   ! du2/dx1 - du1/dx2
         Enstr_tmp = vort(:,1)*vort(:,1) + vort(:,2)*vort(:,2)
     &             + vort(:,3)*vort(:,3)
         Enstr_e   = Enstr_e + 0.5d0 * Enstr_tmp * WdetJ

c....    strain-rate tensor S_ij and its norm S_ij S_ij
         Sij        = zero
         Sij(:,1,1) =        g1yi(:,2)
         Sij(:,2,2) =        g2yi(:,3)
         Sij(:,3,3) =        g3yi(:,4)
         Sij(:,1,2) = 0.5d0*(g2yi(:,2) + g1yi(:,3))
         Sij(:,1,3) = 0.5d0*(g3yi(:,2) + g1yi(:,4))
         Sij(:,2,3) = 0.5d0*(g3yi(:,3) + g2yi(:,4))

         Snorm = sqrt( Sij(:,1,1)*Sij(:,1,1)
     &               + Sij(:,2,2)*Sij(:,2,2)
     &               + Sij(:,3,3)*Sij(:,3,3)
     &               + 2.0d0*( Sij(:,1,2)*Sij(:,1,2)
     &                       + Sij(:,1,3)*Sij(:,1,3)
     &                       + Sij(:,2,3)*Sij(:,2,3) ) )

         Viscdissip_e = Viscdissip_e
     &                + 2.0d0 * nuu * (Snorm*Snorm) * WdetJ

c....    divergence of velocity
         divu_tmp = g1yi(:,2) + g2yi(:,3) + g3yi(:,4)
         divu_e   = divu_e + divu_tmp * WdetJ

c....    metric tensor (gijd)
         gijd(:,1) = dxidx(:,1,1)**2 + dxidx(:,2,1)**2 + dxidx(:,3,1)**2
         gijd(:,2) = dxidx(:,1,2)**2 + dxidx(:,2,2)**2 + dxidx(:,3,2)**2
         gijd(:,3) = dxidx(:,1,3)**2 + dxidx(:,2,3)**2 + dxidx(:,3,3)**2
         gijd(:,4) = dxidx(:,1,1)*dxidx(:,1,2) + dxidx(:,2,1)*dxidx(:,2,2)
     &             + dxidx(:,3,1)*dxidx(:,3,2)
         gijd(:,5) = dxidx(:,1,2)*dxidx(:,1,3) + dxidx(:,2,2)*dxidx(:,2,3)
     &             + dxidx(:,3,2)*dxidx(:,3,3)
         gijd(:,6) = dxidx(:,1,1)*dxidx(:,1,3) + dxidx(:,2,1)*dxidx(:,2,3)
     &             + dxidx(:,3,1)*dxidx(:,3,3)

c....    momentum stabilization parameter (tauM)
         tauM = ( (2.0d0 * dts)**2 
     &        + ( u1 * ( gijd(:,1) * u1 + gijd(:,4) * u2 + gijd(:,6) * u3 )
     &          + u2 * ( gijd(:,4) * u1 + gijd(:,2) * u2 + gijd(:,5) * u3 )
     &          + u3 * ( gijd(:,6) * u1 + gijd(:,5) * u2 + gijd(:,3) * u3 ) ) )
     &        + fff * (rnu**2) * ( gijd(:,1)**2 + gijd(:,2)**2 + gijd(:,3)**2 
     &          + 2.0d0 * ( gijd(:,4)**2 + gijd(:,5)**2 + gijd(:,6)**2 ) )
         
         fact = sqrt(tauM)
         fact = 1.0d0 / fact

c....    continuity stabilization parameter (tauC) and dissipation
         tauC = rho * pt125 * fact / (gijd(:,1)+gijd(:,2)+gijd(:,3)) * ff

         graddivdissip_e = graddivdissip_e
     &                   + tauC * (divu_tmp**2) * WdetJ

      enddo

      return
      end
