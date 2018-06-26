
      subroutine setupdata(nx, ny, nz, cx)
      complex*16 cx(1,1,1)
	real*8 tpi
c
      tpi=8.d0*datan(1.d0)

      do iz = 1, nz 
       z = float(iz-1)*10.0/float(nz)-5.0
       do iy = 1, ny 
        y = float(iy-1)*10.0/float(ny)-5.0
        do ix = 1, nx
          x = float(ix-1)*10.0/float(nx)-5.0
          tmp  = (sin(tpi*x)*sin(tpi*y)+1.)*exp(-z**2)
          cx(ix,iy, iz) = cmplx(tmp, tmp)
        enddo
       enddo
      enddo
      return
      end

