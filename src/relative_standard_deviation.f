c##############################
c  relative standard deviations
c  of potential and field
c  see definition in Ref. in
c  introduction.f file 
c##############################
      subroutine realat_sd_poten(nr_frame, nat, nr_qm
     *, grom2cpmd, charge_activ
     *, max_nr_mm, nr_mm 
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot
     *, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z)
c
      implicit none
      integer nr_frame, nat, nr_qm
     *, grom2cpmd
     *, max_nr_mm, nr_mm
     *, i,j,k, k_c
      double precision charge_activ
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot
     *, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z
     *, sd_v, sd_e, v2, e2
     *, phi, ex, ey, ez
     *, dx, dy, dz, d2, d
     *, normx, normy, normz
      dimension nr_mm(nr_frame)
     *, charge_activ (nat), grom2cpmd(nat)
     *, tr_qm_x(nr_frame, nr_qm)
     *, tr_qm_y(nr_frame, nr_qm)
     *, tr_qm_z(nr_frame, nr_qm)
     *, tr_mm_x(nr_frame, max_nr_mm)
     *, tr_mm_y(nr_frame, max_nr_mm)
     *, tr_mm_z(nr_frame, max_nr_mm)
c
     *, tr_ele_pot(nr_frame, max_nr_mm)
     *, tr_ele_fd_x(nr_frame, max_nr_mm)
     *, tr_ele_fd_y(nr_frame, max_nr_mm)
     *, tr_ele_fd_z(nr_frame, max_nr_mm)
c
      write(*,*)''
      write(*,*)'-------------------------------------'
      write(*,*)'calculating relative standard deviation'
      write(*,*)'see J Chem Theory Comput 2007, 3, 628'
      write(*,*)'for the definition'
      write(*,*)''
c
      sd_v = 0.0d0
      sd_e = 0.0d0
      v2 = 0.0d0
      e2 = 0.0d0
c
      do i=1,nr_frame                !  frame loop
         do j=1, nr_mm(i)            !  MM loop
c
          phi = 0.0d0
          ex = 0.0d0
          ey = 0.0d0
          ez = 0.0d0
          k_c = 0
          do k=1, nat              ! nat loop
            if(grom2cpmd(k) .le. nr_qm)then
               k_c = k_c + 1
               dx=tr_mm_x(i,j) - tr_qm_x(i,k_c)
               dy=tr_mm_y(i,j) - tr_qm_y(i,k_c)
               dz=tr_mm_z(i,j) - tr_qm_z(i,k_c)
               d2 = (dx*dx) + (dy*dy) + (dz*dz)
               d = dsqrt( d2 )
               phi = phi + charge_activ( k ) / d
c
               normx = dx / d
               normy = dy / d
               normz = dz / d
c
               ex = ex + (charge_activ( k )*normx/d2)
               ey = ey + (charge_activ( k )*normy/d2)
               ez = ez + (charge_activ( k )*normz/d2)
            endif
          enddo                      ! nat loop
c
          sd_v = sd_v + (phi -  tr_ele_pot(i, j) )**2.0d0
          v2 = v2 + tr_ele_pot(i, j)**2.0d0
c
          sd_e = sd_e  + ( ex - tr_ele_fd_x(i, j) )**2.0d0
     *   + ( ey - tr_ele_fd_y(i, j) )**2.0d0
     *   + ( ez - tr_ele_fd_z(i, j) )**2.0d0
c
          e2 = e2 + tr_ele_fd_x(i, j)**2.0d0
     *   + tr_ele_fd_y(i, j)**2.0d0
     *   + tr_ele_fd_z(i, j)**2.0d0
c
         enddo                       !  MM loop
      enddo                          !  frame loop
c
      sd_v = dsqrt(sd_v/v2)
      sd_e = dsqrt( sd_e/e2)
c
      write(*,"(a,f15.8)")'SDv = ',sd_v
      write(*,"(a,f15.8)")'SDe = ',sd_e
c
      write(*,*)'calculating relative standard deviation done!'
      write(*,*)' '
c
      return
      end
