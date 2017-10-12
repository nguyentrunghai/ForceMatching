c#################################
c generate target matrix for 
c charge fitting
c this is done for all of the
c QM/MM snapshots
c###############################
      subroutine matr_gen(nr_frame,nr_qm, max_nr_mm
     *,nr_char_tp
     *, nr_li_pot, nr_li_fd, nr_li_hirsh
     *, nr_mm, charg_tp
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z, tr_chj
     *, w_v, w_e, w_h, w_q
     *, tot_qm_charg
     *, tot_matr_size, tot_matr, tot_targ )
c
      implicit none
      integer nr_frame,nr_qm, max_nr_mm,nr_char_tp
     *, nr_li_pot, nr_li_fd, nr_li_hirsh
     *, nr_mm, charg_tp, tot_matr_size
     *, itmp, i, j
      double precision  tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z, tr_chj
     *, tot_matr, tot_targ
     *, w_v, w_e, w_h, w_q, tot_qm_charg
     *, dtmp
       dimension nr_mm(nr_frame), charg_tp(nr_qm)
     *, tr_qm_x(nr_frame,nr_qm)
     *, tr_qm_y(nr_frame,nr_qm)
     *, tr_qm_z(nr_frame,nr_qm)
     *, tr_mm_x(nr_frame, max_nr_mm)
     *, tr_mm_y(nr_frame, max_nr_mm)
     *, tr_mm_z(nr_frame, max_nr_mm)
     *, tr_ele_pot(nr_frame, max_nr_mm)
     *, tr_ele_fd_x(nr_frame, max_nr_mm)
     *, tr_ele_fd_y(nr_frame, max_nr_mm)
     *, tr_ele_fd_z(nr_frame, max_nr_mm)
     *, tr_chj(nr_frame, nr_qm)
     *, tot_matr( tot_matr_size, nr_char_tp)
     *, tot_targ (tot_matr_size)
c
      double precision, dimension (:,:), allocatable :: pot_mat
     *, field_mat,  hirs_mat
      allocate ( pot_mat(nr_li_pot, nr_char_tp + 1) )
      allocate ( field_mat(nr_li_fd, nr_char_tp + 1) )
      allocate ( hirs_mat(nr_li_hirsh, nr_char_tp + 1) )
c
      itmp = nr_li_pot + nr_li_fd + nr_li_hirsh + 1
      if(itmp .ne. tot_matr_size)then
      stop'matr_gen: problem with the total size of the matrix'
      endif
c
      do i=1,tot_matr_size
       do j=1,nr_char_tp
        tot_matr( i, j) = 0.0d0
       enddo
        tot_targ (i) = 0.0d0
      enddo
c
      call pot_matr_gen(nr_frame, nr_qm, max_nr_mm
     *, nr_char_tp, nr_li_pot
     *, nr_mm, charg_tp
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot
     *, pot_mat)
c
      call field_matr_gen(nr_frame, nr_qm, max_nr_mm
     *, nr_char_tp, nr_li_fd
     *, nr_mm, charg_tp
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z 
     *, field_mat) 
c
      call hirsh_matr_gen(nr_frame, nr_qm
     *, nr_char_tp, nr_li_hirsh
     *, charg_tp, tr_chj
     *, hirs_mat)
c pot
      do i=1,nr_li_pot 
         do j=1,nr_char_tp
           tot_matr(i, j) = pot_mat (i,j) * w_v
         enddo
         tot_targ(i) = pot_mat (i, nr_char_tp + 1) * w_v
      enddo
c field
      do i=1,nr_li_fd
         do j=1,nr_char_tp
           tot_matr(i + nr_li_pot, j) = field_mat(i, j) * w_e
         enddo
         tot_targ(i + nr_li_pot) = field_mat(i, nr_char_tp + 1) * w_e
      enddo
c  hirshfeld
      do i=1,nr_li_hirsh
         do j=1,nr_char_tp
           tot_matr(i + nr_li_pot + nr_li_fd, j) = hirs_mat(i, j) * w_h
         enddo
         tot_targ(i + nr_li_pot + nr_li_fd) = 
     *              hirs_mat(i, nr_char_tp + 1) * w_h
      enddo
c  total charge
      do j=1,nr_qm
        if( charg_tp(j) .ne. 0 )then
           tot_matr(tot_matr_size, charg_tp(j) ) = 
     *              tot_matr(tot_matr_size, charg_tp(j) ) + w_q
        endif
      enddo
      tot_targ(tot_matr_size) = tot_qm_charg * w_q 
c######
c eleminate zero rows
c
      do i=1, tot_matr_size
        dtmp = 0.0d0
c
          do j=1,nr_char_tp
           dtmp = dtmp + dabs( tot_matr(i, j) )
          enddo
c
          if(dtmp .eq. 0.0d0)then
            tot_targ( i ) = 0.0d0
        endif
c
      enddo
c
       write(*,*)''
       write(*,*)'---------------------------'
       write(*,*)'colecting all the matrices'
       write(*,*)''
c
      return
      end
c###################################
c    generate target matrix for 
c    electrostatic POTENTIAL
c##################################
      subroutine pot_matr_gen(nr_frame, nr_qm, max_nr_mm
     *, nr_char_tp, nr_li_pot
     *, nr_mm, charg_tp
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot
     *, pot_mat)         !  only pot_mat is the output
c
      implicit none
      integer nr_frame, nr_qm, max_nr_mm
     *,nr_char_tp, nr_li_pot, nr_mm, charg_tp
     *,  li_count ,i, j, k, itmp
      double precision tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot
     *, pot_mat, dx, dy, dz, d
      dimension  nr_mm(nr_frame), charg_tp(nr_qm)
     *, tr_qm_x(nr_frame,nr_qm)
     *, tr_qm_y(nr_frame,nr_qm)
     *, tr_qm_z(nr_frame,nr_qm)
     *, tr_mm_x(nr_frame, max_nr_mm)
     *, tr_mm_y(nr_frame, max_nr_mm)
     *, tr_mm_z(nr_frame, max_nr_mm)
     *, tr_ele_pot(nr_frame, max_nr_mm)
     *, pot_mat(nr_li_pot, nr_char_tp + 1)
c
       write(*,*)''
       write(*,*)'----------------------------------'
       write(*,*)'generating matrix for potential fit'
c  initilize matrix
      do i=1,nr_li_pot
      do j=1,nr_char_tp + 1
       pot_mat(i,j)=0.0d0
      enddo
      enddo
c
      li_count = 0
c
      do i=1,nr_frame              !   frame loop
         do j=1,nr_mm(i)           !   MM loop
         li_count = li_count + 1
           do k=1,nr_qm            !   QM loop
c
             dx=tr_mm_x(i,j) - tr_qm_x(i,k)
             dy=tr_mm_y(i,j) - tr_qm_y(i,k)
             dz=tr_mm_z(i,j) - tr_qm_z(i,k)
             d = dsqrt( (dx*dx) + (dy*dy) + (dz*dz) )
             if( d .eq. 0.0d0 )then
              stop'pot_matr_gen: d .eq. 0.0d0 '
             endif
             itmp = charg_tp (k)
             if(itmp.ne.0)then
             pot_mat( li_count, itmp) =  pot_mat( li_count, itmp) 
     *                                  + ( 1.0d0 / d )
             endif
           enddo                   !   QM loop
c
           pot_mat( li_count, nr_char_tp + 1) = tr_ele_pot(i, j)
c
         enddo                     !   MM loop
      enddo                        !   frame lopp
c
      if( li_count .ne. nr_li_pot) then
      stop'pot_matr_gen: wrong value of nr_li_pot'
      endif
c
       write(*,*)'generating matrix for potential fit done!'
       write(*,*)''
c
      return
      end
c###################################
c    generate target matrix for 
c    electrostatic FIELD
c##################################
      subroutine field_matr_gen(nr_frame, nr_qm, max_nr_mm
     *, nr_char_tp, nr_li_fd
     *, nr_mm, charg_tp
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z 
     *, field_mat)         !  only pot_mat is the output
c
      implicit none
      integer nr_frame, nr_qm, max_nr_mm
     *,nr_char_tp, nr_li_fd, nr_mm, charg_tp
     *,  li_count ,i, j, k, itmp
      double precision tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z
     *, field_mat, dx, dy, dz, d, d2
     *, normx, normy, normz
      dimension  nr_mm(nr_frame), charg_tp(nr_qm)
     *, tr_qm_x(nr_frame,nr_qm)
     *, tr_qm_y(nr_frame,nr_qm)
     *, tr_qm_z(nr_frame,nr_qm)
     *, tr_mm_x(nr_frame, max_nr_mm)
     *, tr_mm_y(nr_frame, max_nr_mm)
     *, tr_mm_z(nr_frame, max_nr_mm)
     *, tr_ele_fd_x(nr_frame, max_nr_mm)
     *, tr_ele_fd_y(nr_frame, max_nr_mm)
     *, tr_ele_fd_z(nr_frame, max_nr_mm)
     *, field_mat(nr_li_fd, nr_char_tp + 1)
c
       write(*,*)''
       write(*,*)'-------------------------------'
       write(*,*)'generating matrix for field fit'
c
      do i=1,nr_li_fd
      do j=1,nr_char_tp + 1
       field_mat(i,j)=0.0d0
      enddo
      enddo
c
      li_count = 0
c
      do i=1,nr_frame              !   frame loop
         do j=1,nr_mm(i)           !   MM loop
c f_x
         li_count = li_count + 1
           do k=1,nr_qm            !   QM loop
c
             dx=tr_mm_x(i,j) - tr_qm_x(i,k)
             dy=tr_mm_y(i,j) - tr_qm_y(i,k)
             dz=tr_mm_z(i,j) - tr_qm_z(i,k)
             d2 = (dx*dx) + (dy*dy) + (dz*dz)
             d = dsqrt(d2)
             if( d .eq. 0.0d0 )then
               stop'field_matr_gen: d .eq. 0.0d0'
             endif
             normx = dx / d                       !  must be x
             itmp = charg_tp (k)
             if(itmp.ne.0)then
             field_mat( li_count, itmp) =  field_mat( li_count, itmp) 
     *                                  + (normx / d2)            !  must be x
             endif
           enddo                   !   QM loop
c
           field_mat( li_count, nr_char_tp + 1) = tr_ele_fd_x (i, j)     !  msut be x
c end f_x   #################
c
c f_y
         li_count = li_count + 1
           do k=1,nr_qm            !   QM loop
c
             dx=tr_mm_x(i,j) - tr_qm_x(i,k)
             dy=tr_mm_y(i,j) - tr_qm_y(i,k)
             dz=tr_mm_z(i,j) - tr_qm_z(i,k)
             d2 = (dx*dx) + (dy*dy) + (dz*dz)
             d = dsqrt(d2)
             normy = dy / d                  !  must be y
             itmp = charg_tp (k)
             if(itmp.ne.0)then
             field_mat( li_count, itmp) =  field_mat( li_count, itmp) 
     *                                  + (normy / d2)       !  msut be y
             endif
           enddo                   !   QM loop
c
           field_mat( li_count, nr_char_tp + 1) = tr_ele_fd_y (i, j)      ! mustbe y
c end f_y   #################
c
c f_z
         li_count = li_count + 1
           do k=1,nr_qm            !   QM loop
c
             dx=tr_mm_x(i,j) - tr_qm_x(i,k)
             dy=tr_mm_y(i,j) - tr_qm_y(i,k)
             dz=tr_mm_z(i,j) - tr_qm_z(i,k)
             d2 = (dx*dx) + (dy*dy) + (dz*dz)
             d = dsqrt(d2)
             normz = dz / d               !  must be z
             itmp = charg_tp (k)
             if(itmp.ne.0)then
             field_mat( li_count, itmp) =  field_mat( li_count, itmp) 
     *                                  + (normz / d2)    ! must be z
             endif
           enddo                   !   QM loop
c
           field_mat( li_count, nr_char_tp + 1) = tr_ele_fd_z (i, j)   !  must be z
c end f_z   #################
c
         enddo                     !   MM loop
      enddo                        !   frame lopp
c
      if( li_count .ne. nr_li_fd) then
      stop'field_matr_gen: wrong value of nr_li_fd'
      endif
c
      write(*,*)'generating matrix for field fit done!'
      write(*,*)''
c
      return
      end
c#########################################
c
c    generate target matrix for Hirshfeld charges
c
c###################################
      subroutine  hirsh_matr_gen(nr_frame, nr_qm
     *, nr_char_tp, nr_li_hirsh
     *, charg_tp, tr_chj
     *, hirs_mat)
c
      implicit none
      integer nr_frame, nr_qm
     *, nr_char_tp, nr_li_hirsh
     *, charg_tp, i, j, li_count, itmp
      double precision tr_chj
     *, hirs_mat
      dimension charg_tp(nr_qm), tr_chj(nr_frame, nr_qm)
     *, hirs_mat(nr_li_hirsh, nr_char_tp + 1)
c
       write(*,*)''
       write(*,*)'----------------------------------'
       write(*,*)'generating matrix for Hirshfeld fit'
c initilize
      do i=1,nr_li_hirsh
         do j=1, nr_char_tp + 1
          hirs_mat(i, j) = 0.0d0
         enddo
      enddo
c
      li_count = 0
      do i=1, nr_frame         !  frame loop
         do j=1, nr_qm         !  QM loop
           li_count = li_count + 1
           itmp = charg_tp (j)
           if(itmp.ne.0)then
             hirs_mat( li_count , itmp ) = 1.0d0
             hirs_mat( li_count , nr_char_tp + 1 ) = tr_chj( i, j)
           endif
         enddo                 !  QM loop
      enddo                    !  frame loop
c
      if( li_count .ne. nr_li_hirsh) then
       stop'hirsh_matr_gen: wrong value of nr_li_hirsh'
      endif
c
      write(*,*)'generating matrix for Hirshfeld fit done!'
      write(*,*)''
c
      return
      end
c####################################
c   generate the size (no. of lines)
c  of the target matrices
c##################################
      subroutine matr_size_gen(nr_frame, nr_qm, nr_mm 
     *, nr_li_pot, nr_li_fd, nr_li_hirsh 
     *, tot_matr_size )
c
      implicit none
      integer nr_frame, nr_qm, nr_mm 
     *, nr_li_pot, nr_li_fd, nr_li_hirsh
     *, tot_matr_size
     *,i, itmp
      dimension nr_mm(nr_frame)
c
      write(*,*)''
      write(*,*)'------------------------'
      write(*,*)'generating matrix size'
c
      itmp=0
      do i=1,nr_frame
      itmp=itmp+nr_mm(i)
      enddo
c
      nr_li_pot=itmp
      nr_li_fd=itmp*3
      nr_li_hirsh=nr_frame*nr_qm
c
      tot_matr_size = nr_li_pot + nr_li_fd + nr_li_hirsh +1
      if(tot_matr_size .le. 0)then
      stop'matr_size_gen: tot_matr_size .le. 0'
      endif
c
      write(*,*)'matrix sizes'
      write(*,"(4I8)")nr_li_pot, nr_li_fd, nr_li_hirsh
     *, tot_matr_size
      write(*,*)'generating matrix size done!'
      write(*,*)''
c
      return
      end
