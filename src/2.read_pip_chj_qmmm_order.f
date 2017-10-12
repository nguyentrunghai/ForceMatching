c#######################################
c  read FM_REF_PIP file
c  this file contains 
c  the electrostatic potential 
c  for nr_frame QM/MM snapshots
c#########################################
      subroutine read_pip( nr_frame, nr_qm, max_nr_mm
     *, nr_mm, frame_ind
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z )
c
      implicit none
c
      integer nr_frame, nr_qm, nr_mm
     *, max_nr_mm
     *, nr_mm_tmp, frame_ind
     *, i, j, iq_count, im_count
c
      character(2) qmmm
c
      double precision tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot
     *, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z
     *, x,y,z,clas_charg, ele_pot, dtmp
     *, fdx, fdy, fdz
c
      dimension nr_mm(nr_frame),frame_ind (nr_frame)
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
c
      write(*,*)' '
      write(*,*)'----------------------'
      write(*,*)'reading FM_REF_PIP'
      open(1,file='FM_REF_PIP',status='old')
      rewind(1)
c
      do i=1, nr_frame     !   frame loop
         read(1,*)nr_mm_tmp, frame_ind (i)
         nr_mm(i)=nr_mm_tmp - nr_qm
         if(max_nr_mm .lt. nr_mm_tmp) then
         stop'read_pip: max_nr_mm .lt. nr_mm_tmp'
         endif
c
          if(nr_mm(i) .le. 0) then
          stop'read_pip: nr_mm(i) .le. 0'
          endif
c
         iq_count=0
         im_count=0
         do j=1,nr_mm_tmp      !  atom loop
           read(1,*)x,y,z,qmmm,clas_charg, ele_pot, dtmp
     *     , dtmp, dtmp, fdx, fdy, fdz
c
            if(qmmm .eq. 'QM') then 
                 iq_count=iq_count+1
                 if(iq_count .gt. nr_qm) then
                   stop'read_pip: iq_count gt nr_qm'
                 endif
c
                 tr_qm_x(i,iq_count)=x
                 tr_qm_y(i,iq_count)=y
                 tr_qm_z(i,iq_count)=z
            else
                if(qmmm .eq. 'MM') then 
                    im_count = im_count+1
                    if(im_count .gt. nr_mm(i) ) then
                    stop'read_pip: im_count gt nr_mm(i)'
                    endif
c
                    tr_mm_x(i, im_count) = x
                    tr_mm_y(i, im_count) = y
                    tr_mm_z(i, im_count) = z
                    tr_ele_pot(i, im_count) = ele_pot
                    tr_ele_fd_x(i, im_count) = fdx
                    tr_ele_fd_y(i, im_count) = fdy
                    tr_ele_fd_z(i, im_count) = fdz
c
                else
                    stop'read_pip: neither MM nor QM found'
                endif
            endif
c
         enddo     !  atom loop
      enddo        !  frame loop
c
      close(1)
      write(*,*)'reading FM_REF_PIP done !'
      write(*,*)''
c
      return
      end
c#############################
c   reading FM_REF_CHJ
c  this file contains 
c  Hirshfeld charges from CPMD
c#############################
      subroutine read_chj(nat, nr_qm, nr_frame, frame_ind
     *,grom2cpmd
     *, tr_chj)
c
      implicit none
c
      integer nat, nr_qm, nr_frame, frame_ind, grom2cpmd
     *, i, j, itmp, icount
      double precision tr_chj_tmp,tr_chj
      dimension frame_ind(nr_frame)
     *, tr_chj_tmp(nr_frame, nr_qm), tr_chj(nr_frame, nr_qm)
     *, grom2cpmd(nat)
c
      write(*,*)''
      write(*,*)'-----------------------'
      write(*,*)'reading FM_REF_CHJ'
      open(1,file='FM_REF_CHJ',status='old')
      rewind(1)
      do i=1,nr_frame    !  frame loop
        read(1,*)itmp
        if( itmp.ne.frame_ind(i) )then
           write(*,*)'frame ',i
           stop'read_chj: frame index does not match!'
        endif
c
         read(1,*)( tr_chj_tmp(i,j), j=1,nr_qm )
c
      enddo     ! frame loop 
      close(1)
c
      write(*,*) 'reordering the chj chagres from cpmd to gromos'
c
      do i=1,nr_frame     ! frame loop
          icount=0
          do j=1,nat      ! gromos ind loop
             itmp =  grom2cpmd (j)
             if(itmp.le.nr_qm)then
                 icount=icount+1
                tr_chj(i,icount)= tr_chj_tmp(i,itmp)
             endif
          enddo           ! gromos ind loop
          if(icount.ne.nr_qm)then
            write(*,*)'frame',i
            stop'read_chj: icount.ne.nr_qm '
          endif
      enddo            !  frame loop
c
      write(*,*)'reading FM_REF_CHJ done!'
      write(*,*)''
c
      return
      end
c#################################
c read QMMM_ORDER
c this file contain conversion rule between CPMD and gromos indexing
c genrate cpmd2grom (cpmd index) --> gromos index
c and     grom2cpmd (gromos index ) --> cpmd index
c#################################
      subroutine read_qmmm_order(nat, nr_qm
     *, cpmd2grom, grom2cpmd)
c
      implicit none 
c
      integer nat, nr_qm
     *, cpmd2grom, grom2cpmd
     *, itmp, i, j
     *, nat_tmp, nr_qm_tmp
c
      character(1)ctmp
c
      dimension cpmd2grom (nat), grom2cpmd(nat)
c
      write(*,*)''
      write(*,*)'------------------------------'
      write(*,*)'reading the rest of QMMM_ORDER'
      open(1,file='QMMM_ORDER',status='old')
      rewind(1)
c
      read(1,*)ctmp
      read(1,*)ctmp,ctmp,ctmp,nat_tmp,ctmp,ctmp,nr_qm_tmp
      read(1,*)ctmp
      read(1,*)ctmp
      if( (nat_tmp.ne.nat) .or. (nr_qm_tmp.ne.nr_qm) )then
         stop'read_qmmm_order: number of atoms is not right'
      endif
c
      do i=1,nat
        read(1,*)cpmd2grom (i),itmp,itmp,itmp,ctmp
      enddo
      close(1)
c
      do i=1,nat
        do j=1,nat
           if(cpmd2grom (j) .eq. i) then
               grom2cpmd(i) = j
           endif
        enddo
      enddo
      write(*,*)'reading the rest of QMMM_ORDER done!'
      write(*,*)''
c
      return
      end
c##################################
c  read top part of 
c QMMM_ORDER to get nat and nr_qm
c##################################
      subroutine read_nat_qm(nat, nr_qm)
c
      implicit none
      integer nat, nr_qm
     *, nat_tmp, nr_qm_tmp
      character(1)ctmp
c
      write(*,*)' '
      write(*,*)'-------------------------'
      write(*,*)'reading top part of QMMM_ORDER'
c
      open(1,file='QMMM_ORDER',status='old')
      rewind(1)
c
      read(1,*)ctmp
      read(1,*)ctmp,ctmp,ctmp,nat_tmp,ctmp,ctmp,nr_qm_tmp
      read(1,*)ctmp
      read(1,*)ctmp
      close(1)
c
      nat = nat_tmp
      nr_qm = nr_qm_tmp
c
      if(nat .eq. 0)then
      stop'read_nat_qm: .eq. 0'
      endif
c
      if(nr_qm .eq. 0)then
      stop'read_nat_qm: nr_qm .eq. 0'
      endif
c
      write(*,"(a,I8)")'total number of atoms ',nat
      write(*,"(a,I8)")'number of QM atoms ',nr_qm
      write(*,*)'reading top part of QMMM_ORDER done!'
      write(*,*)''
c
      return
      end
