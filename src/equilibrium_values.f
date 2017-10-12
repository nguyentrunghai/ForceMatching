c###############################
c  to calculate equlibrium bond lengths
c  in Angstrom
c###############################
      subroutine equ_bond_lengths(nat, nr_frame, max_nr_qmb
     *, nr_qmb, ib_qm, jb_qm
     *, trx, try, trz
     *, req_cal )
c
      implicit none
      integer nat, nr_frame, max_nr_qmb
     *, nr_qmb, ib_qm, jb_qm
     *, i, j, itmp, jtmp
      double precision trx, try, trz
     *, req_cal
     *, bond_dist, dtmp
      dimension ib_qm(max_nr_qmb), jb_qm(max_nr_qmb)
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
     *,  req_cal( max_nr_qmb )
c
      write(*,*)''
      write(*,*)'----------------------------------'
      write(*,*)'calculating averaged bond distance'
c
      do i=1, nr_qmb                         !  qm bond loop
          req_cal(i) = 0.0d0
          itmp = ib_qm(i)
          jtmp = jb_qm(i)
c
          if( (itmp .le. 0) .or. (itmp .gt. nat) )then
        stop'equ_bond_lengths: (itmp .le. 0) .or. (itmp .gt. nat)'
          endif
c
          if( (jtmp .le. 0) .or. (jtmp .gt. nat) )then
        stop'equ_bond_lengths: (jtmp .le. 0) .or. (jtmp .gt. nat)'
          endif
c
          if(nr_frame .eq. 0)then
             stop'equ_bond_lengths: divide by 0'
          endif
c
          do j=1, nr_frame                  !  frame loop
               dtmp = bond_dist( trx(j,itmp), try(j,itmp), trz(j,itmp)
     *,                      trx(j,jtmp), try(j,jtmp), trz(j,jtmp) )
c
               req_cal(i) = req_cal(i) + dtmp
c
          enddo                             !  frame loop
c
               req_cal(i) = req_cal(i)/dble(nr_frame)
c
      enddo
c
      write(*,*)'calculating averaged bond distance done'
      write(*,*)''
c
      return
      end
c###############################
c  to calculate equlibrium angles
c  in rad
c###############################
      subroutine equ_bond_angles(nat, nr_frame, max_nr_qma
     *, nr_qma, ia_qm, ja_qm, ka_qm
     *, trx, try, trz
     *, teq_cal )
c
      implicit none
      integer nat, nr_frame, max_nr_qma
     *, nr_qma, ia_qm, ja_qm, ka_qm
     *, i, j, itmp, jtmp, ktmp
     *, j_test
     *, nr_pi_angl, pi_angl
      character(30)filename
      double precision trx, try, trz
     *, teq_cal, bond_ang, dtmp
     *, pi
      dimension ia_qm(max_nr_qma), ja_qm(max_nr_qma), ka_qm(max_nr_qma)
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
     *, teq_cal(max_nr_qma)
     *, pi_angl( max_nr_qma )
c
      write(*,*)''
      write(*,*)'----------------------------------'
      write(*,*)'calculating averaged bond angles'
c
      pi=dacos(-1.0d0)
c
      do i=1, nr_qma                       !  qm angle loop
          teq_cal(i) = 0.0d0
          itmp = ia_qm( i )
          jtmp = ja_qm( i )
          ktmp = ka_qm( i )
c
          if( (itmp .le. 0) .or. (itmp .gt. nat) )then
        stop'equ_bond_angles: (itmp .le. 0) .or. (itmp .gt. nat)'
          endif
c
          if( (jtmp .le. 0) .or. (jtmp .gt. nat) )then
        stop'equ_bond_angles: (jtmp .le. 0) .or. (jtmp .gt. nat)'
          endif
c
          if( (ktmp .le. 0) .or. (ktmp .gt. nat) )then
        stop'equ_bond_angles: (ktmp .le. 0) .or. (ktmp .gt. nat)'
          endif
c
          if(nr_frame .eq. 0)then
          stop'equ_bond_angles: divide by 0'
          endif
c
          do j=1, nr_frame                 ! frame loop
              dtmp = bond_ang(trx(j,itmp), try(j,itmp), trz(j,itmp)
     *,                       trx(j,jtmp), try(j,jtmp), trz(j,jtmp) 
     *,                      trx(j,ktmp), try(j,ktmp), trz(j,ktmp) )
c
              teq_cal(i) = teq_cal(i) + dtmp
c
          enddo                            ! frame loop
c
           teq_cal(i) = teq_cal(i)/dble(nr_frame)
c
      enddo                                !  qm angle loop
c
c       write(*,*)''
c       write(*,*)'do you have linear 
c     * angles (at1-at2-at3 = 180) in your system'
c       write(*,*)'yes: enter 1'
c       write(*,*)'no: enter 0'
c       read(*,*)j_test
c       if(j_test .eq. 1)then
c          write(*,*)'please provide a file, in which'
c          write(*,*)''
c          write(*,*)'Number of 180-degree angles'
c          write(*,*)'angle 1'
c          write(*,*)'angle 2'
c          write(*,*)'.'
c          write(*,*)'.'
c          write(*,*)'.'
c          write(*,*)''
c          write(*,*)'now, give me the name of that file:'
c          write(*,*)'( 30 characters or less )'
c          read(*,*)filename
c          open(1,file=filename,status='old')
c          rewind(1)
c          read(1,*)nr_pi_angl
c          write(*,"(i8)")nr_pi_angl
c          if( nr_pi_angl .gt. max_nr_qma )then
c           stop'equ_bond_angles: nr_pi_angl .gt. max_nr_qma'
c          endif
c
c          if(nr_pi_angl .gt. 0)then
c            do i=1,nr_pi_angl
c               read(1,*)pi_angl(i)
c               write(*,"(i8)")pi_angl(i)
c            enddo
c          endif
c
c          close(1)
c
c          do i=1, nr_pi_angl
c              do j=1, nr_qma
c
c                if(pi_angl(i) .eq. j)then
c                    teq_cal(j) = pi
c                endif
c
c              enddo
c          enddo
c
c       else
c           if(j_test .ne. 0)then
c              stop'equ_bond_angles: dont know what to do'
c           endif
c
c       endif
c
c
       write(*,*)'calculating averaged bond angles done'
       write(*,*)''
c
      return
      end
c###############################
c  evolution of dihedral angles
c  this is usefule for deciding
c  the number of minima and phase
c###############################
      subroutine evolu_dihedal(nat, nr_frame, max_nr_qmd
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm
     *, pro_imp
     *, trx, try, trz
     *, dih_cal)
c
      implicit none
      integer nat, nr_frame, max_nr_qmd
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm
     *, i,j, itmp, jtmp, ktmp, ltmp
      character(3)pro_imp
      double precision trx, try, trz
     *, dih_cal, dihed_ang, impro_ang, dtmp
      dimension id_qm(max_nr_qmd), jd_qm(max_nr_qmd)
     *, kd_qm(max_nr_qmd), ld_qm(max_nr_qmd)
     *, pro_imp(max_nr_qmd)
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
     *, dih_cal( nr_frame, max_nr_qmd )
c
      write(*,*)''
      write(*,*)'----------------------------------'
      write(*,*)'calculating dihedral angles'
c
      do i=1, nr_qmd                    ! qm dihedral loop
          itmp = id_qm( i )
          jtmp = jd_qm( i )
          ktmp = kd_qm( i )
          ltmp = ld_qm( i )
c
          if( (itmp .le. 0) .or. (itmp .gt. nat) )then
        stop'evolu_dihedal: (itmp .le. 0) .or. (itmp .gt. nat)'
          endif
c
          if( (jtmp .le. 0) .or. (jtmp .gt. nat) )then
        stop'evolu_dihedal: (jtmp .le. 0) .or. (jtmp .gt. nat)'
          endif
c
          if( (ktmp .le. 0) .or. (ktmp .gt. nat) )then
        stop'evolu_dihedal: (ktmp .le. 0) .or. (ktmp .gt. nat)'
          endif
c
          if( (ltmp .le. 0) .or. (ltmp .gt. nat) )then
        stop'evolu_dihedal: (ltmp .le. 0) .or. (ltmp .gt. nat)'
          endif
c
          if(nr_frame .eq. 0)then
          stop'equ_bond_angles: divide by 0'
          endif
c
          do j=1, nr_frame              ! frame loop
             if(pro_imp(i) .eq. 'PRO' )then
               dtmp = dihed_ang(trx(j,itmp), try(j,itmp), trz(j,itmp)
     *,                        trx(j,jtmp), try(j,jtmp), trz(j,jtmp) 
     *,                       trx(j,ktmp), try(j,ktmp), trz(j,ktmp) 
     *,                       trx(j,ltmp), try(j,ltmp), trz(j,ltmp)  )
             else
                 if(pro_imp(i) .eq. 'IMP' )then
                 dtmp = impro_ang(trx(j,ktmp), try(j,ktmp), trz(j,ktmp)
     *,                        trx(j,jtmp), try(j,jtmp), trz(j,jtmp) 
     *,                       trx(j,itmp), try(j,itmp), trz(j,itmp) 
     *,                       trx(j,ltmp), try(j,ltmp), trz(j,ltmp)  )
                 else 
                     stop'evolu_dihedal: wrong PRO IMP'
                 endif
             endif
c
          dih_cal( j, i ) = dtmp
c
          enddo                         ! frame loop
      enddo                             ! qm dihedral loop
c
      write(*,*)'calculating dihedral angles done'
      write(*,*)''
c
      return
      end
c####################################
c  print dihedrals to separate files
c#####################################
      subroutine writ_dihe( nr_frame, max_nr_qmd
     *, nr_qmd
     *, nr_qmd_typ, qmd_typ, pro_imp
     *, dih_cal )
c
      implicit none
      integer nr_frame, max_nr_qmd
     *, nr_qmd, nr_qmd_typ, qmd_typ
     *, i, j, k, icount, typ_count
     *, j_test
     *, nr_arr, itmp
      integer, dimension (:,:), allocatable :: dih_ind
      character(3)pro_imp
      double precision dih_cal
     *, angle_arr, histo, dtmp
      double precision, dimension (:,:,:), allocatable :: dih
      parameter(nr_arr = 181)
      dimension qmd_typ(max_nr_qmd)
     *, pro_imp(max_nr_qmd)
     *, dih_cal( nr_frame, max_nr_qmd )
     *, typ_count( max_nr_qmd )
     *, angle_arr(nr_arr), histo(nr_arr)
c
      allocate( dih_ind( max_nr_qmd, max_nr_qmd ) )
      allocate( dih ( nr_frame, max_nr_qmd, nr_qmd_typ+1 ) )
c
      write(*,*)''
      write(*,*)'---------------------------------'
      write(*,*)'writing dihedral angles'
c      write(*,*)'this info may help you determine phase and 
c     * number of bariers'
c      write(*,*)'please make directory named DIHE_EVOL'
c      write(*,*)'then press enter'
c      read(*,*)
c
      if(nr_qmd_typ .gt. 0) then
      do i=1, nr_qmd_typ           ! type loop
c
          icount = 0
          do j=1, nr_qmd           ! qmd loop
c
            if( qmd_typ(j) .eq. i) then
c
              icount = icount + 1
c
              do k=1,nr_frame
               dih ( k, icount, i ) = dih_cal( k, j )
              enddo
c
              dih_ind( i, icount ) = j
c
            endif
c
          enddo                    ! qmd loop
          typ_count(i) = icount
c
      enddo                        ! type loop
c
      do i=1, nr_qmd_typ
        call indexe_files(i)
c
         write(1,"(a,2i8)")'# type ',i, typ_count(i)
         write(1,"(a,30i15)")'# order,'
     *,             ( dih_ind( i, j ),j=1,typ_count(i) )
         do k=1,nr_frame
           write(1,"(i8,30f15.7)")k
     *, ( dih ( k, icount, i )*180.0d0/3.14159265d0
     *, icount = 1, typ_count(i) )
         enddo
c
        close(1)
      enddo
c
c  histogram
c
      call databin( -180.0d0, 180.0d0
     *, nr_arr, angle_arr )
c
      do i=1,nr_qmd_typ           ! type loop
c
          do j=1,nr_arr
            histo( j ) = 0.0d0
          enddo
c
          icount = 0
          do k=1,nr_frame         ! frame loop
              do j=1,typ_count(i)
                 icount = icount + 1
                 dtmp = dih ( k, j, i )*180.0d0/3.14159265d0
                 call valueloca( angle_arr, nr_arr, dtmp, itmp )
c
                 histo( itmp ) = histo( itmp ) + 1.0d0
c
              enddo
          enddo                   ! frame loop
c
          if(icount .gt. 0)then
          do j=1,nr_arr
            histo( j ) = histo( j )/dble( icount )
          enddo
          endif
c
      call indexe_files1(i)
      write(1,"(a,i8)")'# type ',i
        do j=1,nr_arr
          write(1,"(2f15.7)")angle_arr(j),histo(j)
        enddo
      close(1)
c
      enddo                       ! type loop
c
c      write(*,*)'please check the xxx.dih files'
c      write(*,*)''
c
c      write(*,*)'press 1 to continue or 2 to stop'
c      read(*,*)j_test
c      if(j_test .eq. 2)then
c      stop
c      endif
c
      endif
      write(*,*)'writing dihedral angles done!'
      write(*,*)''
c
      return
      end
c###################################
c  to open an indexed file
c##################################
      subroutine indexe_files(fr_in)
c
      implicit none
c
      integer fr_in
c
      character(20)flname,forma
c
      if (fr_in.lt.10)then
           forma="(a,i1,a)"
      else
          if (fr_in.lt.100)then
          forma="(a,i2,a)"
          else
          forma="(a,i3,a)"
          endif
      endif
c
      write(flname,forma)"evol.dih",fr_in
c
      write(*,"(a,a)")'open file ',flname
c
      OPEN(1,file=flname,status='unknown')
      rewind 1
c
      return
      end
c###################################
c  to open an indexed file
c##################################
      subroutine indexe_files1(fr_in)
c
      implicit none
c
      integer fr_in
c
      character(20)flname,forma
c
      if (fr_in.lt.10)then
           forma="(a,i1,a)"
      else
          if (fr_in.lt.100)then
          forma="(a,i2,a)"
          else
          forma="(a,i3,a)"
          endif
      endif
c
      write(flname,forma)"his.dih",fr_in
c
      write(*,"(a,a)")'open file ',flname
c
      OPEN(1,file=flname,status='unknown')
      rewind 1
c
      return
      end
c##################################
c   to generate bined data array 
c##################################
      subroutine databin(min,max,numb,arr)
      dimension arr(numb)
c
      double precision min,max,arr
     *,valra,bin
c
      integer numb,i
c
      valra=max-min
      bin=valra/dble(numb-1)
c
      arr(1)=min
c
      do i=2,numb
      arr(i)=min+dble(i-1)*bin
      enddo
c
      end
c
c###############################
c    to locate a value in an array
c###############################
      subroutine valueloca(arr,numb,val,no)
      dimension arr(numb)
c
      double precision arr,val
     *,bin,bin2,tmp1,tmp2
c
      integer numb,no,loca1,loca2
c
      bin=(arr(numb)-arr(1))/dble(numb-1)
      bin2=bin/2.0d0
c
      tmp1=val-arr(1)
      tmp2=tmp1/bin
      loca1=int(tmp2)+1
      loca2=loca1+1
c
      if((val-arr(loca1)).le.bin2)then
      no=loca1
      else
      no=loca2
      endif
c
      end
