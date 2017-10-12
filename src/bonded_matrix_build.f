c#################################################
c biuld matrix to fit bonded 
c paremeters
c column ordering: bond, angle, dihedral
c row ordering: cpmd (according to storing of bonded forces)
c###############################################
      subroutine bonded_matrix(nat, nr_qm, nr_frame
     *, max_nr_qmb, max_nr_qma, max_nr_qmd 
     *, grom2cpmd
     *, trx, try, trz
     *, bonded_fx, bonded_fy, bonded_fz
c
     *, nr_qmb, ib_qm, jb_qm, nr_qmb_typ, qmb_typ 
     *, icb_qm, req_cal, NUMBND, RK, REQ 
c
     *, nr_qma, ia_qm, ja_qm, ka_qm, nr_qma_typ, qma_typ
     *, ica_qm, teq_cal, NUMANG, TK, TEQ 
c
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm, nr_qmd_typ, qmd_typ
     *, icd_qm, pro_imp, NPTRA, PK, PN, PHASE 
     *, nr_bonded_row, nr_bonded_col, bonded_matr, bonded_targ )
c
      implicit none
c
      integer nat, nr_qm, nr_frame
     *, max_nr_qmb, max_nr_qma, max_nr_qmd
     *, grom2cpmd
     *, nr_qmb, ib_qm, jb_qm, nr_qmb_typ, qmb_typ
     *, icb_qm, NUMBND
     *, nr_qma, ia_qm, ja_qm, ka_qm, nr_qma_typ, qma_typ
     *, ica_qm, NUMANG
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm, nr_qmd_typ, qmd_typ
     *, icd_qm, NPTRA
     *, nr_bonded_row, nr_bonded_col
     *, nr_col, i, j, k, fr_i
     *, ind1, ind2, ind3, ind4
     *, row1, row2, row3, row4, col
     *, icount
       character(3)pro_imp
       double precision trx, try, trz
     *, bonded_fx, bonded_fy, bonded_fz
     *, req_cal, RK, REQ
     *, teq_cal, TK, TEQ
     *, PK, PN, PHASE
     *, bonded_matr, bonded_targ
     *, k_const
     *, fb1, fb2
     *, fa1, fa2, fa3
     *, fd1, fd2, fd3, fd4
     *, dtmp
      dimension grom2cpmd( nat )
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
     *, bonded_fx(nr_frame, nr_qm)
     *, bonded_fy(nr_frame, nr_qm), bonded_fz(nr_frame, nr_qm)
c
     *, ib_qm(max_nr_qmb), jb_qm(max_nr_qmb), icb_qm(max_nr_qmb)
     *, qmb_typ(max_nr_qmb), req_cal( max_nr_qmb )
     *, RK(NUMBND), REQ(NUMBND)
c
     *, ia_qm(max_nr_qma), ja_qm(max_nr_qma), ka_qm(max_nr_qma)
     *, qma_typ( max_nr_qma ), ica_qm( max_nr_qma )
     *, teq_cal(max_nr_qma), TK(NUMANG), TEQ(NUMANG)
c
     *, id_qm(max_nr_qmd), jd_qm(max_nr_qmd)
     *, kd_qm(max_nr_qmd), ld_qm(max_nr_qmd)
     *, qmd_typ(max_nr_qmd), icd_qm(max_nr_qmd), pro_imp(max_nr_qmd)
     *, PK(NPTRA), PN(NPTRA), PHASE(NPTRA)
c
     *, bonded_matr(nr_bonded_row, nr_bonded_col)
     *, bonded_targ(nr_bonded_row)
c
     *, fb1(3), fb2(3)
     *, fa1(3), fa2(3), fa3(3)
     *, fd1(3), fd2(3), fd3(3), fd4(3)
c
      double precision, dimension (:,:), allocatable :: matr_x
     *, matr_y, matr_z
      double precision, dimension (:), allocatable :: targ_x
     *, targ_y, targ_z
c
      write(*,*)''
      write(*,*)'making matrix for bonded fiting'
c
      nr_col = nr_qmb_typ + nr_qma_typ + nr_qmd_typ
      if(nr_col .ne.  nr_bonded_col)then
        stop'bonded_matrix: nr_col .ne.  nr_bonded_col'
      endif
      allocate( matr_x(nr_qm, nr_col), matr_y(nr_qm, nr_col)
     *, matr_z(nr_qm, nr_col) )
      allocate( targ_x(nr_qm), targ_y(nr_qm), targ_z(nr_qm) )
c
      do i=1, nr_bonded_row
        do j=1,nr_bonded_col
           bonded_matr( i, j ) = 0.0d0
        enddo
        bonded_targ( i ) = 0.0d0
      enddo
c
      icount = 0
      do fr_i = 1, nr_frame                    !  frame loop
c
c       initalize
          do i=1, nr_qm
c
             do j=1, nr_col
                 matr_x(i, j) = 0.0d0
                 matr_y(i, j) = 0.0d0
                 matr_z(i, j) = 0.0d0
             enddo
c
             targ_x( i ) = bonded_fx( fr_i, i )
             targ_y( i ) = bonded_fy( fr_i, i )
             targ_z( i ) = bonded_fz( fr_i, i )
          enddo
c
          k_const = 1.0d0
c################################
          do i=1,nr_qmb              !!!!!!!!!!!!!!!!! bond loop
             ind1 = ib_qm( i )
             ind2 = jb_qm( i )
c
             row1 = grom2cpmd( ind1 )
             row2 = grom2cpmd( ind2 )
c
             do k=1,3
               fb1(k) = 0.0d0
               fb2(k) = 0.0d0
             enddo
c
             if( qmb_typ(i) .eq. 0 )then     !  bond type conditions
               call bon_str_force( 1
     *,            RK( icb_qm(i) ), req_cal(i)  ! REQ( icb_qm(i) ) req_cal(i) 
     *,             trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,             trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,                           fb1, fb2 )
c  atom 1
                if( row1 .le. nr_qm )then
                   targ_x( row1 ) = targ_x( row1 ) - fb1(1)
                   targ_y( row1 ) = targ_y( row1 ) - fb1(2)
                   targ_z( row1 ) = targ_z( row1 ) - fb1(3)
                endif
c  atom 2
                if( row2 .le. nr_qm )then
                   targ_x( row2 ) = targ_x( row2 ) - fb2(1)
                   targ_y( row2 ) = targ_y( row2 ) - fb2(2)
                   targ_z( row2 ) = targ_z( row2 ) - fb2(3)
                endif
c
             else                             ! bond type conditions
               if( (qmb_typ(i) .gt. 0) .and.  ! optimized bonds
     *                 (qmb_typ(i) .le. nr_qmb_typ) )then
                  call bon_str_force( 2
     *,              k_const, req_cal(i)
     *,             trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,             trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,                           fb1, fb2 )
c
                    col = qmb_typ(i)
                    if( col .gt. nr_col )then
                     stop'bonded_matrix: qmb_typ(i) .gt. nr_col'
                    endif
c atom 1
                    if( row1 .le. nr_qm )then
                         matr_x( row1, col ) = 
     *                         matr_x( row1, col ) + fb1(1)
                         matr_y( row1, col ) = 
     *                         matr_y( row1, col ) + fb1(2)
                         matr_z( row1, col ) = 
     *                         matr_z( row1, col ) + fb1(3)
                    endif
c atom 2
                    if( row2 .le. nr_qm )then
                         matr_x( row2, col ) = 
     *                         matr_x( row2, col ) + fb2(1)
                         matr_y( row2, col ) = 
     *                         matr_y( row2, col ) + fb2(2)
                         matr_z( row2, col ) = 
     *                         matr_z( row2, col ) + fb2(3)
                    endif
               else
                   stop'bonded_matrix: wrong bond type'
c
               endif                          ! optimized bonds
c
             endif                            ! bond type conditions
c
          enddo                       !!!!!!!!!!!!!!! bond loop
c##########################################
c
           do i=1,nr_qma              !!!!!!!!!!!!!!! angle loop
              ind1 = ia_qm( i ) 
              ind2 = ja_qm( i )
              ind3 = ka_qm( i )
c
              row1 = grom2cpmd( ind1 )
              row2 = grom2cpmd( ind2 )
              row3 = grom2cpmd( ind3 )
c
              do k=1,3
                fa1(k) = 0.0d0
                fa2(k) = 0.0d0
                fa3(k) = 0.0d0
              enddo
c
              if( qma_typ(i) .eq. 0 )then     !  angle type conditions
                 call ang_ben_for( 1
     *,               TK( ica_qm(i) ),teq_cal(i)   ! TEQ( ica_qm(i) ), teq_cal(i)
     *,               trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,               trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,               trx(fr_i, ind3), try(fr_i, ind3), trz(fr_i, ind3)
     *,                           fa1, fa2, fa3 )
c atom 1
                if( row1 .le. nr_qm )then
                   targ_x( row1 ) = targ_x( row1 ) - fa1(1)
                   targ_y( row1 ) = targ_y( row1 ) - fa1(2)
                   targ_z( row1 ) = targ_z( row1 ) - fa1(3)
                endif
c atom 2
                if( row2 .le. nr_qm )then
                   targ_x( row2 ) = targ_x( row2 ) - fa2(1)
                   targ_y( row2 ) = targ_y( row2 ) - fa2(2)
                   targ_z( row2 ) = targ_z( row2 ) - fa2(3)
                endif
c atom 3
                if( row3 .le. nr_qm )then
                   targ_x( row3 ) = targ_x( row3 ) - fa3(1)
                   targ_y( row3 ) = targ_y( row3 ) - fa3(2)
                   targ_z( row3 ) = targ_z( row3 ) - fa3(3)
                endif
c
              else                            !  angle type conditions
                 if( (qma_typ(i) .gt. 0) .and.
     *                 (qma_typ(i) .le. nr_qma_typ) )then   ! optimized angles
                    call ang_ben_for( 2
     *,               k_const, teq_cal(i)
     *,               trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,               trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,               trx(fr_i, ind3), try(fr_i, ind3), trz(fr_i, ind3)
     *,                           fa1, fa2, fa3 )
c
                    col =  qma_typ(i) + nr_qmb_typ
                   if( col .gt. nr_col )then
                     stop'bonded_matrix: qma_typ(i) + nr_qmb_typ 
     *                      .gt. nr_col'
                    endif
c atom 1
                    if( row1 .le. nr_qm )then
                         matr_x( row1, col ) = 
     *                         matr_x( row1, col ) + fa1(1)
                         matr_y( row1, col ) = 
     *                         matr_y( row1, col ) + fa1(2)
                         matr_z( row1, col ) = 
     *                         matr_z( row1, col ) + fa1(3)
                    endif
c atom 2
                    if( row2 .le. nr_qm )then
                         matr_x( row2, col ) = 
     *                         matr_x( row2, col ) + fa2(1)
                         matr_y( row2, col ) = 
     *                         matr_y( row2, col ) + fa2(2)
                         matr_z( row2, col ) = 
     *                         matr_z( row2, col ) + fa2(3)
                    endif
c atom 3
                    if( row3 .le. nr_qm )then
                         matr_x( row3, col ) = 
     *                         matr_x( row3, col ) + fa3(1)
                         matr_y( row3, col ) = 
     *                         matr_y( row3, col ) + fa3(2)
                         matr_z( row3, col ) = 
     *                         matr_z( row3, col ) + fa3(3)
                    endif
c
                 else
                     stop'bonded_matrix: wrong angle type'
c
                 endif                                      ! optimized angles
c
              endif                           !  angle type conditions
c
           enddo                      !!!!!!!!!!!!!!! angle loop
c#####################################
           do i=1,nr_qmd              !!!!!!!!!!!!!!! dihedral loop
              ind1 = id_qm( i ) 
              ind2 = jd_qm( i )
              ind3 = kd_qm( i ) 
              ind4 = ld_qm( i )
c
              row1 = grom2cpmd( ind1 )
              row2 = grom2cpmd( ind2 )
              row3 = grom2cpmd( ind3 )
              row4 = grom2cpmd( ind4 )
c
              do k=1,3
                 fd1(k) = 0.0d0
                 fd2(k) = 0.0d0
                 fd3(k) = 0.0d0
                 fd4(k) = 0.0d0
              enddo
c
              if( qmd_typ(i) .eq. 0 )then      !  dihedral type conditions
                 if(pro_imp(i) .eq. 'PRO')then   ! PRO - IMP
                   call dihedr_rot_for( 1
     *,        PK( icd_qm(i) ), PN( icd_qm(i) ), PHASE( icd_qm(i) )
     *,           trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,           trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,           trx(fr_i, ind3), try(fr_i, ind3), trz(fr_i, ind3)
     *,           trx(fr_i, ind4), try(fr_i, ind4), trz(fr_i, ind4)
     *,                 fd1, fd2, fd3, fd4 )
c
                 else
                     if(pro_imp(i) .eq. 'IMP')then
                      call impro_for( 1
     *,        PK( icd_qm(i) ), PN( icd_qm(i) ), PHASE( icd_qm(i) )
     *,           trx(fr_i, ind3), try(fr_i, ind3), trz(fr_i, ind3)
     *,           trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,           trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,           trx(fr_i, ind4), try(fr_i, ind4), trz(fr_i, ind4)
     *,                 fd3, fd2, fd1, fd4 )
c
                     else
                         stop'bonded_matrix: wrong PRO IMP'
c
                     endif
                 endif                             ! PRO - IMP
c atom 1
                if( row1 .le. nr_qm )then
                   targ_x( row1 ) = targ_x( row1 ) - fd1(1)
                   targ_y( row1 ) = targ_y( row1 ) - fd1(2)
                   targ_z( row1 ) = targ_z( row1 ) - fd1(3)
                endif
c atom 2
                if( row2 .le. nr_qm )then
                   targ_x( row2 ) = targ_x( row2 ) - fd2(1)
                   targ_y( row2 ) = targ_y( row2 ) - fd2(2)
                   targ_z( row2 ) = targ_z( row2 ) - fd2(3)
                endif
c atom 3
                if( row3 .le. nr_qm )then
                   targ_x( row3 ) = targ_x( row3 ) - fd3(1)
                   targ_y( row3 ) = targ_y( row3 ) - fd3(2)
                   targ_z( row3 ) = targ_z( row3 ) - fd3(3)
                endif
c atom 4
                if( row4 .le. nr_qm )then
                   targ_x( row4 ) = targ_x( row4 ) - fd4(1)
                   targ_y( row4 ) = targ_y( row4 ) - fd4(2)
                   targ_z( row4 ) = targ_z( row4 ) - fd4(3)
                endif
c
              else                             !  dihedral type conditions
                  if( (qmd_typ(i) .gt. 0) .and. 
     *                (qmd_typ(i) .le. nr_qmd_typ) )then  ! optimized dihedrals
c
                      if(pro_imp(i) .eq. 'PRO')then   ! PRO - IMP
                      call dihedr_rot_for( 2
     *,           k_const, PN( icd_qm(i) ), PHASE( icd_qm(i) )
     *,           trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,           trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,           trx(fr_i, ind3), try(fr_i, ind3), trz(fr_i, ind3)
     *,           trx(fr_i, ind4), try(fr_i, ind4), trz(fr_i, ind4)
     *,                 fd1, fd2, fd3, fd4 )
c
                       else
                          if(pro_imp(i) .eq. 'IMP')then
                        call impro_for( 2
     *,           k_const, PN( icd_qm(i) ), PHASE( icd_qm(i) )
     *,           trx(fr_i, ind3), try(fr_i, ind3), trz(fr_i, ind3)
     *,           trx(fr_i, ind2), try(fr_i, ind2), trz(fr_i, ind2)
     *,           trx(fr_i, ind1), try(fr_i, ind1), trz(fr_i, ind1)
     *,           trx(fr_i, ind4), try(fr_i, ind4), trz(fr_i, ind4)
     *,                 fd3, fd2, fd1, fd4 )
c
                        else
                            stop'bonded_matrix: wrong PRO IMP'
c
                         endif
                     endif                             ! PRO - IMP
c
                       col =  qmd_typ(i) + nr_qmb_typ + nr_qma_typ
                     if( col .gt. nr_col )then
                     stop'bonded_matrix:  qmd_typ(i) + nr_qmb_typ 
     *                      + nr_qma_typ .gt. nr_col'
                    endif
c atom 1
                    if( row1 .le. nr_qm )then
                         matr_x( row1, col ) = 
     *                         matr_x( row1, col ) + fd1(1)
                         matr_y( row1, col ) = 
     *                         matr_y( row1, col ) + fd1(2)
                         matr_z( row1, col ) = 
     *                         matr_z( row1, col ) + fd1(3)
                    endif
c atom 2
                    if( row2 .le. nr_qm )then
                         matr_x( row2, col ) = 
     *                         matr_x( row2, col ) + fd2(1)
                         matr_y( row2, col ) = 
     *                         matr_y( row2, col ) + fd2(2)
                         matr_z( row2, col ) = 
     *                         matr_z( row2, col ) + fd2(3)
                    endif
c atom 3
                    if( row3 .le. nr_qm )then
                         matr_x( row3, col ) = 
     *                         matr_x( row3, col ) + fd3(1)
                         matr_y( row3, col ) = 
     *                         matr_y( row3, col ) + fd3(2)
                         matr_z( row3, col ) = 
     *                         matr_z( row3, col ) + fd3(3)
                    endif
c atom 4
                    if( row4 .le. nr_qm )then
                         matr_x( row4, col ) = 
     *                         matr_x( row4, col ) + fd4(1)
                         matr_y( row4, col ) = 
     *                         matr_y( row4, col ) + fd4(2)
                         matr_z( row4, col ) = 
     *                         matr_z( row4, col ) + fd4(3)
                    endif
c
                  else
                      stop'bonded_matrix: wrong dihedral type'
c
                  endif                                   ! optimized dihedrals
c
              endif                            !  dihedral type conditions
c
           enddo                      !!!!!!!!!!!!!!! dihedral loop
c
          do i=1, nr_qm
c x
             icount = icount + 1
             do j=1,nr_bonded_col
               bonded_matr( icount, j ) = matr_x( i, j )
             enddo
             bonded_targ( icount ) = targ_x( i )
c y
             icount = icount + 1
             do j=1,nr_bonded_col
               bonded_matr( icount, j ) = matr_y( i, j )
             enddo
             bonded_targ( icount ) = targ_y( i )
c z
             icount = icount + 1
             do j=1,nr_bonded_col
               bonded_matr( icount, j ) = matr_z( i, j )
             enddo
             bonded_targ( icount ) = targ_z( i )
c
           if(icount .gt. nr_bonded_row)then
             stop'bonded_matrix: icount .gt. nr_bonded_row'
           endif
c
          enddo
c
      enddo                                    !  frame loop
c
      do i=1,nr_bonded_row
       dtmp = 0.0d0
c
          do j=1,nr_bonded_col
             dtmp = dtmp + dabs( bonded_matr( i, j ) )
          enddo
        if(dtmp .eq. 0.0d0)then
            bonded_targ( i ) = 0.0d0
        endif
c
      enddo
c
c      open(1,file='bonded_matrix.dat',status='unknown')
c      rewind(1)
c      do i=1, nr_bonded_row
c      write(1,"(30f15.7)")(bonded_matr( i, j ),j=1,nr_bonded_col)
c     *,  bonded_targ( i )
c      enddo
c      close(1)
c
      write(*,*)'making matrix for bonded fiting done'
      write(*,*)''
c
      return
      end
