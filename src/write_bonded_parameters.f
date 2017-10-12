c###############################
c  to write the optimimized
c  bonded parameters
c###############################
      subroutine write_bonded(nat, nr_qm, nr_frame 
     *, max_nr_qmb, max_nr_qma, max_nr_qmd
     *, ISYMBL, IGRAPH
c
     *, nr_qmb, ib_qm, jb_qm, nr_qmb_typ, qmb_typ
     *, icb_qm, req_cal, NUMBND, RK, REQ
c
     *, nr_qma, ia_qm, ja_qm, ka_qm, nr_qma_typ, qma_typ
     *, ica_qm, teq_cal, NUMANG, TK, TEQ
c
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm, nr_qmd_typ, qmd_typ
     *, icd_qm, pro_imp, NPTRA, PK, PN, PHASE
     *, nr_bonded_col, bonded_solu )
c
      implicit none
      integer nat, nr_qm, nr_frame 
     *, max_nr_qmb, max_nr_qma, max_nr_qmd
c
     *, nr_qmb, ib_qm, jb_qm, nr_qmb_typ, qmb_typ
     *, icb_qm, NUMBND
c
     *, nr_qma, ia_qm, ja_qm, ka_qm, nr_qma_typ, qma_typ
     *, ica_qm, NUMANG
c
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm, nr_qmd_typ, qmd_typ
     *, icd_qm, NPTRA
     *, nr_bonded_col
     *, i, j, icount
      character(4) ISYMBL, IGRAPH
      character(3) pro_imp
      double precision 
     *  req_cal, RK, REQ 
     *, teq_cal, TK, TEQ
     *, PK, PN, PHASE
     *, bonded_solu
     *, bond_av, angle_av
     *, pi
      dimension ISYMBL( nat ), IGRAPH( nat )
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
     *, bonded_solu(nr_bonded_col)
     *, bond_av(nr_qmb_typ), angle_av(nr_qma_typ)
c
      write(*,*)''
      write(*,*)'----------------------------'
      write(*,*)'writing the bonded parameters'
c
      pi=dacos(-1.0d0)
c bonds
      do i=1, nr_qmb_typ              !  bond type loop
        bond_av( i ) = 0.0d0
        icount = 0
c
        do j=1, nr_qmb                !  all bond loop
          if( qmb_typ(j) .eq. i )then
            icount = icount + 1
            bond_av( i ) = bond_av( i ) + req_cal( j )
          endif
        enddo                         !  all bond loop
c
        if(icount .eq. 0)then
          stop'write_bonded: some bond type has no member'
        endif
        bond_av( i ) = bond_av( i ) / dble(icount)
c
      enddo                           !  bond type loop
c
c angles
c
      do i=1, nr_qma_typ              !  angle type loop
         angle_av( i ) = 0.0d0
         icount = 0
c
         do j=1, nr_qma               ! all angle loop
           if( qma_typ(j) .eq. i)then
             icount = icount + 1
             angle_av( i ) = angle_av( i ) + teq_cal( j )
           endif
         enddo                        ! all angle loop
c
         if(icount .eq. 0)then
          stop'write_bonded: some angle type has no member'
        endif
        angle_av( i ) = angle_av( i ) / dble(icount)
c
      enddo                           !  angle type loop
c
      write(*,*)'details of parameters witten to bonded_fit_report.out'
      open(1,file='bonded_fit_report.out',status='unknown')
      rewind(1)
      write(1,*)'Affected bonds'
      write(1,*)''
      write(1,*)'index, at1, at2, name1, name2, type1, type2, fm_type
     *, fm_K, fm_Req, AMBER_K, AMBER_Req'
      write(*,*)''
c
      do i=1,nr_qmb
         if(qmb_typ(i) .gt. 0)then
           write(1,"(3i6,4a6,i6,4f10.4)") i
     *,      ib_qm(i), jb_qm(i)
     *,      IGRAPH( ib_qm(i) ), IGRAPH( jb_qm(i) )
     *,      ISYMBL( ib_qm(i) ), ISYMBL( jb_qm(i) )
     *,      qmb_typ(i)
     *,      bonded_solu( qmb_typ(i) ), bond_av( qmb_typ(i) )
     *,      RK( icb_qm(i) ), REQ( icb_qm(i) )
         endif
      enddo
c########
c
      write(1,*)''
      write(1,*)'---------------------'
      write(1,*)''
      write(1,*)'Affected angles'
      write(1,*)''
      write(1,*)'index, at1, at2, at3, name1, name2, name3
     *, type1, type2, type3, fm_type, fm_KT, fm_Teq
     *, AMBER_KT, AMBER_Teq'
      write(*,*)''
c
      do i=1, nr_qma
         if(qma_typ(i) .gt. 0)then
           write(1,"(4i6,6a6,i6,4f10.4)") i
     *, ia_qm(i), ja_qm(i), ka_qm(i)
     *, IGRAPH( ia_qm(i) ), IGRAPH( ja_qm(i) ), IGRAPH( ka_qm(i) )
     *, ISYMBL( ia_qm(i) ), ISYMBL( ja_qm(i) ), ISYMBL( ka_qm(i) )
     *, qma_typ(i)
     *, bonded_solu( qma_typ(i) + nr_qmb_typ )
     *, angle_av( qma_typ(i) )*180.0d0/pi
     *, TK( ica_qm(i) ), TEQ( ica_qm(i) )*180.0d0/pi
         endif
      enddo
c
      write(1,*)''
      write(1,*)'---------------------'
      write(1,*)''
      write(1,*)'Affected dihedrals'
      write(1,*)''
      write(1,*)'index, at1, at2, at3, at4
     *, name1, name2, name3, name4
     *, type1, type2, type3, type4
     *, fm_type, PRO/IMP
     *, fm_Vn, AMBER_Vn, AMBER_N, AMBER_phase'
c
      do i=1, nr_qmd          !  dihedral loop
         if(qmd_typ(i) .gt. 0)then
          write(1,"(5i6,8a6,i6,a6,4f10.4)") i
     *, id_qm(i), jd_qm(i), kd_qm(i), ld_qm(i)
     *, IGRAPH( id_qm(i) ), IGRAPH( jd_qm(i) )
     *, IGRAPH( kd_qm(i) ), IGRAPH( ld_qm(i) )
     *, ISYMBL( id_qm(i) ), ISYMBL( jd_qm(i) )
     *, ISYMBL( kd_qm(i) ), ISYMBL( ld_qm(i) )
     *, qmd_typ(i), pro_imp(i)  
     *, bonded_solu( qmd_typ(i) + nr_qmb_typ + nr_qma_typ )
     *, PK( icd_qm(i) ), PN( icd_qm(i) )
     *, PHASE( icd_qm(i) )*180.0d0/pi
         endif
      enddo                   !  dihedral loop
c
      close(1)
c
      write(*,*)'the amber modification file written to fm.frcmod'
      open(1,file='fm.frcmod',status='unknown')
      rewind(1)
c
      write(1,"(a13)")'fore matching'
      write(1,"(a4)")'MASS'
      write(1,*)''
      write(1,"(a4)")'BOND'
      do i=1,nr_qmb_typ
        do j=1,nr_qmb
           if(qmb_typ(j) .eq. i)then
               write(1,"(a3,a1,a3,2f15.5)")ISYMBL( ib_qm(j) )
     *,   '-', ISYMBL( jb_qm(j) )
     *,   bonded_solu( qmb_typ(j) ), bond_av( qmb_typ(j) )
           go to 1000
           endif
        enddo
 1000    continue
      enddo
      write(1,*)''
c
      write(1,"(a4)")'ANGL'
      do i=1, nr_qma_typ
         do j=1,nr_qma
           if( qma_typ(j) .eq. i )then
            write(1,"(a3,a1,a3,a1,a3,2f15.5)")ISYMBL( ia_qm(j) ),'-'
     *,   ISYMBL( ja_qm(j) ),'-', ISYMBL( ka_qm(j) )
     *,   bonded_solu( qma_typ(j) + nr_qmb_typ )
     *, angle_av( qma_typ(j) )*180.0d0/pi
           go to 1001
           endif
         enddo
 1001    continue
      enddo
      write(1,*)''
c
      write(1,"(a8)")'DIHEDRAL'
      do i=1,nr_qmd_typ
         do j=1,nr_qmd
c
         if( (qmd_typ(j) .eq. i) .and. (pro_imp(j).eq."PRO") )then
          write(1,"(a3,a1,a3,a1,a3,a1,a3,i6,3f15.5)")
     *  ISYMBL( id_qm(j) ),'-', ISYMBL( jd_qm(j) ),'-'
     *, ISYMBL( kd_qm(j) ),'-', ISYMBL( ld_qm(j) )
     *, 1, bonded_solu( qmd_typ(j) + nr_qmb_typ + nr_qma_typ )
     *, PHASE( icd_qm(j) )*180.0d0/pi
     *, PN( icd_qm(j) )
         go to 1002
         endif
         enddo
 1002    continue
      enddo
      write(1,*)''
c
      write(1,"(a8)")'IMPROPER'
      do i=1,nr_qmd_typ
         do j=1,nr_qmd
c
         if( (qmd_typ(j) .eq. i) .and. (pro_imp(j).eq."IMP") )then
          write(1,"(a3,a1,a3,a1,a3,a1,a3,i6,3f15.5)")
     *  ISYMBL( id_qm(j) ),'-', ISYMBL( jd_qm(j) ),'-'
     *, ISYMBL( kd_qm(j) ),'-', ISYMBL( ld_qm(j) )
     *, 1, bonded_solu( qmd_typ(j) + nr_qmb_typ + nr_qma_typ )
     *, PHASE( icd_qm(j) )*180.0d0/pi
     *, PN( icd_qm(j) )
         go to 1003
         endif
         enddo
 1003    continue
      enddo
      write(1,*)''
c
      write(1,"(a4)")'NONB'
      write(1,*)''
c
      close(1)
c
      write(*,*)'writing the bonded parameters done!'
      write(*,*)''
c
      return
      end
