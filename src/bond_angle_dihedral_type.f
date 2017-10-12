c######################################
c  counting and classifying 
c  unique bond types.
c  note that the bond types here are not 
c  the AMBER index into RK, REQ parameters
c  They include only atoms in QM zone
c  those outside (or dont need to be fitted)
c  will be set to zero.
c  they will be used as index of fitted variables,
c  not to get parameters from AMBER  
c qmb_typ = 0 means dont fit but calculate and put to the right side
c qmb_typ = 1, 2, 3 ... means 'will be fitted'
c####################################
      subroutine bond_types( nat, nr_qm, max_nr_qmb
     *, grom2cpmd
     *, NATOM, ISYMBL, IGRAPH
     *, NBONH, IBH, JBH, ICBH
     *, NBONA, IB, JB, ICB 
     *, nr_qmb, ib_qm, jb_qm
     *, nr_qmb_typ, qmb_typ 
     *, icb_qm )
c
      implicit none
      integer nat, nr_qm, max_nr_qmb
     *, grom2cpmd
     *, NATOM
     *, NBONH, IBH, JBH, ICBH
     *, NBONA, IB, JB, ICB
     *, nr_qmb, nr_qmb_typ
     *, ib_qm, jb_qm, icb_qm, qmb_typ
     *, i, j, icount, jcount, itmp, jtmp
     *, bond_qm_or_not, i_test, j_test
     *, nr_b_add, nr_b_remove
     *, b_add, b_remove
      character(4) ISYMBL, IGRAPH
      character(30)filename,ctmp
      dimension grom2cpmd( nat )
     *, ISYMBL( NATOM ), IGRAPH( NATOM )
     *, IBH(NBONH), JBH(NBONH), ICBH(NBONH)
     *, IB(NBONA), JB(NBONA), ICB(NBONA)
     *, ib_qm(max_nr_qmb), jb_qm(max_nr_qmb), icb_qm(max_nr_qmb)
     *, qmb_typ(max_nr_qmb)
c
     *,  b_add( max_nr_qmb ), b_remove( max_nr_qmb )
c
       write(*,*)' '
       write(*,*)'----------------------'
       write(*,*)'Counting no of bonds involved in the QM zone'
c
c  look for all the bonds (including H) with at least one atom is QM
c  the total number of these will be nr_qmb
c
c   then decide which bonds are surely not optimized e.g. boundary, involve H ... put to 0
c    the remaining will be optimized, and will be initalized to -1
c
c    sort the remmaining bonds (only -1) into ordered qmb_typ = 1, 2, ... nr_qmb_typ
c
      if(NBONH .le. 0)then
        stop'bond_types: NBONH .le. 0'
      endif
c
      icount = 0
      do i=1, NBONH         !   Bond_Inc_Hydro loop
        itmp = abs(IBH(i))/3 + 1
        jtmp = abs(JBH(i))/3 + 1
c
        i_test = bond_qm_or_not( itmp, jtmp
     *, nat, nr_qm, grom2cpmd)
c
        if( i_test .ne. 0 ) then
           icount = icount + 1
           ib_qm( icount ) = itmp
           jb_qm( icount ) = jtmp
           icb_qm( icount ) = ICBH(i)
c
           if(i_test .eq. -1) then
                qmb_typ(icount) = 0
           else
                qmb_typ(icount) = -1
           endif
c
        endif
c
      enddo         !   Bond_Inc_Hydro loop
c#######
c  the decision not to fit bond_inc_H goes here  qmb_typ() = 0
c
      write(*,*)''
      write(*,*)'--------------------'
      write(*,*)'do you want to fit bonds involved Hydrogens ?'
      write(*,*)'no: enter 0'
      write(*,*)'yes: enter 1'
      read(*,*)j_test
c
      if(j_test .eq. 0)then
        write(*,*)'Bonds included H will not be fitted'
        do i=1,icount
           qmb_typ(i) = 0
        enddo
      else
         if(j_test .ne. 1)then
           stop'bond_types: dont know what to do'
         endif
         write(*,*)'Bonds included H will be fitted'
      endif 
c
c#######
c
      if(NBONA .le. 0)then
        stop'bond_types: NBONA .le. 0'
      endif
c
      do i=1, NBONA         !   Bond_Not_H loop
        itmp = abs(IB(i))/3 + 1
        jtmp = abs(JB(i))/3 + 1
c
        i_test = bond_qm_or_not( itmp, jtmp
     *, nat, nr_qm, grom2cpmd)
c
        if( i_test .ne. 0 ) then
           icount = icount + 1
           ib_qm( icount ) = itmp
           jb_qm( icount ) = jtmp
           icb_qm( icount ) = ICB(i)
c
           if(i_test .eq. -1) then
                qmb_typ(icount) = 0
           else
                qmb_typ(icount) = -1
           endif
c
        endif
c
      enddo         !   Bond_Not_H loop
c
      if( icount .gt. max_nr_qmb) then
      stop'bond_types: icount .gt. max_nr_qmb '
      endif
c
      nr_qmb = icount
      if(nr_qmb .eq. 0)then
      stop'bond_types: nr_qmb .eq. 0'
      endif
c
c  print all the bonds in QM zone
c
      write(*,*)''
      write(*,*)'-----------------------------'
      write(*,*)'these are all the bonds in QM zone'
      write(*,*)'-1 in the last column are those affected '
      write(*,*)'0 are those unaffected'
      write(*,*)' '
      write(*,*)'bond_ind; at_ind1, at_ind2, at_name_1, at_name_2;' 
      write(*,*)'at_type_1, at_type_2'
      write(*,*)' '
      do i=1,nr_qmb
          write(*,"(3i6,4a6,i8)")i,ib_qm( i ),jb_qm( i )
     *,    IGRAPH( ib_qm(i) ), IGRAPH( jb_qm(i) )
     *,    ISYMBL( ib_qm(i) ), ISYMBL( jb_qm(i) )
     *, qmb_typ(i)
      enddo
c
c   print all the bonds affected
c
c      write(*,*)''
c      write(*,*)'-----------------------------'
c      write(*,*)'these bonds will be affected'
c      write(*,*)'-1 in the last column are those affected '
c      write(*,*)' '
c      write(*,*)'bond_ind; at_ind1, at_ind2, at_name_1, at_name_2;' 
c      write(*,*)'at_type_1, at_type_2'
c      write(*,*)' '
c      do i=1,nr_qmb
c       if(qmb_typ(i) .eq. -1)then
c          write(*,"(3i6,4a6)")i,ib_qm( i ),jb_qm( i )
c     *,    IGRAPH( ib_qm(i) ), IGRAPH( jb_qm(i) )
c     *,    ISYMBL( ib_qm(i) ), ISYMBL( jb_qm(i) )
c       endif
c      enddo
c
      write(*,*)''
      write(*,*)'---------------------------------------'
      write(*,*)'here, you will have a chance to further '
      write(*,*)'specify which bonds you want or dont want to'
      write(*,*)'include in the fit'
      write(*,*)'if you satisfy with the above list, enter 1'
      write(*,*)'if you want to add and(or) remove'
      write(*,*)'some more bonds, enter 2'
      read(*,*)j_test
c
      nr_b_add = 0
      nr_b_remove = 0
      if(j_test.eq.2)then
          write(*,*)'-----------'
          write(*,*)'now, please prepare a file, in which you specify:'
          write(*,*)'NR_add'
          write(*,*)'bond_ind_1'
          write(*,*)'bond_ind_2'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'bond_ind_NR_add'
          write(*,*)'#'
          write(*,*)'NR_revmove'
          write(*,*)'bond_ind_1'
          write(*,*)'bond_ind_2'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'bond_ind_NR_remove'
          write(*,*)''
          write(*,*)'now, give me the name of that file:'
          write(*,*)'( 30 characters or less )'
          read(*,*)filename
          open(1,file=filename,status='old')
          rewind(1)
c
          read(1,*)nr_b_add
          write(*,"(i8)")nr_b_add
          if(nr_b_add .gt. max_nr_qmb) then
            stop'bond_types: nr_b_add .gt. max_nr_qmb'
          endif
c
          if(nr_b_add .gt. 0)then
             do i=1,nr_b_add
               read(1,*)b_add(i)
               write(*,"(i8)")b_add(i)
             enddo
          endif
c
          read(1,*)ctmp
          write(*,"(a)")ctmp
c
          read(1,*)nr_b_remove
          write(*,"(i8)")nr_b_remove
          if(nr_b_remove .gt. max_nr_qmb) then
            stop'bond_types: nr_b_remove .gt. max_nr_qmb'
          endif
c
          if(nr_b_remove .gt. 0)then
             do i=1,nr_b_remove
               read(1,*)b_remove(i)
               write(*,"(i8)")b_remove(i)
             enddo
          endif
          close(1)
      else
         if(j_test.ne.1)then
            stop'bond_types: dont know what to do'
         endif
      endif
c
c  now add more
c 
      if(nr_b_add .gt. 0)then
         do i=1,nr_b_add
         qmb_typ( b_add(i) ) = -1
         enddo
      endif
c
c now remove some
c
      if(nr_b_remove .gt. 0)then
        do i=1,nr_b_remove
         qmb_typ(b_remove(i)) = 0
        enddo
      endif
c
c  now search among -1 to sort bond types
c
      jcount = 0
      do i=1,nr_qmb
        if( qmb_typ( i ) .eq. -1)then
            jcount = jcount + 1
            qmb_typ( i ) = jcount
c
            do j=i+1,nr_qmb
               if( qmb_typ( j ) .eq. -1)then
                  if( icb_qm( i ) .eq. icb_qm( j ) )then
                   qmb_typ( j ) = jcount
                  endif
                endif
            enddo
        endif
      enddo
c
      nr_qmb_typ = jcount
      if(nr_qmb_typ .eq. 0)then
      write(*,*)'bond_types: WARNING: no bond to fit'
      endif
c
c
      write(*,*)' '
      write(*,*)'--------------'
      write(*,*)'Now, below are new list of bonds'
      write(*,*)'zero in the last column means not fitted'
      write(*,*)' '
      write(*,*)'bond_ind; at_ind1, at_ind2, at_name_1, at_name_2;' 
      write(*,*)'at_type_1, at_type_2, unique_bond_type'
      write(*,*)' '
      do i=1,nr_qmb
          write(*,"(3i6,4a6,i8)")i,ib_qm( i ),jb_qm( i )
     *,    IGRAPH( ib_qm(i) ), IGRAPH( jb_qm(i) )
     *,    ISYMBL( ib_qm(i) ), ISYMBL( jb_qm(i) )
     *,    qmb_typ( i )
      enddo
c
      write(*,*)' '
      write(*,"(a,i8)")'Number of bond types: ',nr_qmb_typ
c
      write(*,*)''
      write(*,*)'finding and classifying bonds done !'
      write(*,*)''
c
      return
      end
c#####################################
c  counting total of angles in QM zone
c and classifying unique angle types
c note that the angle types here are not
c AMBER index into TK, TEQ parameters
c They include only atoms in QM zone
c those outside (or dont need to be fitted)
c will be set to zero.
c they will be used as index of fitted variables,
c not to get parameters from AMBER
c qma_typ = 0 means dont fit but calculate and put to the right side
c qma_typ = 1, 2, 3 ... means 'will be fitted'
c#####################################
      subroutine angle_types(nat, nr_qm, max_nr_qma
     *, grom2cpmd
     *, NATOM, ISYMBL, IGRAPH
     *, NTHETH, ITH, JTH, KTH, ICTH
     *, NTHETA, IT,JT,KT,ICT
     *, nr_qma, ia_qm, ja_qm, ka_qm
     *, nr_qma_typ, qma_typ
     *, ica_qm )
c
      implicit none
      integer nat, nr_qm, max_nr_qma
     *, grom2cpmd
     *, NATOM
     *, NTHETH, ITH, JTH, KTH, ICTH
     *, NTHETA, IT,JT,KT,ICT
     *, nr_qma, ia_qm, ja_qm, ka_qm
     *, nr_qma_typ, qma_typ
     *, ica_qm
     *, i, j, icount, jcount, itmp, jtmp, ktmp
     *, ang_qm_or_not, i_test, j_test
     *, nr_a_add, nr_a_remove
     *, a_add, a_remove
      character(4) ISYMBL, IGRAPH
      character(30)filename,ctmp
      dimension grom2cpmd( nat )
     *, ISYMBL( NATOM ), IGRAPH( NATOM )
     *, ITH(NTHETH), JTH(NTHETH), KTH(NTHETH), ICTH(NTHETH)
     *, IT(NTHETA), JT(NTHETA), KT(NTHETA), ICT(NTHETA)
     *, ia_qm(max_nr_qma), ja_qm(max_nr_qma), ka_qm(max_nr_qma)
     *, qma_typ( max_nr_qma ), ica_qm( max_nr_qma )
c
     *, a_add( max_nr_qma ), a_remove( max_nr_qma )
c
       write(*,*)' '
       write(*,*)'-------------------------------'
       write(*,*)'Counting no of angles involved in the QM zone'
c
c  look for all the angles (including H) with at least one atom is QM
c  the total number of these will be nr_qma
c
c   then decide which angles are surely not optimized e.g. boundary, involve H ... put qma_typ to 0
c    the remaining will be optimized, and will be initalized to -1
c
c    sort the remmaining bonds (only -1) into ordered qma_typ = 1, 2, ... nr_qmb_typ
c
      if(NTHETH .le. 0)then
        stop'angle_types: NTHETH .le. 0'
      endif
c
      icount = 0
      do i=1,NTHETH            !  Angle_Inc_Hydro loop
        itmp = abs( ITH(i) )/3 + 1
        jtmp = abs( JTH(i) )/3 + 1
        ktmp = abs( KTH(i) )/3 + 1
c
        i_test = ang_qm_or_not( itmp, jtmp, ktmp
     *, nat, nr_qm, grom2cpmd )
c
        if(i_test .ne. 0) then
           icount = icount + 1
           ia_qm( icount ) = itmp
           ja_qm( icount ) = jtmp
           ka_qm( icount ) = ktmp
           ica_qm( icount ) = ICTH( i )
c
           if( i_test .eq. -1 )then
               qma_typ( icount ) = 0
           else
               qma_typ( icount ) = -1
           endif
        endif
      enddo                    !  Angle_Inc_Hydro loop
c################
c the decision not to fit angle_inc_H goes here  qma_typ() = 0
c
      write(*,*)''
      write(*,*)'--------------------'
      write(*,*)'do you want to fit angles involved Hydrogens ?'
      write(*,*)'no: enter 0'
      write(*,*)'yes: enter 1'
      read(*,*)j_test
c
      if(j_test .eq. 0)then
        write(*,*)'Angles included H will not be fitted'
        do i=1,icount
           qma_typ(i) = 0
        enddo
      else
         if(j_test .ne. 1)then
           stop'angle_types: dont know what to do'
         endif
      write(*,*)'Angles included H will be fitted'
      endif 
c
c###########
c
      if(NTHETA .le. 0)then
        stop'angle_types: NTHETA .le. 0'
      endif
c
      do i=1,NTHETA       ! angles_Not_H loop
        itmp = abs( IT(i) )/3 + 1
        jtmp = abs( JT(i) )/3 + 1
        ktmp = abs( KT(i) )/3 + 1
c
        i_test = ang_qm_or_not( itmp, jtmp, ktmp
     *, nat, nr_qm, grom2cpmd )
c
        if(i_test .ne. 0) then
           icount = icount + 1
           ia_qm( icount ) = itmp
           ja_qm( icount ) = jtmp
           ka_qm( icount ) = ktmp
           ica_qm( icount ) = ICT( i )
c
           if( i_test .eq. -1 )then
               qma_typ( icount ) = 0
           else
               qma_typ( icount ) = -1
           endif
        endif
        
      enddo               ! angles_Not_H loop
c
      if( icount .gt. max_nr_qma ) then
      stop'angle_types: icount .gt. max_nr_qma '
      endif
c
      nr_qma = icount
      if(nr_qma .eq. 0)then
      stop'angle_types: nr_qma .eq. 0'
      endif
c
c  print all the angles in QM zone 
c 
      write(*,*)''
      write(*,*)'-----------------------------'
      write(*,*)'these are all the angles in QM zone'
      write(*,*)'-1 in the last column are those affected '
      write(*,*)'0 are those unaffected'
      write(*,*)''
      write(*,*)'angle_ind; at_ind1, at_ind2, at_ind3'
      write(*,*)'at_name_1, at_name_2, at_name_3' 
      write(*,*)'at_type_1, at_type_2, at_type_3'
      write(*,*)' '
      do i=1,nr_qma
      write(*,"(4i6,6a6,i8)")i,ia_qm( i ),ja_qm( i ), ka_qm( i )
     *,    IGRAPH( ia_qm(i) ), IGRAPH( ja_qm(i) ), IGRAPH( ka_qm(i) )
     *,    ISYMBL( ia_qm(i) ), ISYMBL( ja_qm(i) ), ISYMBL( ka_qm(i) )
     *,    qma_typ( i )
      enddo
c
c  print all the angles affected
c
c      write(*,*)''
c      write(*,*)'-----------------------------'
c      write(*,*)'these angles will be affected'
c      write(*,*)''
c      write(*,*)'angle_ind; at_ind1, at_ind2, at_ind3'
c      write(*,*)'at_name_1, at_name_2, at_name_3' 
c      write(*,*)'at_type_1, at_type_2, at_type_3'
c      write(*,*)' '
c      do i=1,nr_qma
c         if( qma_typ( i ) .eq. -1 )then
c            write(*,"(4i6,6a6)")i,ia_qm( i ),ja_qm( i ), ka_qm( i )
c     *,     IGRAPH( ia_qm(i) ), IGRAPH( ja_qm(i) ), IGRAPH( ka_qm(i) )
c     *,     ISYMBL( ia_qm(i) ), ISYMBL( ja_qm(i) ), ISYMBL( ka_qm(i) )
c         endif
c      enddo
c
      write(*,*)' '
       write(*,*)'---------------------------------------'
      write(*,*)'here, you will have a chance to further '
      write(*,*)'specify which angles you want or dont want to'
      write(*,*)'include in the fit'
      write(*,*)'if you satisfy with the above list, enter 1'
      write(*,*)'if you want to add and(or) remove'
      write(*,*)'some more angles, enter 2'
      read(*,*)j_test
c
      nr_a_add = 0
      nr_a_remove = 0
      if(j_test.eq.2)then
          write(*,*)'-----------'
          write(*,*)'now, please prepare a file, in which you specify:'
          write(*,*)'NR_add'
          write(*,*)'angle_ind_1'
          write(*,*)'angle_ind_2'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'angle_ind_NR_add'
          write(*,*)'#'
          write(*,*)'NR_revmove'
          write(*,*)'angle_ind_1'
          write(*,*)'angle_ind_2'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'angle_ind_NR_remove'
          write(*,*)''
          write(*,*)'now, give me the name of that file:'
          write(*,*)'( 30 characters or less )'
          read(*,*)filename
          open(1,file=filename,status='old')
          rewind(1)
c
          read(1,*)nr_a_add
          write(*,"(i8)")nr_a_add
          if( nr_a_add .gt. max_nr_qma)then 
            stop'angle_types: nr_a_add .gt.  max_nr_qma'
          endif
c
          if(nr_a_add .gt. 0)then
            do i=1,nr_a_add
              read(1,*)a_add(i)
              write(*,"(i8)")a_add(i)
            enddo
          endif
c
          read(1,*)ctmp
          write(*,"(a)")ctmp
c
          read(1,*)nr_a_remove
          write(*,"(i8)")nr_a_remove
          if( nr_a_remove .gt. max_nr_qma)then 
            stop'angle_types: nr_a_remove .gt.  max_nr_qma'
          endif
c
          if(nr_a_remove .gt. 0)then
             do i=1,nr_a_remove
                read(1,*)a_remove(i)
                write(*,"(i8)")a_remove(i)
             enddo
          endif
c
          close(1)
      else
         if(j_test.ne.1)then
            stop'angle_types: dont know what to do'
         endif
      endif
c
c  now add more
c 
      if(nr_a_add .gt. 0)then
         do i=1,nr_a_add
         qma_typ( a_add(i) ) = -1
         enddo
      endif
c
c now remove some
c
      if(nr_a_remove .gt. 0)then
        do i=1,nr_a_remove
         qma_typ( a_remove(i) ) = 0
        enddo
      endif
c
c  now search among -1 to sort angle types
c
      jcount = 0
      do i=1,nr_qma
        if( qma_typ( i ) .eq. -1)then
            jcount = jcount + 1
            qma_typ( i ) = jcount
c
            do j=i+1,nr_qma
               if( qma_typ( j ) .eq. -1)then
                  if( ica_qm( i ) .eq. ica_qm( j ) )then
                   qma_typ( j ) = jcount
                  endif
                endif
            enddo
        endif
      enddo
c
      nr_qma_typ = jcount
      if(nr_qma_typ .eq. 0)then
         write(*,*)'angle_types: WARNING: no angle to fit'
      endif
c
c
      write(*,*)' '
      write(*,*)'--------------'
      write(*,*)'Now, below are new list of angles'
      write(*,*)'zero in the last column means not fitted'
      write(*,*)' '
      write(*,*)'angle_ind; at_ind1, at_ind2, at_ind3'
      write(*,*)'at_name_1, at_name_2, at_name_3' 
      write(*,*)'at_type_1, at_type_2, at_type_3, unique_angle_type'
      write(*,*)' '
      do i=1,nr_qma
            write(*,"(4i6,6a6,i8)")i,ia_qm( i ),ja_qm( i ), ka_qm( i )
     *,     IGRAPH( ia_qm(i) ), IGRAPH( ja_qm(i) ), IGRAPH( ka_qm(i) )
     *,     ISYMBL( ia_qm(i) ), ISYMBL( ja_qm(i) ), ISYMBL( ka_qm(i) )
     *,      qma_typ( i )
      enddo
c
      write(*,*)' '
      write(*,"(a,i8)")'No of angle types: ',nr_qma_typ
c
      write(*,*)''
      write(*,*)'finding and classifying angles done !'
      write(*,*)'--------------------------------'
      write(*,*)''
      return
      end
c#######################################
c  counting and classifying 
c   unique dihedral types 
c  note that the dihedral types here are not 
c  the AMBER index into PK, PN, PHASE parameters
c  They include only atoms in QM zone
c  those outside (or dont need to be fitted)
c  will be set to zero.
c  they will be used as index of fitted variables,
c  not to get parameters from AMBER  
c qmd_typ = 0 means dont fit but calculate and put to the right side
c qmd_typ = 1, 2, 3 ... means 'will be fitted'
c##########################################
      subroutine dihedral_types(nat, nr_qm, max_nr_qmd
     *, grom2cpmd
     *, NATOM, ISYMBL, IGRAPH
     *, NPHIH, IPH,JPH,KPH,LPH,ICPH
     *, NPHIA, IP,JP,KP,LP,ICP 
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm
     *, nr_qmd_typ, qmd_typ
     *, icd_qm, pro_imp )
c
      implicit none
      integer nat, nr_qm, max_nr_qmd
     *, grom2cpmd
     *, NATOM
     *, NPHIH, IPH,JPH,KPH,LPH,ICPH
     *, NPHIA, IP,JP,KP,LP,ICP
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm
     *, nr_qmd_typ, qmd_typ
     *, icd_qm
     *, i, j, icount, jcount
     *, itmp, jtmp, ktmp, ltmp
     *, dih_qm_or_not, i_test, j_test
     *, nr_d_add, d_add
      character(3)pro_imp
      character(4) ISYMBL, IGRAPH
      character(30)filename
      dimension grom2cpmd( nat )
     *, ISYMBL( NATOM ), IGRAPH( NATOM )
     *, IPH(NPHIH), JPH(NPHIH), KPH(NPHIH), LPH(NPHIH), ICPH(NPHIH)
     *, IP(NPHIA), JP(NPHIA), KP(NPHIA), LP(NPHIA), ICP(NPHIA)
     *, id_qm(max_nr_qmd), jd_qm(max_nr_qmd)
     *, kd_qm(max_nr_qmd), ld_qm(max_nr_qmd)
     *, icd_qm(max_nr_qmd), qmd_typ(max_nr_qmd)
     *, pro_imp(max_nr_qmd)
c
     *, d_add( max_nr_qmd )
c
       write(*,*)' '
       write(*,*)'----------------------'
       write(*,*)'Counting no of diherals involved in the QM zone'
c
c  look for all the dihedrals (including H) with at least one atom is QM
c  the total number of these will be nr_qmd
c
c   then decide which dihedrals will be optimized qmd_typ .ne. -1
c
c    sort the optimized dihedrals into types 1, 2, 3 ... nr_qmd_typ
c
      if( NPHIH .le. 0 )then
        stop'dihedral_types: NPHIH .le. 0'
      endif
c
      icount = 0
      do i=1, NPHIH                ! DIHEDRALS_INC_HYDROGEN
         itmp = abs(IPH(i))/3 + 1
         jtmp = abs(JPH(i))/3 + 1
         ktmp = abs(KPH(i))/3 + 1
         ltmp = abs(LPH(i))/3 + 1
c
         i_test = dih_qm_or_not(itmp,jtmp,ktmp,ltmp
     *,    nat, nr_qm, grom2cpmd)
c
         if ( i_test .ne. 0 ) then
              icount = icount + 1
c
              id_qm( icount ) = itmp
              jd_qm( icount ) = jtmp
              kd_qm( icount ) = ktmp
              ld_qm( icount ) = ltmp
              icd_qm( icount ) = ICPH(i)
c
              qmd_typ(icount) = 0
c
              if( LPH(i) .ge. 0)then
                pro_imp( icount ) ='PRO'
              else
                pro_imp( icount ) ='IMP'
              endif
c
         endif
c
      enddo                        ! DIHEDRALS_INC_HYDROGEN
c
c
      if( NPHIA .le. 0 )then
        stop'dihedral_types: NPHIA .le. 0'
      endif
c
      do i=1,NPHIA                 ! DIHEDRALS_NOT_HYDROGEN
         itmp = abs(IP(i))/3 + 1
         jtmp = abs(JP(i))/3 + 1
         ktmp = abs(KP(i))/3 + 1
         ltmp = abs(LP(i))/3 + 1
c
         i_test = dih_qm_or_not(itmp,jtmp,ktmp,ltmp
     *,    nat, nr_qm, grom2cpmd)
c
          if ( i_test .ne. 0 ) then
              icount = icount + 1
c
              id_qm( icount ) = itmp
              jd_qm( icount ) = jtmp
              kd_qm( icount ) = ktmp
              ld_qm( icount ) = ltmp
              icd_qm( icount ) = ICP(i)
c
              qmd_typ(icount) = 0
c
              if( LP(i) .ge. 0)then
                pro_imp( icount ) ='PRO'
              else
                pro_imp( icount ) ='IMP'
              endif
c
          endif
c
      enddo                        ! DIHEDRALS_NOT_HYDROGEN
c
      if( icount .gt. max_nr_qmd) then
      stop'dihedral_types: icount .gt. max_nr_qmd '
      endif
c
       nr_qmd = icount
       if(nr_qmd .eq. 0)then
        stop'dihedral_types: nr_qmd .eq. 0'
       endif
c
c
c  print all the dihedrals in QM zone
c
      write(*,*)''
      write(*,*)'-----------------------------'
      write(*,*)'these are all the dihedrals in QM zone'
      write(*,*)' '
      write(*,*)'dihe_ind; at_ind1, at_ind2, at_ind3, at_ind4'
      write(*,*)'at_name1, at_name2, at_name3, at_name4' 
      write(*,*)'at_type1, at_type2, at_type3, at_type4'
      write(*,*)'PRO or IMP, index to amber, type'
      write(*,*)' '
      do i=1,nr_qmd
          write(*,"(5i6,8a6,a10,2i6)")i
     *,   id_qm( i ),jd_qm( i ), kd_qm( i ), ld_qm( i )
c 
     *,   IGRAPH( id_qm(i) ), IGRAPH( jd_qm(i) )
     *,   IGRAPH( kd_qm(i)), IGRAPH( ld_qm(i) )
c
     *,   ISYMBL( id_qm(i) ), ISYMBL( jd_qm(i) )
     *,   ISYMBL( kd_qm(i) ), ISYMBL( ld_qm(i) )
     *,   pro_imp(i), icd_qm(i),  qmd_typ(i)
      enddo
c
      write(*,*)''
      write(*,*)'-------------------'
      write(*,*)'by default, the dihedrals will not be fitted'
      write(*,*)'so the user has to explicitly select'
      write(*,*)'the ones he/she wants to fit'
      write(*,*)' '
      write(*,*)'if you dont want to fit any dihedral: enter 1'
      write(*,*)'if you want to include some dihedrals into the LS fit:'
      write(*,*)'enter 2'
      read(*,*)j_test
      nr_d_add = 0
      if(j_test.eq.2)then
          write(*,*)'-----------'
          write(*,*)'now, please prepare a file, in which you specify:'
          write(*,*)'NR_add'
          write(*,*)'dihe_ind_1'
          write(*,*)'dihe_ind_2'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'.'
          write(*,*)'dihe_ind_NR_add'
          write(*,*)' '
c
          write(*,*)'now, give me the name of that file:'
          write(*,*)'( 30 characters or less )'
          read(*,*)filename
          open(1,file=filename,status='old')
          rewind(1)
c
          read(1,*)nr_d_add
          write(*,"(i8)")nr_d_add
          if(nr_d_add .gt. max_nr_qmd) then
            stop'dihedral_types: nr_d_add .gt. max_nr_qmd'
          endif
c
          if(nr_d_add .gt. 0)then
            do i=1,nr_d_add
              read(1,*)d_add(i)
              write(*,"(i8)")d_add(i)
            enddo
          endif
c
          close(1)
      else
         if(j_test.ne.1)then
            stop'dihedral_types: dont know what to do'
         endif
      endif
c
c  now add more
c 
      if(nr_d_add .gt. 0)then
         do i=1,nr_d_add
         qmd_typ( d_add(i) ) = -1
         enddo
      endif
c
c
c  now search among -1 to sort bond types
c
      jcount = 0
      do i=1,nr_qmd
        if( qmd_typ( i ) .eq. -1)then
            jcount = jcount + 1
            qmd_typ( i ) = jcount
c
            do j=i+1,nr_qmd
               if( qmd_typ( j ) .eq. -1)then
                  if( icd_qm( i ) .eq. icd_qm( j ) )then
                   qmd_typ( j ) = jcount
                  endif
                endif
            enddo
        endif
      enddo
c
      nr_qmd_typ = jcount
c
      write(*,*)' '
      write(*,*)'--------------'
      write(*,*)'Now, below are new list of dihedrals'
      write(*,*)'zero in the last column means not fitted'
      write(*,*)' '
      write(*,*)'dihe_ind; at_ind1, at_ind2, at_ind3, at_ind4'
      write(*,*)'at_name1, at_name2, at_name3, at_name4' 
      write(*,*)'at_type1, at_type2, at_type3, at_type4'
      write(*,*)'PRO or IMP,index to amber, type'
      write(*,*)' '
      do i=1,nr_qmd
          write(*,"(5i6,8a6,a8,2i6)")i
     *,   id_qm( i ),jd_qm( i ), kd_qm( i ), ld_qm( i )
c 
     *,   IGRAPH( id_qm(i) ), IGRAPH( jd_qm(i) )
     *,   IGRAPH( kd_qm(i)), IGRAPH( ld_qm(i) )
c
     *,   ISYMBL( id_qm(i) ), ISYMBL( jd_qm(i) )
     *,   ISYMBL( kd_qm(i) ), ISYMBL( ld_qm(i) )
     *,   pro_imp(i), icd_qm(i), qmd_typ(i)
      enddo
c
      write(*,*)'Counting no of diherals involved in the QM zone DONE'
      write(*,*)''
c
      return
      end
c#################################
c   to see whether a bond is 
c   completely not QM:  0
c   partly QM:         -1
c   completely QM:      1
c###
c   Psedo_CAP - MM  : 0 (among -1 , if DUM found -> 0)
c   Psedo_CAP - QM  : -1 
c  ########################
c in case DUM is used : if DUM found - > 0
c##################################
      integer function bond_qm_or_not(ib,jb
     *, nat, nr_qm, grom2cpmd)
c
      implicit none
      integer ib,jb
     *, nat, nr_qm, grom2cpmd
     *, itmp,icount,i
      dimension grom2cpmd(nat)
     *,itmp(2)
c
      itmp(1) = grom2cpmd( ib )
      itmp(2) = grom2cpmd( jb )
c
      icount = 0
      do i=1,2
        if( itmp(i) .le. nr_qm) then
            icount = icount + 1
        endif
      enddo
c
      if(icount .eq. 0) bond_qm_or_not =  0
      if(icount .eq. 1) bond_qm_or_not = -1   
      if(icount .eq. 2) bond_qm_or_not =  1
c
      return
      end
c#################################
c   to see whether a angle is 
c   completely not QM:  0
c   partly QM:         -1
c   completely QM:      1
c##################################
      integer function ang_qm_or_not(it,jt,kt
     *, nat, nr_qm, grom2cpmd)
c
      implicit none
c
      integer it,jt,kt
     *, nat, nr_qm, grom2cpmd
     *, itmp,icount,i
      dimension grom2cpmd(nat)
     *,itmp(3)
c
      itmp(1) = grom2cpmd( it )
      itmp(2) = grom2cpmd( jt )
      itmp(3) = grom2cpmd( kt )
c
      icount = 0
      do i=1,3
       if( itmp(i) .le. nr_qm) then
            icount = icount + 1
        endif
      enddo
c
      if(icount .eq. 0) ang_qm_or_not =  0
      if(icount .eq. 1) ang_qm_or_not = -1
      if(icount .eq. 2) ang_qm_or_not = -1
      if(icount .eq. 3) ang_qm_or_not =  1
c
      return
      end
c#################################
c   to see whether a dihedral is 
c   completely not QM:  0
c   partly QM:         -1
c   completely QM:      1
c##################################
      integer function dih_qm_or_not(ip,jp,kp,lp
     *, nat, nr_qm, grom2cpmd)
c
      implicit none
c
      integer ip,jp,kp,lp
     *, nat, nr_qm, grom2cpmd
     *, itmp,icount,i
      dimension grom2cpmd(nat)
     *,itmp(4)
c
      itmp(1) = grom2cpmd( ip )
      itmp(2) = grom2cpmd( jp )
      itmp(3) = grom2cpmd( kp )
      itmp(4) = grom2cpmd( lp )
c
      icount = 0
      do i=1,4
       if( itmp(i) .le. nr_qm) then
            icount = icount + 1
        endif
      enddo
c
      if(icount .eq. 0) dih_qm_or_not =  0
      if(icount .eq. 1) dih_qm_or_not = -1
      if(icount .eq. 2) dih_qm_or_not = -1
      if(icount .eq. 3) dih_qm_or_not = -1
      if(icount .eq. 4) dih_qm_or_not =  1
c
      return
      end
