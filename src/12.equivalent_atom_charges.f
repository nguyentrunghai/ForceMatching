c##################################
c to read equivalent 
c charges and set the 
c charge type
c#################################
      subroutine equivalent_charg(nat, nr_qm
     *, cpmd2grom, IGRAPH
     *, nr_char_tp
     *, charg_tp)
c
      implicit none
      integer nat, nr_qm
     *, cpmd2grom
     *, nr_char_tp, charg_tp
     *, i, j, k, icount, jcount
     *, j_test, itmp
     *, max_d
     *, nr_not_fit
     *, not_fit_at
     *, nr_equi_list, nr_at_in_list
     *, at_ind
     *, char_typ_tmp, grom_ind
      parameter(max_d=1000)
      character(30)filename,ctmp
      character(4) IGRAPH
      dimension cpmd2grom(nat)
     *, IGRAPH( nat )
     *, charg_tp(nr_qm)
     *, not_fit_at(max_d)
     *, nr_at_in_list(max_d)
     *, at_ind(max_d, max_d)
     *, char_typ_tmp( nr_qm )
     *, grom_ind( nr_qm )
c
      write(*,*)''
      write(*,*)'-------------------------------'
      write(*,*)'classifying QM atoms into charge types'
      write(*,*)'0, 1, 2, ... nr_char_tp'
      write(*,*)'(0 means not included in the least square fit)'
      write(*,*)'please provide me a file, in which'
      write(*,*)''
      write(*,*)'number of atoms not fitted'
      write(*,*)'atom atom atom ...'
      write(*,*)'#'
      write(*,*)'number of equivalent lists'
      write(*,*)'number of atoms in the 1st list'
      write(*,*)'atom atom atom ...'
      write(*,*)'number of atoms in the 2ndst list'
      write(*,*)'atom atom atom ...'
      write(*,*)'.'
      write(*,*)'.'
      write(*,*)'.'
      write(*,*)''
      write(*,*)'where atom is gromos index'
c
      write(*,*)'now give me the name of that file'
      write(*,*)'( 30 characters or less )'
      read(*,*)filename
      open(1,file=filename,status='old')
      rewind(1)
c
      write(*,*)''
      read(1,*)nr_not_fit
      write(*,"(i5)")nr_not_fit
      if((nr_not_fit .gt. 0) .and. (nr_not_fit .lt. max_d) )then
         read(1,*)(not_fit_at(i), i=1,nr_not_fit )
         write(*,"(10i5)")(not_fit_at(i), i=1,nr_not_fit )
      else
        if(nr_not_fit .ne. 0)then
           stop'equivalent_charg: wrong nr_not_fit'
        endif
      endif
c
      read(1,*)ctmp
      write(*,"(a)")ctmp
c
      read(1,*)nr_equi_list
      write(*,"(i5)")nr_equi_list
      if(nr_equi_list .gt. 0)then
         do i=1,nr_equi_list
              read(1,*)nr_at_in_list(i)
               write(*,"(i5)")nr_at_in_list(i)
              if((nr_at_in_list(i) .gt. 0) .and. 
     *                  (nr_at_in_list(i) .lt. max_d))then
c
                read(1,*)( at_ind(i,j), j=1,nr_at_in_list(i) )
                write(*,"(10I5)")( at_ind(i,j), j=1,nr_at_in_list(i) )
              else
                  stop'equivalent_charg: wrong nr_at_in_list(i)'
c
              endif
         enddo
      endif
      close(1)
c
      do j=1, nr_qm
         char_typ_tmp( j ) = -1
      enddo
c not fit  ####
      if( nr_not_fit .gt. 0 )then
        do i=1,nr_not_fit
           do j=1,nr_qm
c
           if( cpmd2grom(j) .eq. not_fit_at(i) )then
              if(char_typ_tmp( j ) .eq. -1)then
                  char_typ_tmp( j ) = 0
              endif
           endif
c
           enddo
        enddo
      endif
c
c fited atom read from equivalency file
      if( nr_equi_list .gt. 0)then
      do k=1,nr_equi_list              ! number of lines loop
         do i=1,nr_at_in_list(k)         ! atom loop
            do j=1,nr_qm
c
               if( cpmd2grom(j) .eq. at_ind( k, i ) )then
                  if(char_typ_tmp( j ) .eq. -1)then
                      char_typ_tmp( j ) = k
                   endif
               endif
c
            enddo
         enddo                           ! atom loop
      enddo                            ! number of lines loop
      endif
c
c set the rest
      icount = nr_equi_list
      do j=1,nr_qm
        if(char_typ_tmp( j ) .eq. -1)then
           icount = icount + 1
           char_typ_tmp( j ) = icount
        endif
      enddo
c
      nr_char_tp = icount
      if(nr_char_tp .le. 0)then
      stop'equivalent_charg: nr_char_tp .le. 0'
      endif
c
c   change to gromos order
c
      jcount = 0
      do i=1,nat                ! gromos loop
        do j=1, nr_qm           ! cpmd loop
c
            if(cpmd2grom(j) .eq. i)then
                jcount = jcount + 1
                charg_tp( jcount ) = char_typ_tmp( j )
                grom_ind(jcount) = i
            endif
c
        enddo                   ! cpmd loop
      enddo                     ! gromos loop
c######
c write
      write(*,*)''
      write(*,*)'below is the list of qm atom and charge types'
      write(*,*)''
      do j=1,nr_qm
      write(*,"(2i5,a6,i5)")j,grom_ind(j)
     *,IGRAPH( grom_ind(j) ),charg_tp( j )
      enddo
      write(*,"(i8)")nr_char_tp
      write(*,*)''
c
      write(*,*)'if you satisfy with the above list'
      write(*,*)'enter 1'
      write(*,*)'if you want to input an explicit charge type file'
      write(*,*)'enter 2'
      write(*,*)'to stop, enter any other number'
      read(*,*)j_test
c
      if(j_test .eq. 2)then
         write(*,*)'please provide me a file, in which'
         write(*,*)''
         write(*,*)'number-of-qm-atoms   number-of-charge-types'
         write(*,*)'gromos-index   charge-type'
         write(*,*)'gromos-index   charge-type'
         write(*,*)'gromos-index   charge-type'
         write(*,*)'.'
         write(*,*)'.'
         write(*,*)'.'
c
         write(*,*)''
         write(*,*)'charge-type should be 0, 1, 2
     * ... number-of-charge-types'
         write(*,*)'now give me the name of that file'
         write(*,*)'( 30 characters or less )'
         read(*,*)filename
         open(1,file=filename,status='old')
         rewind(1)
c
          read(1,*)itmp,nr_char_tp
         if(itmp .ne. nr_qm)then
         stop'read_charg_typ: number of QM is not correct'
         endif
c
         do i=1,nr_qm
         read(1,*)itmp, charg_tp(i)
         enddo
         close(1)
c
         write(*,*)''
         write(*,*)'below is new charge types, please check'
         write(*,*)''
         do j=1,nr_qm
         write(*,"(2i5,a6,i5)")j,grom_ind(j)
     *,       IGRAPH( grom_ind(j) ),charg_tp( j )
         enddo
         write(*,"(i8)")nr_char_tp
         write(*,*)''
c
      write(*,*)''
      write(*,*)'enter 1 to continune'
      write(*,*)'or any other number of stop'
      read(*,*)j_test
      if(j_test .ne. 1)then
      stop'equivalent_charg: stopped'
      endif
c
      else
           if(j_test .ne. 1)then
             stop'equivalent_charg: dont know what to do'
           endif
c
      endif
c
      write(*,*)'classifying QM atoms into charge types done!'
      write(*,*)''
c
      return
      end
c##################################
c to estimate total charge of
c the fitted atoms, from force field
c##################################
      subroutine esti_tota_charg(nat, nr_qm
     *, grom2cpmd, CHARGE
     *, charg_tp
     *, tot_qm_charg )
c
      implicit none
      integer nat, nr_qm
     *, grom2cpmd
     *, charg_tp
     *, i, icount, j_test
      double precision CHARGE, tot_qm_charg
     *, dtmp
      dimension grom2cpmd(nat)
     *, CHARGE(nat)
     *, charg_tp(nr_qm)
c
      write(*,*)''
      write(*,*)'----------------'
      write(*,*)'estimating total charge'
c
      dtmp = 0.0d0
      icount = 0
      do i=1,nat
        if(grom2cpmd(i) .le. nr_qm)then
            icount = icount + 1
            if( charg_tp(icount) .ne. 0 )then
             dtmp = dtmp + ( CHARGE(i)/18.22230d0 )
            endif
        endif
      enddo
      tot_qm_charg = dtmp
c
      write(*,"(a,f15.7)")'Q = ',tot_qm_charg
      write(*,*)''
      write(*,*)'enter 1 to continue'
      write(*,*)'enter 2 to change the above value'
      read(*,*)j_test
      if(j_test .eq. 2)then
        write(*,*)'so you want to enter anather number'
        write(*,*)'Q = '
        read(*,*)tot_qm_charg
        write(*,"(a,f15.7)")'new value',tot_qm_charg
      else
          if(j_test .ne. 1)then
             stop'esti_tota_charg: dont know what to do'
          endif
      endif
c
       write(*,*)'estimating total charge done!'
       write(*,*)''
c
      return
      end
c##################################
c  read charg_type.dat
c################################
      subroutine read_charg_typ(nr_qm, nr_char_tp
     *, charg_tp)
c
      implicit none
      integer nr_qm, nr_char_tp
     *, itmp,i, charg_tp
      dimension charg_tp(nr_qm)
c
      write(*,*)''
      write(*,*)'reading charg_type.dat'
      open(1,file='charg_type.dat',status='old')
      rewind(1)
      read(1,*)itmp,nr_char_tp
      if(itmp .ne. nr_qm)then
      stop'read_charg_typ: number of QM is not correct'
      endif
c
      do i=1,nr_qm
      read(1,*)itmp, charg_tp(i)
      enddo
      close(1)
c
      write(*,*)'reading charg_type.dat done !'
      write(*,*)''      
c
      return
      end

