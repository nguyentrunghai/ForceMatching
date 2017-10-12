c################################
c  store ab initio bonded forces
c################################
      subroutine bonded_force_store(nr_qm, nr_frame
     *,felec_x, felec_y, felec_z
     *,fvdw_x, fvdw_y, fvdw_z
     *,ref_fx, ref_fy, ref_fz
     *,bonded_fx, bonded_fy, bonded_fz)
c
      implicit none
      integer nr_qm, nr_frame, i, j
      double precision felec_x, felec_y, felec_z
     *, fvdw_x, fvdw_y, fvdw_z
     *, ref_fx, ref_fy, ref_fz
     *, bonded_fx, bonded_fy, bonded_fz
     *, const
      dimension felec_x(nr_frame, nr_qm), felec_y(nr_frame, nr_qm)
     *, felec_z(nr_frame, nr_qm)
c
     *, fvdw_x(nr_frame, nr_qm), fvdw_y(nr_frame, nr_qm)
     *, fvdw_z(nr_frame, nr_qm)
c
     *, ref_fx(nr_frame, nr_qm), ref_fy(nr_frame, nr_qm)
     *, ref_fz(nr_frame, nr_qm)
c
     *, bonded_fx(nr_frame, nr_qm), bonded_fy(nr_frame, nr_qm)
     *, bonded_fz(nr_frame, nr_qm)
c
      write(*,*)''
      write(*,*)'--------------------------------'
      write(*,*)'storing ab initio bonded forces'
      const = 1185.82105d0
c
      do i=1, nr_frame           !  frame loop
        do j=1, nr_qm            !  qm atom loop
c
         bonded_fx(i, j) = (ref_fx(i, j)*const) - felec_x(i, j) 
     *                                          - fvdw_x(i, j)
         bonded_fy(i, j) = (ref_fy(i, j)*const) - felec_y(i, j) 
     *                                          - fvdw_y(i, j)
         bonded_fz(i, j) = (ref_fz(i, j)*const) - felec_z(i, j) 
     *                                          - fvdw_z(i, j)
c
        enddo                    !  qm atom loop
      enddo                      !  frame loop
c
      write(*,*)'storing ab initio bonded forces done!'
      write(*,*)'in AMBER unit'
      write(*,*)' the ordering is CPMD'
      write(*,*)''
      return
      end
c#############################
c  read FM_REF_FORCES file
c  this file contatins ab initio forces
c  from CPMD
c##############################
      subroutine read_fm_forces(nr_qm, nr_frame, frame_ind
     *,ref_fx, ref_fy, ref_fz)
c
      implicit none
      integer nr_qm, nr_frame, frame_ind
     *, i, j, itmp, itmp1,itmp2
      double precision ref_fx, ref_fy, ref_fz
     *, dtmp
      dimension frame_ind(nr_frame)
     *, ref_fx( nr_frame, nr_qm ), ref_fy( nr_frame, nr_qm )
     *, ref_fz( nr_frame, nr_qm )
c
      write(*,*)''
      write(*,*)'------------------------------'
      write(*,*)'reading FM_REF_FORCES file'
      open(1,file='FM_REF_FORCES',status='old')
      rewind(1)
c
      do i=1, nr_frame       !   frame loop
        read(1,*)itmp1,itmp2
c
        if(itmp1 .ne. frame_ind(i) )then
          write(*,"(a,2I8)")'frame ',i,itmp1
          stop'read_fm_forces: wrong frame index'
        endif

        if(itmp2 .ne. nr_qm)then
           write(*,"(a,2I8)")'frame ',i,itmp1
           stop'read_fm_forces: wrong No. QM atoms'
        endif
c
          do j=1,nr_qm
          read(1,*)itmp,dtmp,dtmp,dtmp
     *,  ref_fx( i, j ), ref_fy( i, j ), ref_fz( i, j )
c
          enddo
c
      enddo                  !   frame loop
      close(1)
c
      write(*,*)'reading FM_REF_FORCES file done!'
      write(*,*)'the unit of forces stored in FM_REF_FORCES is au' 
      write(*,*)' the ordering is CPMD'
      write(*,*)''
      return
      end
c##############################
c  calculate vdW forces
c  vdW parameters read from Amber topology
c#################################
      subroutine vdw_forces(nat, nr_qm, nr_frame
     *, cpmd2grom
     *, NATOM, NTYPES, IAC, ICO, CN1, CN2
     *, excl_fact
     *, trx, try, trz 
     *, fvdw_x, fvdw_y, fvdw_z)
c
      implicit none
      integer nat, nr_qm, nr_frame, cpmd2grom
     *, NATOM, NTYPES, IAC, ICO
     *, i, j, k,  g_ind, vdw_ind
      double precision CN1, CN2
     *, excl_fact, trx, try, trz
     *,  fvdw_x, fvdw_y, fvdw_z
     *, dx, dy, dz, d2, d, n_x, n_y, n_z
     *, a12, b6, fcom
      dimension cpmd2grom(nat)
     *, IAC (NATOM), ICO (NTYPES*NTYPES)
     *, CN1( NTYPES*(NTYPES+1)/2 ), CN2( NTYPES*(NTYPES+1)/2 )
     *, excl_fact( nat, nat )
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
     *, fvdw_x(nr_frame, nr_qm), fvdw_y(nr_frame, nr_qm)
     *, fvdw_z(nr_frame, nr_qm)
c
      write(*,*)''
      write(*,*)'--------------------------------'
      write(*,*)'calculating van der Waals forces'
      if(nat .ne. NATOM)then
      stop'vdw_forces: wrong nat or NATOM'
      endif
c
      do i=1, nr_frame
        do j=1,nr_qm
          fvdw_x(i, j) = 0.0d0
          fvdw_y(i, j) = 0.0d0
          fvdw_z(i, j) = 0.0d0
        enddo
      enddo
c
c###
      do i=1, nr_frame            !  frame loop
         do j=1, nr_qm            ! qm lopp
            g_ind = cpmd2grom( j )
c
            do k=1,nat            !  all atom loop
            vdw_ind = ICO( NTYPES*( IAC(g_ind) -1 ) + IAC(k) )
            if(vdw_ind .gt. 0)then  !  vdw_ind must positive
            a12 = CN1( vdw_ind )
            b6  = CN2( vdw_ind )
c
              dx = trx(i,g_ind) - trx(i,k)
              dy = try(i,g_ind) - try(i,k)
              dz = trz(i,g_ind) - trz(i,k)
              d2 = (dx*dx) + (dy*dy) + (dz*dz)
              d = dsqrt(d2)
              if(d2 .ne. 0.0d0)then          ! d2 must be non-zero
                n_x = dx/d
                n_y = dy/d
                n_z = dz/d
                fcom = ( (12.0d0*a12)/(d**13) - (6.0d0*b6)/(d**7) )
                fcom = fcom * excl_fact( g_ind ,k )
c
                fvdw_x(i, j) = fvdw_x(i, j) + ( fcom * n_x )
                fvdw_y(i, j) = fvdw_y(i, j) + ( fcom * n_y )
                fvdw_z(i, j) = fvdw_z(i, j) + ( fcom * n_z )
              endif                          ! d2 must be non-zero
c
            endif                   !  vdw_ind must positive
            enddo                 !  all atom loop
         enddo                    ! qm loop
      enddo                       !  frame lopp
c
      write(*,*)'calculating van der Waals forces, done!'
      write(*,*)'vdW forces have been stored with 
     * the CMPD ordering, they are in AMBER unit'
      write(*,*)''
c
      return
      end
C##############################
c  calculating electrostatic forces
c  using the newly fitted point charges 
C##############################
      subroutine elect_forces(nat, nr_qm, nr_frame
     *, cpmd2grom
     *, charge_activ, excl_fact
     *, trx, try, trz
     *, felec_x, felec_y, felec_z)
c
      implicit none
      integer nat, nr_qm, nr_frame, cpmd2grom
     *, i, j, k, g_ind
      double precision  charge_activ, excl_fact
     *, trx, try, trz
     *, felec_x, felec_y, felec_z
     *, charge, dx ,dy, dz, d2, d
     *, n_x, n_y, n_z, fcom
      dimension cpmd2grom(nat)
     *, charge_activ (nat), excl_fact( nat, nat )
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
     *, felec_x(nr_frame, nr_qm), felec_y(nr_frame, nr_qm)
     *, felec_z(nr_frame, nr_qm)
     *, charge(nat)
c
      write(*,*)''
      write(*,*)'-------------------------------'
      write(*,*)'calculating elctrostatic forces'
c
      do i=1,nat
      charge( i ) = charge_activ( i )*18.22230d0
      enddo
c
      do i=1,nr_frame
        do j=1,nr_qm
            felec_x(i, j) = 0.0d0
            felec_y(i, j) = 0.0d0
            felec_z(i, j) = 0.0d0
        enddo
      enddo
c###############
      do i=1,nr_frame             ! frame loop
c
         do j=1,nr_qm             ! qm loop
           g_ind = cpmd2grom( j )
c
           do k=1,nat             ! all atom loop
c
            dx = trx(i,g_ind) - trx(i,k)
            dy = try(i,g_ind) - try(i,k)
            dz = trz(i,g_ind) - trz(i,k)
            d2 = (dx*dx) + (dy*dy) + (dz*dz)
            d = dsqrt(d2)
            if(d2 .ne. 0.0d0)then
              n_x = dx/d
              n_y = dy/d
              n_z = dz/d
              fcom = charge(g_ind)*charge(k)*excl_fact( g_ind ,k )/d2
c
              felec_x(i, j) = felec_x(i, j) + ( fcom * n_x )
              felec_y(i, j) = felec_y(i, j) + ( fcom * n_y )
              felec_z(i, j) = felec_z(i, j) + ( fcom * n_z )
            endif
c 
           enddo                  ! all atom loop
c
         enddo                    ! qm loop
      enddo                       ! frame loop
c
      write(*,*)'calculating elctrostatic forces, done!'
      write(*,*)'electrostatice forces have been stored with 
     * the CMPD ordering, they are in AMBER unit'
      write(*,*)''
c
      return
      end
c##################################
c  build the exclude factor matrix
c  for electrostatic forces calculation
c#################################
      subroutine nonb_exclude_fact(IBH,JBH,NBONH,   IB,JB,NBONA
     *, ITH,JTH,KTH,NTHETH,   IT,JT,KT,NTHETA
     *, IPH,JPH,KPH,LPH,NPHIH,   IP,JP,KP,LP,NPHIA 
     *, nat, excl_fact)
      implicit none
      integer IBH,JBH,NBONH,   IB,JB,NBONA
     *, ITH,JTH,KTH,NTHETH,   IT,JT,KT,NTHETA
     *, IPH,JPH,KPH,LPH,NPHIH,   IP,JP,KP,LP,NPHIA
     *, nat,i,j, itmp, jtmp
      double precision excl_fact, nonb_14_fac
      dimension IBH(NBONH),JBH(NBONH)
     *, IB(NBONA),JB(NBONA)
     *, ITH(NTHETH),JTH(NTHETH),KTH(NTHETH)
     *, IT(NTHETA),JT(NTHETA),KT(NTHETA)
     *, IPH(NPHIH),JPH(NPHIH),KPH(NPHIH),LPH(NPHIH)
     *, IP(NPHIA),JP(NPHIA),KP(NPHIA),LP(NPHIA)
     *, excl_fact(nat,nat)
c
      nonb_14_fac = 0.830d0
c
      write(*,*)''
      write(*,*)'--------------------------------'
      write(*,*)'generating nonb excluding matrix'
      write(*,"(f15.8)")0.830d0
c
      do i=1,nat
         do j=1,nat
c
         excl_fact(i,j) = 1.0d0
c
         enddo
      enddo
c
      do i=1,nat
      excl_fact(i,i) = 0.0d0
      enddo
c bonds
      do i=1,NBONH
      itmp = int( abs( IBH(i) )/3 ) + 1
      jtmp = int( abs( JBH(i) )/3 ) + 1
      excl_fact(itmp, jtmp) = 0.0d0
      excl_fact(jtmp, itmp) = 0.0d0
      enddo
      do i=1,NBONA
      itmp = int( abs( IB(i) )/3 ) + 1
      jtmp = int( abs( JB(i) )/3 ) + 1
      excl_fact(itmp, jtmp) = 0.0d0
      excl_fact(jtmp, itmp) = 0.0d0
      enddo
c angles
      do i=1,NTHETH
      itmp = int( abs( ITH(i) )/3 ) + 1
      jtmp = int( abs( KTH(i) )/3 ) + 1
      excl_fact(itmp, jtmp) = 0.0d0
      excl_fact(jtmp, itmp) = 0.0d0
      enddo
      do i=1,NTHETA
      itmp = int( abs( IT(i) )/3 ) + 1
      jtmp = int( abs( KT(i) )/3 ) + 1
      excl_fact(itmp, jtmp) = 0.0d0
      excl_fact(jtmp, itmp) = 0.0d0
      enddo
c dihedrals
      do i=1,NPHIH
      itmp = int( abs( IPH(i) )/3 ) + 1
      jtmp = int( abs( LPH(i) )/3 ) + 1
      excl_fact(itmp, jtmp) = excl_fact(itmp, jtmp)*nonb_14_fac
      excl_fact(jtmp, itmp) = excl_fact(jtmp, itmp)*nonb_14_fac
      enddo
      do i=1,NPHIA
      itmp = int( abs( IP(i) )/3 ) + 1
      jtmp = int( abs( LP(i) )/3 ) + 1
      excl_fact(itmp, jtmp) = excl_fact(itmp, jtmp)*nonb_14_fac
      excl_fact(jtmp, itmp) = excl_fact(jtmp, itmp)*nonb_14_fac
      enddo
c
      write(*,*)'generating nonb excluding matrix done!'
      write(*,*)''
      return
      end
c#####################################
c  read QM/MM trajectory
c#####################################
      subroutine read_traje(nr_frame, frame_ind, nat
     *, cpmd2grom
     *, trx, try, trz)
c
      implicit none 
      integer nr_frame, frame_ind, nat
     *,itmp, no_new_dat,new_dat, max_d
     *, nln,nfr, i,j,check,menb_check, icount
     *, cpmd2grom, cent_ind
     *, j_test
      parameter(max_d = 100000)
      character(10)ctmp
      double precision  trx, try, trz
     *, lx,ly,lz
     *, x, y, z, dtmp, const
     *, cent_x, cent_y, cent_z
     *, xtmp, ytmp, ztmp, dx, dy, dz
      dimension frame_ind (nr_frame)
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
     *, new_dat(max_d)
     *, x(nat) ,y(nat), z(nat)
     *, cpmd2grom(nat)
c
      write(*,*)''
      write(*,*)'--------------------------------------------'
      write(*,*)'reading TRAJECTORY'
      write(*,*)''
      write(*,*)'to wrap PBC, please tell me the size of the box'
      write(*,*)'this is the box of the total QM/MM system'
      write(*,*)'not just QM'
      write(*,*)'you can find it in gromos.inp'
      write(*,*)'now enter lx ly lz in ANGSTROM'
      write(*,*)''
c
      write(*,*)'lx = '
      read(*,*)lx
c
      write(*,*)'ly = '
      read(*,*)ly
c
      write(*,*)'lz = '
      read(*,*)lz
c
      write(*,*)''
      write(*,*)'here is what you input (in ANGSTROM)'
      write(*,"(3f15.8)")lx, ly, lz
      write(*,*)''
      write(*,*)'enter 1 to continue'
      write(*,*)'or any other number to stop'
      read(*,*)j_test
      if(j_test .ne. 1)then
        stop
      endif
c
      write(*,*)''
      write(*,*)'the reading is time consuming 
     * if you have  a big TRAJECTORY file'
      write(*,*)''
c
      call no_atoms(itmp)
      if(itmp .ne. nat) then
      stop'read_traje: wrong nat'
      endif
c
      call no_lines(nln, no_new_dat, new_dat, max_d, nat)
c      write(*,"('No of lines in TRAJECTORY:',i10)")nln
      nfr = nln
      write(*,"('No of frames in TRAJECTORY:',i10)")nfr
      write(*,"('But I only read ',i10,a)")nr_frame,' frames'
c
      open(1,file='TRAJECTORY',status='old')
      rewind(1)
c
      write(*,*)'reading coordinates'
c
      const=0.529177208590d0
c
      icount=1
      do i=1,nfr      ! frame loop
      if(icount .le. nr_frame) then
          check = menb_check(i,no_new_dat,new_dat, max_d)
          if(check .eq. 1)then
             read(1,*)ctmp
             write(*,*)ctmp
          endif
c
          do j=1,nat
            read(1,*)itmp,x(j),y(j),z(j),dtmp,dtmp,dtmp
          enddo
c
          if(itmp .eq. frame_ind(icount) )then
            write(*,"(a,i8)")'loading frame ',itmp
             do j=1,nat
               trx(icount,j)=x(j)*const
               try(icount,j)=y(j)*const
               trz(icount,j)=z(j)*const
             enddo
             icount=icount+1    
           endif
      endif
      enddo           ! frame loop
c
      if( icount.ne.(nr_frame+1) )then
      stop'read_traje: wrong icount'
      endif
      close(1)
c
      write(*,*)'Trajecroy has been stored in angstrom'
c
      write(*,*)'PBC wraping'
c
      cent_ind = cpmd2grom(1)
c
      do i=1,nr_frame      ! frame loop
         cent_x = trx( i, cent_ind )
         cent_y = try( i, cent_ind )
         cent_z = trz( i, cent_ind )
c
         do j=1,nat        ! atom loop
c##########################################
c  x
            xtmp = trx( i, j )
            dx =  xtmp - cent_x
c
 203        if( dx .lt. (-lx/2.0d0) )then
               xtmp = xtmp + lx
               dx = xtmp - cent_x
               goto 203
            endif
c
 204        if( dx .gt. (lx/2.0d0) )then
            xtmp = xtmp - lx
            dx= xtmp - cent_x
            goto 204
            endif
            trx( i, j ) = xtmp
c################################
c  y
            ytmp = try( i, j )
            dy =  ytmp - cent_y
c
 205        if( dy .lt. (-ly/2.0d0) )then
               ytmp = ytmp + ly
               dy = ytmp - cent_y
               goto 205
            endif
c
 206        if( dy .gt. (ly/2.0d0) )then
            ytmp = ytmp - ly
            dy = ytmp - cent_y
            goto 206
            endif
            try( i, j ) = ytmp
c####################################
c  z
            ztmp = trz( i, j )
            dz =  ztmp - cent_z
c
 207        if( dz .lt. (-lz/2.0d0) )then
               ztmp = ztmp + lz
               dz = ztmp - cent_z
               goto 207
            endif
c
 208        if( dz .gt. (lz/2.0d0) )then
            ztmp = ztmp - lz
            dz = ztmp - cent_z
            goto 208
            endif
            trz( i, j ) = ztmp
c####################################
c
         enddo             ! atom loop
      enddo                ! frame loop
c
      write(*,*)'reading TRAJECTORY done!'
      return
      end
c###########################################
c  read no of atoms in file TRAJECTRORY
c##########################################
      subroutine no_atoms(nat)
c
      implicit none
c
      integer nat
     *,icount, itmp1,itmp2
c
      write(*,*)'reading number of atoms in TRAJECTORY'
      open(1,file='TRAJECTORY',status='old')
      rewind(1)
c
      read(1,*)itmp1
      icount=1
 200  continue
      read(1,*)itmp2
      if(itmp2.eq.itmp1)then
      icount=icount+1
      itmp1=itmp2
      goto 200
      endif
c
      close (1)
c
      nat=icount
c
      return
      end
c###########################
c   read number of frames in TRAJECTORY
c###########################
      subroutine no_lines(nfr, no_new_dat, new_dat, max_d
     *,nat)
c
      implicit none
c
      integer nfr
     *, itmp
     *, io, new_dat
     *, no_new_dat
     *, i, nat, max_d
      character(10)ctmp
c
      double precision  icount
c
      dimension new_dat(max_d)
c
      write(*,*)'reading number of lines in TRAJECTORY'
      open(1,file='TRAJECTORY',status='old')
      rewind(1)
c
      icount=0.0d0
      itmp=0
 201  continue
      read(1,*,IOSTAT=io)ctmp
      if(ctmp .eq. "<<<<<<")then
      itmp=itmp+1
      write(*,"(i8,xxx,a10)")itmp,ctmp
c
      new_dat(itmp)=int(icount)/nat+1
c
      endif
c
      if( io.eq.0 )then
      if (ctmp .ne. "<<<<<<") then
      icount = icount + 1.0d0
      endif
      goto 201
      endif
c
      nfr=int(icount)/nat
c
       no_new_dat = itmp
c
       do i=1,itmp
       write(*,"(i8)")new_dat(i)
       enddo
c
       close(1)
c
      return
      end
c#############################
c   check member of an array
c#############################
      integer function menb_check(ind,no_new_dat,new_dat, max_d)
c
      implicit none
c
      integer ind,no_new_dat,new_dat
     *,i, max_d
c
      dimension new_dat(max_d)
c
      menb_check=0
c
      do i=1,no_new_dat
      if(new_dat(i).eq.ind)then
        menb_check=1
      endif
      enddo
c
      return
      end
