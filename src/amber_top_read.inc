c##################################
c  to read AMBER_TOP
c  this is a bit messy, need to improve
c##################################
      write(*,*)''
      write(*,*)'------------------'
      write(*,*)'reading AMBER_TOP'
      open(1,file='AMBER_TOP',status='old')
c
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG POINTERS")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( amb_pt(i_ab), i_ab = 1, amb_nr_pt )
      write(*,"(10I8)")( amb_pt(i_ab), i_ab = 1, amb_nr_pt )
      if(amb_pt(1) .ne. nat) then
      stop'read POINTERS: wrong nat!'
      endif
c

      NATOM = amb_pt (1)
      NTYPES = amb_pt (2)
      NBONH = amb_pt (3)
      MBONA = amb_pt (4)
      NTHETH = amb_pt (5)
      MTHETA = amb_pt (6)
      NPHIH = amb_pt (7)
      MPHIA = amb_pt (8)
      NHPARM =  amb_pt (9)
      NPARM = amb_pt (10)
      NNB =  amb_pt (11)
      NRES = amb_pt (12)
      NBONA = amb_pt (13)
      NTHETA = amb_pt (14)
      NPHIA = amb_pt (15)
      NUMBND = amb_pt (16)
      NUMANG = amb_pt (17)
      NPTRA =  amb_pt (18)
      NATYP = amb_pt (19)
      NPHB = amb_pt (20)
      IFPERT = amb_pt (21)
      NBPER = amb_pt (22)
      NGPER = amb_pt (23)
      NDPER = amb_pt (24)
      MBPER = amb_pt (25)
      MGPER = amb_pt (26)
      MDPER = amb_pt (27)
      IFBOX = amb_pt (28)
      NMXRS = amb_pt (29)
      IFCAP = amb_pt (30)
      NUMEXTRA = amb_pt (31)
c
      allocate( IAC (NATOM) )
      allocate( NUMEX (NATOM) )
      allocate( ICO (NTYPES*NTYPES) )
      allocate( IPRES (NRES) )
      allocate( IBH(NBONH), JBH(NBONH), ICBH(NBONH) )
      allocate( IB(NBONA), JB(NBONA), ICB(NBONA) )
      allocate( ITH(NTHETH), JTH(NTHETH), KTH(NTHETH), ICTH(NTHETH) )
      allocate( IT(NTHETA), JT(NTHETA), KT(NTHETA), ICT(NTHETA) )
      allocate( IPH(NPHIH), JPH(NPHIH), KPH(NPHIH), LPH(NPHIH)
     *         , ICPH(NPHIH) )
      allocate( IP(NPHIA), JP(NPHIA), KP(NPHIA), LP(NPHIA), ICP(NPHIA) )
      allocate( INB(NNB) )
c      allocate( NSP(NSPM) )
      allocate( IGRAPH(NATOM) )
      allocate( LBRES(NRES) )
      allocate( ISYMBL(NATOM) )
      allocate( CHARGE(NATOM) )
      allocate( AMASS(NATOM) )
      allocate( RK(NUMBND) )
      allocate( REQ(NUMBND) )
      allocate( TK(NUMANG) )
      allocate( TEQ(NUMANG) )
      allocate( PK(NPTRA) )
      allocate( PN(NPTRA) )
      allocate( PHASE(NPTRA) )
      allocate( SOLTY(NATYP) )
      allocate( CN1( NTYPES*(NTYPES+1)/2 ) )
      allocate( CN2( NTYPES*(NTYPES+1)/2 ) )
c
c ATOM_NAME 
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG ATOM_NAME")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(20a4)")( IGRAPH(i_ab), i_ab = 1,NATOM )
c      write(*,"(20a4)")( IGRAPH(i_ab), i_ab = 1,NATOM )
c CHARGE 
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG CHARGE")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( CHARGE(i_ab), i_ab = 1,NATOM )
c      write(*,"(5E16.8)")( CHARGE(i_ab), i_ab = 1,NATOM )
c MASS 
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG MASS")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( AMASS(i_ab), i_ab = 1,NATOM )
c      write(*,"(5E16.8)")( AMASS(i_ab), i_ab = 1,NATOM )
c ATOM_TYPE_INDEX 
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG ATOM_TYPE_INDEX")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( IAC(i_ab), i_ab = 1,NATOM )
c      write(*,"(10I8)")( IAC(i_ab), i_ab = 1,NATOM )
c NUMBER_EXCLUDED_ATOMS
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG NUMBER_EXCLUDED_ATOMS")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( NUMEX(i_ab), i_ab = 1,NATOM )
c      write(*,"(10I8)")( NUMEX(i_ab), i_ab = 1,NATOM )
c NONBONDED_PARM_INDEX
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG NONBONDED_PARM_INDEX")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( ICO(i_ab), i_ab = 1,NTYPES*NTYPES )
c      write(*,"(10I8)")( ICO(i_ab), i_ab = 1,NTYPES*NTYPES )
c RESIDUE_LABEL
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG RESIDUE_LABEL")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(20A4)")( LBRES(i_ab), i_ab = 1,NRES )
c      write(*,"(20A4)")( LBRES(i_ab), i_ab = 1,NRES )
c RESIDUE_POINTER
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG RESIDUE_POINTER")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( IPRES(i_ab), i_ab = 1,NRES )
c      write(*,"(10I8)")( IPRES(i_ab), i_ab = 1,NRES )
c BOND_FORCE_CONSTANT
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG BOND_FORCE_CONSTANT")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( RK(i_ab), i_ab = 1,NUMBND )
c      write(*,"(5E16.8)")( RK(i_ab), i_ab = 1,NUMBND )
c BOND_EQUIL_VALUE
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG BOND_EQUIL_VALUE")
          read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( REQ(i_ab), i_ab = 1,NUMBND )
c      write(*,"(5E16.8)")( REQ(i_ab), i_ab = 1,NUMBND )
c ANGLE_FORCE_CONSTANT
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG ANGLE_FORCE_CONSTANT")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( TK(i_ab), i_ab = 1,NUMANG )
c      write(*,"(5E16.8)")( TK(i_ab), i_ab = 1,NUMANG )
c ANGLE_EQUIL_VALUE
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG ANGLE_EQUIL_VALUE")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( TEQ(i_ab), i_ab = 1,NUMANG )
c      write(*,"(5E16.8)")( TEQ(i_ab), i_ab = 1,NUMANG )
c DIHEDRAL_FORCE_CONSTANT
      rewind(1)
      flag='a'
      do while (flag .ne."%FLAG DIHEDRAL_FORCE_CONSTANT")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( PK(i_ab), i_ab = 1,NPTRA )
c      write(*,"(5E16.8)")( PK(i_ab), i_ab = 1,NPTRA )
c DIHEDRAL_PERIODICITY
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG DIHEDRAL_PERIODICITY")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( PN(i_ab), i_ab = 1,NPTRA )
c      write(*,"(5E16.8)")( PN(i_ab), i_ab = 1,NPTRA )
c DIHEDRAL_PHASE
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG DIHEDRAL_PHASE")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( PHASE(i_ab), i_ab = 1,NPTRA )
c      write(*,"(5E16.8)")( PHASE(i_ab), i_ab = 1,NPTRA )
c LENNARD_JONES_ACOEF
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG LENNARD_JONES_ACOEF")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( CN1(i_ab), i_ab = 1,NTYPES*(NTYPES+1)/2 )
c      write(*,"(5E16.8)")( CN1(i_ab), i_ab = 1,NTYPES*(NTYPES+1)/2 )
c LENNARD_JONES_BCOEF
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG LENNARD_JONES_BCOEF")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(5E16.8)")( CN2(i_ab), i_ab = 1,NTYPES*(NTYPES+1)/2 )
c      write(*,"(5E16.8)")( CN2(i_ab), i_ab = 1,NTYPES*(NTYPES+1)/2 )
c BONDS_INC_HYDROGEN
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG BONDS_INC_HYDROGEN")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( IBH(i_ab),JBH(i_ab),ICBH(i_ab), i_ab = 1,NBONH )
c      write(*,"(10I8)")(IBH(i_ab),JBH(i_ab),ICBH(i_ab), i_ab = 1,NBONH)
c BONDS_WITHOUT_HYDROGEN
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG BONDS_WITHOUT_HYDROGEN")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")(IB(i_ab),JB(i_ab),ICB(i_ab), i_ab = 1,NBONA)
c      write(*,"(10I8)")(IB(i_ab),JB(i_ab),ICB(i_ab), i_ab=1,NBONA)
c ANGLES_INC_HYDROGEN
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG ANGLES_INC_HYDROGEN")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")(ITH(i_ab),JTH(i_ab),KTH(i_ab),ICTH(i_ab)
     *                  , i_ab=1,NTHETH)
c      write(*,"(10I8)")(ITH(i_ab),JTH(i_ab),KTH(i_ab),ICTH(i_ab)
c     *                    , i_ab=1,NTHETH)
c ANGLES_WITHOUT_HYDROGEN
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG ANGLES_WITHOUT_HYDROGEN")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( IT(i_ab),JT(i_ab),KT(i_ab),ICT(i_ab)
     *        , i_ab = 1,NTHETA )
c      write(*,"(10I8)")( IT(i_ab),JT(i_ab),KT(i_ab),ICT(i_ab)
c     *        , i_ab = 1,NTHETA )
c DIHEDRALS_INC_HYDROGEN
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG DIHEDRALS_INC_HYDROGEN")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( IPH(i_ab),JPH(i_ab),KPH(i_ab),LPH(i_ab)
     *                  ,ICPH(i_ab), i_ab=1,NPHIH )
c      write(*,"(10I8)")( IPH(i_ab),JPH(i_ab),KPH(i_ab),LPH(i_ab)
c     *                  ,ICPH(i_ab), i_ab=1,NPHIH )
c DIHEDRALS_WITHOUT_HYDROGEN
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG DIHEDRALS_WITHOUT_HYDROGEN")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( IP(i_ab),JP(i_ab),KP(i_ab),LP(i_ab)
     *                 ,ICP(i_ab), i_ab=1,NPHIA )
c      write(*,"(10I8)")( IP(i_ab),JP(i_ab),KP(i_ab),LP(i_ab)
c     *                  ,ICP(i_ab), i_ab=1,NPHIA )
c EXCLUDED_ATOMS_LIST
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG EXCLUDED_ATOMS_LIST")
         read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(10I8)")( INB(i_ab), i_ab = 1,NNB )
c      write(*,"(10I8)")( INB(i_ab), i_ab = 1,NNB )
c AMBER_ATOM_TYPE
      rewind(1)
      flag='a'
      do while (flag.ne."%FLAG AMBER_ATOM_TYPE")
        read(1,"(a50)")flag
      enddo
c
      write(*,"(a50)")flag
      read(1,"(a50)")flag
      read(1,"(20A4)")( ISYMBL(i_ab), i_ab=1,NATOM )
c      write(*,"(20A4)")( ISYMBL(i_ab), i_ab=1,NATOM )
c
      write(*,*)'reading AMBER_TOP done!'
      write(*,*)''
      close(1)
