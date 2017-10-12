c
      program main_FM
c
      implicit none
c
      integer nr_frame, nr_qm, max_nr_mm
     *, nat
     *, i, j, j_test
     *, nr_char_tp
     *, nr_li_pot, nr_li_fd, nr_li_hirsh
     *, tot_matr_size
     *, max_nr_qmb, nr_qmb
     *, nr_qmb_typ
     *, max_nr_qma, nr_qma
     *, nr_qma_typ
     *, max_nr_qmd, nr_qmd
     *, nr_qmd_typ
     *, nr_bonded_row, nr_bonded_col
c
      integer NATOM, NTYPES, NBONH, MBONA, NTHETH, MTHETA
     *, i_ab
     *, NPHIH, MPHIA, NHPARM, NPARM, NNB, NRES
     *, NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA
     *, NATYP, NPHB, IFPERT, NBPER, NGPER, NDPER
     *, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP
     *, NUMEXTRA, NCOPY
     *, amb_nr_pt
c
      parameter(
     *  max_nr_mm = 5000
     *, max_nr_qmb = 5000
     *, max_nr_qma = 5000
     *, max_nr_qmd = 5000 )
      parameter(amb_nr_pt = 31)
c#############
      integer, dimension (:), allocatable ::
     *  IAC        ! FLAG ATOM_TYPE_INDEX (NATOM) FORMAT(1OI8)
     *, NUMEX      ! FLAG NUMBER_EXCLUDED_ATOMS (NATOM) FORMAT(10I8)
     *, ICO        ! FLAG NONBONDED_PARM_INDEX  (NTYPES*NTYPES) FORMAT(10I8)
     *, IPRES      ! FLAG RESIDUE_POINTER  (NRES) FORMAT(10I8) 
     *, IBH,JBH,ICBH ! FLAG BONDS_INC_HYDROGEN ( NBONH ) FORMAT(10I8)
     *, IB,JB,ICB  ! FLAG BONDS_WITHOUT_HYDROGEN (NBONA) FORMAT(10I8)
     *, ITH,JTH,KTH,ICTH ! FLAG ANGLES_INC_HYDROGEN (NTHETH) FORMAT(10I8)
     *, IT,JT,KT,ICT ! FLAG ANGLES_WITHOUT_HYDROGEN (NTHETA) FORMAT(10I8)
     *, IPH,JPH,KPH,LPH,ICPH ! FLAG DIHEDRALS_INC_HYDROGEN (NPHIH) FORMAT(10I8)
     *, IP,JP,KP,LP,ICP ! FLAG DIHEDRALS_WITHOUT_HYDROGEN (NPHIA) FORMAT(10I8)
     *, INB        ! FLAG EXCLUDED_ATOMS_LIST (NNB) FORMAT(10I8) 
     *, NSP        ! FLAG ATOMS_PER_MOLECULE (NSPM) FORMAT(10I8)
c#############
      integer, dimension (:), allocatable :: nr_mm
      integer, dimension (:), allocatable :: frame_ind
      integer, dimension (:), allocatable :: cpmd2grom
      integer, dimension (:), allocatable :: grom2cpmd
      integer, dimension (:), allocatable :: charg_tp
c
      integer, dimension (:), allocatable :: ib_qm
      integer, dimension (:), allocatable :: jb_qm
      integer, dimension (:), allocatable :: qmb_typ
      integer, dimension (:), allocatable :: icb_qm
c
      integer, dimension (:), allocatable :: ia_qm
      integer, dimension (:), allocatable :: ja_qm
      integer, dimension (:), allocatable :: ka_qm
      integer, dimension (:), allocatable :: qma_typ
      integer, dimension (:), allocatable :: ica_qm
c
      integer, dimension (:), allocatable :: id_qm
      integer, dimension (:), allocatable :: jd_qm
      integer, dimension (:), allocatable :: kd_qm
      integer, dimension (:), allocatable :: ld_qm
      integer, dimension (:), allocatable :: qmd_typ
      integer, dimension (:), allocatable :: icd_qm
c#####
      integer, dimension (:), allocatable :: amb_pt
c#####
      character(50) flag
c 
      character(4), dimension (:), allocatable ::
     *  IGRAPH     ! FLAG ATOM_NAME (NATOM) FORMAT(20a4)
     *, LBRES      ! FLAG RESIDUE_LABEL (NRES) FORMAT(20A4)
     *, ISYMBL     ! FLAG AMBER_ATOM_TYPE, (NATOM) FORMAT(20A4)
c#########
      character(3), dimension (:), allocatable :: pro_imp
c
      double precision 
     *  w_v, w_e, w_h, w_q, tot_qm_charg
     *, dtmp
c
      double precision, dimension (:,:), allocatable :: tr_qm_x
      double precision, dimension (:,:), allocatable :: tr_qm_y
      double precision, dimension (:,:), allocatable :: tr_qm_z
c
      double precision, dimension (:,:), allocatable :: tr_mm_x
      double precision, dimension (:,:), allocatable :: tr_mm_y
      double precision, dimension (:,:), allocatable :: tr_mm_z
c
      double precision, dimension (:,:), allocatable :: tr_ele_pot
c
      double precision, dimension (:,:), allocatable :: tr_ele_fd_x
      double precision, dimension (:,:), allocatable :: tr_ele_fd_y
      double precision, dimension (:,:), allocatable :: tr_ele_fd_z
c
      double precision, dimension (:,:), allocatable :: tr_chj
c
      double precision, dimension (:), allocatable :: charge_activ
c
      double precision, dimension (:), allocatable :: req_cal
      double precision, dimension (:), allocatable :: teq_cal
      double precision, dimension (:,:), allocatable :: dih_cal
c#####
      double precision, dimension (:), allocatable ::
     *  CHARGE     ! FLAG CHARGE (NATOM) FORMAT(5E16.8)
     *, AMASS      ! FLAG MASS  (NATOM) FORMAT(5E16.8)
     *, RK         ! FLAG BOND_FORCE_CONSTANT (NUMBND) FORMAT(5E16.8)
     *, REQ        ! FLAG BOND_EQUIL_VALUE  (NUMBND) FORMAT(5E16.8)
     *, TK         ! FLAG ANGLE_FORCE_CONSTANT (NUMANG) FORMAT(5E16.8) 
     *, TEQ        ! FLAG ANGLE_EQUIL_VALUE  (NUMANG)  FORMAT(5E16.8) 
     *, PK         ! FLAG DIHEDRAL_FORCE_CONSTANT (NPTRA) FORMAT(5E16.8)
     *, PN         ! FLAG DIHEDRAL_PERIODICITY  (NPTRA) FORMAT(5E16.8)
     *, PHASE      ! FLAG DIHEDRAL_PHASE  (NPTRA) FORMAT(5E16.8)
     *, SOLTY      ! FLAG SOLTY  (NATYP) FORMAT(5E16.8)
     *, CN1        ! FLAG LENNARD_JONES_ACOEF ( NTYPES*(NTYPES+1)/2 ) FORMAT(5E16.8)
     *, CN2        ! FLAG LENNARD_JONES_BCOEF  ( NTYPES*(NTYPES+1)/2 ) 
c#####
      double precision, dimension (:,:), allocatable :: tot_matr
      double precision, dimension (:), allocatable :: tot_targ, solut
      double precision, dimension (:,:), allocatable :: excl_fact
      double precision, dimension (:,:), allocatable :: trx, try, trz
      double precision, dimension (:,:), allocatable :: felec_x
     *, felec_y, felec_z
      double precision, dimension (:,:), allocatable :: fvdw_x
     *, fvdw_y, fvdw_z
      double precision, dimension (:,:), allocatable :: ref_fx
     *, ref_fy, ref_fz
      double precision, dimension (:,:), allocatable :: bonded_fx
     *, bonded_fy, bonded_fz
      double precision, dimension (:,:), allocatable :: bonded_matr
      double precision, dimension (:), allocatable :: bonded_targ
      double precision, dimension (:), allocatable :: bonded_solu
c
      call introduce()
c
      call read_nat_qm(nat, nr_qm)
      allocate( cpmd2grom(nat) )
      allocate( grom2cpmd(nat) )
      call read_qmmm_order(nat, nr_qm, cpmd2grom, grom2cpmd)
c
      write(*,*)''
      write(*,*)'-------------------------------------'
      write(*,*)'please enter number of frames (integer)
     * stored in the FM_REF_XXX files'
      read(*,*)nr_frame
      if(nr_frame .le. 0)then
       stop'main: nothing to do'
      endif
c
      allocate( nr_mm(nr_frame) )
      allocate( frame_ind(nr_frame) )
      allocate( charg_tp(nr_qm) )
c
      allocate( ib_qm(max_nr_qmb) )
      allocate( jb_qm(max_nr_qmb) )
      allocate( qmb_typ(max_nr_qmb) )
      allocate( icb_qm(max_nr_qmb) )
      allocate( req_cal( max_nr_qmb ) )
c
      allocate( ia_qm(max_nr_qma) )
      allocate( ja_qm(max_nr_qma) )
      allocate( ka_qm(max_nr_qma) )
      allocate( qma_typ( max_nr_qma ) )
      allocate( ica_qm( max_nr_qma ) )
      allocate( teq_cal(max_nr_qma) )
c
      allocate( id_qm(max_nr_qmd) )
      allocate( jd_qm(max_nr_qmd) )
      allocate( kd_qm(max_nr_qmd) )
      allocate( ld_qm(max_nr_qmd) )
      allocate( qmd_typ(max_nr_qmd) )
      allocate( icd_qm(max_nr_qmd) )
      allocate(  pro_imp(max_nr_qmd) )
      allocate( dih_cal( nr_frame, max_nr_qmd ) )
c
      allocate( tr_qm_x(nr_frame, nr_qm) )
      allocate( tr_qm_y(nr_frame, nr_qm) )
      allocate( tr_qm_z(nr_frame, nr_qm) )
c
      allocate( tr_mm_x(nr_frame, max_nr_mm) )
      allocate( tr_mm_y(nr_frame, max_nr_mm) )
      allocate( tr_mm_z(nr_frame, max_nr_mm) )
c
      allocate( tr_ele_pot(nr_frame, max_nr_mm) )
      allocate( tr_ele_fd_x(nr_frame, max_nr_mm) )
      allocate( tr_ele_fd_y(nr_frame, max_nr_mm) )
      allocate( tr_ele_fd_z(nr_frame, max_nr_mm) )
c
      allocate( tr_chj(nr_frame, nr_qm) )
c
      allocate( charge_activ (nat) )
c
c#######
      allocate( amb_pt(amb_nr_pt) )
      include 'amber_top_read.inc'
c#######
c
      write(*,*)' '
      write(*,*)'-------------------------------'
      write(*,*)'preparing some parameters for charge fiting'
c
      w_v = 1.0d0
      w_e = 0.20d0
      w_h = 0.0020d0
      w_q = 10000000.0d0
c
      write(*,*)'below are default values of the weights'
      write(*,"(a20,f20.5)")'potential w_v: ',w_v
      write(*,"(a20,f20.5)")'field w_e: ',w_e
      write(*,"(a20,f20.5)")'Hirshfeld w_h: ',w_h
      write(*,"(a20,e20.5)")'total charge w_q: ',w_q
c
      write(*,*)' '
      write(*,*)'now, if you want to keep them: enter 1'
      write(*,*)'or to change them: enter 2'
      read(*,*)j_test
      if(j_test .eq. 2)then
         write(*,*)'w_v = '
         read(*,*)w_v
c
         write(*,*)'w_e = '
         read(*,*)w_e
c
         write(*,*)'w_h = '
         read(*,*)w_h
c
         write(*,*)'w_q = '
         read(*,*)w_q
c
         write(*,*)''
         write(*,*)'new weights'
         write(*,"(a20,f20.5)")'potential w_v: ',w_v
         write(*,"(a20,f20.5)")'field w_e: ',w_e
         write(*,"(a20,f20.5)")'Hirshfeld w_h: ',w_h
         write(*,"(a20,e20.5)")'total charge w_q: ',w_q
         write(*,*)''
c
      else
          if(j_test .ne. 1)then
          stop'main: dont know what to do'
          endif
c
      endif
c
      w_v = dsqrt( w_v )
      w_e = dsqrt( w_e )
      w_h = dsqrt( w_h )
      w_q = dsqrt( w_q )
c
      call equivalent_charg(nat, nr_qm
     *, cpmd2grom, IGRAPH
     *, nr_char_tp
     *, charg_tp)
c
      call esti_tota_charg(nat, nr_qm
     *, grom2cpmd, CHARGE
     *, charg_tp
     *, tot_qm_charg )
c
      call read_pip( nr_frame, nr_qm, max_nr_mm
     *, nr_mm, frame_ind
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z )
c
      call  read_chj(nat, nr_qm, nr_frame, frame_ind, grom2cpmd
     *, tr_chj)
c
      call matr_size_gen(nr_frame, nr_qm, nr_mm 
     *, nr_li_pot, nr_li_fd, nr_li_hirsh
     *, tot_matr_size )
c
      allocate ( tot_matr( tot_matr_size, nr_char_tp) )
      allocate ( tot_targ (tot_matr_size) )
      allocate ( solut (nr_char_tp) )
c
      call matr_gen(nr_frame,nr_qm, max_nr_mm,nr_char_tp
     *, nr_li_pot, nr_li_fd, nr_li_hirsh
     *, nr_mm, charg_tp
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot, tr_ele_fd_x, tr_ele_fd_y,tr_ele_fd_z,tr_chj
     *, w_v, w_e, w_h, w_q
     *, tot_qm_charg
     *, tot_matr_size, tot_matr, tot_targ )
c
      call least_quare(tot_matr_size, nr_char_tp
     *, tot_matr, tot_targ
     *, solut)
c
      call write_charge(NATOM,nat, nr_qm
     *, nr_char_tp, solut
     *, charg_tp, grom2cpmd
     *, IGRAPH , CHARGE
     *, charge_activ)
c
      call realat_sd_poten(nr_frame, nat, nr_qm
     *, grom2cpmd, charge_activ
     *, max_nr_mm, nr_mm 
     *, tr_qm_x, tr_qm_y, tr_qm_z
     *, tr_mm_x, tr_mm_y, tr_mm_z
     *, tr_ele_pot
     *, tr_ele_fd_x, tr_ele_fd_y, tr_ele_fd_z)
c
      write(*,*)''
      write(*,*)'-----------------------------------'
      write(*,*)'Now, I am done with the chagre fit !'
      write(*,*)'enter 1 to continue'
      write(*,*)'or any other number to stop'
      read(*,*)j_test
      if(j_test .ne. 1)then
        stop'main: stopped'
      endif
      write(*,*)''
c
c##################################
      write(*,*)'#################################'
      write(*,*)'####### START BONDED FIT ########'
      write(*,*)'#################################'
c
      call bond_types(nat, nr_qm, max_nr_qmb
     *, grom2cpmd
     *, NATOM, ISYMBL, IGRAPH
     *, NBONH, IBH, JBH, ICBH
     *, NBONA, IB, JB, ICB 
     *, nr_qmb, ib_qm, jb_qm
     *, nr_qmb_typ, qmb_typ 
     *, icb_qm )
c
       call angle_types(nat, nr_qm, max_nr_qma
     *, grom2cpmd
     *, NATOM, ISYMBL, IGRAPH
     *, NTHETH, ITH, JTH, KTH, ICTH
     *, NTHETA, IT,JT,KT,ICT
     *, nr_qma, ia_qm, ja_qm, ka_qm
     *, nr_qma_typ, qma_typ
     *, ica_qm )
c
      call dihedral_types(nat, nr_qm, max_nr_qmd
     *, grom2cpmd
     *, NATOM, ISYMBL, IGRAPH
     *, NPHIH, IPH,JPH,KPH,LPH,ICPH
     *, NPHIA, IP,JP,KP,LP,ICP
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm
     *, nr_qmd_typ, qmd_typ
     *, icd_qm, pro_imp)
c
      allocate( excl_fact( nat, nat ) )
      call nonb_exclude_fact(IBH,JBH,NBONH,   IB,JB,NBONA
     *, ITH,JTH,KTH,NTHETH,   IT,JT,KT,NTHETA
     *, IPH,JPH,KPH,LPH,NPHIH,   IP,JP,KP,LP,NPHIA 
     *, nat, excl_fact)
c
      allocate(trx(nr_frame,nat), try(nr_frame,nat)
     *, trz(nr_frame,nat) )
      call read_traje( nr_frame, frame_ind, nat
     *, cpmd2grom
     *, trx, try, trz)
c
      call write_pdb(nat, nr_frame, NRES
     *, IPRES, LBRES
     *, IGRAPH
     *, trx, try, trz)
c
      allocate( felec_x(nr_frame, nr_qm)
     *, felec_y(nr_frame, nr_qm), felec_z(nr_frame, nr_qm) )
      call elect_forces( nat, nr_qm, nr_frame
     *, cpmd2grom
     *, charge_activ, excl_fact
     *, trx, try, trz
     *, felec_x, felec_y, felec_z )
c
      allocate( fvdw_x(nr_frame, nr_qm)
     *, fvdw_y(nr_frame, nr_qm), fvdw_z(nr_frame, nr_qm) )
      call vdw_forces(nat, nr_qm, nr_frame
     *, cpmd2grom
     *, NATOM, NTYPES, IAC, ICO, CN1, CN2
     *, excl_fact
     *, trx, try, trz 
     *, fvdw_x, fvdw_y, fvdw_z)
c
      allocate( ref_fx( nr_frame, nr_qm )
     *, ref_fy( nr_frame, nr_qm ), ref_fz( nr_frame, nr_qm ) )
      call read_fm_forces( nr_qm, nr_frame, frame_ind
     *,ref_fx, ref_fy, ref_fz )
c
      allocate( bonded_fx(nr_frame, nr_qm)
     *, bonded_fy(nr_frame, nr_qm), bonded_fz(nr_frame, nr_qm) )
      call bonded_force_store(nr_qm, nr_frame
     *,felec_x, felec_y, felec_z
     *,fvdw_x, fvdw_y, fvdw_z
     *,ref_fx, ref_fy, ref_fz
     *,bonded_fx, bonded_fy, bonded_fz)
c
c#####################
      call equ_bond_lengths(nat, nr_frame, max_nr_qmb
     *, nr_qmb, ib_qm, jb_qm
     *, trx, try, trz 
     *, req_cal )
c
      call equ_bond_angles(nat, nr_frame, max_nr_qma
     *, nr_qma, ia_qm, ja_qm, ka_qm
     *, trx, try, trz
     *, teq_cal )
c
      call evolu_dihedal(nat, nr_frame, max_nr_qmd
     *, nr_qmd, id_qm, jd_qm, kd_qm, ld_qm
     *, pro_imp
     *, trx, try, trz
     *, dih_cal)
c
      call writ_dihe( nr_frame, max_nr_qmd
     *, nr_qmd
     *, nr_qmd_typ, qmd_typ, pro_imp
     *, dih_cal )
c
c  bonded fit
c
      nr_bonded_row = nr_frame * nr_qm * 3
      nr_bonded_col = nr_qmb_typ + nr_qma_typ + nr_qmd_typ
c
      allocate(bonded_matr(nr_bonded_row, nr_bonded_col))
      allocate(bonded_targ(nr_bonded_row))
      allocate(bonded_solu(nr_bonded_col))
c
      call bonded_matrix(nat, nr_qm, nr_frame
     *, max_nr_qmb, max_nr_qma, max_nr_qmd 
     *, grom2cpmd
     *, trx, try, trz
     *, bonded_fx, bonded_fy, bonded_fz
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
      call least_quare(nr_bonded_row, nr_bonded_col
     *, bonded_matr, bonded_targ
     *, bonded_solu)
c
      call write_bonded( nat, nr_qm, nr_frame 
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
      write(*,*)''
      write(*,*)'thank you for using me, bye'
      write(*,*)''
c
      stop
      end
c      include 'introduction.f'
c      include 'read_pip_chj_qmmm_order.f'
c      include 'electro_matrix.f'
c      include 'least_square.f'
c      include 'write_charge.f'
c      include 'nonbonded_forces_cal.f'
c      include 'bond_angle_dihedral_type.f'
c      include 'equilibrium_values.f'
c      include 'bond_angle_dihedr_function.f'
c      include 'force_bon_ang_dih_imp.f'
c      include 'bonded_matrix_build.f'
c      include 'equivalent_atom_charges.f'
c      include 'relative_standard_deviation.f'
c      include 'write_bonded_parameters.f'


