c#########################
c  write the optimzed partial 
c  charges for qm atoms
c########################
      subroutine write_charge(NATOM,nat, nr_qm
     *, nr_char_tp, solut
     *, charg_tp, grom2cpmd
     *, IGRAPH, CHARGE
     *, charge_activ)
c
      implicit none
      integer nr_char_tp, charg_tp, grom2cpmd
     *, nat, nr_qm, NATOM
     *,i,icount,itmp, jtmp
      character(4) IGRAPH
      double precision solut, CHARGE, charge_activ
     *, dtmp
      dimension solut(nr_char_tp)
     *, charg_tp(nr_qm), grom2cpmd(nat)
     *, IGRAPH(NATOM), CHARGE(NATOM)
     *, charge_activ (nat)
c
      write(*,*)''
      write(*,*)'-------------------------'
      write(*,*)'writing optimized charges'
      if(nat .ne. NATOM) then
      stop'write_charge: wrong either NATOM or nat'
      endif
c
      do i=1,NATOM
      charge_activ(i)=CHARGE(i)/18.22230d0
      enddo
c
      open(1,file='optimized_charges.out',status='unknown')
      rewind(1)
c
      icount=0
      dtmp = 0.0d0
c
      do i=1,NATOM
      itmp = grom2cpmd(i)
      if(itmp .le. nr_qm) then
           icount = icount + 1
           jtmp = charg_tp (icount)
           if(jtmp .ne. 0)then
             write(1,"(2I8,a8,f15.6,f15.5)")icount, i, IGRAPH(i)
     *                             ,solut(jtmp),solut(jtmp)
c
           charge_activ(i) = solut(jtmp)
c
           dtmp = dtmp + solut(jtmp)
c
           endif
      endif
      enddo
c
      write(1,*)' '
      write(1,"(a,f15.8)")'Q = ', dtmp
c
      close(1)
c
      write(*,*)'the whole set of chagres is updated'
      write(*,*)'they are in gromos ordering
     * and in the unit of e'
c
      write(*,*)'optimized charges have been written to
     * optimized_charges.out'
      write(*,*)'writing optimized charges done!'
      write(*,*)''
c
c      open(1,file='active_charges.dat',status='unknown')
c      rewind(1)
c      do i=1,NATOM
c      write(1,"(I8,a8,f15.6)")i,IGRAPH(i), charge_activ(i)
c      enddo
c      close(1)
c
      return
      end
c#######################################
c   write pdb for the first frame
c#######################################
      subroutine write_pdb(nat, nr_frame, NRES
     *, IPRES, LBRES
     *, IGRAPH
     *, trx, try, trz)
c
      implicit none
      integer nat, nr_frame
     *, NRES, IPRES
     *, i,j,icount
      character(4) IGRAPH, LBRES
      double precision trx, try, trz
      dimension IPRES (NRES), LBRES(NRES)
     *, IGRAPH(nat)
     *, trx(nr_frame,nat), try(nr_frame,nat), trz(nr_frame,nat)
c
      write(*,*)''
      write(*,*)'--------------------------------'
      write(*,*)'writing the first frame to pdb'
      write(*,*)'please check file first_frame.pdb'
c
      open(1,file='first_frame.pdb',status='unknown')
      rewind(1)
c
      icount = 0
      do i=1,NRES-1
         do j= IPRES(i), IPRES(i+1)-1
         icount = icount + 1
           write(1,"(a4,i7,x,a3,xx,a3,i6,xxxx,3f8.3,2f6.2)")'ATOM'
     *,      j,IGRAPH(j)
     *,            LBRES(i),i
     *,   trx(1,j), try(1,j), trz(1,j), 1.00d0, 0.00d0
         enddo
      enddo
c
      do i = icount+1, nat
        write(1,"(a4,i7,x,a3,xx,a3,i6,xxxx,3f8.3,2f6.2)")'ATOM'
     *,      i,IGRAPH(i)
     *,            LBRES(NRES),NRES
     *,   trx(1,i), try(1,i), trz(1,i), 1.00d0, 0.00d0
      enddo
c
      close(1)
c
      write(*,*)'writing the first frame to pdb done'
      write(*,*)''
c
      return
      end
