c##############################
c some info about the code
c#############################
      subroutine introduce()
c
      implicit none
c
      write(*,*)''
      write(*,*)'--------------------------------------------------'
      write(*,"(20x,a)")'FORCE MATCHING'
      write(*,"(10x,a)")'Trung Hai Nguyen'
      write(*,"(10x,a)")'Email: nguyentrungh@gmail.com'
      write(*,*)'--------------------------------------------------'
      write(*,*)''
      write(*,*)'--------------------------------------------------'
      write(*,"(2x,a)")'Reference:'
      write(*,"(4x,a)")'Patrick Maurer, Alessandro Laio, Hakan W.
     * Hugosson,'
      write(*,"(4x,a)")'Maria Carola Colombo and Ursula Rothlisberger'
      write(*,"(4x,a)")'Automated Parametrization of Biomolecular'
      write(*,"(4x,a)")'Force Fields from Quantum Mechanics/Molecular'
      write(*,"(4x,a)")'Mechanics (QM/MM) Simulations through 
     *Force Matching'
      write(*,"(4x,a)")'J. Chem. Theory Comput. 2007, 3, 628.'
      write(*,*)'--------------------------------------------------'
      write(*,*)''
      write(*,*)'--------------------------------------------------'
      write(*,"(2x,a)")'Required inputs:'
      write(*,"(5x,a)")'AMBER_TOP (amber topology, prepared by Leap)'
      write(*,"(5x,a)")'FM_REF_PIP (electrostatic grid from CPMD)'
      write(*,"(5x,a)")'FM_REF_CHJ (Hirshfeld charges from CPMD)'
      write(*,"(5x,a)")'FM_REF_FORCES (ab initio forces from CPMD)'
      write(*,"(5x,a)")'QMMM_ORDER (from CPMD)'
      write(*,"(5x,a)")'TRAJECTORY (from CPMD)'
      write(*,*)'--------------------------------------------------'
      write(*,*)''
      write(*,*)'--------------------------------------------------'
      write(*,"(2x,a)")'Outputs'
      write(*,"(5x,a)")'Many output flies, but most important:'
      write(*,"(5x,a)")'optimized_charges.out (final optimized charges)'
      write(*,"(5x,a)")'fm.frcmod (final bonded parameters)'
      write(*,*)'--------------------------------------------------'
      write(*,*)''
      write(*,*)'--------------------------------------------------'
      write(*,"(8x,a)")'IMPORTANT NOTES:'
      write(*,"(2x,a)")'The code can only provide new charges, 
     *force constants and equilibrium length values.'
      write(*,"(2x,a)")'It cannot change the TOPOLOGY 
     *of the molecule.'
      write(*,*)'Therefore, the user needs to take care of 
     * the correct definition of e.g. atom types, bond, 
     * angle types ... in 
     * the AMBER_TOP file, using LEAP program from 
     * AMBERtool.'
      write(*,*)'For dihedrals, the code can only provide new 
     *force constants (barrier heights).'
      write(*,*)'The user needs to correctly define the number of 
     *minima and phase in advance in the AMBER_TOP file.'
      write(*,*)'--------------------------------------------------'
      write(*,*)''
c
      write(*,*)'Press enter to continue'
      read(*,*)
c
      return
      end
