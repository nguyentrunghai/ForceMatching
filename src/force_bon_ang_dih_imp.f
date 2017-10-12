c##############################################
c  force on two particles due to
c  bond stretching
c
c   flag = 1 calculate forces
c
c   flag = 2 calculate derivative wrt k
c
c############################################
      subroutine bon_str_force(flag
     *,k_const,r_eq
     *,x1,y1,z1
     *,x2,y2,z2
     *,f1,f2)
c
      implicit none
c
      integer flag
c
      double precision k_const,r_eq
     *,x1,y1,z1
     *,x2,y2,z2
     *,f1,f2
     *,d,dx,dy,dz
     *,norx,nory,norz
     *,dtmp
c
      dimension f1(3),f2(3)
c
      dx=x2-x1
      dy=y2-y1 
      dz=z2-z1
c
      d=(dx*dx)+(dy*dy)+(dz*dz)
      d=dsqrt(d)
      if(d.eq.0.0d0)then
      stop" bon_str_force: dist = 0"
      endif
c
      norx=dx/d
      nory=dy/d
      norz=dz/d
c
      if(flag.eq.1)then
      dtmp=2.0d0*k_const*(d-r_eq)
      else
          if(flag.eq.2)then
            dtmp=2.0d0*(d-r_eq)
          else
             stop"bon_str_force: wrong flag"
          endif
      endif
c
      f1(1)=dtmp*norx
      f1(2)=dtmp*nory
      f1(3)=dtmp*norz
c
      f2=-1.0d0*f1
c
      return
      end
c#####################################
c  force of angle bend
c  W. F. van Gunsteren et al. Biomolecular
c simulation: the gromos96 manual and user
c guide
c###################################
      subroutine ang_ben_for(flag
     *,k_theta,thet_eq
     *,xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,fi,fj,fk)
c
      implicit none
c
      integer flag
c
      double precision k_theta,thet_eq
     *,xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,fi,fj,fk
c
      double precision rij_x, rij_y, rij_z, dij
     *,rkj_x, rkj_y, rkj_z, dkj
     *,theta,bond_ang,pi
     *,derv,dconst1,dconst2
c
      dimension fi(3),fj(3),fk(3)  
c
      rij_x = xi - xj
      rij_y = yi - yj
      rij_z = zi - zj
      dij = rij_x*rij_x + rij_y*rij_y + rij_z*rij_z
      dij = dsqrt( dij )
c
      rkj_x = xk - xj
      rkj_y = yk - yj
      rkj_z = zk - zj
      dkj = rkj_x*rkj_x + rkj_y*rkj_y + rkj_z*rkj_z
      dkj = dsqrt( dkj )
c
      theta = bond_ang( xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk )
c
      pi=dacos(-1.0d0)
c
      if((theta.eq.pi).or.(theta.eq.0.0d0))then
      stop'ang_ben_for: theta is exactly pi or zero'
      endif
c
      if(flag.eq.1)then
      derv = 2.0d0*k_theta*( theta - thet_eq )
      else
      if(flag.eq.2)then
        derv = 2.0d0*( theta - thet_eq )
      else
      stop'ang_ben_for: flag is not correct'
      endif
      endif
c
      dconst1 = -derv/dsin(theta)/dij
      fi(1) = dconst1*( rij_x*dcos(theta)/dij - rkj_x/dkj )
      fi(2) = dconst1*( rij_y*dcos(theta)/dij - rkj_y/dkj )
      fi(3) = dconst1*( rij_z*dcos(theta)/dij - rkj_z/dkj )
c
      dconst2 = -derv/dsin(theta)/dkj
      fk(1) = dconst2*( rkj_x*dcos(theta)/dkj - rij_x/dij )
      fk(2) = dconst2*( rkj_y*dcos(theta)/dkj - rij_y/dij )
      fk(3) = dconst2*( rkj_z*dcos(theta)/dkj - rij_z/dij)
c
      fj(1) = -fi(1) -fk(1)
      fj(2) = -fi(2) -fk(2)
      fj(3) = -fi(3) -fk(3)
c
      return
      end
c######################################
c  to calculate the angle bending forces
c  acing on three particles
c  particles 1 and 2 at the two ends
c  particle 3 at the center
c
c  flag = 1  force
c
c   flag = 2 derivative wrt k_theta
c J. Comp Phys. 100 17-24  "An algorithm for
c    calculating angle-dpendent forces on vector
c       Computers "
c######################################
      subroutine ang_ben_for_1(flag
     *,k_theta,thet_eq
     *,x1,y1,z1
     *,x3,y3,z3
     *,x2,y2,z2
     *,f1,f3,f2)
c
      implicit none
c
      integer flag
c
      double precision k_theta,thet_eq
     *,x1,y1,z1
     *,x3,y3,z3
     *,x2,y2,z2
     *,f1,f3,f2
     *,f1_x,f1_y,f1_z
     *,f2_x,f2_y,f2_z
     *,f3_x,f3_y,f3_z
c
      double precision r31_x,r31_y,r31_z,d31
     *,r32_x,r32_y,r32_z,d32
     *,r21_x,r21_y,r21_z,d21
     *,p31,p32
     *,theta,bond_ang,dconst
     *,ft_x,ft_y,ft_z
     *,pi
c
      dimension p31(3,3),p32(3,3)
     *,f1(3),f3(3),f2(3)
c 31
      r31_x=x1-x3
      r31_y=y1-y3
      r31_z=z1-z3
      d31=(r31_x*r31_x)+(r31_y*r31_y)+(r31_z*r31_z)
      d31=dsqrt(d31)
      call p_matrix(r31_x,r31_y,r31_z,p31)
c 32
      r32_x=x2-x3
      r32_y=y2-y3
      r32_z=z2-z3
      d32=(r32_x*r32_x)+(r32_y*r32_y)+(r32_z*r32_z)
      d32=dsqrt(d32)
      call p_matrix(r32_x,r32_y,r32_z,p32)
c 21
      r21_x=x1-x2
      r21_y=y1-y2
      r21_z=z1-z2
      d21=(r21_x*r21_x)+(r21_y*r21_y)+(r21_z*r21_z)
      d21=dsqrt(d21)
c
      theta=bond_ang(x1,y1,z1
     *,x3,y3,z3
     *,x2,y2,z2)
c
      pi=dacos(-1.0d0)
c
      if((theta.eq.pi).or.(theta.eq.0.0d0))then
      stop'ang_ben_for ER, theta is exactly pi or zero'
      endif
c
      if(flag.eq.1)then
      dconst=-2.0d0*k_theta*(theta-thet_eq)/dsin(theta)
      else
          if(flag.eq.2)then
          dconst=-2.0d0*(theta-thet_eq)/dsin(theta)
          else
          stop'ang_ben_for ER, flag is not correct'
          endif
      endif
c
      if ( (d31.eq.0.0d0).or. (d32.eq.0.0d0) )then
      stop'ang_ben_for ER d31 or d32 exactly zero'
      endif
c
      dconst=dconst/d31/d32
      ft_x=dconst*r21_x
      ft_y=dconst*r21_y
      ft_z=dconst*r21_z
c######################################
c f1
      call matr_vect(p31,ft_x,ft_y,ft_z
     *,f1_x,f1_y,f1_z)
c
      f1(1)=f1_x
      f1(2)=f1_y
      f1(3)=f1_z
c
c f2
      call matr_vect(p32,ft_x,ft_y,ft_z
     *,f2_x,f2_y,f2_z)
       f2_x=-1.0d0*f2_x
       f2_y=-1.0d0*f2_y
       f2_z=-1.0d0*f2_z
c
       f2(1)=f2_x
       f2(2)=f2_y
       f2(3)=f2_z
c
c f3
      f3_x=-1.0d0*(f1_x+f2_x)
      f3_y=-1.0d0*(f1_y+f2_y)
      f3_z=-1.0d0*(f1_z+f2_z)
c
      f3(1)=f3_x
      f3(2)=f3_y
      f3(3)=f3_z
c
      return
      end
c###########################################
c dihedral rotation force acting on 4 particles
c  i,j,k,l
c  flag = 1  force
c  flag = 2  derivative wrt force constant
c  H. Bekker et al. Force and Virial
c of the Torsional-angle-dependent potentials
c  J. Comp. Chem. 1995, 16, 527-533
c###########################################
      subroutine dihedr_rot_for(flag
     *,vn_const,n_const,gama
     *,xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl
     *,fi,fj,fk,fl)
c
      implicit none
c
      integer flag
c
      double precision vn_const,n_const,gama
     *,xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl
     *,fi,fj,fk,fl
     *,fi_x,fi_y,fi_z
     *,fj_x,fj_y,fj_z
     *,fk_x,fk_y,fk_z
     *,fl_x,fl_y,fl_z
c
      double precision rij_x,rij_y,rij_z,d2ij
     *,rkj_x,rkj_y,rkj_z,d2kj
     *,rkl_x,rkl_y,rkl_z,d2kl
     *,mx,my,mz,m2
     *,nx,ny,nz,n2
     *,phi,dihed_ang
     *,const1,const2,const3
     *,dotproduct,pi
c
      dimension fi(3),fj(3),fk(3),fl(3)
c ij
      rij_x=xi-xj
      rij_y=yi-yj
      rij_z=zi-zj
      d2ij=(rij_x*rij_x)+(rij_y*rij_y)+(rij_z*rij_z)
c kj
      rkj_x=xk-xj
      rkj_y=yk-yj
      rkj_z=zk-zj
      d2kj=(rkj_x*rkj_x)+(rkj_y*rkj_y)+(rkj_z*rkj_z)
c kl
      rkl_x=xk-xl
      rkl_y=yk-yl
      rkl_z=zk-zl
      d2kl=(rkl_x*rkl_x)+(rkl_y*rkl_y)+(rkl_z*rkl_z)
c
      call cross_product(rij_x,rij_y,rij_z
     *,rkj_x,rkj_y,rkj_z
     *,mx,my,mz)
      m2=(mx*mx)+(my*my)+(mz*mz)
c
      call cross_product(rkj_x,rkj_y,rkj_z
     *,rkl_x,rkl_y,rkl_z
     *,nx,ny,nz)
      n2=(nx*nx)+(ny*ny)+(nz*nz)
c
      phi=dihed_ang(xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl)
c
      if( d2kj.eq.0.0d0  )then
       stop'dihedr_for ER, d2kj zero'
      endif
c
      if ( m2 .eq. 0.0d0 )then
      stop'dihedr_for ER, m2 zero'
      endif
c
      if ( n2 .eq. 0.0d0 )then
      stop'dihedr_for ER, n2 zero'
      endif
c
      pi=dacos(-1.0d0)
      if((phi.eq.pi).or.(phi.eq.0.0d0))then
      stop'dihedr_for ER, phi is exactly pi or zero'
      endif
c
      if ( flag .eq. 1 ) then
      const1= -1.0d0 * vn_const * n_const*dsin(n_const*phi-gama)
      else
         if ( flag .eq. 2 ) then
           const1= -1.0d0 * n_const*dsin(n_const*phi-gama)
          else
           stop'dihedr_ro_for ER: flag is not correct'
         endif
      endif
      const1=const1*dsqrt(d2kj)
c fi
      fi_x=-const1*mx/m2
      fi_y=-const1*my/m2
      fi_z=-const1*mz/m2
c fl
      fl_x=const1*nx/n2
      fl_y=const1*ny/n2
      fl_z=const1*nz/n2
c     
      const2=dotproduct(rij_x,rij_y,rij_z
     *,rkj_x,rkj_y,rkj_z)
      const2=const2/d2kj
c
      const3=dotproduct(rkl_x,rkl_y,rkl_z
     *,rkj_x,rkj_y,rkj_z)
      const3=const3/d2kj
c fj
      fj_x=-1.0d0*fi_x+const2*fi_x-const3*fl_x
      fj_y=-1.0d0*fi_y+const2*fi_y-const3*fl_y
      fj_z=-1.0d0*fi_z+const2*fi_z-const3*fl_z
c fk
      fk_x=-1.0d0*fl_x-const2*fi_x+const3*fl_x
      fk_y=-1.0d0*fl_y-const2*fi_y+const3*fl_y
      fk_z=-1.0d0*fl_z-const2*fi_z+const3*fl_z
c
c  output
c
      fi(1)=fi_x
      fi(2)=fi_y
      fi(3)=fi_z
c
      fl(1)=fl_x
      fl(2)=fl_y
      fl(3)=fl_z
c
      fj(1)=fj_x
      fj(2)=fj_y
      fj(3)=fj_z
c
      fk(1)=fk_x
      fk(2)=fk_y
      fk(3)=fk_z
c
      return
      end
c#########################################
c  calculate force for the 
c  improper dihedrals with the same
c  energy function as that of dihedral
c W. F. van Gunsteren et al. Biomolecular
c simulation: the gromos96 manual and user
c guide
c#######################################
      subroutine impro_for(flag
     *,vn_cons,n_const,gama
     *,xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl
     *,fi,fj,fk,fl)
c
      implicit none
c
      integer flag
c
      double precision vn_cons,n_const,gama
     *,xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl
     *,fi,fj,fk,fl
c
      double precision zeta,impro_ang
     *,derv
     *,rij,rkj,rkl
     *,rmj,dmj,rnk,dnk
     *,dotproduct
     *,dkj,dconst,dconst1,dconst2
c
      dimension fi(3),fj(3),fk(3),fl(3)
     *,rij(3),rkj(3),rkl(3)
     *,rmj(3),rnk(3)
c
      zeta=impro_ang(
     * xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl)
c
      if( flag.eq.1 )then
      derv = -1.0d0 * vn_cons * n_const*dsin( n_const*zeta - gama )
      else
        if ( flag.eq.2 )then
          derv = -1.0d0 * n_const * dsin( n_const*zeta - gama )
        else
        stop'impro_for: flag is not correct'
        endif
      endif
c
c rij
      rij(1)=xi-xj
      rij(2)=yi-yj
      rij(3)=zi-zj
c rkj
      rkj(1)=xk-xj
      rkj(2)=yk-yj
      rkj(3)=zk-zj
      dkj = rkj(1)*rkj(1) + rkj(2)*rkj(2) + rkj(3)*rkj(3)
      dkj = dsqrt(dkj)
c rkl
      rkl(1)=xk-xl
      rkl(2)=yk-yl
      rkl(3)=zk-zl
c rmj
      call cross_product( rij(1), rij(2), rij(3)
     *, rkj(1), rkj(2), rkj(3)
     *, rmj(1), rmj(2), rmj(3) )
       dmj = rmj(1)*rmj(1) + rmj(2)*rmj(2) + rmj(3)*rmj(3)
       dmj = dsqrt(dmj)
c rnk
      call cross_product( rkj(1), rkj(2), rkj(3)
     *, rkl(1), rkl(2), rkl(3)
     *, rnk(1), rnk(2), rnk(3) )
      dnk = rnk(1)*rnk(1) + rnk(2)*rnk(2) + rnk(3)*rnk(3)
      dnk = dsqrt(dnk)
c fi!!!!!
      if( dmj .eq. 0.0d0 )then
        stop'impro_for: dmj .eq. 0.0d0'
      endif
          dconst = -1.0d0*derv*dkj/dmj/dmj
      fi=dconst*rmj
c fl !!!
      if( dnk .eq. 0.0d0 )then
        stop'impro_for: dnk .eq. 0.0d0'
      endif
         dconst = derv*dkj/dnk/dnk
      fl = dconst*rnk
c fj!!!!!
      dconst1 = dotproduct ( rij(1), rij(2), rij(3)
     *, rkj(1), rkj(2), rkj(3) )
      if( dkj .eq. 0.0d0 )then
        stop'impro_for: dkj .eq. 0.0d0'
      endif
        dconst1 = dconst1/dkj/dkj - 1.0d0
c 
      dconst2 = dotproduct ( rkl(1), rkl(2), rkl(3)
     *, rkj(1), rkj(2), rkj(3) )
      if( dkj .eq. 0.0d0 )then
        stop'impro_for: dkj .eq. 0.0d0'
      endif
         dconst2 = dconst2/dkj/dkj
c
      fj = dconst1*fi - dconst2*fl
c fk
      fk = -fi-fj-fl
c
      return
      end
c######################################
c  to calculate the P matrix for a vector
c  
c######################################
      subroutine p_matrix(x,y,z
     *,p)
c
      implicit none
c
      integer i,j
c
      double precision x,y,z,p
     *,pt,d2
c
      dimension p(3,3)
     *,pt(3,3)
c
      d2=(x*x)+(y*y)+(z*z)
c
      pt(1,1)=x*x
      pt(2,2)=y*y
      pt(3,3)=z*z
c
      pt(1,2)=x*y
      pt(1,3)=x*z
c
      pt(2,1)=x*y
      pt(2,3)=y*z
c
      pt(3,1)=x*z
      pt(3,2)=y*z
c
      do i=1,3
      do j=1,3
c
      pt(i,j)=pt(i,j)/d2
c
      enddo
      enddo
c
      do i=1,3
      do j=1,3
c
      if(i.eq.j)then
      p(i,j)=1.0d0-pt(i,j)
      else
      p(i,j)=-1.0d0*pt(i,j)
      endif
c
      enddo
      enddo
c
      return
      end
c######################################
c to apply a matrix to a vector,
c  resulting in a vector
c#####################################
      subroutine matr_vect(p,x,y,z
     *,x1,y1,z1)
c
      implicit none
c
      double precision p,x,y,z
     *,x1,y1,z1
c
      dimension p(3,3)
c
      x1=p(1,1)*x+p(1,2)*y+p(1,3)*z
      y1=p(2,1)*x+p(2,2)*y+p(2,3)*z
      z1=p(3,1)*x+p(3,2)*y+p(3,3)*z
c
c
      return
      end
c
