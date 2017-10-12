c##########################
c  calculate bond distance
c  between x1,y1,z1 and
c  x2,y2,z2
c##########################
      double precision function bond_dist(x1,y1,z1
     *,x2,y2,z2)
c
      implicit none
c
      double precision x1,y1,z1
     *,x2,y2,z2
     *,d,dx,dy,dz
c
      dx=x2-x1
      dy=y2-y1
      dz=z2-z1 
c
      d=(dx*dx)+(dy*dy)+(dz*dz)
      d=dsqrt(d)
c
      bond_dist=d
c
      return
      end
c#####################################
c   calculate the bond angle (radian)
c   formed by x1,y1,z1
c             x2,y2,z2
c   and       x3,y3,z3
c######################################
      double precision function bond_ang(x1,y1,z1
     *,x2,y2,z2
     *,x3,y3,z3)
c
      implicit none
c
      double precision x1,y1,z1
     *,x2,y2,z2
     *,x3,y3,z3
     *,ax,ay,az
     *,bx,by,bz
c
      double precision ra,rb,dotpro,cosin
     *,angtem
c
      ax=x1-x2
      ay=y1-y2
      az=z1-z2
c
      bx=x3-x2
      by=y3-y2
      bz=z3-z2
c
c  modulus
c
      ra=(ax*ax)+(ay*ay)+(az*az)
      ra=dsqrt(ra)
c
      rb=(bx*bx)+(by*by)+(bz*bz)
      rb=dsqrt(rb)
c
c  dot product     
c
      dotpro=(ax*bx)+(ay*by)+(az*bz)
c
      cosin=dotpro/ra/rb
c
      angtem=dacos(cosin)
c
      bond_ang=angtem
c
      return
      end
c##############################
c   calculate the dihedral
c   angles between x1,y1,z1
c                  x2,y2,z2
c                  x3,y3,z3
c                  x4,y4,z4
c   rotate around 2-3 bond
c#############################
      double precision function dihed_ang(x1,y1,z1
     *,x2,y2,z2
     *,x3,y3,z3
     *,x4,y4,z4)
c
      implicit none
c
      double precision x1,y1,z1
     *,x2,y2,z2
     *,x3,y3,z3
     *,x4,y4,z4
c
      double precision r12x,r12y,r12z
     *,r43x,r43y,r43z
     *,r32x,r32y,r32z,d32
     *,r12r32,r43r32
     *,rx,ry,rz,dr
     *,sx,sy,sz,ds
     *,rs,posdihed,pi
     *,r32xs_x,r32xs_y,r32xs_z
     *,mixpro,dihed
c
c      pi=3.14159265358979311599796346854418516d0
c 
c r12
      r12x=x1-x2
      r12y=y1-y2
      r12z=z1-z2
c
c r43
      r43x=x4-x3
      r43y=y4-y3
      r43z=z4-z3
c
c r32
      r32x=x3-x2
      r32y=y3-y2
      r32z=z3-z2
c
      d32=(r32x*r32x)+(r32y*r32y)+(r32z*r32z)
c
      d32=dsqrt(d32)
c
      r32x=r32x/d32
      r32y=r32y/d32
      r32z=r32z/d32
c
c  dot products
c
      r12r32=(r12x*r32x)+(r12y*r32y)+(r12z*r32z)
c
      r43r32=(r43x*r32x)+(r43y*r32y)+(r43z*r32z)
c
c  R vector
c
      rx=r12x-r12r32*r32x
      ry=r12y-r12r32*r32y
      rz=r12z-r12r32*r32z
c
      dr=(rx*rx)+(ry*ry)+(rz*rz)
      dr=dsqrt(dr)
c
      rx=rx/dr
      ry=ry/dr
      rz=rz/dr
c
c  S vector
c
      sx=r43x-r43r32*r32x
      sy=r43y-r43r32*r32y
      sz=r43z-r43r32*r32z
c
      ds=(sx*sx)+(sy*sy)+(sz*sz)
      ds=dsqrt(ds)
c
c   be careful of dividing for zero
c
      sx=sx/ds
      sy=sy/ds
      sz=sz/ds
c
      rs=(rx*sx)+(ry*sy)+(rz*sz)
c
      posdihed=dacos(rs)
c
c  sign of the dihedral angle
c
      r32xs_x=(sy*r32z)-(sz*r32y)
      r32xs_y=(sz*r32x)-(sx*r32z)
      r32xs_z=(sx*r32y)-(sy*r32x)
c
      mixpro=(rx*r32xs_x)+(ry*r32xs_y)+(rz*r32xs_z)
c
      if(mixpro.ge.0.0d0)then
      dihed=posdihed
      else
      dihed=-1.0d0*posdihed
      endif
c
      dihed_ang=dihed
c
      return
      end
c########################################
c calculate improper angle
c i,j,k,l 
c where i is the central atom
c W. F. van Gunsteren et al. Biomolecular
c simulation: the gromos96 manual and user
c guide 
c  the atom order is i,j,k,l
c but in AMBER k,j,i,l
c########################################
      double precision function impro_ang(
     * xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl)
c
      implicit none
c
      double precision xi,yi,zi
     *,xj,yj,zj
     *,xk,yk,zk
     *,xl,yl,zl
     *,rij,rkj,rkl
     *,rmj,dmj,rnk,dnk
     *,ang_tmp,dotproduct
     *,sign_tmp
c
      dimension rij(3),rkj(3),rkl(3)
     *,rmj(3),rnk(3)
c rij
      rij(1)=xi-xj
      rij(2)=yi-yj
      rij(3)=zi-zj
c rkj
      rkj(1)=xk-xj
      rkj(2)=yk-yj
      rkj(3)=zk-zj 
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
c
c  be careful of dividing for zero
c
      ang_tmp = dotproduct ( rmj(1), rmj(2), rmj(3)
     *, rnk(1), rnk(2), rnk(3) )
c
      if( (dmj .ne. 0.0d0) .and. (dnk .ne. 0) )then
        ang_tmp = ang_tmp/dmj/dnk
      endif
      ang_tmp = dacos ( ang_tmp )
c
c
      sign_tmp = dotproduct ( rij(1), rij(2), rij(3)
     *,rnk(1), rnk(2), rnk(3) )
c
      if( sign_tmp .ne. 0.0d0 )then
      sign_tmp = sign_tmp/dabs(sign_tmp)
      endif
c
       impro_ang = sign_tmp * ang_tmp
c
      return
      end
c#####################################
c  calculate the cross product
c######################################
      subroutine cross_product(x1,y1,z1
     *,x2,y2,z2
     *,x3,y3,z3)
c
      implicit none
c
      double precision x1,y1,z1
     *,x2,y2,z2
     *,x3,y3,z3
c
      x3=(y1*z2)-(y2*z1)
c
      y3=(x2*z1)-(x1*z2)
c
      z3=(x1*y2)-(x2*y1)
c
      return
      end
c############################
c  dot product
c############################
      double precision function dotproduct(
     * x1,y1,z1
     *,x2,y2,z2)
c
      implicit none
c
      double precision x1,y1,z1
     *,x2,y2,z2
c
      dotproduct=(x1*x2)+(y1*y2)+(z1*z2)
c
      return
      end
c

