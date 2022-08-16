# testing the re2/map/par readers
#

      using Revise
      using PyCall
      using PyPlot 

#      include("JNek_PARALLEL.jl")         # JNek_PARALLEL
      
      include("/home/prabal/workstation/git/julia_scripts/nek_scripts/JNek_IO.jl")               # JNek_IO

#      import JNek_IO
      
      using MPI

      nid0  = 0

      if (MPI.Initialized() == false)
        MPI.Init()
      end  
        
      comm = MPI.COMM_WORLD
      rank = MPI.Comm_rank(comm)

      rea   = "box.rea"
      re2   = "taylor.re2"
      map   = "taylor.ma2"
      fld1   = "prctaylor0.f00001"
      fld2   = "prctaylor0.f00002"
      fld3   = "prctaylor0.f00003"
      fld4   = "prctaylor0.f00004"
      fld5   = "prctaylor0.f00005"
      fld6   = "prctaylor0.f00006"
      fld7   = "prctaylor0.f00007"
      fld8   = "exttaylor0.f00001"
      fld9   = "tmptaylor0.f00001"
      fld10  = "tmptaylor0.f00002"
      fld11  = "tmptaylor0.f00003"
      fld12  = "tmptaylor0.f00004"
      fld13  = "tmptaylor0.f00005"

      wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl = JNek_IO.read_re2(re2,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x1,y1,z1,u1,v1,w1,p1,t1 = JNek_IO.read_fld(fld1,MPI,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x2,y2,z2,u2,v2,w2,p2,t2 = JNek_IO.read_fld(fld2,MPI,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x3,y3,z3,u3,v3,w3,p3,t3 = JNek_IO.read_fld(fld3,MPI,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x4,y4,z4,u4,v4,w4,p4,t4 = JNek_IO.read_fld(fld4,MPI,nid0)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x5,y5,z5,u5,v5,w5,p5,t5 = JNek_IO.read_fld(fld5,MPI,nid0)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x6,y6,z6,u6,v6,w6,p6,t6 = JNek_IO.read_fld(fld6,MPI,nid0)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x7,y7,z7,u7,v7,w7,p7,t7 = JNek_IO.read_fld(fld7,MPI,nid0)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x8,y8,z8,u8,v8,w8,p8,t8 = JNek_IO.read_fld(fld8,MPI,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x9,y9,z9,u9,v9,w9,p9,t9 = JNek_IO.read_fld(fld9,MPI,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x10,y10,z10,u10,v10,w10,p10,t10 = JNek_IO.read_fld(fld10,MPI,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x11,y11,z11,u11,v11,w11,p11,t11 = JNek_IO.read_fld(fld11,MPI,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x12,y12,z12,u12,v12,w12,p12,t12 = JNek_IO.read_fld(fld12,MPI,nid0)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x13,y13,z13,u13,v13,w13,p13,t13 = JNek_IO.read_fld(fld13,MPI,nid0)
     
#       pp = p3[:,1,1,1:3];
#       pvec = pp[:];
#
#       xx = x3[:,1,1,1:3];
#       xvec = xx[:];

       pp2 = p2[4,:,1,[1 4 7]];
       pvec2 = pp2[:];

       yy2 = y2[4,:,1,[1 4 7]];
       yvec = yy2[:];

#       pp5 = p5[4,:,1,[1 4 7]];
#       pvec5 = pp5[:];
#
#       yy5 = y5[4,:,1,[1 4 7]];
#       yvec = yy5[:];
#
#       pp6 = p6[4,:,1,[1 4 7]];
#       pvec6 = pp6[:];
#
#       yy6 = y6[4,:,1,[1 4 7]];
#       yvec = yy6[:];
#
#       pp7 = p7[4,:,1,[1 4 7]];
#       pvec7 = pp7[:];
#
#       yy7 = y7[4,:,1,[1 4 7]];
#       yvec = yy7[:];
#
#
#       plot(yvec,pvec5)
#       plot(yvec,pvec6)
#       plot(yvec,pvec7)

      if rank == 0
        println("Done")
      end  

#      MPI.Finalize()
