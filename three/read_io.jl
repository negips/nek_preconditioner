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
      re2   = "three.re2"
      map   = "three.ma2"
      fld1   = "prcthree0.f00001"
      fld2   = "prcthree0.f00002"
      fld3   = "prcthree0.f00003"
      fld4   = "prcthree0.f00004"
      fld5   = "prcthree0.f00005"
      fld6   = "prcthree0.f00006"
      fld7   = "prcthree0.f00007"
      fld8   = "extthree0.f00001"
      fld9   = "tmpthree0.f00001"
      fld10  = "tmpthree0.f00002"
      fld11  = "tmpthree0.f00003"
      fld12  = "tmpthree0.f00004"
      fld13  = "tmpthree0.f00005"

#      wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl = JNek_IO.read_re2(re2,nid0)
      re2fld1 = JNek_IO.read_re2_struct(re2,nid0)

      fld1 = JNek_IO.read_fld_struct(fld1,MPI,nid0)

      fld2 = JNek_IO.read_fld_struct(fld2,MPI,nid0)

      fld3 = JNek_IO.read_fld_struct(fld3,MPI,nid0)

      fld4 = JNek_IO.read_fld_struct(fld4,MPI,nid0)

#      fld5 = JNek_IO.read_fld_struct(fld5,MPI,nid0)
#
#      fld6 = JNek_IO.read_fld_struct(fld6,MPI,nid0)
#
#      fld7 = JNek_IO.read_fld_struct(fld7,MPI,nid0)
#
#      fld8 = JNek_IO.read_fld_struct(fld8,MPI,nid0)

      fld9 = JNek_IO.read_fld_struct(fld9,MPI,nid0)

      fld10 = JNek_IO.read_fld_struct(fld10,MPI,nid0)

      fld11 = JNek_IO.read_fld_struct(fld11,MPI,nid0)

      fld12 = JNek_IO.read_fld_struct(fld12,MPI,nid0)

      fld13 = JNek_IO.read_fld_struct(fld13,MPI,nid0)
     
#       pp7 = fld7.p[4,:,1,[1 4 7]];
#       pvec7 = pp7[:];
#
#       yy7 = fld7.y[4,:,1,[1 4 7]];
#       yvec = yy7[:];
#
#       plot(yvec,pvec7)

      if rank == 0
        println("Done")
      end  

#      MPI.Finalize()
