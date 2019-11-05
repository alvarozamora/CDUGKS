import "regent"

-- Helper modules to handle PNG files and command line arguments
local Config = require("config")
--local Dump = require("dump") -- TODO

-- Some C APIs
local c     = regentlib.c
local sqrt  = regentlib.sqrt(double)
local cmath = terralib.includec("math.h")
local PI = cmath.M_PI

-- Field space for Simulation Parameters
fspace params{

  -- Simulation Parameters
  N  : int32[3],
  NV : int32[3],
  Nc : int64,
  Nv : int64,
  R : double,
  K : double,
  g : double,
  w : double,
  Cv : double,
  ur : double,
  Tr : double,
  Pr : double,
  effD : int32,
  BCs : int32[3],
  Vmin : double[3],
  Vmax : double[3],

}

-- Field space for reduced distribution functions and their gradients
fspace grid{

  -- Reduced Distribution Functions  
  g : double,
  b : double,

}

-- Field space for rectangular mesh
fspace mesh{

  -- Cell Position
  x    : double,
  y    : double,
  z    : double,
  
  -- Cell Size
  dx   : double,
  dy   : double,
  dz   : double,

}

-- Field space for conserved variables at cell center and interfaces
fspace W{

    -- Conserved Variables
    rho  : double,
    rhov : double[3],
    rhoE : double

}

-- Field space for velocity mesh
fspace vmesh
{

  -- Velocity and Weight
  v : double,
  w : double

}


task TestProblem(r_params : region(ispace(int1d), params), testProblem : int32)
where
  reads writes(r_params)
do
  for e in r_params do

	-- Sod Shock
        if testProblem == 1 then

                --Dimensionality
                r_params[e].effD = 1

                --Spatial Resolution
                r_params[e].N[0]  = 256
		r_params[e].N[1]  = 1
 		r_params[e].N[2]  = 1
                
		--Velocity Resolution
		r_params[e].NV[0] = 256
		r_params[e].NV[1] = 1
		r_params[e].NV[2] = 1
                
		r_params[e].Vmin[0] = -10
		r_params[e].Vmin[1] = 0
		r_params[e].Vmin[2] = 0
                
		r_params[e].Vmax[0] = 10
		r_params[e].Vmax[1] = 0
		r_params[e].Vmax[2] = 0
                
		-- Number of Cells
		r_params[e].Nc = r_params[e].N[0]*r_params[e].N[1]*r_params[e].N[2]
                r_params[e].Nv  = r_params[e].NV[0]*r_params[e].NV[1]*r_params[e].NV[2]

                -- Boundary Conditions
                r_params[e].BCs[0] = 1 
		r_params[e].BCs[1] = 0
		r_params[e].BCs[2] = 0


                -- Physical Parameters
                r_params[e].R   = 0.5           			-- Gas Constant
                r_params[e].K   = 2.0          				-- Internal DOF
                r_params[e].Cv  = (3+r_params[e].K)*r_params[e].R/2     -- Specific Heat
                r_params[e].g   = (r_params[e].K+5)/(r_params[e].K+3)   -- gamma -- variable name taken
                r_params[e].w   = 0.5           			-- Viscosity exponent
                r_params[e].ur  = 1e-4          			-- Reference Visc
                r_params[e].Tr  = 1.0           			-- Reference Temp
                r_params[e].Pr  = 2/3.          			-- Prandtl Number

 
        -- Kelvin-Helmholtz
        elseif testProblem == 2 then 

                --Dimensionality
                r_params[e].effD = 2

		var num : int32 = 128
                --Spatial Resolution
                r_params[e].N[0]  = num
		r_params[e].N[1]  = num
 		r_params[e].N[2]  = 1
                
		--Velocity Resolution
		r_params[e].NV[0] = num
		r_params[e].NV[1] = num
		r_params[e].NV[2] = 1
                
		r_params[e].Vmin[0] = -10
		r_params[e].Vmin[1] = -10
		r_params[e].Vmin[2] = 0
                
		r_params[e].Vmax[0] = 10
		r_params[e].Vmax[1] = 10
		r_params[e].Vmax[2] = 0
                
		-- Number of Cells
		r_params[e].Nc = r_params[e].N[0]*r_params[e].N[1]*r_params[e].N[2]
                r_params[e].Nv  = r_params[e].NV[0]*r_params[e].NV[1]*r_params[e].NV[2]


                -- Boundary Conditions
                r_params[e].BCs[0] = 0
		r_params[e].BCs[1] = 0
		r_params[e].BCs[2] = 0


                -- Physical Parameters
                r_params[e].R   = 0.5          				-- Gas Constant
                r_params[e].K   = 2.0           			-- Internal DOF
                r_params[e].Cv  = (3+r_params[e].K)*r_params[e].R/2     -- Specific Heat
                r_params[e].g   = (r_params[e].K+5)/(r_params[e].K+3)   -- gamma -- variable name taken
                r_params[e].w   = 0.5           			-- Viscosity exponent
                r_params[e].ur  = 1e-4          			-- Reference Visc
                r_params[e].Tr  = 1.0           			-- Reference Temp
                r_params[e].Pr  = 2/3.          			-- Prandtl Number
	end
  end
end

terra dumpdouble(f : &c.FILE, val : double)
  var a : double[1]
  a[0] = val
  c.fwrite(&a, 8, 1, f)
end

terra dumpint32(f : &c.FILE, val : int32)
  var a : int32[1]
  a[0] = val
  c.fwrite(&a, 4, 1, f)
end

terra dumpbool(f : &c.FILE, val : bool)
  var a : int32[1]
  if val==true then
    a[0] = 1
  else 
    a[0] = 0
  end
  c.fwrite(&a, 4, 1, f)
end

task factorize2d(parallelism : int) : int2d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y = i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int2d { size_x, size_y }
end

task factorize3d(parallelism : int) : int2d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y = i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int2d { size_x, size_y }
end


terra wait_for(x : int) return 1 end
task block_task(r_image : region(ispace(int1d), particle))
where
  reads writes(r_image)
do
  return 1
end


task Dump(r_W : region(ispace(int1d), W), iter : int32)
where
  reads (r_W)
do
  var filename : int8[1000]
  c.sprintf([&int8](filename), './Data/rho_%4d',iter)
  var g = c.fopen(filename,'wb')

  for e in r_W do
    dumpdouble(g, r_W[e].rho)
  end
  __fence(__execution, __block)
  return 1
end

task toplevel()
  var config : Config
  config:initialize_from_command()

  -- Simulation Parameters
  var testProblem : int32 = 1
  var r_params = region(ispace(int1d, 1), params)
  TestProblem(r_params, testProblem)
  var N  : int32[3] = r_params[0].N
  var NV : int32[3] = r_params[0].NV
  var Nc : int64 = r_params[0].Nc
  var Nv : int64 = r_params[0].Nv
  var R : double = r_params[0].R
  var K : double = r_params[0].K
  var g : double = r_params[0].g
  var w : double = r_params[0].w
  var Cv : double = r_params[0].Cv
  var ur : double = r_params[0].ur
  var Tr : double = r_params[0].Tr
  var Pr : double = r_params[0].Pr
  var effD : int32 = r_params[0].effD
  var BCs : int32[3] = r_params[0].BCs
  var Vmin : double[3] = r_params[0].Vmin
  var Vmax : double[3] = r_params[0].Vmax

  __fence(__execution, __block) 
  c.printf("Simulation Parameters\n")
  if testProblem > 0 then c.printf("testProblem = %d\n", testProblem) end
  c.printf("N = {%d, %d, %d}, NV = {%d, %d, %d}, effD = %d\n", N[0], N[1], N[2], NV[0], NV[1], NV[2], effD)
  c.printf("BCs = {%d, %d, %d}, Vmin = {%f, %f, %f}, Vmax = {%f, %f, %f}\n", BCs[0], BCs[1], BCs[2], Vmin[0], Vmin[1], Vmin[2], Vmax[0], Vmax[1], Vmax[2])
  c.printf("R = %f, K = %f, g = %f, Cv = %f\n", R, K, g, Cv)
  c.printf("w = %f, ur = %f, Tr = %f, Pr = %f\n", w, ur, Tr, Pr)

end

regentlib.start(toplevel)
