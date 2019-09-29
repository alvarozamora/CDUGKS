import "regent"

-- Helper modules to handle PNG files and command line arguments
local Config = require("config")
--local Dump = require("dump") -- TODO
local coloring   = require("coloring_util")

-- Some C APIs
local c     = regentlib.c
local sqrt  = regentlib.sqrt(double)
local cmath = terralib.includec("math.h")
local PI = cmath.M_PI
local isnan = regentlib.isnan(double)

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
  Tf : double,
  dtdump: double

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

		var num : int32 = 8192
                --Spatial Resolution
                r_params[e].N[0]  = num
		r_params[e].N[1]  = 1
 		r_params[e].N[2]  = 1
                
		--Velocity Resolution
		r_params[e].NV[0] = num
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
                r_params[e].Nv = r_params[e].NV[0]*r_params[e].NV[1]*r_params[e].NV[2]

                -- Boundary Conditions
                r_params[e].BCs[0] = 1 
		r_params[e].BCs[1] = 0
		r_params[e].BCs[2] = 0


                -- Physical Parameters
                r_params[e].R   = 0.5           			-- Gas Constant
                r_params[e].K   = 2.0          				-- Internal DOF
                r_params[e].Cv  = (3+r_params[e].K)*r_params[e].R/2.0   -- Specific Heat
                r_params[e].g   = (r_params[e].K+5)/(r_params[e].K+3.0) -- gamma -- variable name taken
                r_params[e].w   = 0.5           			-- Viscosity exponent
                r_params[e].ur  = 1e-4          			-- Reference Visc
                r_params[e].Tr  = 1.0           			-- Reference Temp
                r_params[e].Pr  = 2.0/3.0          			-- Prandtl Number

		-- Simulation Parameters
		r_params[e].Tf = 0.15					-- Stop Time
		r_params[e].dtdump = r_params[e].Tf/200			-- Time Between Dumps
 
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
                r_params[e].Nv = r_params[e].NV[0]*r_params[e].NV[1]*r_params[e].NV[2]


                -- Boundary Conditions
                r_params[e].BCs[0] = 0
		r_params[e].BCs[1] = 0
		r_params[e].BCs[2] = 0


                -- Physical Parameters
                r_params[e].R   = 0.5          				-- Gas Constant
                r_params[e].K   = 2.0           			-- Internal DOF
                r_params[e].Cv  = (3+r_params[e].K)*r_params[e].R/2.0   -- Specific Heat
                r_params[e].g   = (r_params[e].K+5)/(r_params[e].K+3.0) -- gamma -- variable name taken
                r_params[e].w   = 0.5           			-- Viscosity exponent
                r_params[e].ur  = 1e-4          			-- Reference Visc
                r_params[e].Tr  = 1.0           			-- Reference Temp
                r_params[e].Pr  = 2.0/3.0          			-- Prandtl Number

		-- Simulation Parameters
		r_params[e].Tf = 1.2					-- Stop Time
		r_params[e].dtdump = r_params[e].Tf/200			-- Time Between Dumps
	end
  end
end

task NewtonCotes(vxmesh : region(ispace(int1d), vmesh), vymesh : region(ispace(int1d), vmesh),
                 vzmesh : region(ispace(int1d), vmesh),
                 NV : int32[3], Vmin : double[3], Vmax : double[3])
where
  reads writes(vxmesh, vymesh, vzmesh)
do
  -- Newton-Cotes Weights
  var k : int32
  var n : int32
  var b : double
  var a : double
  var dh : double

  --c.assert(NV[0]%4 == 0 and (NV[1]%4 == 0 or NV[1] == 1) and (NV[2]%4 == 0 or NV[2] == 1)) -- Using np = 4 Newton Cotes


  -- First Dimension
  a = Vmin[0]
  b = Vmax[0]
  n = (NV[0])/4  -- no longer subtracting 1
  dh = (b-a)/(NV[0]-1)*4 -- mine (Al)

  for kx = 0, NV[0] do
    vxmesh[kx].v = Vmin[0] + kx*dh/4
  end 

  for kx = 0, n do

    --vxmesh[4*kx].w   = 14.0
    --vxmesh[4*kx+1].w = 32.0
    --vxmesh[4*kx+2].w = 12.0
    --vxmesh[4*kx+3].w = 32.0

    --printf("Changing %d\n", 4*kx)
    --printf("Changing %d\n", 4*kx+1)
    --printf("Changing %d\n", 4*kx+2)
    --printf("Changing %d\n", 4*kx+3)

    -- TODO Trying Lower Order Integration
    vxmesh[4*kx].w   = 1.0
    vxmesh[4*kx+1].w = 1.0
    vxmesh[4*kx+2].w = 1.0
    vxmesh[4*kx+3].w = 1.0

  end
  
  --vxmesh[0].w = 7.0
  --vxmesh[NV[0]-1].w = 7.0


  -- Check Weights
  --for kx = 0, NV[0] do
  --  c.printf("Main.hh vxmesh[%d].v = %f\n", kx, vxmesh[kx].v)
  --end
  --for kx = 0, NV[0] do
  --  c.printf("Main.hh vxmesh[%d].w = %f\n", kx, vxmesh[kx].w)
  --end
  
   
  for kx = 0, NV[0] do
    vxmesh[kx].w = vxmesh[kx].w*dh/4. -- TODO Trying Lower Order Integration old factor was /90.}
  end

  -- Second Dimension
  if NV[1] >= 4 then
    a = Vmin[1]
    b = Vmax[1]
    n = (NV[1])/4  -- no longer subtracting 1
    -- dh=(b-a)/n -- old from dugks
    dh = (b-a)/(NV[1]-1)*4 -- mine (Al)

    for ky = 0, NV[1] do
      vymesh[ky].v = Vmin[1] + ky*dh/4
    end 

    for ky = 0, n do

      --vymesh[4*ky].w   = 14.0
      --vymesh[4*ky+1].w = 32.0
      --vymesh[4*ky+2].w = 12.0
      --vymesh[4*ky+3].w = 32.0

      -- TODO Trying Lower Order Integration
      vymesh[4*ky].w   = 1.0
      vymesh[4*ky+1].w = 1.0
      vymesh[4*ky+2].w = 1.0
      vymesh[4*ky+3].w = 1.0
    end
  
    --vymesh[0].w = 7.0
    --vymesh[NV[1]-1].w = 7.0

    for ky = 0, NV[1] do
      vymesh[ky].w = vymesh[ky].w*dh/4. -- TODO Trying Lower Order Integration old factor was /90.}
    end
  elseif NV[1] == 1 then
    vymesh[0].v = 0
    vymesh[0].w = 1
  else
    c.printf("Error with Newton Cotes (number of vy points)\n")
  end 

  -- Check Weights
  --for ky = 0, NV[1] do
  --  c.printf("Main.hh vymesh[%d].v = %f\n", ky, vymesh[ky].v)
  --end
  --for ky = 0, NV[1] do
  --  c.printf("Main.hh vymesh[%d].w = %f\n", ky, vymesh[ky].w)
  --end

  -- Third Dimension
  if NV[2] >= 4 then
    a = Vmin[2]
    b = Vmax[2]
    n = (NV[2])/4  -- no longer subtracting 1
    -- dh=(b-a)/n -- old from dugks
    dh = (b-a)/(NV[2]-1)*4 -- mine (Al)

    for kz = 0, NV[2] do
      vzmesh[kz].v = Vmin[2] + k*dh/4
    end 

    for kz = 0, n do
      --vzmesh[4*kz].w   = 14.0
      --vzmesh[4*kz+1].w = 32.0
      --vzmesh[4*kz+2].w = 12.0
      --vzmesh[4*kz+3].w = 32.0

      -- TODO Trying Lower Order Integration
      vzmesh[4*kz].w   = 1.0
      vzmesh[4*kz+1].w = 1.0
      vzmesh[4*kz+2].w = 1.0
      vzmesh[4*kz+3].w = 1.0
    end

    --vzmesh[0].w = 7.0
    --vzmesh[NV[2]-1].w = 7.0


    for kz = 0, NV[2] do
      vzmesh[kz].w = vzmesh[kz].w*dh/4. -- TODO Trying Lower Order Integration old factor was /90.}
    end
  elseif NV[2] == 1 then
    vzmesh[0].v = 0
    vzmesh[0].w = 1
  else
    c.printf("Error with Newton Cotes (number of vz points)\n")
  end 

  -- Check Weights
  --for kz = 0, NV[2] do
  --  c.printf("Main.hh vzmesh[%d].v = %f\n", kz, vzmesh[kz].v)
  --end
  --for kz = 0, NV[2] do
  --  c.printf("Main.hh vzmesh[%d].w = %f\n", kz, vzmesh[kz].w)
  --end

end

task InitializeMesh(r_mesh : region(ispace(int3d), mesh), N : int32[3], MeshType : int32)
where
  reads writes(r_mesh)
do
  -- User-Specified Mesh
  if MeshType == 0 then
    -- TODO
  
  -- Rectangular Mesh
  elseif MeshType == 1 then
    
    for e in r_mesh do

      var Nx : int32 = N[0]
      var Ny : int32 = N[1]
      var Nz : int32 = N[2]

      var dX : double = 1.0/Nx
      var dY : double = 1.0/Ny
      var dZ : double = 1.0/Nz

      r_mesh[e].x = dX*(e.x + 1.0/2.0)
      r_mesh[e].y = dY*(e.y + 1.0/2.0)
      r_mesh[e].z = dZ*(e.z + 1.0/2.0)

      r_mesh[e].dx = dX
      r_mesh[e].dy = dY
      r_mesh[e].dz = dZ

    end
  -- Nested Rectangular 
  elseif MeshType == 2 then 
    -- TODO  
  end
end

task InitializeW(r_W : region(ispace(int3d), W), 
                 r_mesh : region(ispace(int3d), mesh),
                 N : int32[3], NV : int32[3], testProblem : int32, R : double, Cv : double)
where
  reads writes(r_W),
  reads(r_mesh)
do
  if testProblem == 1 then
    
    var pL : double = 1.0
    var pR : double = 1.0/8.0
    var PL : double = 1.0
    var PR : double = 1.0/10.0

    for e in r_W do
      if r_mesh[e].x <= 0.5 then
        r_W[e].rho = pL
	r_W[e].rhov[0] = 0
	r_W[e].rhov[1] = 0
	r_W[e].rhov[2] = 0
        r_W[e].rhoE = Cv*PL/R + 0.5*r_W[e].rhov[0]*r_W[e].rhov[0]/r_W[e].rho
      else
        r_W[e].rho = pR
	r_W[e].rhov[0] = 0
	r_W[e].rhov[1] = 0
	r_W[e].rhov[2] = 0
        r_W[e].rhoE = Cv*PR/R + 0.5*r_W[e].rhov[0]*r_W[e].rhov[0]/r_W[e].rho
      end
      --c.printf("W[%d] = {%f, %f, %f}\n", e.x, r_W[e].rho, r_W[e].rhov[1], r_W[e].rhoE)
    end
  elseif testProblem == 2 then
    
    var p1 : double = 2.0
    var p2 : double = 1.0
    var P1 : double = 2.0
    var P2 : double = 2.0
    var vrel : double = 2.0
    var amp : double = 0.04

    for e in r_W do
      if (0.25 <= r_mesh[e].y and r_mesh[e].y <= 0.75) then
        r_W[e].rho = p1
       	r_W[e].rhov[0] = vrel/2
       	r_W[e].rhov[1] = amp*cmath.sin(2*PI*r_mesh[e].x)*r_W[e].rho
        r_W[e].rhov[2] = 0
        r_W[e].rhoE = Cv*P1/R + 0.5*(r_W[e].rhov[0]*r_W[e].rhov[0] + r_W[e].rhov[1]*r_W[e].rhov[1])/r_W[e].rho
      else
	r_W[e].rho = p2
        r_W[e].rhov[0] = -vrel/2
        r_W[e].rhov[1] = amp*cmath.sin(2*PI*r_mesh[e].x)*r_W[e].rho
        r_W[e].rhov[2] = 0
        r_W[e].rhoE = Cv*P2/R + 0.5*(r_W[e].rhov[0]*r_W[e].rhov[0] + r_W[e].rhov[1]*r_W[e].rhov[1])/r_W[e].rho
      end      
    end
  end
end

terra Temperature(E : double, u : double, g : double, R : double)
  return (g - 1)/R*(E - 0.5*u*u)
end

terra geq(c2 : double, rho : double, T : double, R : double, effD : int32)
  
  var x : double = rho*cmath.exp(-c2/(2*R*T))*cmath.pow(2*PI*R*T, -double(effD)/2.0)

  return x
end

terra visc(T : double, ur : double, Tr : double, w : double)
  return ur*cmath.pow(T/Tr,w)
end

terra sgn(x : double)
  var sx : double = double(0 < x) - double(x < 0)
  return sx
end

terra abs(x : double)
  return x*sgn(x)
end

terra VanLeer(L : double, C : double, R : double, xL : double, xC : double, xR : double)
  
  if xR < xC then xR = xR + 1.0 end
  if xL > xC then xL = xL - 1.0 end

  var s1 : double = (C - L)/(xC - xL)
  var s2 : double = (R - C)/(xR - xC)

  if (C == L and C == R) then
    return 0.
  elseif (xC == xL or xC == xR) then
    return 0.
  else 
    return (sgn(s1) + sgn(s2))*(abs(s1) * abs(s2))/(abs(s1) + abs(s2))
  end  
end


terra BC(i : int32, j : int32, k : int32, Dim : int32, BCs : int32[3], N : int32[3])

  var IL : int32
  var JL : int32
  var KL : int32
  var IR : int32
  var JR : int32
  var KR : int32

  -- Periodic Boundary Conditions
  if (Dim == 0 and BCs[0] == 0) then
    IL = (i - 1 + N[0])%N[0]
    IR = (i + 1)%N[0]
    JL = j
    JR = j
    KL = k
    KR = k
  elseif (Dim == 1 and BCs[1] == 0) then
    IL = i
    IR = i
    JL = (j - 1 + N[1])%N[1]
    JR = (j + 1)%N[1]
    KL = k
    KR = k
  elseif (Dim == 2 and BCs[2] == 0) then
    IL = i
    IR = i
    JL = j
    JR = j
    KL = (k - 1 + N[2])%N[2]
    KR = (k + 1)%N[2]
  -- Dirichlet Boundary Conditions
  elseif (Dim == 0 and BCs[0] == 1) then
    IL = i - 1
    IR = i + 1
    if IL <  0 then IL = 0 end
    if IR == N[0] then IR = N[0] - 1 end
    JL = j
    JR = j
    KL = k
    KR = k
  elseif (Dim == 1 and BCs[1] == 1) then
    IL = i
    IR = i
    JL = j - 1
    JR = j + 1
    if JL <  0 then JL = 0 end
    if JR == N[1] then JR = N[1] - 1 end
    KL = k
    KR = k
  elseif (Dim == 2 and BCs[2] == 1) then
    IL = i
    IR = i
    JL = j
    JR = j
    KL = k - 1
    KR = k + 1
    if KL <  0 then KL = 0 end
    if KR == N[2] then KR = N[2] - 1 end
  -- Neumann Boundary Conditions
  elseif (Dim == 0 and BCs[0] == 2) then
    IL = i - 1
    IR = i + 1
    if IL <  0 then IL = 0 end
    if IR == N[0] then IR = N[0] - 1 end
    JL = j
    JR = j
    KL = k
    KR = k
  elseif (Dim == 1 and BCs[1] == 2) then
    IL = i
    IR = i
    JL = j - 1
    JR = j + 1
    if JL <  0 then JL = 0 end
    if JR == N[1] then JR = N[1] - 1 end
    KL = k
    KR = k
  elseif (Dim == 2 and BCs[2] == 2) then
    IL = i
    IR = i
    JL = j
    JR = j
    KL = k - 1
    KR = k + 1
    if KL <  0 then KL = 0 end
    if KR == N[2] then KR = N[2] - 1 end
  end 
  
  var LR : int32[6] 
  LR[0] = IL
  LR[1] = JL
  LR[2] = KL
  LR[3] = IR
  LR[4] = JR
  LR[5] = KR

  return LR
end

terra TimeStep(calcdt : double, dumptime : double, tend : double)
  var timestep : double = cmath.fmin(cmath.fmin(calcdt, dumptime), tend)
  return timestep
end 



-- Step 1: Phibar at interface
-- Step 1a: Phibar at Cell Center.
task Step1a(r_grid : region(ispace(int6d), grid),
            r_gridbarp : region(ispace(int6d), grid),
            r_S : region(ispace(int6d), grid),
            r_W : region(ispace(int3d), W),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            dt : double, R : double, K : double, Cv : double,
            g : double, w : double, ur : double, Tr : double,
            Pr : double, effD : int32)
where
  reads(r_grid, r_W, vxmesh, vymesh, vzmesh), 
  reads writes(r_S),
  reads writes(r_gridbarp)
do
  var tg : double
  var tb : double    
  var c2 : double
  var u : double
  var Xi : double[3]
  var T : double

  var slo : int3d = {r_grid.bounds.lo.x, r_grid.bounds.lo.y, r_grid.bounds.lo.z}
  var shi : int3d = {r_grid.bounds.hi.x, r_grid.bounds.hi.y, r_grid.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)
  
  var rhotest : double = 0
  for s in s3 do

    u = 0
    for d = 0, effD do
      u += r_W[s].rhov[d]/r_W[s].rho*r_W[s].rhov[d]/r_W[s].rho
    end
    u = sqrt(u)

    T = Temperature(r_W[s].rhoE/r_W[s].rho, u, g, R)
    if T < 0 then
      c.printf("T < 0, r_W[s].rhoE = %f, r_W[s].rho = %f, u = %f, g = %f, R = %f\n", r_W[s].rhoE, r_W[s].rho, u, g, R)
      regentlib.assert(T >= 0, "Negative Temperature\n")
    end

    tg = visc(T, ur, Tr, w)/r_W[s].rho/R/T
    tb = tg/Pr

    --if s.x == r_gridbarp.bounds.lo.x + 1 then c.printf("Step1a[%d]: u = %f, T = %f\n", s.x, u, T) end
    for v in v3 do
      var e : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}

      -- For Now...
      r_S[e].g = 0.
      r_S[e].b = 0.

      c2 = 0
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v     
      for d = 0, effD do
        c2 += (Xi[d]-r_W[s].rhov[d]/r_W[s].rho)*(Xi[d]-r_W[s].rhov[d]/r_W[s].rho)
      end

      var g_eq : double = geq(c2, r_W[s].rho, T, R, effD)
      var b_eq : double = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0

      r_gridbarp[e].g = (tg - dt/4.)/tg*r_grid[e].g + dt/(4.*tg)*g_eq + dt/4.*r_S[e].g
      r_gridbarp[e].b = (tb - dt/4.)/tb*r_grid[e].b + dt/(4.*tb)*b_eq + dt/4.*r_S[e].b

      if s.x == r_gridbarp.bounds.lo.x + 1 then rhotest += r_gridbarp[e].g*vxmesh[v.x].w end
      --if v.x == 128 then c.printf("phi to phibarplus[%d]: {%f to %f}, {%f to %f}\n", s.x, r_grid[e].g, r_gridbarp[e].g, r_grid[e].b, r_gridbarp[e].b) end
    
      if (isnan(r_gridbarp[e].g) == 1 or isnan(r_gridbarp[e].b) == 1) then

        c.printf("Step 1a: gp = %.12f, bp = %.12f, g = %.12f, b = %.12f, g_eq = %.12f, Sg = %.12f, Sb = %.12f, taus = {%.12f, %.12f}\n", r_gridbarp[e].g, r_gridbarp[e].b, r_grid[e].g, r_grid[e].b, g_eq, r_S[e].g, r_S[e].b, tg, tb)

        regentlib.assert(not [bool](isnan(r_gridbarp[e].g)), "Step 1a\n")
        regentlib.assert(not [bool](isnan(r_gridbarp[e].b)), "Step 1a\n")
    
      end

    end
  end
  --c.printf("rhotest = %f\n", rhotest)
end

--Step 1b: compute gradient of phibar to compute phibar at interface. compute phibar at interface.
task Step1b_sigx(r_gridbarp : region(ispace(int6d), grid),
            r_sig : region(ispace(int7d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plx_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            plx_gridbarp : region(ispace(int6d), grid),
            prx_gridbarp : region(ispace(int6d), grid),
            vxmesh : region(ispace(int1d), vmesh),            
            vymesh : region(ispace(int1d), vmesh),            
            vzmesh : region(ispace(int1d), vmesh),            
            BCs : int32[3], N : int32[3], effD : int32)
where 
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_gridbarp, prx_gridbarp),
  reads writes(r_sig)
do
  var Dim : int32 = 0 -- change when copy

  -- Cell Centers 
  var xC : double
  var xL : double
  var xR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  
  for s in s3 do
    xC = r_mesh[s].x -- change when copy
  
    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3] 
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do
      
      var e : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}
      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var eL : int6d = {IL, JL, KL, v.x, v.y, v.z}
      var eR : int6d = {IR, JR, KR, v.x, v.y, v.z}

      -- Computing phisigma, at cell 
      var gbpL : double
      var gbpR : double
      var bbpL : double
      var bbpR : double

      if r_gridbarp.bounds.lo.x == s.x then -- change when copy
        gbpL = plx_gridbarp[eL].g 
        bbpL = plx_gridbarp[eL].b 
        xL = plx_mesh[eL3].x
      else
        gbpL = r_gridbarp[eL].g
        bbpL = r_gridbarp[eL].b
        xL = r_mesh[eL3].x
      end

      if r_gridbarp.bounds.hi.x == s.x then -- change when copy
        gbpR = prx_gridbarp[eR].g 
        bbpR = prx_gridbarp[eR].b 
        xR = prx_mesh[eR3].x
      else
        gbpR = r_gridbarp[eR].g
        bbpR = r_gridbarp[eR].b
        xR = r_mesh[eR3].x
      end


      r_sig[e7].g = VanLeer(gbpL, r_gridbarp[e].g, gbpR, xL, xC, xR)
      r_sig[e7].b = VanLeer(bbpL, r_gridbarp[e].b, bbpR, xL, xC, xR)

      if s.x == 128 and v.x == 128 then
        --c.printf("r_sig[%d] = {%f, %f}\n", e7.x, r_sig[e7].g, r_sig[e7].b)
      end

      -- NAN checker
      if (isnan(r_sig[e7].g) == 1 or isnan(r_sig[e7].b) == 1) then

        c.printf("Step 1b_sigx: r_sig.g = %f, r_sig.b = %f\n", r_sig[e7].g, r_sig[e7].b)

        regentlib.assert(not [bool](isnan(r_sig[e7].g)), "Step 1b_sigx\n")
        regentlib.assert(not [bool](isnan(r_sig[e7].b)), "Step 1b_sigx\n")

      end
  
    end
  end
end

--Step 1b: compute gradient of phibar to compute phibar at interface. compute phibar at interface.
task Step1b_sigy(r_gridbarp : region(ispace(int6d), grid),
            r_sig : region(ispace(int7d), grid),
            r_mesh : region(ispace(int3d), mesh),
            ply_mesh : region(ispace(int3d), mesh),
            pry_mesh : region(ispace(int3d), mesh),
            ply_gridbarp : region(ispace(int6d), grid),
            pry_gridbarp : region(ispace(int6d), grid),
            vxmesh : region(ispace(int1d), vmesh),            
            vymesh : region(ispace(int1d), vmesh),            
            vzmesh : region(ispace(int1d), vmesh),            
            BCs : int32[3], N : int32[3], effD : int32)
where 
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, ply_mesh, pry_mesh, ply_gridbarp, pry_gridbarp),
  reads writes(r_sig)
do
  var Dim : int32 = 1 -- change when copy

  -- Cell Centers 
  var yC : double
  var yL : double
  var yR : double

  
  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)
  
  for s in s3 do
  
    yC = r_mesh[s].y -- change when copy

    var i : int32 = s.x
    var j : int32 = s.y
    var k : int32 = s.z

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3] 
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}
  
    for v in v3 do
      
      var e : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}
      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var eL : int6d = {IL, JL, KL, v.x, v.y, v.z}
      var eR : int6d = {IR, JR, KR, v.x, v.y, v.z}

      -- Computing phisigma, at cell 
      var gbpL : double
      var gbpR : double
      var bbpL : double
      var bbpR : double

      if r_gridbarp.bounds.lo.y == s.y then -- change when copy
        gbpL = ply_gridbarp[eL].g 
        bbpL = ply_gridbarp[eL].b 
        yL = ply_mesh[eL3].y
      else
        gbpL = r_gridbarp[eL].g
        bbpL = r_gridbarp[eL].b
        yL = r_mesh[eL3].y
      end

      if r_gridbarp.bounds.hi.y == s.y then -- change when copy
        gbpR = pry_gridbarp[eR].g 
        bbpR = pry_gridbarp[eR].b 
        yR = pry_mesh[eR3].y
      else
        gbpR = r_gridbarp[eR].g
        bbpR = r_gridbarp[eR].b
        yR = r_mesh[eR3].y
      end

      r_sig[e7].g = VanLeer(gbpL, r_gridbarp[e].g, gbpR, yL, yC, yR)
      r_sig[e7].b = VanLeer(bbpL, r_gridbarp[e].b, bbpR, yL, yC, yR)

      -- NAN checker
      if (isnan(r_sig[e7].g) == 1 or isnan(r_sig[e7].b) == 1) then

        c.printf("Step 1b_sigy: r_sig.g = %f, r_sig.b = %f\n", r_sig[e7].g, r_sig[e7].b)

        regentlib.assert(not [bool](isnan(r_sig[e7].g)), "Step 1b_sigy\n")
        regentlib.assert(not [bool](isnan(r_sig[e7].b)), "Step 1b_sigy\n")

      end
  
    end
  end
end

--Step 1b: compute gradient of phibar to compute phibar at interface. compute phibar at interface.
task Step1b_sigz(r_gridbarp : region(ispace(int6d), grid),
            r_sig : region(ispace(int7d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_gridbarp : region(ispace(int6d), grid),
            prz_gridbarp : region(ispace(int6d), grid),
            vxmesh : region(ispace(int1d), vmesh),            
            vymesh : region(ispace(int1d), vmesh),            
            vzmesh : region(ispace(int1d), vmesh),            
            BCs : int32[3], N : int32[3], effD : int32)
where 
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_gridbarp, prz_gridbarp),
  reads writes(r_sig)
do
  var Dim : int32 = 2 -- change when copy

  -- Cell Centers 
  var zC : double
  var zL : double
  var zR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)
  
  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    var i : int32 = s.x
    var j : int32 = s.y
    var k : int32 = s.z

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3] 
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}
  
    for v in v3 do
      
      var e : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}
      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var eL : int6d = {IL, JL, KL, v.x, v.y, v.z}
      var eR : int6d = {IR, JR, KR, v.x, v.y, v.z}

      -- Computing phisigma, at cell 
      var gbpL : double
      var gbpR : double
      var bbpL : double
      var bbpR : double

      if r_gridbarp.bounds.lo.z == s.z then -- change when copy
        gbpL = plz_gridbarp[eL].g 
        bbpL = plz_gridbarp[eL].b 
        zL = plz_mesh[eL3].z
      else
        gbpL = r_gridbarp[eL].g
        bbpL = r_gridbarp[eL].b
        zL = r_mesh[eL3].z
      end

      r_sig[e7].g = VanLeer(gbpL, r_gridbarp[e].g, gbpR, zL, zC, zR)
      r_sig[e7].b = VanLeer(bbpL, r_gridbarp[e].b, bbpR, zL, zC, zR)
     

      -- NAN checker
      if (isnan(r_sig[e7].g) == 1 or isnan(r_sig[e7].b) == 1) then

        c.printf("Step 1b_sigz: r_sig.g = %f, r_sig.b = %f\n", r_sig[e7].g, r_sig[e7].b)

        regentlib.assert(not [bool](isnan(r_sig[e7].g)), "Step 1b_sigz\n")
        regentlib.assert(not [bool](isnan(r_sig[e7].b)), "Step 1b_sigz\n")

      end
  
    end
  end
end

      
task Step1b_sigx_x(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plx_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            plx_sig : region(ispace(int7d), grid),
            prx_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_sig, prx_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    xC = r_mesh[s].x -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.x == s.x then
        gsigL = plx_sig[eL7].g 
        bsigL = plx_sig[eL7].b 
        xL = plx_mesh[eL3].x
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        xL = r_mesh[eL3].x
      end

      if r_sig.bounds.hi.x == s.x then
        gsigR = prx_sig[eR7].g
        bsigR = prx_sig[eR7].b
        xR = prx_mesh[eR3].x
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        xR = r_mesh[eR3].x
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, xL, xC, xR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, xL, xC, xR)
    end
  end
end

task Step1b_sigx_x2(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            prx_sig : region(ispace(int7d), grid),
            prx_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, prx_mesh, prx_sig, prx_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 0

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    var sC : double = r_mesh[s].dx -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
      var eR8 : int8d = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      var gsig : double
      var bsig : double
      var gsig2 : double
      var bsig2 : double
      var swap : double = 1.0
      if (vxmesh[v.x].v < 0 and Dim == 0) then
        if s.x == r_sigb.bounds.hi.x then
          gsig = prx_sig[eR7].g
          bsig = prx_sig[eR7].b
          gsig2 = prx_sig2[eR8].g
          bsig2 = prx_sig2[eR8].b
          sC = prx_mesh[eR3].dx
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dx
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[s].dx
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.x == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.x, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigy_x(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plx_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            plx_sig : region(ispace(int7d), grid),
            prx_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_sig, prx_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    xC = r_mesh[s].x -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
    
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.x == s.x then
        gsigL = plx_sig[eL7].g 
        bsigL = plx_sig[eL7].b 
        xL = plx_mesh[eL3].x
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        xL = r_mesh[eL3].x
      end

      if r_sig.bounds.hi.x == s.x then
        gsigR = prx_sig[eR7].g
        bsigR = prx_sig[eR7].b
        xR = prx_mesh[eR3].x
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        xR = r_mesh[eR3].x
      end

      r_sig2[e8].g = r_sig[e7].g + (xC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, xL, xC, xR)
      r_sig2[e8].b = r_sig[e7].b + (xC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, xL, xC, xR)
    end
  end
end

task Step1b_sigz_x(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plx_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            plx_sig : region(ispace(int7d), grid),
            prx_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_sig, prx_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    xC = r_mesh[s].x -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.x == s.x then
        gsigL = plx_sig[eL7].g 
        bsigL = plx_sig[eL7].b 
        xL = plx_mesh[eL3].x
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        xL = r_mesh[eL3].x
      end

      if r_sig.bounds.hi.x == s.x then
        gsigR = prx_sig[eR7].g
        bsigR = prx_sig[eR7].b
        xR = prx_mesh[eR3].x
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        xR = r_mesh[eR3].x
      end

      r_sig2[e8].g = r_sig[e7].g + (xC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, xL, xC, xR)
      r_sig2[e8].b = r_sig[e7].b + (xC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, xL, xC, xR)
    end
  end
end

task Step1b_sigx_y(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            ply_mesh : region(ispace(int3d), mesh),
            pry_mesh : region(ispace(int3d), mesh),
            ply_sig : region(ispace(int7d), grid),
            pry_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, ply_mesh, pry_mesh, ply_sig, pry_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    yC = r_mesh[s].y -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.y == s.y then
        gsigL = ply_sig[eL7].g 
        bsigL = ply_sig[eL7].b 
        yL = ply_mesh[eL3].y
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        yL = r_mesh[eL3].y
      end

      if r_sig.bounds.hi.y == s.y then
        gsigR = pry_sig[eR7].g
        bsigR = pry_sig[eR7].b
        yR = pry_mesh[eR3].y
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        yR = r_mesh[eR3].y
      end

      r_sig2[e8].g = r_sig[e7].g + (yC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, yL, yC, yR)
      r_sig2[e8].b = r_sig[e7].b + (yC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, yL, yC, yR)
    end
  end
end

task Step1b_sigy_y(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            ply_mesh : region(ispace(int3d), mesh),
            pry_mesh : region(ispace(int3d), mesh),
            ply_sig : region(ispace(int7d), grid),
            pry_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, ply_mesh, pry_mesh, ply_sig, pry_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    yC = r_mesh[s].y -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.y == s.y then
        gsigL = ply_sig[eL7].g 
        bsigL = ply_sig[eL7].b 
        yL = ply_mesh[eL3].y
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        yL = r_mesh[eL3].y
      end

      if r_sig.bounds.hi.y == s.y then
        gsigR = pry_sig[eR7].g
        bsigR = pry_sig[eR7].b
        yR = pry_mesh[eR3].y
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        yR = r_mesh[eR3].y
      end

      r_sig2[e8].g = r_sig[e7].g + (yC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, yL, yC, yR)
      r_sig2[e8].b = r_sig[e7].b + (yC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, yL, yC, yR)
    end
  end
end

task Step1b_sigz_y(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            ply_mesh : region(ispace(int3d), mesh),
            pry_mesh : region(ispace(int3d), mesh),
            ply_sig : region(ispace(int7d), grid),
            pry_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, ply_mesh, pry_mesh, ply_sig, pry_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    yC = r_mesh[s].y -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}

     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.y == s.y then
        gsigL = ply_sig[eL7].g 
        bsigL = ply_sig[eL7].b 
        yL = ply_mesh[eL3].y
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        yL = r_mesh[eL3].y
      end

      if r_sig.bounds.hi.y == s.y then
        gsigR = pry_sig[eR7].g
        bsigR = pry_sig[eR7].b
        yR = pry_mesh[eR3].y
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        yR = r_mesh[eR3].y
      end

      r_sig2[e8].g = r_sig[e7].g + (yC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, yL, yC, yR)
      r_sig2[e8].b = r_sig[e7].b + (yC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, yL, yC, yR)
    end
  end
end

task Step1b_sigx_z(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_sig, prz_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.z == s.z then
        gsigL = plz_sig[eL7].g 
        bsigL = plz_sig[eL7].b 
        zL = plz_mesh[eL3].z
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        zL = r_mesh[eL3].z
      end

      if r_sig.bounds.hi.z == s.z then
        gsigR = prz_sig[eR7].g
        bsigR = prz_sig[eR7].b
        zR = prz_mesh[eR3].z
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        zR = r_mesh[eR3].z
      end

      r_sig2[e8].g = r_sig[e7].g + (zC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, zL, zC, zR)
      r_sig2[e8].b = r_sig[e7].b + (zC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, zL, zC, zR)
    end
  end
end

task Step1b_sigy_z(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_sig, prz_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.z == s.z then
        gsigL = plz_sig[eL7].g 
        bsigL = plz_sig[eL7].b 
        zL = plz_mesh[eL3].z
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        zL = r_mesh[eL3].z
      end

      if r_sig.bounds.hi.z == s.z then
        gsigR = prz_sig[eR7].g
        bsigR = prz_sig[eR7].b
        zR = prz_mesh[eR3].z
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        zR = r_mesh[eR3].z
      end

      r_sig2[e8].g = r_sig[e7].g + (zC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, zL, zC, zR)
      r_sig2[e8].b = r_sig[e7].b + (zC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, zL, zC, zR)
    end
  end
end

task Step1b_sigz_z(r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_sig, prz_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(s.x, s.y, s.z, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[1]
    var KL : int32 = bc[2]
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

      var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
      var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}
      var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}
     
      var gsigL : double
      var gsigR : double
      var bsigL : double
      var bsigR : double

      if r_sig.bounds.lo.z == s.z then
        gsigL = plz_sig[eL7].g 
        bsigL = plz_sig[eL7].b 
        zL = plz_mesh[eL3].z
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        zL = r_mesh[eL3].z
      end

      if r_sig.bounds.hi.z == s.z then
        gsigR = prz_sig[eR7].g
        bsigR = prz_sig[eR7].b
        zR = prz_mesh[eR3].z
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        zR = r_mesh[eR3].z
      end

      r_sig2[e8].g = r_sig[e7].g + (zC/2.0)*VanLeer(gsigL, r_sig[e7].g, gsigR, zL, zC, zR)
      r_sig2[e8].b = r_sig[e7].b + (zC/2.0)*VanLeer(bsigL, r_sig[e7].b, bsigR, zL, zC, zR)
    end
  end
end


task Step1b_b(r_sig : region(ispace(int7d), grid),
              r_mesh : region(ispace(int3d), mesh),
              r_gridbarp : region(ispace(int6d), grid),
              r_gridbarpb : region(ispace(int7d), grid),
              prx_gridbarp : region(ispace(int6d), grid),
              pry_gridbarp : region(ispace(int6d), grid),
              prz_gridbarp : region(ispace(int6d), grid),
              prx_sig : region(ispace(int7d), grid),
              pry_sig : region(ispace(int7d), grid),
              prz_sig : region(ispace(int7d), grid),
	      prx_mesh : region(ispace(int3d), mesh),
	      pry_mesh : region(ispace(int3d), mesh),
	      prz_mesh : region(ispace(int3d), mesh),
              vxmesh : region(ispace(int1d), vmesh),
              vymesh : region(ispace(int1d), vmesh),
              vzmesh : region(ispace(int1d), vmesh),
              BCs : int32[3], N : int32[3], effD : int32)
where
  reads writes(r_gridbarpb),
  reads(r_sig, r_gridbarp, r_mesh.{dx,dy,dz}, vxmesh, vymesh, vzmesh),
  reads(prx_gridbarp, pry_gridbarp, prz_gridbarp, prx_sig, pry_sig, prz_sig, prx_mesh.dx, pry_mesh.dy, prz_mesh.dz)
do

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var rhotest : double = 0

  --c.printf("slo = {%d, %d, %d}, shi = {%d, %d, %d}\n", slo.x, slo.y, slo.z, shi.x, shi.y, shi.z)
  for s in s3 do

    var e3 : int3d = {s.x, s.y, s.z}

    for Dim = 0, effD do

      -- Gathering Right Indices
      var bc : int32[6] = BC(s.x, s.y, s.z, Dim, BCs, N)
      var IR : int32 = bc[3]
      var JR : int32 = bc[4]
      var KR : int32 = bc[5]
      var eR3 : int3d = {IR, JR, KR}

      for v in v3 do      

        var e6 : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}
        var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}

        -- Gathering Right Indices 
        var eR6 : int6d = {IR, JR, KR, v.x, v.y, v.z}
        var eR7 : int7d = {IR, JR, KR, Dim, v.x, v.y, v.z}

        -- Dot Product is just a single product when using rectangular mesh
        var swap : double = 1.0
        var gb : double = r_gridbarp[e6].g
        var bb : double = r_gridbarp[e6].b
        var gsig : double = r_sig[e7].g
        var bsig : double = r_sig[e7].b
      
        var sC : double[3] 
        sC[0] = r_mesh[e3].dx
        sC[1] = r_mesh[e3].dy
        sC[2] = r_mesh[e3].dz

        if (vxmesh[v.x].v < 0 and Dim == 0) then
          swap = -1
          if s.x == r_mesh.bounds.hi.x then
            --c.printf("eR7 = {%d, %d, %d, %d, %d, %d, %d}\n", s.x, eR7.y, eR7.z, eR7.w, eR7.v, eR7.u, eR7.t)
            gsig = prx_sig[eR7].g
            bsig = prx_sig[eR7].b
            gb = prx_gridbarp[eR6].g
            bb = prx_gridbarp[eR6].b
            sC[Dim] = prx_mesh[eR3].dx
           else
            gsig = r_sig[eR7].g
            bsig = r_sig[eR7].b
            gb = r_gridbarp[eR6].g
            bb = r_gridbarp[eR6].b
            sC[Dim] = r_mesh[eR3].dx
           end
        elseif (vymesh[v.y].v < 0 and Dim == 1) then
          gsig = pry_sig[eR7].g
          bsig = pry_sig[eR7].b
          gb = pry_gridbarp[eR6].g
          bb = pry_gridbarp[eR6].b
          swap = -1
          sC[Dim] = pry_mesh[eR3].dy
        elseif (vzmesh[v.z].v < 0 and Dim == 2) then
          gsig = prz_sig[eR7].g
          bsig = prz_sig[eR7].b
          gb = prz_gridbarp[eR6].g
          bb = prz_gridbarp[eR6].b
          swap = -1
          sC[Dim] = prz_mesh[eR3].dz
        end
     
        if s.x == 128 and v.x == 128 then
          --c.printf("gb = %f, dgb = %f\n", gb, swap*sC[Dim]/2.0*gsig)
        end

        -- TODO need to change sC to sR when swap, doesnt currently matter for sod/KHI/RTI bc dx_i = dx_0
        r_gridbarpb[e7].g = gb + swap*sC[Dim]/2.0*gsig
        r_gridbarpb[e7].b = bb + swap*sC[Dim]/2.0*bsig
        --c.printf("Interpolating g[%d][%d] {%f, %f}\n", s.x, v.x, gb, r_gridbarpb[e7].g)

        if s.x == 129 then
          rhotest += vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*r_gridbarpb[e7].g
        end

        -- NAN checker
        if (isnan(r_gridbarpb[e7].g) == 1 or isnan(r_gridbarpb[e7].b) == 1) then
 
          c.printf("Step 1b: r_gridbarp.g = %f, r_gridbarp.b = %f, r_sig.g = %f, r_sig.b = %f\n", gb, bb, gsig, bsig)

          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].g)), "Step 1b\n")
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].b)), "Step 1b\n")
    
        end

      end
    end 
  end
  c.printf("Step1b rhotest[129] = %f\n", rhotest)
end

-- Step 1c: Compute phibar at interface by interpolating w/ phisigma2, x-Xi*dt/2
task Step1c(r_gridbarpb : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),        
            vymesh : region(ispace(int1d), vmesh),        
            vzmesh : region(ispace(int1d), vmesh),        
            r_sigb : region(ispace(int8d), grid),
            dt : double, BCs : int32[3],  N : int32[3], effD : int32)
where
  reads(vxmesh, vymesh, vzmesh, r_sigb), 
  reads writes(r_gridbarpb)
do     
  -- Compute gbar/bbar @ t=n+1/2  with interface sigma
  var Xi : double[3]
  
  var slo : int3d = {r_gridbarpb.bounds.lo.x, r_gridbarpb.bounds.lo.y, r_gridbarpb.bounds.lo.z}
  var shi : int3d = {r_gridbarpb.bounds.hi.x, r_gridbarpb.bounds.hi.y, r_gridbarpb.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var rhotest : double = 0
  for v in v3 do
    Xi[0] = vxmesh[v.x].v
    Xi[1] = vymesh[v.y].v
    Xi[2] = vzmesh[v.z].v
    for s in s3 do
      for Dim2 = 0, effD do
        for Dim = 0, effD do
          -- Gather Indices
          var e8 : int8d = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
          var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}

          r_gridbarpb[e7].g = r_gridbarpb[e7].g - dt/2.0*Xi[Dim]*r_sigb[e8].g
          r_gridbarpb[e7].b = r_gridbarpb[e7].b - dt/2.0*Xi[Dim]*r_sigb[e8].b
    
          if s.x == 129 then rhotest += r_gridbarpb[e7].g*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w end
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].g)), "Step 1c\n")
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].b)), "Step 1c\n")
        
        end 
      end 
    end
  end 
  --c.printf("Step1c rhotest[129] = %f\n", rhotest)
end

--Step 2: Microflux
--Step 2a: Compute W at interface.
task Step2a(r_gridbarpb : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            r_Wb   : region(ispace(int4d), W),
            dt : double, effD : int32)
where
  reads(r_gridbarpb, vxmesh, vymesh, vzmesh),
  reads writes(r_Wb)
do    
  var slo : int3d = {r_Wb.bounds.lo.x, r_Wb.bounds.lo.y, r_Wb.bounds.lo.z}
  var shi : int3d = {r_Wb.bounds.hi.x, r_Wb.bounds.hi.y, r_Wb.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Reset field space
  fill(r_Wb.rho, 0)

  -- First do density at boundary, density is needed for others.
  for s in s3 do
    for Dim = 0, effD do
      var e4 : int4d = {s.x, s.y, s.z, Dim}
      for v in v3 do
        var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
        r_Wb[e4].rho = r_Wb[e4].rho + vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*r_gridbarpb[e7].g
      end
    end
  end
      
  -- Then do momentum and energy
  -- Initialize, then iterate over contributions in velocity space
  var U : double[3]
  for s in s3 do
    for Dim = 0, effD do

      var e4 : int4d = {s.x, s.y, s.z, Dim}  
      for d = 0, effD do
        r_Wb[e4].rhov[d] = dt/2.0*r_Wb[e4].rho*0 -- TODO: In future, replace 0 with acceleration field 
      end
      r_Wb[e4].rhoE = dt/2.*r_Wb[e4].rho*0 -- TODO: In future replace 0 with u.dot(a), vel dot acc

      for v in v3 do

        U[0] = vxmesh[v.x].v
        U[1] = vymesh[v.y].v
        U[2] = vzmesh[v.z].v
 

        var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}    

        for d = 0, effD do
          r_Wb[e4].rhov[d] = r_Wb[e4].rhov[d] + vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*U[d]*r_gridbarpb[e7].g 
        end

        for d = effD, 3 do
          r_Wb[e4].rhov[d] = 0 -- Bug Preventer
        end
          
        r_Wb[e4].rhoE = r_Wb[e4].rhoE + vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*r_gridbarpb[e7].b
      end
    end
  end

  -- NAN checker
  for e in r_Wb do
    
    --c.printf("r_Wb[%d] = {%f, %f, %f}\n", e.x, r_Wb[e].rho, r_Wb[e].rhov[0], r_Wb[e].rhoE)
    regentlib.assert(not [bool](isnan(r_Wb[e].rho)), "Step 2a rho\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[0])), "Step 2a rhov0\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[1])), "Step 2a rhov1\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[2])), "Step 2a rhov2\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhoE)), "Step 2a rhoE\n")
    
  end

end

-- Step 2b: compute original phi at interface using gbar, W at interface
-- Memory Recycling: phibar @ interface is used to store phi @ interface.
task Step2b(r_gridbarpb : region(ispace(int7d), grid), 
            r_Wb      : region(ispace(int4d), W),
            vxmesh    : region(ispace(int1d), vmesh),
            vymesh    : region(ispace(int1d), vmesh),
            vzmesh    : region(ispace(int1d), vmesh),
            dt : double, R : double, K : double, Cv : double, g : double,
            w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads writes(r_gridbarpb),
  reads(r_Wb, vxmesh, vymesh, vzmesh)
do
  var slo : int3d = {r_Wb.bounds.lo.x, r_Wb.bounds.lo.y, r_Wb.bounds.lo.z}
  var shi : int3d = {r_Wb.bounds.hi.x, r_Wb.bounds.hi.y, r_Wb.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)
  

  for s in s3 do

    var tg : double
    var tb : double
    var u : double
    var T : double
    var Xi : double[3]
    var e3 : int3d = {s.x, s.y, s.z}

    for Dim = 0, effD do

      var e4 : int4d = {s.x, s.y, s.z, Dim}

      u = 0 
      for d = 0, effD do
        u = u + r_Wb[e4].rhov[d]/r_Wb[e4].rho*r_Wb[e4].rhov[d]/r_Wb[e4].rho
      end
      u = sqrt(u)
      
      T = Temperature(r_Wb[e4].rhoE/r_Wb[e4].rho, u, g, R)
      if T < 0 then
        c.printf("T < 0, r_Wb[{%d, %d, %d, %d}].rhoE = %f, r_Wb[e4].rho = %f, u = %f, g = %f, R = %f\n", e4.x, e4.y, e4.z, e4.w, r_Wb[e4].rhoE, r_Wb[e4].rho, u, g, R)
        regentlib.assert(T >= 0, "Negative Temperature\n")
      end

      tg = visc(T, ur, Tr, w)/r_Wb[e4].rho/R/T
      tb = tg/Pr
      --c.printf("tg = %f, tb = %f\n", tg, tb)

      for v in v3 do

        var e : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
        var c2 : double = 0
        var g_eq : double
        var b_eq : double

        Xi[0] = vxmesh[v.x].v
        Xi[1] = vymesh[v.y].v
        Xi[2] = vzmesh[v.z].v

        for d = 0, effD do
          c2 =  c2 + (Xi[d] - r_Wb[e4].rhov[d]/r_Wb[e4].rho)*(Xi[d] - r_Wb[e4].rhov[d]/r_Wb[e4].rho)
        end

        g_eq = geq(c2, r_Wb[e4].rho, T, R, effD)
        b_eq = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0

        if (isnan(g_eq) == 1 or isnan(b_eq) == 1) then

          c.printf("c2 = %f, r_Wb[e4].rho = %f, T = %f, R = %f, effD = %d\n", c2, r_Wb[e4].rho, T, R, effD) 
      
          regentlib.assert(not  [bool](isnan(g_eq)), "Step 2b\n")
          regentlib.assert(not  [bool](isnan(b_eq)), "Step 2b\n")
    
        end

        -- this is actually the original distribution function, recycling memory from gbar
        r_gridbarpb[e].g = 2*tg/(2*tg + dt/2.)*r_gridbarpb[e].g + dt/(4*tg + dt)*g_eq + dt*tg/(4*tg + dt)*0 -- TODO replace this last *0 with source term 
        r_gridbarpb[e].b = 2*tb/(2*tb + dt/2.)*r_gridbarpb[e].b + dt/(4*tb + dt)*b_eq + dt*tb/(4*tb + dt)*0 -- TODO replace this last *0 with source term

        if (isnan(r_gridbarpb[e].g) == 1 or isnan(r_gridbarpb[e].b) == 1) then

          c.printf("gbar = %.10f, bbar = %.10f, g_eq = %.10f, tg = %.10f, tb = %.10f\n", r_gridbarpb[e].g, r_gridbarpb[e].b, g_eq, tg, tb)
      
          regentlib.assert(not [bool](isnan(r_gridbarpb[e].g)), "Step 2b\n")
          regentlib.assert(not [bool](isnan(r_gridbarpb[e].b)), "Step 2b\n")
    
        end

      end
    end 
  end
end

-- Step 2c: Compute Microflux F at interface at half timestep using W/phi at interface.
task Step2c(r_gridbarpb : region(ispace(int7d), grid),
            r_F       : region(ispace(int6d), grid),
            r_mesh    : region(ispace(int3d), mesh),
            vxmesh    : region(ispace(int1d), vmesh),
            vymesh    : region(ispace(int1d), vmesh),
            vzmesh    : region(ispace(int1d), vmesh),
            plx_gridbarpb : region(ispace(int7d), grid),
            ply_gridbarpb : region(ispace(int7d), grid),
            plz_gridbarpb : region(ispace(int7d), grid),
            BCs : int32[3], R : double, K : double, Cv : double, g : double,
            w : double, Tr : double, Pr : double, effD : int32, N : int32[3])
where
  reads(r_gridbarpb, vxmesh, vymesh, vzmesh, r_mesh, plx_gridbarpb, ply_gridbarpb, plz_gridbarpb),
  reads writes(r_F)
do    
  var A : double[3]
  var Xi : double[3]

  fill(r_F.g, 0)
  fill(r_F.b, 0)
  __fence(__execution, __block)

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    var e3 : int3d = {s.x, s.y, s.z}

    A[0] = r_mesh[e3].dy*r_mesh[e3].dz
    A[1] = r_mesh[e3].dx*r_mesh[e3].dz
    A[2] = r_mesh[e3].dx*r_mesh[e3].dy


    for Dim = 0, effD do
      var bc : int32[6] =  BC(s.x, s.y, s.z, Dim, BCs, N)
      var IL : int32 = bc[0]
      var IR : int32 = bc[1]
      var JL : int32 = bc[2]
      var JR : int32 = bc[3]
      var KL : int32 = bc[4]
      var KR : int32 = bc[5]

      -- Boundary Conditions
      var right : double = 1.0 
      var left : double = 1.0

      -- Dirichlet Boundary Conditions
      if (Dim == 0 and BCs[0] == 1 and (s.x == 0 or s.x == N[0] - 1)) then
        left = 0
        right = 0
      end 
      if (Dim == 1 and BCs[1] == 1 and (s.y == 0 or s.y == N[1] - 1)) then
        left = 0
        right = 0
      end
      if (Dim == 2 and BCs[2] == 1 and (s.z == 0 or s.z == N[2] - 1)) then
        left = 0
        right = 0
      end

      -- Neumann Boundary Conditions
      if (Dim == 0 and BCs[0] == 1) then
        if s.x == 0 then 
          left = 0
        elseif s.x == N[0] - 1 then 
          right = 0
        end
      end 
      if (Dim == 1 and BCs[1] == 1) then
        if s.y == 0 then 
          left = 0
        elseif s.y == N[1] - 1 then 
          right = 0
        end
      end
      if (Dim == 2 and BCs[2] == 1) then
        if s.z == 0 then 
          left = 0
        elseif s.z == N[2] - 1 then 
          right = 0
        end
      end

      var eL3 : int3d = {IL, JL, KL}
      
      for v in v3 do 
        -- Gather Left Indices
        var e7 : int7d = {s.x, s.y, s.z, Dim, v.x, v.y, v.z}
        var e6 : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}
        var eL : int6d = {IL, JL, KL, v.x, v.y, v.z}
        var eL7 : int7d = {IL, JL, KL, Dim, v.x, v.y, v.z}

        Xi[0] = vxmesh[v.x].v
        Xi[1] = vymesh[v.y].v
        Xi[2] = vzmesh[v.z].v
  
        var gL : double
        var bL : double

        if Dim == 0 then
          if s.x == s3.bounds.lo.x then
            gL = plx_gridbarpb[eL7].g
            bL = plx_gridbarpb[eL7].b
          else
            gL = r_gridbarpb[eL7].g
            bL = r_gridbarpb[eL7].b
          end
        end

        if Dim == 1 then
          if s.y == s3.bounds.lo.y then
            gL = ply_gridbarpb[eL7].g
            bL = ply_gridbarpb[eL7].b
          else
            gL = r_gridbarpb[eL7].g
            bL = r_gridbarpb[eL7].b
          end
        end

        if Dim == 2 then
          if s.z == s3.bounds.lo.z then
            gL = plz_gridbarpb[eL7].g
            bL = plz_gridbarpb[eL7].b
          else
            gL = r_gridbarpb[eL7].g
            bL = r_gridbarpb[eL7].b
          end
        end
           
        r_F[e6].g = r_F[e6].g + Xi[Dim]*A[Dim]*(right*r_gridbarpb[e7].g - left*gL)
        r_F[e6].b = r_F[e6].b + Xi[Dim]*A[Dim]*(right*r_gridbarpb[e7].b - left*bL)

        if (v.x == 128) then
          --c.printf("r_F[%d]= {%f, %f), right/left = {%f, %f}, gL/gR = {%f, %f}, bL/bR = {%f, %f}\n", s.x, r_F[e6].g, r_F[e6].b, right, left, gL, r_gridbarpb[e7].g, bL, r_gridbarpb[e7].b)
        end
        regentlib.assert(not [bool](isnan(r_F[e6].g)), "Step 2c\n")
        regentlib.assert(not [bool](isnan(r_F[e6].b)), "Step 2c\n")
      end 
    end
  end 
end

--Step 3: Source Terms
task Step3()
  -- 
end

--Step 4: Update Conservative Variables W at cell center at next timestep
--Step 5: Update Phi at cell center at next time step
task Step4and5(r_grid : region(ispace(int6d), grid),
               r_W    : region(ispace(int3d), W),
               r_mesh : region(ispace(int3d), mesh),
               r_F    : region(ispace(int6d), grid),
               vxmesh : region(ispace(int1d), vmesh),
               vymesh : region(ispace(int1d), vmesh),
               vzmesh : region(ispace(int1d), vmesh),
               dt : double, BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
               g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(vxmesh, vymesh, vzmesh, r_mesh, r_F),
  reads writes(r_W, r_grid)
do
  var V : double      -- Volume of Cell
  var Xi : double[3]  -- Discrete Velocity 
  var uo : double     -- Old Flow Velocity
  var To : double     -- Old Temp
  var tgo : double    -- Old visc
  var tbo : double    -- Old visc
  var u : double      -- Flow Velocity
  var T : double      -- Temperature 
  var tg : double     -- visc
  var tb : double     -- visc 
  var c2 : double     -- peculiar velocity sq
  var g_eq : double   -- Equils
  var b_eq : double
  var g_eqo : double
  var b_eqo : double

  var slo : int3d = {r_F.bounds.lo.x, r_F.bounds.lo.y, r_F.bounds.lo.z}
  var shi : int3d = {r_F.bounds.hi.x, r_F.bounds.hi.y, r_F.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  

  var gold : double[256]
  var gmid : double[256]
  var gnew : double[256]
  var bold : double[256]
  var bmid : double[256]
  var bnew : double[256]
  
  for s in s3 do

    var e3 : int3d = {s.x, s.y, s.z}
    
    V = r_mesh[e3].dx*r_mesh[e3].dy*r_mesh[e3].dz

    -- Compute old flow velocity
    uo = 0 -- old flow velocity
    for d = 0, effD do
      uo += r_W[e3].rhov[d]/r_W[e3].rho*r_W[e3].rhov[d]/r_W[e3].rho
    end
    uo = sqrt(uo)
    regentlib.assert(bool(uo>=0), "uo")

    -- Compute old temperature
    To = Temperature(r_W[e3].rhoE/r_W[e3].rho, uo, g, R)
    regentlib.assert(bool(To>=0), "To")
  
    -- Compute old taus
    tgo = visc(To, ur, Tr, w)/r_W[e3].rho/R/To
    tbo = tgo/Pr 
  
    for v in v3 do 
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v
      var e6 : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}

      -- Compute old eq's
      c2 = 0
      for d = 0, effD do
        c2 += (Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)*(Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)
      end
      g_eqo = geq(c2, r_W[e3].rho, To, R, effD)
      b_eqo = g_eqo*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3-effD+K)*R*To)/2

      -- First Update phi
      var i : int32 = s.x
      var j : int32 = s.y
      var k : int32 = s.z

      if v.x == 128 then 
        gold[s.x] = r_grid[e6].g
        bold[s.x] = r_grid[e6].b
      end
      -- Update First Step (terms involving oldW)
      if ((BCs[0] == 1 and i == 0) or (BCs[0] == 1 and i == N[0] - 1) or
          (BCs[1] == 1 and j == 0 and effD > 1) or (BCs[1] == 1 and j == N[1] - 1 and effD > 1) or
          (BCs[2] == 1 and k == 0 and effD > 2) or (BCs[2] == 1 and k == N[2] - 1 and effD > 2)) then
       
        r_grid[e6].g = r_grid[e6].g
        r_grid[e6].b = r_grid[e6].b

        --c.printf("Not updating : s = {%d, %d, %d}\n", s.x, s.y, s.z)
      else
        r_grid[e6].g = r_grid[e6].g + dt/2.0*(g_eqo-r_grid[e6].g)/tgo - dt/V*r_F[e6].g + dt*0 -- TODO replace 0 with source term
        r_grid[e6].b = r_grid[e6].b + dt/2.0*(b_eqo-r_grid[e6].b)/tbo - dt/V*r_F[e6].b + dt*0 -- TODO replace 0 with source term
      end

      if v.x == 128 then 
        gmid[s.x] = r_grid[e6].g
        bmid[s.x] = r_grid[e6].b
      end

    end
  end
  __fence(__execution, __block)


  for s in s3 do
    var e3 : int3d = {s.x, s.y, s.z}
    var drho : double = 0
    var drhov : double = 0
    var dE : double = 0

    for v in v3 do 
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v

      var e6 : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}
      -- Step 4: Update W at cell center
      drho  += -dt*(r_F[e6].g/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w
      drhov += -dt*(r_F[e6].g/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*Xi[0]
      dE    += -dt*(r_F[e6].b/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w
      r_W[e3].rho = r_W[e3].rho - dt*(r_F[e6].g/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w -- TODO replace 0 with source term.
      for d = 0, effD do
        r_W[e3].rhov[d] = r_W[e3].rhov[d] - dt*Xi[d]*(r_F[e6].g/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w -- TODO replace 0 with source term
      end
  
      r_W[e3].rhoE = r_W[e3].rhoE - dt*(r_F[e6].b/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w -- TODO replace 0 with source term      
    end
    --c.printf("updated W[%d] = {%f, %f, %f}, dW = {%f, %f, %f}\n", s.x, r_W[e3].rho, r_W[e3].rhov[0], r_W[e3].rhoE, drho, drhov, dE)
  end     

  -- Second Update Phi at cell center using new tau/W
  -- Compute flow velocity u

  for s in s3 do
    var e3 : int3d = {s.x, s.y, s.z}

    u = 0 
    for d = 0, effD do 
      u += r_W[e3].rhov[d]/r_W[e3].rho*r_W[e3].rhov[d]/r_W[e3].rho
    end
    u = sqrt(u)
    regentlib.assert(bool(u>=0), "u")

    -- Compute T
    T = Temperature(r_W[e3].rhoE/r_W[e3].rho, u, g, R)
    if T < 0 then
      c.printf("r_W[%d].rhoE = %f, r_W[e3].rho = %f, u = %f\n", e3.x, r_W[e3].rhoE, r_W[e3].rho, u)
    end
    regentlib.assert(bool(T>=0), "T")
  
    -- Compute new taus
    tg = visc(T, ur, Tr, w)/r_W[e3].rho/R/T 
    tb = tg/Pr 
    -- printf("taus are = {%f, %f}\n", tg, tb)

    for v in v3 do
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v

      -- Compute new eq's
      c2 = 0  -- reset c2 from before
      for d = 0, effD do
        c2 += (Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)*(Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)
      end
      g_eq = geq(c2, r_W[e3].rho, T, R, effD)
      b_eq = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3-effD+K)*R*T)/2


      -- Second Update phi
      var i : int32 = s.x
      var j : int32 = s.y
      var k : int32 = s.z
      var e6 : int6d = {s.x, s.y, s.z, v.x, v.y, v.z} 

      if ((BCs[0] == 1 and i == 0) or (BCs[0] == 1 and i == N[0] - 1) or
          (BCs[1] == 1 and j == 0 and effD > 1) or (BCs[1] == 1 and j == N[1] - 1 and effD > 1) or
          (BCs[2] == 1 and k == 0 and effD > 2) or (BCs[2] == 1 and k == N[2] - 1 and effD > 2)) then
                
        r_grid[e6].g = r_grid[e6].g
        r_grid[e6].b = r_grid[e6].b
        --c.printf("Not updating : s = {%d, %d, %d}\n", s.x, s.y, s.z)
      else
        r_grid[e6].g = (r_grid[e6].g + dt/2.0*g_eq/tg)/(1+dt/2.0/tg)
        r_grid[e6].b = (r_grid[e6].b + dt/2.0*b_eq/tb)/(1+dt/2.0/tb)

        --r_grid[e6].g = r_grid[e6].g/(1+dt/2.0/tg)
        --r_grid[e6].b = r_grid[e6].b/(1+dt/2.0/tb)
      end

      if v.x == 128 then 
        gnew[s.x] = r_grid[e6].g
        bnew[s.x] = r_grid[e6].b
        --c.printf("phi old/mid/new[%d]: g = {%f, %f, %f}, b = {%f, %f, %f}, phieq = {%f, %f}\n", s.x, gold[s.x], gmid[s.x], gnew[s.x], bold[s.x], bmid[s.x], bnew[s.x], g_eq, b_eq)
      end

      if isnan(r_grid[e6].g) == 1 then
        c.printf("Step4and5: g_eq = %f, tg = %f, tgo = %f, r_F[e6].g = %f\n", g_eq, tg, tgo, r_F[e6].g)
      end 
      if isnan(r_grid[e6].b) == 1 then
        c.printf("Step4and5: b_eq = %f, tb = %f, tbo = %f, r_F[e6].b = %f\n", b_eq, tb, tbo, r_F[e6].b)
      end 
      regentlib.assert(not [bool](isnan(r_grid[e6].g)), "Step4and5\n")
      regentlib.assert(not [bool](isnan(r_grid[e6].b)), "Step4and5\n")
    end
  end --TODO Indentation
end



task MaxwellianInitialization(r_grid  : region(ispace(int6d), grid),
         r_mesh : region(ispace(int3d), mesh),
         r_W    : region(ispace(int3d), W),
         vxmesh : region(ispace(int1d), vmesh),
         vymesh : region(ispace(int1d), vmesh),
         vzmesh : region(ispace(int1d), vmesh),
         testProblem: int32, R : double, K : double, Cv : double, 
         g : double, w : double, ur : double, Tr : double, Pr : double,
         N : int32[3], NV : int32[3], effD : int32)
where 
  reads writes(r_grid, r_mesh, r_W, vxmesh, vymesh, vzmesh)
do
  var T : double
  var u : double 

  var rhotest : double = 0
  var Etest : double = 0

  var slo : int3d = {r_grid.bounds.lo.x, r_grid.bounds.lo.y, r_grid.bounds.lo.z}
  var shi : int3d = {r_grid.bounds.hi.x, r_grid.bounds.hi.y, r_grid.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do

    u = 0
    for dim = 0, effD do
      u += r_W[s].rhov[dim]/r_W[s].rho*r_W[s].rhov[dim]/r_W[s].rho
    end
    u = sqrt(u)
      
    T = Temperature(r_W[s].rhoE/r_W[s].rho, u, g, R)
    
    for v in v3 do
    
      var e : int6d = {s.x, s.y, s.z, v.x, v.y, v.z}

      var c2 : double = 0
      var Xi : double[3]  
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v

      for d = 0, effD do
        c2 += (Xi[d] - r_W[s].rhov[d]/r_W[s].rho)*(Xi[d] - r_W[s].rhov[d]/r_W[s].rho) --TODO could be bugged
      end

      r_grid[e].g = geq(c2, r_W[s].rho, T, R, effD)
      r_grid[e].b = r_grid[e].g*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0
  
      if s.x == 0 and s.y == 0 and s.z == 0 then
        rhotest += r_grid[e].g*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w
        Etest += r_grid[e].b*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w
      end
    end
  end
  c.printf("rhotest = %f, Etest = %f\n", rhotest, Etest)
end

task InitializeGrid(r_grid  : region(ispace(int6d), grid),
                r_mesh : region(ispace(int3d), mesh),
		r_W    : region(ispace(int3d), W),
		vxmesh : region(ispace(int1d), vmesh),
		vymesh : region(ispace(int1d), vmesh),
		vzmesh : region(ispace(int1d), vmesh),
                testProblem: int32, R : double, K : double, Cv : double, g : double, w : double, ur : double, Tr : double, Pr : double,
		N : int32[3], NV : int32[3], effD : int32)
where
  reads writes(r_grid, r_mesh, r_W, vxmesh, vymesh, vzmesh)
do

  -- InitializeTestProblem
  if testProblem == 0 then
    -- TODO : User Defined
  elseif testProblem > 0 or testProblem == -1 then
    c.printf("Maxwellian Initialization\n")
    MaxwellianInitialization(r_grid, r_mesh, r_W, vxmesh, vymesh, vzmesh, testProblem, R, K, Cv, g, w, ur, Tr, Pr, N, NV, effD)
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

task factorize1d(parallelism : int) : int3d
  var sizex = parallelism
  return int3d {sizex, 1, 1} 
end

task factorize2d(parallelism : int) : int3d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x : int32 = 1
  var size_y : int32 = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y =  i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int3d { size_x, size_y, 1 }
end

task factorize3d(parallelism : int) : int3d
  var limit = [int](cmath.pow([double](parallelism), 1/3))
  var size_x = 1
  var size_y = 1
  var size_z = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_z = i, parallelism / i
      if size_x > size_z then
        size_x, size_z = size_z, size_x
      end
    end

    if parallelism % i*i == 0 then
      size_x, size_y, size_z = i, i, parallelism/i*i
    end
  end
  return int3d { size_x, size_y, size_z }
end

task factorize(parallelism: int, effD : int32)

  var f6 : int6d = {1, 1, 1, 1, 1, 1}
  if effD == 1 then
    var f3 = factorize1d(parallelism)
    f6.x, f6.y, f6.z = f3.x, f3.y, f3.z    
  elseif effD == 2 then
    var f3 = factorize2d(parallelism)
    f6.x, f6.y, f6.z = f3.x, f3.y, f3.z    
  elseif effD == 3 then
    var f3 = factorize3d(parallelism)
    f6.x, f6.y, f6.z = f3.x, f3.y, f3.z    
  end

  return f6
end

terra wait_for(x : int) return 1 end
task block_task(r_image : region(ispace(int1d), particle))
where
  reads writes(r_image)
do
  return 1
end


task Dump(r_W : region(ispace(int3d), W), iter : int32)
where
  reads (r_W)
do
  var filename : int8[1000]
  c.sprintf([&int8](filename), './Data/rho_%04d',iter)
  var g = c.fopen(filename,'wb')

  for e in r_W do
    dumpdouble(g, r_W[e].rho)
  end
  __fence(__execution, __block)
  return 1
end



task toplevel()
  var config : Config
  config:initialize_from_command() -- TODO : CPUs, output bool

  -- Simulation Parameters
  var testProblem : int32 = config.testproblem
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
  var Tf : double = r_params[0].Tf
  var dtdump : double = r_params[0].dtdump

  __fence(__execution, __block) 
  c.printf("Simulation Parameters\n")
  if testProblem > 0 then c.printf("testProblem = %d\n", testProblem) end
  c.printf("N = {%d, %d, %d}, NV = {%d, %d, %d}, effD = %d\n", N[0], N[1], N[2], NV[0], NV[1], NV[2], effD)
  c.printf("BCs = {%d, %d, %d}, Vmin = {%f, %f, %f}, Vmax = {%f, %f, %f}\n", BCs[0], BCs[1], BCs[2], Vmin[0], Vmin[1], Vmin[2], Vmax[0], Vmax[1], Vmax[2])
  c.printf("R = %f, K = %f, g = %f, Cv = %f\n", R, K, g, Cv)
  c.printf("w = %f, ur = %f, Tr = %f, Pr = %f\n", w, ur, Tr, Pr)

  -- Create regions for distribution functions and gradients
  var r_grid      = region(ispace(int6d, {N[0], N[1], N[2], NV[0], NV[1], NV[2]}), grid)
  var r_gridbarp  = region(ispace(int6d, {N[0], N[1], N[2], NV[0], NV[1], NV[2]}), grid)
  var r_gridbarpb = region(ispace(int7d, {N[0], N[1], N[2], effD, NV[0], NV[1], NV[2]}), grid)
  var r_sig       = region(ispace(int7d, {N[0], N[1], N[2], effD, NV[0], NV[1], NV[2]}), grid)
  var r_sig2      = region(ispace(int8d, {N[0], N[1], N[2], effD, effD, NV[0], NV[1], NV[2]}), grid)
  var r_sigb      = region(ispace(int8d, {N[0], N[1], N[2], effD, effD, NV[0], NV[1], NV[2]}), grid)
 
  -- Create regions for mesh and conserved variables (cell center and interface)
  var r_mesh = region(ispace(int3d, {N[0], N[1], N[2]}), mesh)
  var r_W    = region(ispace(int3d, {N[0], N[1], N[2]}), W)
  var r_Wb   = region(ispace(int4d, {N[0], N[1], N[2], effD}), W)
 
  -- Create regions for velocity space and initialize
  var vxmesh = region(ispace(int1d, NV[0]), vmesh) 
  var vymesh = region(ispace(int1d, NV[1]), vmesh) 
  var vzmesh = region(ispace(int1d, NV[2]), vmesh) 
  NewtonCotes(vxmesh, vymesh, vzmesh, NV, Vmin, Vmax)

  -- Create regions for source terms and flux
  var r_S = region(ispace(int6d, {N[0], N[1], N[2], NV[0], NV[1], NV[2]}), grid)
  var r_F = region(ispace(int6d, {N[0], N[1], N[2], NV[0], NV[1], NV[2]}), grid)


  -- Create partitions for regions
  var f6 : int6d = factorize(config.cpus, effD)
  
  var f3 : int3d = {f6.x, f6.y, f6.z}
  var f4 : int4d = {f6.x, f6.y, f6.z, 1}
  var f7 : int7d = {f6.x, f6.y, f6.z, 1, f6.w, f6.v, f6.u}
  var f8 : int8d = {f6.x, f6.y, f6.z, 1, 1, f6.w, f6.v, f6.u}
  var p6 = ispace(int6d, f6)
  var p7 = ispace(int7d, f7)
  var p8 = ispace(int8d, f8)
  var p3 = ispace(int3d, f3)
  var p4 = ispace(int4d, f4)
  var p_grid = partition(equal, r_grid, p6)
  var p_gridbarp = partition(equal, r_gridbarp, p6)
  var p_gridbarpb = partition(equal, r_gridbarpb, p7)
  var p_sig = partition(equal, r_sig, p7)
  var p_sig2 = partition(equal, r_sig2, p8)
  var p_sigb = partition(equal, r_sigb, p8)
  var p_mesh = partition(equal, r_mesh, p3)
  var p_W = partition(equal, r_W, p3)
  var p_Wb = partition(equal, r_Wb, p4)
  var p_S = partition(equal, r_S, p6)
  var p_F = partition(equal, r_F, p6)
  __fence(__execution, __block)
  c.printf("Equal Partitions Done\n")

  -- Create coloring for partitions for left/right ghost regions
  var c3Lx = coloring.create()
  var c6Lx = coloring.create()
  var c7Lx = coloring.create()
  var c8Lx = coloring.create()
  var c3Ly = coloring.create()
  var c6Ly = coloring.create()
  var c7Ly = coloring.create()
  var c8Ly = coloring.create()
  var c3Lz = coloring.create()
  var c6Lz = coloring.create()
  var c7Lz = coloring.create()
  var c8Lz = coloring.create()

  var c3Rx = coloring.create()
  var c6Rx = coloring.create()
  var c7Rx = coloring.create()
  var c8Rx = coloring.create()
  var c3Ry = coloring.create()
  var c6Ry = coloring.create()
  var c7Ry = coloring.create()
  var c8Ry = coloring.create()
  var c3Rz = coloring.create()
  var c6Rz = coloring.create()
  var c7Rz = coloring.create()
  var c8Rz = coloring.create()
 
  -- Create Rects for colorings for partitions
  for col7 in p_gridbarpb.colors do
    var bounds = p_gridbarpb[col7].bounds
    
    -- Leftmost and Rightmost indices
    var il : int32 = bounds.lo.x
    var jl : int32 = bounds.lo.y
    var kl : int32 = bounds.lo.z
    var ir : int32 = bounds.hi.x
    var jr : int32 = bounds.hi.y
    var kr : int32 = bounds.hi.z

    var lx = BC(il, jl, kl, 0, BCs, N)[0]
    var rx = BC(ir, jr, kr, 0, BCs, N)[3]
    var ly = BC(il, jl, kl, 1, BCs, N)[1]
    var ry = BC(ir, jr, kr, 1, BCs, N)[4]
    var lz = BC(il, jl, kl, 2, BCs, N)[2]
    var rz = BC(ir, jr, kr, 2, BCs, N)[5]

    var col3 : int3d = {col7.x, col7.y, col7.z}
    var col6 : int6d = {col7.x, col7.y, col7.z, col7.v, col7.u, col7.t}
    var col8 : int8d = {col7.x, col7.y, col7.z, col7.w, 0, col7.v, col7.u, col7.t}

    -- for reference -- terra BC(i : int32, j : int32, k : int32, Dim : int32, BCs : int32[3], N : int32[3])

    var rleftx3 : rect3d = { {lx, bounds.lo.y, bounds.lo.z}, 
                          {lx, bounds.hi.y, bounds.hi.z}}
    var rlefty3 : rect3d = { {bounds.lo.x, ly, bounds.lo.z}, 
                          {bounds.hi.x, ly, bounds.hi.z}}
    var rleftz3 : rect3d = { {bounds.lo.x, bounds.lo.y, lz}, 
                          {bounds.hi.x, bounds.hi.y, lz}}

    var rrightx3 : rect3d = { {rx, bounds.lo.y, bounds.lo.z}, 
                          {rx, bounds.hi.y, bounds.hi.z}}
    var rrighty3 : rect3d = { {bounds.lo.x, ry, bounds.lo.z}, 
                          {bounds.hi.x, ry, bounds.hi.z}}
    var rrightz3 : rect3d = { {bounds.lo.x, bounds.lo.y, rz}, 
                          {bounds.hi.x, bounds.hi.y, rz}}
    __fence(__execution, __block)
    c.printf("rect3d Done\n")

    var rleftx6 : rect6d = { {lx, bounds.lo.y, bounds.lo.z, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {lx, bounds.hi.y, bounds.hi.z, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rlefty6 : rect6d = { {bounds.lo.x, ly, bounds.lo.z, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                        {bounds.hi.x,  ly, bounds.hi.z, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rleftz6 : rect6d = { {bounds.lo.x, bounds.lo.y, lz, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {bounds.hi.x, bounds.hi.y, lz, bounds.hi.v, bounds.hi.u, bounds.hi.t}}

    var rrightx6 : rect6d = { {rx, bounds.lo.y, bounds.lo.z, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {rx, bounds.hi.y, bounds.hi.z, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rrighty6 : rect6d = { {bounds.lo.x, ry, bounds.lo.z, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, ry, bounds.hi.z, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rrightz6 : rect6d = { {bounds.lo.x, bounds.lo.y, rz, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, rz, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    __fence(__execution, __block)
    c.printf("rect6d Done\n")

    var rleftx7 : rect7d = { {lx, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {lx, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rlefty7 : rect7d = { {bounds.lo.x, ly, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {bounds.hi.x, ly, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rleftz7 : rect7d = { {bounds.lo.x, bounds.lo.y, lz, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {bounds.hi.x, bounds.hi.y, lz, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}

    c.printf("Making rrightx7 : rx = %d\n", rx)
    var rrightx7 : rect7d = { {rx, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {rx, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rrighty7 : rect7d = { {bounds.lo.x, ry, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, ry, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rrightz7 : rect7d = { {bounds.lo.x, bounds.lo.y, rz, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, rz, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}

    __fence(__execution, __block)
    c.printf("rect7d Done\n")

    var rleftx8 : rect8d = { {lx, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {lx, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rlefty8 : rect8d = { {bounds.lo.x, ly, bounds.lo.z, bounds.lo.w, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {bounds.hi.x, ly, bounds.hi.z, bounds.hi.w, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rleftz8 : rect8d = { {bounds.lo.x, bounds.lo.y, lz, bounds.lo.w, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                         {bounds.hi.x, bounds.hi.y, lz, bounds.hi.w, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}

    var rrightx8 : rect8d = { {rx, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {rx, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rrighty8 : rect8d = { {bounds.lo.x, ry, bounds.lo.z, bounds.lo.w, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, ry, bounds.hi.z, bounds.hi.w, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}
    var rrightz8 : rect8d = { {bounds.lo.x, bounds.lo.y, rz, bounds.lo.w, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, rz, bounds.hi.w, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t}}


    __fence(__execution, __block)
    c.printf("rect8d Done\n")

    -- Color in left strips
    coloring.color_domain(c3Lx, col3, rleftx3)
    coloring.color_domain(c3Ly, col3, rlefty3)
    coloring.color_domain(c3Lz, col3, rleftz3)
    __fence(__execution, __block)
    c.printf("3d Coloring Done\n")

    coloring.color_domain(c6Lx, col6, rleftx6)
    coloring.color_domain(c6Ly, col6, rlefty6)
    coloring.color_domain(c6Lz, col6, rleftz6)
    __fence(__execution, __block)
    c.printf("6d Coloring Done\n")

    coloring.color_domain(c7Lx, col7, rleftx7)
    coloring.color_domain(c7Ly, col7, rlefty7)
    coloring.color_domain(c7Lz, col7, rleftz7)
    __fence(__execution, __block)
    c.printf("7d Coloring Done\n")

    coloring.color_domain(c8Lx, col8, rleftx8)
    coloring.color_domain(c8Ly, col8, rlefty8)
    coloring.color_domain(c8Lz, col8, rleftz8)
    __fence(__execution, __block)
    c.printf("8d Coloring Done\n")


    -- Color in right strips
    coloring.color_domain(c3Rx, col3, rrightx3)
    coloring.color_domain(c3Ry, col3, rrighty3)
    coloring.color_domain(c3Rz, col3, rrightz3)

    coloring.color_domain(c6Rx, col6, rrightx6)
    coloring.color_domain(c6Ry, col6, rrighty6)
    coloring.color_domain(c6Rz, col6, rrightz6)

    coloring.color_domain(c7Rx, col7, rrightx7)
    coloring.color_domain(c7Ry, col7, rrighty7)
    coloring.color_domain(c7Rz, col7, rrightz7)

    coloring.color_domain(c8Rx, col8, rrightx8)
    coloring.color_domain(c8Ry, col8, rrighty8)
    coloring.color_domain(c8Rz, col8, rrightz8)
    __fence(__execution, __block)
    c.printf("Coloring Done\n")
  end

  -- Create Partitions
  var plx_mesh = partition(disjoint, r_mesh, c3Lx, p3)
  var ply_mesh = partition(disjoint, r_mesh, c3Ly, p3)
  var plz_mesh = partition(disjoint, r_mesh, c3Lz, p3)
  var prx_mesh = partition(disjoint, r_mesh, c3Rx, p3)
  var pry_mesh = partition(disjoint, r_mesh, c3Ry, p3)
  var prz_mesh = partition(disjoint, r_mesh, c3Rz, p3)
  __fence(__execution, __block)
  c.printf("Mesh Strips Done\n")
  
  var plx_gridbarp = partition(disjoint, r_gridbarp, c6Lx, p6)
  __fence(__execution, __block)
  c.printf("one gridbarp Strips Done\n")
  var ply_gridbarp = partition(disjoint, r_gridbarp, c6Ly, p6)
  var plz_gridbarp = partition(disjoint, r_gridbarp, c6Lz, p6)
  var prx_gridbarp = partition(disjoint, r_gridbarp, c6Rx, p6)
  var pry_gridbarp = partition(disjoint, r_gridbarp, c6Ry, p6)
  var prz_gridbarp = partition(disjoint, r_gridbarp, c6Rz, p6)
  __fence(__execution, __block)
  c.printf("gridbarp Strips Done\n")

  var plx_gridbarpb = partition(disjoint, r_gridbarpb, c7Lx, p7)
  __fence(__execution, __block)
  c.printf("one gridbarpb Strips Done\n")
  var ply_gridbarpb = partition(disjoint, r_gridbarpb, c7Ly, p7)
  var plz_gridbarpb = partition(disjoint, r_gridbarpb, c7Lz, p7)
  var prx_gridbarpb = partition(disjoint, r_gridbarpb, c7Rx, p7)
  var pry_gridbarpb = partition(disjoint, r_gridbarpb, c7Ry, p7)
  var prz_gridbarpb = partition(disjoint, r_gridbarpb, c7Rz, p7)
  __fence(__execution, __block)
  c.printf("gridbarpbb Strips Done\n")

  var plx_sig = partition(disjoint, r_sig, c7Lx, p7)
  var ply_sig = partition(disjoint, r_sig, c7Ly, p7)
  var plz_sig = partition(disjoint, r_sig, c7Lz, p7)
  var prx_sig = partition(disjoint, r_sig, c7Rx, p7)
  var pry_sig = partition(disjoint, r_sig, c7Ry, p7)
  var prz_sig = partition(disjoint, r_sig, c7Rz, p7)
  __fence(__execution, __block)
  c.printf("sig Strips Done\n")


  var plx_sig2 = partition(disjoint, r_sig2, c8Lx, p8)
  var ply_sig2 = partition(disjoint, r_sig2, c8Ly, p8)
  var plz_sig2 = partition(disjoint, r_sig2, c8Lz, p8)
  var prx_sig2 = partition(disjoint, r_sig2, c8Rx, p8)
  var pry_sig2 = partition(disjoint, r_sig2, c8Ry, p8)
  var prz_sig2 = partition(disjoint, r_sig2, c8Rz, p8)
  __fence(__execution, __block)
  c.printf("sig2 Strips Done\n")
  c.printf("Left/Right Strip Partitioning Done\n")

  var plx_sigb = partition(disjoint, r_sigb, c8Lx, p8)
  var ply_sigb = partition(disjoint, r_sigb, c8Ly, p8)
  var plz_sigb = partition(disjoint, r_sigb, c8Lz, p8)
  var prx_sigb = partition(disjoint, r_sigb, c8Rx, p8)
  var pry_sigb = partition(disjoint, r_sigb, c8Ry, p8)
  var prz_sigb = partition(disjoint, r_sigb, c8Rz, p8)
  __fence(__execution, __block)
  c.printf("sigb Strips Done\n")
  c.printf("Left/Right Strip Partitioning Done\n")

  --Initialize r_mesh
  var MeshType : int32 = 1
  InitializeMesh(r_mesh, N, MeshType) --TODO Needs more input for nested, user-def etc.
  __fence(__execution, __block)
  c.printf("Mesh Initialized\n")

  --Initialize r_W
  __demand(__parallel)
  for col3 in p_W.colors do
    InitializeW(p_W[col3], p_mesh[col3], N, NV, testProblem, R, Cv)
  end
  __fence(__execution, __block)
  c.printf("W Initialized\n")
  
  --Initialize r_grid
  __demand(__parallel)
  for col6 in p_grid.colors do
    var col3 : int3d = {col6.x, col6.y, col6.z}
    InitializeGrid(p_grid[col6], p_mesh[col3], p_W[col3], vxmesh, vymesh, vzmesh, testProblem, R, K, Cv, g, w, ur, Tr, Pr, N, NV, effD)
  end
  __fence(__execution, __block)
  c.printf("Grid Initialized\n")

  --Timestep
  var CFL : double = 0.95 -- Safety Factor
  var dxmin : double = 1.0/cmath.fmax(cmath.fmax(N[0],N[1]),N[2]) -- Smallest Cell Width (TODO : Non-Uniform Meshes)
  var umax : double  = 4.0 -- Estimated maximum flow velocity, TODO calculate at each iteration for stronger problems
  var calcdt : double = CFL*dxmin/(umax + sqrt(Vmax[0]*Vmax[0] + Vmax[1]*Vmax[1] + Vmax[2]*Vmax[2]))
  
  var Tsim : double = 0.0  -- Sim time
  var Tdump : double = 0.0 -- Time since last dump 
  
  var iter : int32 = 0
  var dumpiter : int32 = 0
  if testProblem > 0 then 
    Dump(r_W, dumpiter) -- Initial Conditions
    __fence(__execution, __block)
    c.printf("Dump %d\n", dumpiter)
  end
  while Tsim < Tf do --and iter < 3 do
    iter += 1

    var dt = TimeStep(calcdt, dtdump-Tdump, Tf-Tsim)

    __fence(__execution, __block)
    --c.printf("Starting Step1a\n")
    -- Step 1a
    __demand(__parallel)
    for col6 in p_grid.colors do 
      var col3 : int3d = {col6.x, col6.y, col6.z}
      Step1a(p_grid[col6], p_gridbarp[col6], p_S[col6], p_W[col3], vxmesh, vymesh, vzmesh, dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    end
    __fence(__execution, __block)
    --c.printf("Step1a Complete\n")

    -- Step 1b: Compute Gradient Sigma
    __demand(__parallel)
    for col6 in p_grid.colors do 
      var col3 : int3d = {col6.x, col6.y, col6.z}
      var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
      Step1b_sigx(p_gridbarp[col6], p_sig[col7], p_mesh[col3], plx_mesh[col3], prx_mesh[col3], plx_gridbarp[col6], prx_gridbarp[col6], vxmesh, vymesh, vzmesh, BCs, N, effD)
    end
    if effD > 1 then
      __demand(__parallel)
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        Step1b_sigy(p_gridbarp[col6], p_sig[col7], p_mesh[col3], ply_mesh[col3], pry_mesh[col3], ply_gridbarp[col6], pry_gridbarp[col6], vxmesh, vymesh, vzmesh, BCs, N, effD)
      end  
    end
    if effD > 2 then
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        Step1b_sigz(p_gridbarp[col6], p_sig[col7], p_mesh[col3], plz_mesh[col3], prz_mesh[col3], plz_gridbarp[col6], prz_gridbarp[col6], vxmesh, vymesh, vzmesh, BCs, N, effD)
      end  
    end
    __fence(__execution, __block)
    --c.printf("Sig Complete\n")

    -- Step 1b: Compute Gradient of Gradient Sigma, Sigma2
    for col6 in p_grid.colors do 
      var col3 : int3d = {col6.x, col6.y, col6.z}
      var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
      var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
      Step1b_sigx_x(p_sig[col7], p_sig2[col8], p_mesh[col3], plx_mesh[col3], prx_mesh[col3], plx_sig[col7], prx_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
    end  
    for col6 in p_grid.colors do 
      var col3 : int3d = {col6.x, col6.y, col6.z}
      var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
      var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
      Step1b_sigx_x2(p_sig[col7], p_sig2[col8], p_sigb[col8], p_mesh[col3], prx_mesh[col3], prx_sig[col7], prx_sig2[col8], vxmesh, vymesh, vzmesh, BCs, N)
    end  
    if effD > 1 then
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigx_y(p_sig[col7], p_sig2[col8], p_mesh[col3], ply_mesh[col3], pry_mesh[col3], ply_sig[col7], pry_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigy_y(p_sig[col7], p_sig2[col8], p_mesh[col3], ply_mesh[col3], pry_mesh[col3], ply_sig[col7], pry_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigy_x(p_sig[col7], p_sig2[col8], p_mesh[col3], plx_mesh[col3], prx_mesh[col3], plx_sig[col7], prx_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
    end
    if effD > 2 then
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.z, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigz_x(p_sig[col7], p_sig2[col8], p_mesh[col3], plx_mesh[col3], prx_mesh[col3], plx_sig[col7], prx_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigz_y(p_sig[col7], p_sig2[col8], p_mesh[col3], ply_mesh[col3], pry_mesh[col3], ply_sig[col7], pry_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigz_z(p_sig[col7], p_sig2[col8], p_mesh[col3], plz_mesh[col3], prz_mesh[col3], plz_sig[col7], prz_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigx_z(p_sig[col7], p_sig2[col8], p_mesh[col3], plz_mesh[col3], prz_mesh[col3], plz_sig[col7], prz_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
      for col6 in p_grid.colors do 
        var col3 : int3d = {col6.x, col6.y, col6.z}
        var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
        var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
        Step1b_sigy_z(p_sig[col7], p_sig2[col8], p_mesh[col3], plz_mesh[col3], prz_mesh[col3], plz_sig[col7], prz_sig[col7], vxmesh, vymesh, vzmesh, BCs, N)
      end  
    end
    __fence(__execution, __block)
    --c.printf("Sig2 Complete\n")

    -- Step 1b_b: Interpolate to velocity dependent time in past
    __demand(__parallel)
    for col7 in p_gridbarpb.colors do
      var col3 : int3d = {col7.x, col7.y, col7.z}
      var col6 : int6d = {col7.x, col7.y, col7.z, col7.v, col7.u, col7.t}
      var col8 : int8d = {col7.x, col7.y, col7.z, col7.w, 0, col7.v, col7.u, col7.t}
      Step1b_b(p_sig[col7], p_mesh[col3], p_gridbarp[col6], p_gridbarpb[col7], prx_gridbarp[col6], pry_gridbarp[col6], prz_gridbarp[col6], prx_sig[col7], pry_sig[col7], prz_sig[col7], prx_mesh[col3], pry_mesh[col3], prz_mesh[col3],vxmesh, vymesh, vzmesh, BCs, N, effD)
    end
    __fence(__execution, __block)
    --c.printf("Step1b_b Complete\n")


    -- Step 1c
    __demand(__parallel)
    for col6 in p_grid.colors do
      var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
      var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
    
      Step1c(p_gridbarpb[col7], vxmesh, vymesh, vzmesh, p_sigb[col8], dt, BCs, N, effD)
    end
    __fence(__execution, __block)
    --c.printf("Step1c Complete\n")

    -- Step 2a
    __demand(__parallel)
    for col7 in p_gridbarpb.colors do
      var col4 : int4d = {col7.x, col7.y, col7.z, 0}
      Step2a(p_gridbarpb[col7], vxmesh, vymesh, vzmesh, p_Wb[col4], dt, effD)
    end

    -- Step 2b
    __demand(__parallel)
    for col7 in p_gridbarpb.colors do
      var col4 : int4d = {col7.x, col7.y, col7.z, 0}
      Step2b(p_gridbarpb[col7], p_Wb[col4], vxmesh, vymesh, vzmesh, dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    end

    -- Step 2c
    __demand(__parallel)
    for col7 in p_gridbarpb.colors do
      var col3 : int3d = {col7.x, col7.y, col7.z}
      var col6 : int6d = {col7.x, col7.y, col7.z, col7.v, col7.u, col7.t}
      Step2c(p_gridbarpb[col7], p_F[col6], p_mesh[col3], vxmesh, vymesh, vzmesh, plx_gridbarpb[col7], ply_gridbarpb[col7], plz_gridbarpb[col7], BCs, R, K, Cv, g, w, Tr, Pr, effD, N)
    end

    -- Step 3
    Step3() -- TODO
  
    -- Step 4 and 5
    __demand(__parallel)
    for col6 in p_grid.colors do
      var col3 : int3d = {col6.x, col6.y, col6.z}
      Step4and5(p_grid[col6], p_W[col3], p_mesh[col3], p_F[col6], vxmesh, vymesh, vzmesh, dt, BCs, R, K, Cv, N, g, w, ur, Tr, Pr, effD)
    end
    __fence(__execution, __block)

    if dt < calcdt then
      dumpiter += 1
      Tdump = 0
      Dump(r_W, dumpiter)
      c.printf("Dump %d\n", dumpiter)
    else
      Tdump += dt
    end

    Tsim += dt

    c.printf("Iteration = %d, Tsim = %f\n", iter, Tsim)

  end
end

regentlib.start(toplevel)
