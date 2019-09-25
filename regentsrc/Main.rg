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

		var num : int32 = 256
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
                r_params[e].Nv  = r_params[e].NV[0]*r_params[e].NV[1]*r_params[e].NV[2]

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
                r_params[e].Nv  = r_params[e].NV[0]*r_params[e].NV[1]*r_params[e].NV[2]


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

 
      for i = 0, Nx do
        for j = 0, Ny do
          for k = 0, Nz do

            var idx : int3d = {i, j, k}

            r_mesh[idx].x = dX*(i + 1.0/2.0)
            r_mesh[idx].y = dY*(j + 1.0/2.0)
            r_mesh[idx].z = dZ*(k + 1.0/2.0)

            r_mesh[idx].dx = dX
            r_mesh[idx].dy = dY
            r_mesh[idx].dz = dZ

          end
        end
      end
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
  end
  if (Dim == 1 and BCs[1] == 0) then
    IL = i
	IR = i
	JL = (j - 1 + N[1])%N[1]
    JR = (j + 1)%N[1]
    KL = k
	KR = k
  end
  if (Dim == 2 and BCs[2] == 0) then
    IL = i
	IR = i
	JL = j
	JR = j
	KL = (k - 1 + N[2])%N[2]
    KR = (k + 1)%N[2]
  end


  -- Dirichlet Boundary Conditions
  if (Dim == 0 and BCs[0] == 1) then
    IL = i - 1
    IR = i + 1
    if IL <  0 then IL = 0 end
    if IR == N[0] then IR = N[0] - 1 end
    JL = j
	JR = j
	KL = k
	KR = k
  end
  if (Dim == 1 and BCs[1] == 1) then
    IL = i
	IR = i
	JL = j - 1
    JR = j + 1
    if JL <  0 then JL = 0 end
    if JR == N[1] then JR = N[1] - 1 end
    KL = k
	KR = k
  end
  if (Dim == 2 and BCs[2] == 1) then
    IL = i
	IR = i
	JL = j
	JR = j
	KL = k - 1
    KR = k + 1
    if KL <  0 then KL = 0 end
    if KR == N[2] then KR = N[2] - 1 end
  end


  -- Neumann Boundary Conditions
  if (Dim == 0 and BCs[0] == 2) then
    IL = i - 1
    IR = i + 1
    if IL <  0 then IL = 0 end
    if IR == N[0] then IR = N[0] - 1 end
    JL = j
	JR = j
	KL = k
	KR = k
  end
  if (Dim == 1 and BCs[1] == 2) then
    IL = i
	IR = i
	JL = j - 1
    JR = j + 1
    if JL <  0 then JL = 0 end
    if JR == N[1] then JR = N[1] - 1 end
    KL = k
	KR = k
  end
  if (Dim == 2 and BCs[2] == 2) then
    IL = i
	IR = i
	JL = j
    JR = j
   	KL = k - 1
   	KR = k + 1
    if KL <  0 then KL = 0 end
    if KR == N[2] then KR = N[2] - 1 end

    var LR : int32[6] 
    LR[0] = IL
    LR[1] = JL
    LR[2] = KL
    LR[3] = IR
    LR[4] = JR
    LR[5] = KR

    return LR
  end
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

  for e in r_grid do
    var sidx : int3d = {e.w, e.v, e.u}
 
    u = 0
    for dim = 0, effD do
      u += r_W[sidx].rhov[dim]/r_W[sidx].rho*r_W[sidx].rhov[dim]/r_W[sidx].rho
    end
    u = sqrt(u)

    T = Temperature(r_W[sidx].rhoE/r_W[sidx].rho, u, g, R)

    tg = visc(T, ur, Tr, w)/r_W[sidx].rho/R/T
    tb = tg/Pr

    -- For Now...
    r_S[e].g = 0.
    r_S[e].b = 0.

    c2 = 0
    Xi[0] = vxmesh[e.x].v
    Xi[1] = vymesh[e.y].v
    Xi[2] = vzmesh[e.z].v
    for dim = 0, effD do
      c2 += (Xi[dim]-r_W[sidx].rhov[dim]/r_W[sidx].rho)*(Xi[dim]-r_W[sidx].rhov[dim]/r_W[sidx].rho)
    end

    var g_eq : double = geq(c2, r_W[sidx].rho, T, R, effD)
    var b_eq : double = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0

    r_gridbarp[e].g = (2*tg - dt/2.)/(2.*tg)*r_grid[e].g + dt/(4.*tg)*g_eq + dt/4.*r_S[e].g
    r_gridbarp[e].b = (2*tb - dt/2.)/(2.*tb)*r_grid[e].b + dt/(4.*tb)*b_eq + dt/4.*r_S[e].b

    if (isnan(r_gridbarp[e].g) == 1 or isnan(r_gridbarp[e].b) == 1) then

      c.printf("Step 1a: gp = %.12f, bp = %.12f, g = %.12f, b = %.12f, g_eq = %.12f, Sg = %.12f, Sb = %.12f, taus = {%.12f, %.12f}\n", r_gridbarp[e].g, r_gridbarp[e].b, r_grid[e].g, r_grid[e].b, g_eq, r_S[e].g, r_S[e].b, tg, tb)
    
    end
    

    regentlib.assert(not [bool](isnan(r_gridbarp[e].g)), "Step 1a\n")
    regentlib.assert(not [bool](isnan(r_gridbarp[e].b)), "Step 1a\n")

  end
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
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3], 
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where 
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_gridbarp, prx_gridbarp),
  writes(r_sig)
do
  var Dim : int32 = 0 -- change when copy

  -- Cell Centers 
  var xC : double
  var xL : double
  var xR : double

  var olds : int3d = {-1, -1, -1}  
  
  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})
  
  for s in s3 do
    xC = r_mesh[s].x -- change when copy
  
    var i : int32 = s.x
    var j : int32 = s.y
    var k : int32 = s.z

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim, BCs, N)
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
      r_sig[e7].b = VanLeer(bbpL, r_gridbarp[e].b, gbpR, xL, xC, xR)

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
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3], 
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where 
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, ply_mesh, pry_mesh, ply_gridbarp, pry_gridbarp),
  writes(r_sig)
do
  var Dim : int32 = 1 -- change when copy

  -- Cell Centers 
  var yC : double
  var yL : double
  var yR : double

  var olds : int3d = {-1, -1, -1}  
  
  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})
  
  for s in s3 do
  
    yC = r_mesh[s].y -- change when copy

    var i : int32 = s.x
    var j : int32 = s.y
    var k : int32 = s.z

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim, BCs, N)
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
      r_sig[e7].b = VanLeer(bbpL, r_gridbarp[e].b, gbpR, yL, yC, yR)

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
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_gridbarp : region(ispace(int6d), grid),
            prz_gridbarp : region(ispace(int6d), grid),
            vxmesh : region(ispace(int1d), vmesh),            
            vymesh : region(ispace(int1d), vmesh),            
            vzmesh : region(ispace(int1d), vmesh),            
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3], 
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where 
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_gridbarp, prz_gridbarp),
  writes(r_sig)
do
  var Dim : int32 = 2 -- change when copy

  -- Cell Centers 
  var zC : double
  var zL : double
  var zR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})
  
  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    var i : int32 = s.x
    var j : int32 = s.y
    var k : int32 = s.z

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim, BCs, N)
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
      r_sig[e7].b = VanLeer(bbpL, r_gridbarp[e].b, gbpR, zL, zC, zR)

      -- NAN checker
      if (isnan(r_sig[e7].g) == 1 or isnan(r_sig[e7].b) == 1) then

        c.printf("Step 1b_sigz: r_sig.g = %f, r_sig.b = %f\n", r_sig[e7].g, r_sig[e7].b)

        regentlib.assert(not [bool](isnan(r_sig[e7].g)), "Step 1b_sigz\n")
        regentlib.assert(not [bool](isnan(r_sig[e7].b)), "Step 1b_sigz\n")

      end
  
    end
  end
end

      
task Step1b_sigx_x(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plx_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_gridbarp, prx_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    xC = r_mesh[s].x -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
end

task Step1b_sigy_x(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plx_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_gridbarp, prx_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    xC = r_mesh[s].x -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
end

task Step1b_sigz_x(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plx_mesh : region(ispace(int3d), mesh),
            prx_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plx_mesh, prx_mesh, plx_gridbarp, prx_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    xC = r_mesh[s].x -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
end

task Step1b_sigy_y(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            ply_mesh : region(ispace(int3d), mesh),
            pry_mesh : region(ispace(int3d), mesh),
            ply_sig : region(ispace(int7d), grid),
            pry_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, ply_mesh, pry_mesh, ply_gridbarp, pry_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    yC = r_mesh[s].y -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
end

task Step1b_sigz_y(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            ply_mesh : region(ispace(int3d), mesh),
            pry_mesh : region(ispace(int3d), mesh),
            ply_sig : region(ispace(int7d), grid),
            pry_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, ply_mesh, pry_mesh, ply_gridbarp, pry_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    yC = r_mesh[s].y -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
end

task Step1b_sigx_z(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_gridbarp, prz_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
end

task Step1b_sigy_z(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_gridbarp, prz_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
end

task Step1b_sigx_z(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            plz_mesh : region(ispace(int3d), mesh),
            prz_mesh : region(ispace(int3d), mesh),
            plz_sig : region(ispace(int7d), grid),
            prz_sig : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh, plz_mesh, prz_mesh, plz_gridbarp, prz_gridbarp),
  writes(r_sig)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var s3 = ispace(int3d, r_mesh.bounds.lo, r_mesh_bounds.hi)
  var v3 = ispace(int3d, {vxmesh.bounds.lo.x, vymesh.bounds.lo.x, vzmesh.bounds.lo.x}, {vxmesh.bounds.hi.x, vymesh.bounds.hi.x, vzmesh.bounds.hi.x})

  for s in s3 do
    zC = r_mesh[s].z -- change when copy

    -- Gather Left and Right Indices
    var bc : int32[6] = BC(i, j, k, Dim2, BCs, N)
    var IL : int32 = bc[0]
    var JL : int32 = bc[0]
    var KL : int32 = bc[0]
    var IR : int32 = bc[0]
    var JR : int32 = bc[0]
    var KR : int32 = bc[0]

    var eL3 : int3d = {IL, JL, KL}
    var eR3 : int3d = {IR, JR, KR}

    for v in v3 do

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

        var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
        var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}
     
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
              BCs : int32[3], N : int32[3])
where
  reads writes(r_gridbarpb),
  reads(r_sig, r_gridbarp, r_mesh.{dx,dy,dz}, vxmesh.v, vymesh.v, vzmesh.v),
  reads(prx_gridbarp, pry_gridbarp, prz_gridbarp, prx_sig, pry_sig, prz_sig, prx_mesh.dx, pry_mesh.dy, prz_mesh.dz)
do
  for e in r_gridbarpb do
    var e3 : int3d = {e.v, e.u, e.t}
    var e6 : int6d = {e.x, e.y, e.z, e.v, e.u, e.t}

    -- Gathering Right Indices -- (Note: e.w = Dim)
    var bc : int32[6] = BC(e.v, e.u, e.t, e.w, BCs, N)
    var IR : int32 = bc[3]
    var JR : int32 = bc[4]
    var KR : int32 = bc[5]
    var eR : int6d = {e.x, e.y, e.z, IR, JR, KR}
    var eR3 : int3d = {IR, JR, KR}
    var eR6 : int6d = {e.x, e.y, e.z, IR, JR, KR}
    var eR7 : int7d = {e.x, e.y, e.z, e.w, IR, JR, KR}

    -- Dot Product is just a single product when using rectangular mesh
    var swap : double = 1.0
    var gb : double = r_gridbarp[e6].g
    var bb : double = r_gridbarp[e6].b
    var gsig : double = r_sig[e].g
    var bsig : double = r_sig[e].b
      
    var sC : double[3] 
    sC[0] = r_mesh[e3].dx
    sC[1] = r_mesh[e3].dy
    sC[2] = r_mesh[e3].dz
    -- Note : e.w = Dim
    if (vxmesh[e.x].v < 0 and e.w == 0) then
      gsig = prx_sig[eR7].g
      bsig = prx_sig[eR7].b
      gb = prx_gridbarp[eR6].g
      bb = prx_gridbarp[eR6].b
      swap = -1
      sC[e.w] = prx_mesh[eR3].dx
    elseif (vymesh[e.y].v < 0 and e.w == 1) then
      gsig = pry_sig[eR7].g
      bsig = pry_sig[eR7].b
      gb = pry_gridbarp[eR6].g
      bb = pry_gridbarp[eR6].b
      swap = -1
      sC[e.w] = pry_mesh[eR3].dy
    elseif (vzmesh[e.z].v < 0 and e.w == 2) then
      gsig = prz_sig[eR7].g
      bsig = prz_sig[eR7].b
      gb = prz_gridbarp[eR6].g
      bb = prz_gridbarp[eR6].b
      swap = -1
      sC[e.w] = prz_mesh[eR3].dz
    end

     
    -- TODO need to change sC to sR when swap, doesnt currently matter for sod/KHI/RTI bc dx_i = dx_0
    r_gridbarpb[e].g = gb + swap*sC[e.w]/2.0*gsig
    r_gridbarpb[e].b = bb + swap*sC[e.w]/2.0*bsig

    -- NAN checker
    if (isnan(r_gridbarpb[e].g) == 1 or isnan(r_gridbarpb[e].b) == 1) then

      c.printf("Step 1b: r_gridbarp.g = %f, r_gridbarp.b = %f, r_sig.g = %f, r_sig.b = %f\n", gb, bb, gsig, bsig)

      regentlib.assert(not [bool](isnan(r_gridbarpb[e].g)), "Step 1b\n")
      regentlib.assert(not [bool](isnan(r_gridbarpb[e].b)), "Step 1b\n")
    
    end

  end
end

-- Step 1c: Compute phibar at interface by interpolating w/ phisigma2, x-Xi*dt/2
task Step1c(r_gridbar : region(ispace(int7d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),        
            vymesh : region(ispace(int1d), vmesh),        
            vzmesh : region(ispace(int1d), vmesh),        
            r_sig2 : region(ispace(int8d), grid),
            dt : double, BCs : int32[3], R : double, K : double, Cv : double, N : int32[3], 
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(vxmesh, vymesh, vzmesh, r_sig2, r_gridbarpb),
  reads writes(r_gridbar)
do     
  -- Compute gbar/bbar @ t=n+1/2  with interface sigma
  var Xi : double[3] 

  for e in r_gridbar do
    Xi[0] = vxmesh[e.x].v
    Xi[1] = vxmesh[e.y].v
    Xi[2] = vxmesh[e.z].v
    

    -- Phibar at Interface, at t = n+1/2
    r_gridbar[e].g = r_gridbarpb[e].g
    r_gridbar[e].b = r_gridbarpb[e].b

    var dg : double = 0
    var db : double = 0

    for Dim2 = 0, effD do
      var IL2 : int32 
      var IR2 : int32
      var JL2 : int32
      var JR2 : int32
      var KL2 : int32
      var KR2 : int32
   
      var i : int32 = e.v
      var j : int32 = e.u
      var k : int32 = e.t

      -- Periodic Boundary Conditions
      if (Dim2 == 0 and BCs[0] == 0) then 
        IL2 = (i - 1 + N[0])%N[0] 
        IR2 = (i + 1)%N[0] 
        JL2 = j
        JR2 = j
        KL2 = k
        KR2 = k
      end 
      if (Dim2 == 1 and BCs[1] == 0) then
        IL2 = i
        IR2 = i
        JL2 = (j - 1 + N[1])%N[1]
        JR2 = (j + 1)%N[1] 
        KL2 = k
        KR2 = k
      end
      if (Dim2 == 2 and BCs[2] == 0) then
        IL2 = i
        IR2 = i
        JL2 = j
        JR2 = j
        KL2 = (k - 1 + N[2])%N[2] 
        KR2 = (k + 1)%N[2] 
      end


      -- Dirichlet Boundary Conditions
      if (Dim2 == 0 and BCs[0] == 1) then
        IL2 = i - 1 
        IR2 = i + 1 
        if IL2 <  0 then IL2 = 0 end 
        if IR2 == N[0] then IR2 = N[0] - 1 end
        JL2 = j
        JR2 = j
        KL2 = k
        KR2 = k
      end 
      if (Dim2 == 1 and BCs[1] == 1) then
        IL2 = i
        IR2 = i
        JL2 = j - 1 
        JR2 = j + 1
        if JL2 <  0 then JL2 = 0 end 
        if JR2 == N[1] then JR2 = N[1] - 1 end
        KL2 = k
        KR2 = k
      end
      if (Dim2 == 2 and BCs[2] == 1) then
        IL2 = i
        IR2 = i
        JL2 = j
        JR2 = j
        KL2 = k - 1 
        KR2 = k + 1 
        if KL2 <  0 then KL2 = 0 end 
        if KR2 == N[2] then KR2 = N[2] - 1 end
      end
    
    
      -- Neumann Boundary Conditions
      if (Dim2 == 0 and BCs[0] == 2) then
        IL2 = i - 1 
        IR2 = i + 1 
        if IL2 <  0 then IL2 = 0 end 
        if IR2 == N[0] then IR2 = N[0] - 1 end
        JL2 = j
        JR2 = j
        KL2 = k
        KR2 = k
      end 
      if (Dim2 == 1 and BCs[1] == 2) then
        IL2 = i
        IR2 = i
        JL2 = j - 1 
        JR2 = j + 1
        if JL2 <  0 then JL2 = 0 end 
        if JR2 == N[1] then JR2 = N[1] - 1 end
        KL2 = k
        KR2 = k
      end
      if (Dim2 == 2 and BCs[2] == 2) then
        IL2 = i
        IR2 = i
        JL2 = j
        JR2 = j
        KL2 = k - 1 
        KR2 = k + 1 
        if KL2 <  0 then KL2 = 0 end 
        if KR2 == N[2] then KR2 = N[2] - 1 end
      end

      -- Gather Left and Right Indices
      var eL2 : int6d = {e.x, e.y, e.z, IL2, JL2, KL2}
      var eR2 : int6d = {e.x, e.y, e.z, IR2, JR2, KR2}
      var eL2_3 : int3d = {IL2, JL2, KL2}
      var eR2_3 : int3d = {IR2, JR2, KR2}
      var eL2_7 : int7d = {e.x, e.y, e.z, Dim2, IL2, JL2, KL2}
      var eR2_7 : int7d = {e.x, e.y, e.z, Dim2, IR2, JR2, KR2}



      var interpidx : int7d = e
      var swap : double = 1.0

      if (vxmesh[e.x].v < 0 and Dim2 == 0) then
        interpidx = eR2_7
        swap = -1
      elseif (vymesh[e.y].v < 0 and Dim2 == 1) then
        interpidx = eR2_7 
        swap = -1
      elseif (vzmesh[e.z].v < 0 and Dim2 == 2) then
        interpidx = eR2_7 
        swap = -1
      end


      var interpol8 : int8d = {interpidx.x, interpidx.y, interpidx.z, interpidx.w, Dim2, interpidx.v, interpidx.u, interpidx.t} 
         
      r_gridbar[e].g = r_gridbar[e].g - dt/2.0*Xi[Dim2]*r_sig2[interpol8].g
      r_gridbar[e].b = r_gridbar[e].b - dt/2.0*Xi[Dim2]*r_sig2[interpol8].b
    
      regentlib.assert(not [bool](isnan(r_gridbar[e].g)), "Step 1c\n")
      regentlib.assert(not [bool](isnan(r_gridbar[e].b)), "Step 1c\n")

    end
  end 
end

--Step 2: Microflux
--Step 2a: Interpolate W to interface.
task Step2a(r_gridbar : region(ispace(int7d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            r_Wb   : region(ispace(int4d), W),
            dt : double, R : double, K : double, Cv : double, g : double,
            w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(r_gridbar, vxmesh, vymesh, vzmesh),
  reads writes(r_Wb)
do    
  -- Compute conserved variables W at t+1/2
  
  -- Reset field space
  fill(r_Wb.rho, 0)
 
  -- First do density at boundary, density is needed for others.
  for e in r_gridbar do
    var e4 : int4d = {e.w, e.v, e.u, e.t}  
    r_Wb[e4].rho = r_Wb[e4].rho + vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w*r_gridbar[e].g
  end
      
  -- Then do momentum and energy
  -- First initialize w/ source terms
  for e in r_Wb do
    for v = 0, effD do
      r_Wb[e].rhov[v] = dt/2.0*r_Wb[e].rho*0 -- TODO: In future, replace 0 with acceleration field 
    end
    r_Wb[e].rhoE = dt/2.*r_Wb[e].rho*0 -- TODO: In future replace 0 with u.dot(a), vel dot acc
  end

  -- Then iterate over contributions in velocity space
  var U : double[3]
  for e in r_gridbar do
   
    U[0] = vxmesh[e.x].v
    U[1] = vymesh[e.y].v
    U[2] = vzmesh[e.z].v

    var e4 : int4d = {e.w, e.v, e.u, e.t}  

    for v = 0, effD do
      r_Wb[e4].rhov[v] = r_Wb[e4].rhov[v] + vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w*U[v]*r_gridbar[e].g 
    end
          
    r_Wb[e4].rhoE = r_Wb[e4].rhoE + vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w*r_gridbar[e].b


  end
  for e in r_Wb do

    if e.x == 0 then
      --c.printf("rho[%d] = %f\n", e.y, r_Wb[e].rho)
    end

    regentlib.assert(not [bool](isnan(r_Wb[e].rho)), "Step 2a\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[0])), "Step 2a\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[1])), "Step 2a\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[2])), "Step 2a\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhoE)), "Step 2a\n")
    
  end

end

-- Step 2b: compute original phi at interface using gbar, W at interface
-- Memory Recycling: phibar @ interface is used to store phi @ interface.
task Step2b(r_gridbar : region(ispace(int7d), grid), 
            r_Wb      : region(ispace(int4d), W),
            vxmesh    : region(ispace(int1d), vmesh),
            vymesh    : region(ispace(int1d), vmesh),
            vzmesh    : region(ispace(int1d), vmesh),
            dt : double, R : double, K : double, Cv : double, g : double,
            w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads writes(r_gridbar),
  reads(r_Wb, vxmesh, vymesh, vzmesh)
do
  for e in r_gridbar do
    var tg : double
    var tb : double
    var u : double
    var T : double
    var g_eq : double
    var b_eq : double
    var c2 : double
    var Xi : double[3]


    var e3 : int3d = {e.v, e.u, e.t}
    var e4 : int4d = {e.w, e.v, e.u, e.t}


    u = 0 
    for v = 0, effD do
      u = u + r_Wb[e4].rhov[v]/r_Wb[e4].rho*r_Wb[e4].rhov[v]/r_Wb[e4].rho
    end
    u = sqrt(u)
      
    T = Temperature(r_Wb[e4].rhoE/r_Wb[e4].rho, u, g, R)

    tg = visc(T, ur, Tr, w)/r_Wb[e4].rho/R/T
    tb = tg/Pr
    --c.printf("tg = %f, tb = %f\n", tg, tb)


    c2 = 0
    Xi[0] = vxmesh[e.x].v
    Xi[1] = vymesh[e.y].v
    Xi[2] = vzmesh[e.z].v

    for v = 0, effD do
      c2 =  c2 + (Xi[v] - r_Wb[e4].rhov[v]/r_Wb[e4].rho)*(Xi[v] - r_Wb[e4].rhov[v]/r_Wb[e4].rho)
    end

    g_eq = geq(c2, r_Wb[e4].rho, T, R, effD)
    b_eq = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0

    -- this is actually the original distribution function, recycling memory from gbar
    r_gridbar[e].g = 2*tg/(2*tg + dt/2.)*r_gridbar[e].g + dt/(4*tg + dt)*g_eq + dt*tg/(4*tg + dt)*0 -- TODO replace this last *0 with source term 
    r_gridbar[e].b = 2*tb/(2*tb + dt/2.)*r_gridbar[e].b + dt/(4*tb + dt)*b_eq + dt*tb/(4*tb + dt)*0 -- TODO replace this last *0 with source term

    if (isnan(r_gridbar[e].g) == 1 or isnan(r_gridbar[e].b) == 1) then

      c.printf("gbar = %.10f, bbar = %.10f, g_eq = %.10f, tg = %.10f, tb = %.10f\n", r_gridbar[e].g, r_gridbar[e].b, g_eq, tg, tb)
      
      regentlib.assert(not [bool](isnan(r_gridbar[e].g)), "Step 2b\n")
      regentlib.assert(not [bool](isnan(r_gridbar[e].b)), "Step 2b\n")
    
    end 
  end
end

-- Step 2c: Compute Microflux F at interface at half timestep using W/phi at interface.
task Step2c(r_gridbar : region(ispace(int7d), grid),
            r_F       : region(ispace(int6d), grid),
            r_mesh    : region(ispace(int3d), mesh),
            vxmesh    : region(ispace(int1d), vmesh),
            vymesh    : region(ispace(int1d), vmesh),
            vzmesh    : region(ispace(int1d), vmesh),
            BCs : int32[3], R : double, K : double, Cv : double, g : double,
            w : double, Tr : double, Pr : double, effD : int32, N : int32[3])
where
  reads(r_gridbar, vxmesh, vymesh, vzmesh, r_mesh),
  reads writes(r_F)
do    
  var sold : int3d = {-1, -1, -1}
  var A : double[3]
  var Xi : double[3]

  fill(r_F.g, 0)
  fill(r_F.b, 0)

  for e in r_F do
    var e3 : int3d = {e.w, e.v, e.u}
    if (not e3.x == sold.x and not e3.y == sold.y and not e3.z == sold.z) then
      A[0] = r_mesh[e3].dy*r_mesh[e3].dz
      A[1] = r_mesh[e3].dx*r_mesh[e3].dz
      A[2] = r_mesh[e3].dx*r_mesh[e3].dy
    end

    for Dim = 0, effD do
      var IL : int32 
      var IR : int32
      var JL : int32
      var JR : int32
      var KL : int32
      var KR : int32

      var e7 : int7d = {e.x, e.y, e.z, Dim, e.w, e.v, e.u}

      var i : int32 = e.w
      var j : int32 = e.v
      var k : int32 = e.u

      -- Boundary Conditions
      var right : double = 1.0 
      var left : double = 1.0

      -- Periodic Boundary Conditions
      if (Dim == 0 and BCs[0] == 0) then 
        IL = (i - 1 + N[0])%N[0] 
        IR = (i + 1)%N[0] 
        JL = j
        JR = j
        KL = k
        KR = k
      end 
      if (Dim == 1 and BCs[1] == 0) then
        IL = i
        IR = i
        JL = (j - 1 + N[1])%N[1]
        JR = (j + 1)%N[1] 
        KL = k
        KR = k
      end
      if (Dim == 2 and BCs[2] == 0) then
        IL = i
        IR = i
        JL = j
        JR = j
        KL = (k - 1 + N[2])%N[2] 
        KR = (k + 1)%N[2] 
      end


      -- Dirichlet Boundary Conditions
      if (Dim == 0 and BCs[0] == 1) then
        IL = i - 1 
        IR = i + 1 
        if IL <  0 then 
          IL = 0  
          left = 0
          right = 0
        end
        if IR == N[0] then 
          IR = N[0] - 1 
          left = 0
          right = 0
        end
        JL = j
        JR = j
        KL = k
        KR = k
      end 
      if (Dim == 1 and BCs[1] == 1) then
        IL = i
        IR = i
        JL = j - 1 
        JR = j + 1
        if JL <  0 then 
          JL = 0  
          left = 0
          right = 0
        end
        if JR == N[1] then 
          JR = N[1] - 1 
          left = 0
          right = 0
        end
        KL = k
        KR = k
      end
      if (Dim == 2 and BCs[2] == 1) then
        IL = i
        IR = i
        JL = j
        JR = j
        KL = k - 1 
        KR = k + 1 
        if KL <  0 then 
          KL = 0  
          left = 0
          right = 0
        end
        if KR == N[2] then 
          KR = N[2] - 1 
          left = 0
          right = 0
        end
      end
      
      
      -- Neumann Boundary Conditions
      if (Dim == 0 and BCs[0] == 2) then
        IL = i - 1 
        IR = i + 1 
        if IL <  0 then 
          IL = 0  
          left = 0
        end
        if IR == N[0] then 
          IR = N[0] - 1 
          right = 0
        end
        JL = j
        JR = j
        KL = k
        KR = k
      end 
      if (Dim == 1 and BCs[1] == 2) then
        IL = i
        IR = i
        JL = j - 1 
        JR = j + 1
        if JL <  0 then 
          JL = 0  
          left = 0
        end
        if JR == N[1] then 
          JR = N[1] - 1 
          right = 0
        end
        KL = k
        KR = k
      end
      if (Dim == 2 and BCs[2] == 2) then
        IL = i
        IR = i
        JL = j
        JR = j
        KL = k - 1 
        KR = k + 1 
        if KL <  0 then 
          KL = 0  
          left = 0
        end
        if KR == N[2] then 
          KR = N[2] - 1 
          right = 0
        end
      end
      
      -- Gather Left and Right Indices
      var eL : int6d = {e.x, e.y, e.z, IL, JL, KL}
      var eR : int6d = {e.x, e.y, e.z, IR, JR, KR}
      var eL3 : int3d = {IL, JL, KL}
      var eR3 : int3d = {IR, JR, KR}
      var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
      var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}

      Xi[0] = vxmesh[e.x].v
      Xi[1] = vymesh[e.y].v
      Xi[2] = vzmesh[e.z].v
  
      r_F[e].g = r_F[e].g + Xi[Dim]*A[Dim]*(right*r_gridbar[e7].g - left*r_gridbar[eL7].g)
      r_F[e].b = r_F[e].b + Xi[Dim]*A[Dim]*(right*r_gridbar[e7].b - left*r_gridbar[eL7].b)


      regentlib.assert(not [bool](isnan(r_F[e].g)), "Step 2c\n")
      regentlib.assert(not [bool](isnan(r_F[e].b)), "Step 2c\n")

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
  var V : double 
  var Xi : double[3]
  var uo : double
  var To : double
  var tgo : double
  var tbo : double
  var u : double
  var T : double
  var tg : double
  var tb : double
  var c2 : double
  var g_eq : double
  var b_eq : double
  var g_eqo : double
  var b_eqo : double

  for e in r_grid do
    var e3 : int3d = {e.w, e.v, e.u}
    
    V = r_mesh[e3].dx*r_mesh[e3].dy*r_mesh[e3].dz

  

  Xi[0] = vxmesh[e.x].v
  Xi[1] = vymesh[e.y].v
  Xi[2] = vzmesh[e.z].v

  -- Compute old flow velocity
  uo = 0 -- old flow velocity
  for v = 0, effD do
    uo += r_W[e3].rhov[v]/r_W[e3].rho*r_W[e3].rhov[v]/r_W[e3].rho
  end
  uo = sqrt(uo)
 
  regentlib.assert(bool(uo>=0), "uo")
  -- Compute old temperature
  To = Temperature(r_W[e3].rhoE/r_W[e3].rho, uo, g, R)
  regentlib.assert(bool(To>=0), "To")
  
  -- Compute old taus
  tgo = visc(To, ur, Tr, w)/r_W[e3].rho/R/To
  tbo = tgo/Pr 

  -- Compute old eq's
  c2 = 0
  for v = 0, effD do
    c2 += (Xi[v]-r_W[e3].rhov[v]/r_W[e3].rho)*(Xi[v]-r_W[e3].rhov[v]/r_W[e3].rho)
  end
  g_eqo = geq(c2, r_W[e3].rho, To, R, effD)
  b_eqo = g_eqo*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3-effD+K)*R*To)/2


  -- Step 4: Update W at cell center
  r_W[e3].rho = r_W[e3].rho - (dt/V*r_F[e].g + dt*0)*vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w -- TODO replace 0 with source term.
  for v = 0, effD do
      r_W[e3].rhov[v] = r_W[e3].rhov[v] - dt/V*r_F[e].g*Xi[v]*vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w
  end
  
  r_W[e3].rhoE = r_W[e3].rhoE - dt/V*r_F[e].b*vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w

  -- Step 5: Update Phi at cell center, need new eq and tau
  -- Compute flow velocity u
  u = 0 
  for v = 0, effD do 
    u += r_W[e3].rhov[v]/r_W[e3].rho*r_W[e3].rhov[v]/r_W[e3].rho
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

  -- Compute new eq's
  c2 = 0  -- reset c2 from before
  for v = 0, effD do
    c2 += (Xi[v]-r_W[e3].rhov[v]/r_W[e3].rho)*(Xi[v]-r_W[e3].rhov[v]/r_W[e3].rho)
  end
  g_eq = geq(c2, r_W[e3].rho, T, R, effD)
  b_eq = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3-effD+K)*R*T)/2


  -- Update phi
  var i : int32 = e.w
  var j : int32 = e.v
  var k : int32 = e.u

  if ((BCs[0] == 1 and i == 0) or (BCs[0] == 1 and i == N[0] - 1) or
      (BCs[1] == 1 and j == 0) or (BCs[1] == 1 and j == N[1] - 1) or
      (BCs[2] == 1 and k == 0) or (BCs[2] == 1 and k == N[2] - 1)) then
    r_grid[e].g = r_grid[e].g
    r_grid[e].b = r_grid[e].b
  else
    r_grid[e].g = (r_grid[e].g + dt/2.0*(g_eq/tg + (g_eqo-r_grid[e].g)/tgo) - dt/V*r_F[e].g + dt*0)/(1+dt/2.0/tg) -- TODO replace 0 with source term
    r_grid[e].b = (r_grid[e].b + dt/2.0*(b_eq/tb + (b_eqo-r_grid[e].b)/tbo) - dt/V*r_F[e].b + dt*0)/(1+dt/2.0/tb) -- TODO replace 0 with source term
  end

  if isnan(r_grid[e].g) == 1 then
    c.printf("Step4and5: g_eq = %f, tg = %f, tgo = %f, r_F[e].g = %f\n", g_eq, tg, tgo, r_F[e].g)
  end 
  if isnan(r_grid[e].b) == 1 then
    c.printf("Step4and5: b_eq = %f, tb = %f, tbo = %f, r_F[e].b = %f\n", b_eq, tb, tbo, r_F[e].b)
  end 
  regentlib.assert(not [bool](isnan(r_grid[e].g)), "Step4and5\n")
  regentlib.assert(not [bool](isnan(r_grid[e].b)), "Step4and5\n")

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

  for e in r_grid do 
    var s : int3d = {e.w, e.v, e.u}

    u = 0
    for dim = 0, effD do
      u += r_W[s].rhov[dim]/r_W[s].rho*r_W[s].rhov[dim]/r_W[s].rho
    end
    u = sqrt(u)
      
    T = Temperature(r_W[s].rhoE/r_W[s].rho, u, g, R)

    var c2 : double = 0
    var Xi : double[3]  
    Xi[0] = vxmesh[e.x].v
    Xi[1] = vymesh[e.y].v
    Xi[2] = vzmesh[e.z].v

    for dim = 0, effD do
      c2 += (Xi[dim] - r_W[s].rhov[dim]/r_W[s].rho)*(Xi[dim] - r_W[s].rhov[dim]/r_W[s].rho) --TODO could be bugged
    end

    r_grid[e].g = geq(c2, r_W[s].rho, T, R, effD)
    r_grid[e].b = r_grid[e].g*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0
  
    if s.x == 0 and s.y == 0 and s.z == 0 then
      rhotest += r_grid[e].g*vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w
      Etest += r_grid[e].b*vxmesh[e.x].w*vymesh[e.y].w*vzmesh[e.z].w
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
      size_x, size_y =  i, parallelism / i, 1
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
    f6.w, f6.v, f6.u = f3.x, f3.y, f3.x    
  elseif effD == 2 then
    var f3 = factorize2d(parallelism)
    f6.w, f6.v, f6.u = f3.x, f3.y, f3.x    
  elseif effD == 3 then
    var f3 = factorize3d(parallelism)
    f6.w, f6.v, f6.u = f3.x, f3.y, f3.x    
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
  var r_grid      = region(ispace(int6d, {NV[0], NV[1], NV[2], N[0], N[1], N[2]}), grid)
  var r_gridbarp  = region(ispace(int6d, {NV[0], NV[1], NV[2], N[0], N[1], N[2]}), grid)
  var r_gridbarpb = region(ispace(int7d, {NV[0], NV[1], NV[2], effD, N[0], N[1], N[2]}), grid)
  var r_gridbar   = region(ispace(int7d, {NV[0], NV[1], NV[2], effD, N[0], N[1], N[2]}), grid)
  var r_sig       = region(ispace(int7d, {NV[0], NV[1], NV[2], effD, N[0], N[1], N[2]}), grid)
  var r_sig2      = region(ispace(int8d, {NV[0], NV[1], NV[2], effD, effD, N[0], N[1], N[2]}), grid)
 
  -- Create regions for mesh and conserved variables (cell center and interface)
  var r_mesh = region(ispace(int3d, {N[0], N[1], N[2]}), mesh)
  var r_W    = region(ispace(int3d, {N[0], N[1], N[2]}), W)
  var r_Wb   = region(ispace(int4d, {effD, N[0], N[1], N[2]}), W)
 
  -- Create regions for velocity space and initialize
  var vxmesh = region(ispace(int1d, NV[0]), vmesh) 
  var vymesh = region(ispace(int1d, NV[1]), vmesh) 
  var vzmesh = region(ispace(int1d, NV[2]), vmesh) 
  NewtonCotes(vxmesh, vymesh, vzmesh, NV, Vmin, Vmax)

  -- Create regions for source terms and flux
  var r_S = region(ispace(int6d, {NV[0], NV[1], NV[2], N[0], N[1], N[2]}), grid)
  var r_F = region(ispace(int6d, {NV[0], NV[1], NV[2], N[0], N[1], N[2]}), grid)


  -- Create partitions for regions
  var f6 : int6d = factorize(config.cpus, effD)
  var f3 : int3d = {f6.w, f6.v, f6.u}
  var f4 : int4d = {1, f6.w, f6.v, f6.u}
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
  var p_gridbar = partition(equal, r_gridbar, p7)
  var p_sig = partition(equal, r_sig, p7)
  var p_sig2 = partition(equal, r_sig2, p8)
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
    var il : int32 = bounds.lo.v
    var jl : int32 = bounds.lo.u
    var kl : int32 = bounds.lo.t
    var ir : int32 = bounds.hi.v
    var jr : int32 = bounds.hi.u
    var kr : int32 = bounds.hi.t

    var col3 : int3d = {col7.v, col7.u, col7.t}
    var col6 : int6d = {col7.x, col7.y, col7.z, col7.v, col7.u, col7.t}
    var col8 : int8d = {col7.x, col7.y, col7.z, 1, col7.w, col7.v, col7.u, col7.t}

    -- for reference -- terra BC(i : int32, j : int32, k : int32, Dim : int32, BCs : int32[3], N : int32[3])

    var rleftx3 : rect3d = { {BC(il, jl, kl, 0, BCs, N)[0], bounds.lo.u, bounds.lo.t}, 
                          {BC(il, jl, kl, 0, BCs, N)[0] + 1, bounds.hi.u, bounds.hi.t}}
    var rlefty3 : rect3d = { {bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1], bounds.lo.t}, 
                          {bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1] + 1, bounds.hi.t}}
    var rleftz3 : rect3d = { {bounds.lo.v, bounds.lo.u, BC(il, jl, kl, 2, BCs, N)[2]}, 
                          {bounds.lo.v, bounds.hi.u, BC(il, jl, kl, 2, BCs, N)[2] + 1}}

    var rrightx3 : rect3d = { {BC(ir, jr, kr, 0, BCs, N)[3], bounds.lo.u, bounds.lo.t}, 
                          {BC(ir, jr, kr, 0, BCs, N)[3] + 1, bounds.hi.u, bounds.hi.t}}
    var rrighty3 : rect3d = { {bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4], bounds.lo.t}, 
                          {bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4] + 1, bounds.hi.t}}
    var rrightz3 : rect3d = { {bounds.lo.v, bounds.lo.u, BC(ir, jr, kr, 2, BCs, N)[5]}, 
                          {bounds.lo.v, bounds.hi.u, BC(ir, jr, kr, 2, BCs, N)[5] + 1}}
    __fence(__execution, __block)
    c.printf("rect3d Done\n")

    var rleftx6 : rect6d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, BC(il, jl, kl, 0, BCs, N)[0], bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, BC(il, jl, kl, 0, BCs, N)[0] + 1, bounds.hi.u, bounds.hi.t}}
    var rlefty6 : rect6d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1], bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1] + 1, bounds.hi.t}}
    var rleftz6 : rect6d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.v, bounds.lo.u, BC(il, jl, kl, 2, BCs, N)[2]}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.lo.v, bounds.hi.u, BC(il, jl, kl, 2, BCs, N)[2] + 1}}

    var rrightx6 : rect6d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, BC(ir, jr, kr, 0, BCs, N)[3], bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, BC(ir, jr, kr, 0, BCs, N)[3] + 1, bounds.hi.u, bounds.hi.t}}
    var rrighty6 : rect6d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4], bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4] + 1, bounds.hi.t}}
    var rrightz6 : rect6d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.v, bounds.lo.u, BC(ir, jr, kr, 2, BCs, N)[5]}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.lo.v, bounds.hi.u, BC(ir, jr, kr, 2, BCs, N)[5] + 1}}
    __fence(__execution, __block)
    c.printf("rect6d Done\n")

    var rleftx7 : rect7d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, BC(il, jl, kl, 0, BCs, N)[0], bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, BC(il, jl, kl, 0, BCs, N)[0] + 1, bounds.hi.u, bounds.hi.t}}
    var rlefty7 : rect7d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1], bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1] + 1, bounds.hi.t}}
    var rleftz7 : rect7d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, BC(il, jl, kl, 2, BCs, N)[2]}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.lo.v, bounds.hi.u, BC(il, jl, kl, 2, BCs, N)[2] + 1}}

    var rrightx7 : rect7d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, BC(ir, jr, kr, 0, BCs, N)[3], bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, BC(ir, jr, kr, 0, BCs, N)[3] + 1, bounds.hi.u, bounds.hi.t}}
    var rrighty7 : rect7d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4], bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4] + 1, bounds.hi.t}}
    var rrightz7 : rect7d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, BC(ir, jr, kr, 2, BCs, N)[5]}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.lo.v, bounds.hi.u, BC(ir, jr, kr, 2, BCs, N)[5] + 1}}
    __fence(__execution, __block)
    c.printf("rect7d Done\n")

    var rleftx8 : rect8d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, BC(il, jl, kl, 0, BCs, N)[0], bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, 1, BC(il, jl, kl, 0, BCs, N)[0] + 1, bounds.hi.u, bounds.hi.t}}
    var rlefty8 : rect8d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1], bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, 1, bounds.lo.v, BC(il, jl, kl, 1, BCs, N)[1] + 1, bounds.hi.t}}
    var rleftz8 : rect8d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, bounds.lo.v, bounds.lo.u, BC(il, jl, kl, 2, BCs, N)[2]}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, 1, bounds.lo.v, bounds.hi.u, BC(il, jl, kl, 2, BCs, N)[2] + 1}}

    var rrightx8 : rect8d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, BC(ir, jr, kr, 0, BCs, N)[3], bounds.lo.u, bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, 1, BC(ir, jr, kr, 0, BCs, N)[3] + 1, bounds.hi.u, bounds.hi.t}}
    var rrighty8 : rect8d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4], bounds.lo.t}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, 1, bounds.lo.v, BC(ir, jr, kr, 1, BCs, N)[4] + 1, bounds.hi.t}}
    var rrightz8 : rect8d = { {bounds.lo.x, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, bounds.lo.v, bounds.lo.u, BC(ir, jr, kr, 2, BCs, N)[5]}, 
                          {bounds.hi.x, bounds.hi.y, bounds.hi.z, bounds.hi.w, 1, bounds.lo.v, bounds.hi.u, BC(ir, jr, kr, 2, BCs, N)[5] + 1}}
    __fence(__execution, __block)
    c.printf("rect8d Done\n")

    -- Color in left strips
    coloring.color_domain(c3Lx, col3, rleftx3)
    coloring.color_domain(c3Ly, col3, rlefty3)
    coloring.color_domain(c3Lz, col3, rleftz3)

    coloring.color_domain(c6Lx, col6, rleftx6)
    coloring.color_domain(c6Ly, col6, rlefty6)
    coloring.color_domain(c6Lz, col6, rleftz6)

    coloring.color_domain(c7Lx, col7, rleftx7)
    coloring.color_domain(c7Ly, col7, rlefty7)
    coloring.color_domain(c7Lz, col7, rleftz7)

    coloring.color_domain(c8Lx, col8, rleftx8)
    coloring.color_domain(c8Ly, col8, rlefty8)
    coloring.color_domain(c8Lz, col8, rleftz8)

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

  var plx_sig = partition(disjoint, r_sig, c7Lx, p7)
  var ply_sig = partition(disjoint, r_sig, c7Ly, p7)
  var plz_sig = partition(disjoint, r_sig, c7Lz, p7)
  var prx_sig = partition(disjoint, r_sig, c7Rx, p7)
  var pry_sig = partition(disjoint, r_sig, c7Ry, p7)
  var prz_sig = partition(disjoint, r_sig, c7Rz, p7)
  __fence(__execution, __block)
  c.printf("sig Strips Done\n")

  var plx_gridbar = partition(disjoint, r_gridbar, c7Lx, p7)
  var ply_gridbar = partition(disjoint, r_gridbar, c7Ly, p7)
  var plz_gridbar = partition(disjoint, r_gridbar, c7Lz, p7)
  var prx_gridbar = partition(disjoint, r_gridbar, c7Rx, p7)
  var pry_gridbar = partition(disjoint, r_gridbar, c7Ry, p7)
  var prz_gridbar = partition(disjoint, r_gridbar, c7Rz, p7)
  __fence(__execution, __block)
  c.printf("gridbar Strips Done\n")

  var plx_sig2 = partition(disjoint, r_sig2, c8Lx, p8)
  var ply_sig2 = partition(disjoint, r_sig2, c8Ly, p8)
  var plz_sig2 = partition(disjoint, r_sig2, c8Lz, p8)
  var prx_sig2 = partition(disjoint, r_sig2, c8Rx, p8)
  var pry_sig2 = partition(disjoint, r_sig2, c8Ry, p8)
  var prz_sig2 = partition(disjoint, r_sig2, c8Rz, p8)
  __fence(__execution, __block)
  c.printf("sig2 Strips Done\n")
  c.printf("Left/Right Strip Partitioning Done\n")

  --Initialize r_mesh
  var MeshType : int32 = 1
  InitializeMesh(r_mesh, N, MeshType) --TODO Needs more input for nested, user-def etc.

  --Initialize r_W
  InitializeW(r_W, r_mesh, N, NV, testProblem, R, Cv)
  
  --Initialize r_grid
  InitializeGrid(r_grid, r_mesh, r_W, vxmesh, vymesh, vzmesh, testProblem, R, K, Cv, g, w, ur, Tr, Pr, N, NV, effD)

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
    c.printf("Dump %d\n", dumpiter)
  end
  while Tsim < Tf do -- and iter < 1 do
    iter += 1

    var dt = TimeStep(calcdt, dtdump-Tdump, Tf-Tsim)

    for col6 in p_grid.colors do 
      var col3 : int3d = {col6.w, col6.v, col6.u}
      Step1a(p_grid[col6], p_gridbarp[col6], p_S[col6], p_W[col3], vxmesh, vymesh, vzmesh, dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    end

    for col6 in p_grid.colors do 
      var col3 : int3d = {col6.w, col6.v, col6.u}
      var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
      var col8 : int8d = {col6.x, col6.y, col6.z, 0, 0, col6.w, col6.v, col6.u}
      Step1b_a(p_gridbarp[col6], p_gridbarpb[col7], p_sig[col7], p_sig2[col8], p_mesh[col3], plx_mesh[col3], ply_mesh[col3], plz_mesh[col3], prx_mesh[col3], pry_mesh[col3], prz_mesh[col3], plx_gridbarp[col6], ply_gridbarp[col6], plz_gridbarp[col6], prx_gridbarp[col6], pry_gridbarp[col6], prz_gridbarp[col6], plx_sig[col7], ply_sig[col7], plz_sig[col7], prx_sig[col7], pry_sig[col7], prz_sig[col7], vxmesh, vymesh, vzmesh, BCs, R, K, Cv, N, g, w, ur, Tr, Pr, effD)
    end

    for col6 in p_grid.colors do 
      var col3 : int3d = {col6.w, col6.v, col6.u}
      var col7 : int7d = {col6.x, col6.y, col6.z, 0, col6.w, col6.v, col6.u}
      
      Step1b_b(p_sig[col7], p_mesh[col3], p_gridbarp[col6], p_gridbarpb[col7], prx_gridbarp[col6], pry_gridbarp[col6], prz_gridbarp[col6], prx_sig[col7], pry_sig[col7], prz_sig[col7], prx_mesh[col3], pry_mesh[col3], prz_mesh[col3], vxmesh, vymesh, vzmesh, BCs, N)
    end

    Step1c(r_gridbar, r_gridbarpb, vxmesh, vymesh, vzmesh, r_sig2, dt, BCs, R, K, Cv, N, g, w, ur, Tr, Pr, effD)
    Step2a(r_gridbar, vxmesh, vymesh, vzmesh, r_Wb, dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    Step2b(r_gridbar, r_Wb, vxmesh, vymesh, vzmesh, dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    Step2c(r_gridbar, r_F, r_mesh, vxmesh, vymesh, vzmesh, BCs, R, K, Cv, g, w, Tr, Pr, effD, N)
    Step3() -- TODO
    Step4and5(r_grid, r_W, r_mesh, r_F, vxmesh, vymesh, vzmesh, dt, BCs, R, K, Cv, N, g, w, ur, Tr, Pr, effD)
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
