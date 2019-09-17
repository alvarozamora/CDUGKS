import "regent"

-- Helper modules to handle PNG files and command line arguments
local Config = require("config")
--local Dump = require("dump") -- TODO

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
task Step1b(r_gridbarp : region(ispace(int6d), grid),
            r_gridbarpb : region(ispace(int7d), grid),
            r_sig : region(ispace(int7d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int3d), mesh),
            vxmesh : region(ispace(int1d), vmesh),            
            vymesh : region(ispace(int1d), vmesh),            
            vzmesh : region(ispace(int1d), vmesh),            
            BCs : int32[3], R : double, K : double, Cv : double, N : int32[3], 
            g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where 
  reads(r_gridbarp, r_mesh, vxmesh, vymesh, vzmesh),
  reads writes(r_sig),
  writes(r_sig2),
  reads writes(r_gridbarpb)
do
  var xC : double[3]
  var sC : double[3]
  var xL : double[3]
  var sL : double[3]
  var xR : double[3]
  var sR : double[3]
  var xL2 : double[3]
  var sL2 : double[3]
  var xR2 : double[3]
  var sR2 : double[3]

  var olds : int3d = {-1, -1, -1}  
  for e in r_gridbarp do
    var sidx : int3d = {e.w, e.v, e.u}
    if (not sidx.x == olds.x and not sidx.y == olds.y and not sidx.z == olds.z) then
      xC[0] = r_mesh[sidx].x
      xC[1] = r_mesh[sidx].y
      xC[2] = r_mesh[sidx].z
  
      sC[0] = r_mesh[sidx].dx
      sC[1] = r_mesh[sidx].dy
      sC[2] = r_mesh[sidx].dz
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
      end
      
      -- Gather Left and Right Indices
      var eL : int6d = {e.x, e.y, e.z, IL, JL, KL}
      var eR : int6d = {e.x, e.y, e.z, IR, JR, KR}
      var eL3 : int3d = {IL, JL, KL}
      var eR3 : int3d = {IR, JR, KR}
      var eL7 : int7d = {e.x, e.y, e.z, Dim, IL, JL, KL}
      var eR7 : int7d = {e.x, e.y, e.z, Dim, IR, JR, KR}

      

      -- Gather position of left/right cell centers
      xL[0] = r_mesh[eL3].x
      xL[1] = r_mesh[eL3].y
      xL[2] = r_mesh[eL3].z
      xR[0] = r_mesh[eR3].x
      xR[1] = r_mesh[eR3].y
      xR[2] = r_mesh[eR3].z
      
      -- Gather cell size of left/right cells
      sL[0] = r_mesh[eL3].dx
      sL[1] = r_mesh[eL3].dy
      sL[2] = r_mesh[eL3].dz
      sR[0] = r_mesh[eR3].dx
      sR[1] = r_mesh[eR3].dy
      sR[2] = r_mesh[eR3].dz
      

      -- Computing phisigma, at cell 
      r_sig[e7].g = VanLeer(r_gridbarp[eL].g, r_gridbarp[e].g, r_gridbarp[eR].g, xL[Dim], xC[Dim], xR[Dim])
      r_sig[e7].b = VanLeer(r_gridbarp[eR].g, r_gridbarp[e].b, r_gridbarp[eR].b, xL[Dim], xC[Dim], xR[Dim])

      --printf("Checking bsigma input[%d][%d]: bbarp[idxL] = %f, bbarp[idx] = %f, bbarp[idxR] = %f, xL[Dim] = %f, xC[Dim] = %f, xR[Dim] = %f\n", sidx, vx, bbarp[idxL], bbarp[idx], bbarp[idxR], xL[Dim], xC[Dim], xR[Dim])

      

      -- Computing phisigma, at interface
      for Dim2 = 0, effD do
        var IL2 : int32 
        var IR2 : int32
        var JL2 : int32
        var JR2 : int32
        var KL2 : int32
        var KR2 : int32

        var e8 : int8d = {e.x, e.y, e.z, Dim, Dim2, e.w, e.v, e.u}

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


        -- Gather position of left/right cell centers -- TODO Maybe Optimize by only updating Dim2
        xL2[0] = r_mesh[eL2_3].x
        xL2[1] = r_mesh[eL2_3].y
        xL2[2] = r_mesh[eL2_3].z
        xR2[0] = r_mesh[eR2_3].x
        xR2[1] = r_mesh[eR2_3].y
        xR2[2] = r_mesh[eR2_3].z
      
        -- Gather cell size of left/right cells
        sL2[0] = r_mesh[eL2_3].dx
        sL2[1] = r_mesh[eL2_3].dy
        sL2[2] = r_mesh[eL2_3].dz
        sR2[0] = r_mesh[eR2_3].dx
        sR2[1] = r_mesh[eR2_3].dy
        sR2[2] = r_mesh[eR2_3].dz


        -- Dim  is vector component that is being interpolated.
        -- Dim2 is direction of interpolation.

        r_sig2[e8].g = r_sig[e7].g + (sC[Dim2]/2.0)*VanLeer(r_sig[eL2_7].g, r_sig[e7].g, r_sig[eR2_7].g, xL2[Dim2], xC[Dim2], xR2[Dim2])
        r_sig2[e8].b = r_sig[e7].b + (sC[Dim2]/2.0)*VanLeer(r_sig[eL2_7].b, r_sig[e7].b, r_sig[eR2_7].b, xL2[Dim2], xC[Dim2], xR2[Dim2])
      end
  


      -- Dot Product is just a single product when using rectangular mesh
      var interpidx : int7d = e7
      var swap : double = 1.0

      if (vxmesh[e.x].v < 0 and Dim == 0) then
        interpidx = eR7
        swap = -1
      elseif (vymesh[e.y].v < 0 and Dim == 1) then
        interpidx = eR7 
        swap = -1
      elseif (vzmesh[e.z].v < 0 and Dim == 2) then
        interpidx = eR7 
        swap = -1
      end

      var interpid  : int6d = {interpidx.x, interpidx.y, interpidx.z, interpidx.v, interpidx.u, interpidx.t}
      
      -- TODO need to change sC to sR when swap, doesnt currently matter for sod/KHI/RTI bc dx_i = dx_0
      r_gridbarpb[e7].g = r_gridbarp[interpid].g + swap*sC[Dim]/2.0*r_sig[interpidx].g
      r_gridbarpb[e7].b = r_gridbarp[interpid].b + swap*sC[Dim]/2.0*r_sig[interpidx].b


    if (isnan(r_gridbarpb[e7].g) == 1 or isnan(r_gridbarpb[e7].b) == 1) then

      c.printf("Step 1b: r_gridbarp.g = %f, r_gridbarp.b = %f, r_sig.g = %f, r_sig.b = %f\n", r_gridbarp[interpid].g, r_gridbarp[interpid].b, r_sig[interpidx].g, r_sig[interpidx].b)

      regentlib.assert(not [bool](isnan(r_gridbarpb[e7].g)), "Step 1b\n")
      regentlib.assert(not [bool](isnan(r_gridbarpb[e7].b)), "Step 1b\n")
    
    end


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
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y size_z =  i, 1, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int2d { size_x, size_y, 1 }
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

task factorize(parallelism: int, effD : int32) : wild
  if effD == 2 then
    

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
 
  -- Create partitions for regions
  -- Create an equal partition of the grid
  var p_grid_colors = ispace(int2d, factorize2d(config.parallelism))
  var p_grid = partition(equal, r_grid, p_grid_colors)
  
  -- Create partitions for left/right ghost regions
  var plx = region(ispace(int7d, {NV[0], NV[1], NV[2], 1, 1, N[1], N[2]}), ghost)
  var prx = region(ispace(int7d, {NV[0], NV[1], NV[2], 1, 1, N[1], N[2]}), ghost)
  var ply = region(ispace(int7d, {NV[0], NV[1], NV[2], 1, N[0], 1, N[2]}), ghost)
  var pry = region(ispace(int7d, {NV[0], NV[1], NV[2], 1, N[0], 1, N[2]}), ghost)
  var plz = region(ispace(int7d, {NV[0], NV[1], NV[2], 1, N[0], N[1], 1}), ghost)
  var prz = region(ispace(int7d, {NV[0], NV[1], NV[2], 1, N[0], N[1], 1}), ghost)

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

    Step1a(r_grid, r_gridbarp, r_S, r_W, vxmesh, vymesh, vzmesh, dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    Step1b(r_gridbarp, r_gridbarpb, r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, BCs, R, K, Cv, N, g, w, ur, Tr, Pr, effD)
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
