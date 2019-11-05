import "regent"

-- Helper modules to handle PNG files and command line arguments
local EdgeConfig = require("edge_config")
local coloring   = require("coloring_util")

-- Some C APIs
local c     = regentlib.c
local sqrt  = regentlib.sqrt(double)
local cmath = terralib.includec("math.h")
local PI = cmath.M_PI

-- Field space for Grid
fspace Grid2D
{
   s : int8,
   bx: int8,
   by: int8,
   l : int32,
   f : bool
}

terra uniform(data : &c.drand48_data)
  var flip : double[1]
  c.srand48_r(c.legion_get_current_time_in_nanos(), data)
  c.drand48_r(data, [&double](flip))
  return flip[0]
end


task initialize(r_grid : region(ispace(int2d), Grid2D),
                r_rng : region(ispace(int2d), c.drand48_data[1]),
                p : double)
where
  reads writes(r_grid),
  reads writes(r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  for e in r_grid do 
    var flip : double = uniform(data)
    if flip > p/2 then
      r_grid[e].s = 1
    else
      r_grid[e].s = -1
    end

    r_grid[e].bx = 0
    r_grid[e].by = 0
  end
  return 1
end


task Bondsx(r_grid : region(ispace(int2d), Grid2D),
           r_rng : region(ispace(int2d), c.drand48_data[1]),
           p : double)
where
  reads(r_grid.s), writes(r_grid.bx),
  reads writes(r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  var bounds = r_grid.bounds
  for e in r_grid do
    var flip : double = uniform(data)
    if flip > p and e.x <  bounds.hi.x then
      if r_grid[e].s==r_grid[e+{1,0}].s then
        r_grid[e].bx = 1
      else
        r_grid[e].bx = 0
      end
    end
  end
end

task Bondsy(r_grid : region(ispace(int2d), Grid2D),
           r_rng : region(ispace(int2d), c.drand48_data[1]),
           p : double)
where
  reads(r_grid.s), writes(r_grid.by),
  reads writes(r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  var bounds = r_grid.bounds
  for e in r_grid do
    var flip2 : double = uniform(data)
    if flip2 > p and e.y < bounds.hi.y then
      if r_grid[e].s==r_grid[e+{0,1}].s then
        r_grid[e].by = 1
      else
        r_grid[e].by = 0
      end
    end
  end
end

task Boundsx(r_bound : region(ispace(int2d), Grid2D),
            r_rng : region(ispace(int2d), c.drand48_data[1]),
            p : double, L : int32)
where
 reads(r_bound.s), writes(r_bound.bx),
 reads writes(r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  var bounds = r_bound.bounds
  for e in r_bound do
    if e.x == bounds.lo.x then
      if r_bound[e].s == r_bound[e+{1,0}].s then
        var flip : double = uniform(data)
        if flip > p then
          r_bound[e].bx = 1
          --c.printf("e = {%d, %d}\n", e.x, e.y)
        end
      end
    end
  end
end

task Boundsy(r_bound : region(ispace(int2d), Grid2D),
            r_rng : region(ispace(int2d), c.drand48_data[1]),
            p : double, L : int32)
where
 reads(r_bound.s), writes(r_bound.by),
 reads writes(r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  var bounds = r_bound.bounds
  for e in r_bound do
    if e.y == bounds.lo.y then
      if r_bound[e].s == r_bound[e+{0,1}].s then
        var flip : double = uniform(data)
        if flip > p then
          r_bound[e].by = 1
        end
      end
    end
  end
end


task Globalx(r_bound : region(ispace(int2d), Grid2D))
where
 reads(r_bound.bx), reads writes(r_bound.{l,f})
do
  var check : int32 = 0
  var bounds = r_bound.bounds
  for e in r_bound do
    if e.x == bounds.lo.x then
      if r_bound[e].bx == 1 then
        if r_bound[e].l < r_bound[e+{1,0}].l then
          r_bound[e+{1,0}].l = r_bound[e].l
          r_bound[e+{1,0}].f = r_bound[e].f
          check += 1     
        elseif r_bound[e].l > r_bound[e+{1,0}].l then  
          r_bound[e].l = r_bound[e+{1,0}].l
          r_bound[e].f = r_bound[e+{1,0}].f
          check += 1     
        end
      end
    end
  end
  return check
end

task Globaly(r_bound : region(ispace(int2d), Grid2D))
where
 reads(r_bound.by), reads writes(r_bound.{l,f})
do
  var check : int32 = 0
  var bounds = r_bound.bounds
  for e in r_bound do
    if e.y == bounds.lo.y then
      if r_bound[e].by == 1 then
        if r_bound[e].l < r_bound[e+{0,1}].l then
          r_bound[e+{0,1}].l = r_bound[e].l
          r_bound[e+{0,1}].f = r_bound[e].f
          check += 1     
        elseif r_bound[e].l > r_bound[e+{0,1}].l then  
          r_bound[e].l = r_bound[e+{0,1}].l
          r_bound[e].f = r_bound[e+{0,1}].f
          check += 1     
        end
      end
    end
  end
  return check
end


task SubInitialize(r_grid : region(ispace(int2d), Grid2D),
                   r_rng : region(ispace(int2d), c.drand48_data[1]),		   
                   L : int32, p : double)
where
  writes(r_grid.{l,f}),
  reads writes (r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  for e in r_grid do
    r_grid[e].l = e.x + e.y*L

    var flip : double = uniform(data)
    if flip > p then
      r_grid[e].f = true
    else
      r_grid[e].f = false
    end
  end
end


task Local(r_grid : region(ispace(int2d), Grid2D))
where
  reads writes (r_grid.{l,f}),
  reads (r_grid.{bx,by})
do
  var check : int32 = 0
  var bounds = r_grid.bounds
  for e in r_grid do
    if e.x < bounds.hi.x and r_grid[e].bx==1 then 
      if r_grid[e].l > r_grid[e+{1,0}].l then
        check += 1
        r_grid[e].l = r_grid[e+{1,0}].l
        r_grid[e].f = r_grid[e+{1,0}].f
      elseif r_grid[e].l < r_grid[e+{1,0}].l then
        check += 1
        r_grid[e+{1,0}].l = r_grid[e].l
        r_grid[e+{1,0}].f = r_grid[e].f
      end
    end    

    if e.y < bounds.hi.y and r_grid[e].by==1 then 
      if r_grid[e].l > r_grid[e+{0,1}].l then
        check += 1
        r_grid[e].l = r_grid[e+{0,1}].l
        r_grid[e].f = r_grid[e+{0,1}].f
      elseif r_grid[e].l < r_grid[e+{0,1}].l then
        check += 1
        r_grid[e+{0,1}].l = r_grid[e].l
        r_grid[e+{0,1}].f = r_grid[e].f
      end
    end    
  end
  return check
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

terra wait_for(x : int) return 1 end
task block_task(r_image : region(ispace(int2d), Grid2D))
where
  reads writes(r_image)
do
  return 1
end


task Dump(r_grid : region(ispace(int2d), Grid2D))
where 
  reads (r_grid)
do
  var f = c.fopen('Grid','wb')
    
  for e in r_grid do
    dumpint32(f,e.s)
  end
  
  for e in r_grid do
    dumpint32(f,e.l)
  end

  for e in r_grid do
    dumpint32(f,e.bx)
  end
    
  for e in r_grid do
    dumpint32(f,e.by)
  end
end

task Flip(r_grid : region(ispace(int2d), Grid2D))
where 
  reads(r_grid.f),
  reads writes (r_grid.s)
do
  for e in r_grid do
    if r_grid[e].f == true then
      r_grid[e].s = -r_grid[e].s
    end
  end
end


task Edgex(r_left : region(ispace(int2d), Grid2D),
           r_right : region(ispace(int2d), Grid2D),
           L : int32)
where
  reads(r_left.bx), reads writes(r_left.{l,f}),
  reads writes(r_right.{l,f})
do
  var check : int32 = 0
  for e in r_left do
    if r_left[e].bx == 1 then
      if r_left[e].l < r_right[e-{L-1,0}].l then
        r_right[e-{L-1,0}].l = r_left[e].l
        r_right[e-{L-1,0}].f = r_left[e].f
        check += 1
      elseif r_left[e].l > r_right[e-{L-1,0}].l then
        r_left[e].l = r_right[e-{L-1,0}].l
        r_left[e].f = r_right[e-{L-1,0}].f
        check += 1
      end
    end
  end
  return check
end

task Edgey(r_left : region(ispace(int2d), Grid2D),
           r_right : region(ispace(int2d), Grid2D),
           L : int32)
where
  reads(r_left.by), reads writes(r_left.{l,f}),
  reads writes(r_right.{l,f})
do
  var check : int32 = 0
  for e in r_left do
    if r_left[e].by == 1 then
      if r_left[e].l < r_right[e-{0,L-1}].l then
        r_right[e-{0,L-1}].l = r_left[e].l
        r_right[e-{0,L-1}].f = r_left[e].f
        check += 1
      elseif r_left[e].l > r_right[e-{0,L-1}].l then
        r_left[e].l = r_right[e-{0,L-1}].l
        r_left[e].f = r_right[e-{0,L-1}].f
        check += 1
      end
    end
  end
  return check
end


task BondEdgex(r_left : region(ispace(int2d), Grid2D),
           r_right : region(ispace(int2d), Grid2D),
           r_rng : region(ispace(int2d), c.drand48_data[1]),
           p : double, L : int32)
where
  reads(r_left.s), writes(r_left.bx),
  reads(r_right.s),
  reads writes(r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  for e in r_left do
    if r_left[e].s == r_right[e-{L-1,0}].s then
      var flip : double = uniform(data)
      if flip > p then
        r_left[e].bx = 1
      end
    end
  end
end

task BondEdgey(r_left : region(ispace(int2d), Grid2D),
           r_right : region(ispace(int2d), Grid2D),
           r_rng : region(ispace(int2d), c.drand48_data[1]),
           p : double, L : int32)
where
  reads(r_left.s), writes(r_left.by),
  reads(r_right.s),
  reads writes(r_rng)
do
  var data : &c.drand48_data
  for e in r_rng do
    data = [&c.drand48_data](@e)
  end

  for e in r_left do
    if r_left[e].s == r_right[e-{0,L-1}].s then
      var flip : double = uniform(data)
      if flip > p then
        r_left[e].by = 1
      end
    end
  end
end



task Mag2D(r_grid : region(ispace(int2d), Grid2D), L : int32)
where
  reads(r_grid.s)
do
  var N : int32 = L*L
  var mag : double 
  for e in r_grid do
    mag += r_grid[e].s
  end

  return mag/N
end


task toplevel()
  var config : EdgeConfig
  config:initialize_from_command()

  -- Create a logical region for Grid
  var L : int32 = 9
  var dim : int32 = 2
  var size : int2d = {L,L}
  var r_grid = region(ispace(int2d, size), Grid2D) 
  -- Create an equal partition of the grid
  var p_grid_colors = ispace(int2d, factorize2d(config.parallelism))
  var p_grid = partition(equal, r_grid, p_grid_colors)
  var r_rng = region(p_grid_colors, c.drand48_data[1])
  var p_rng = partition(equal, r_rng, p_grid_colors)

  -- Create partition for boundaries (For interior regions)
  var c_bound1 = coloring.create()
  var c_bound2 = coloring.create()

  for color in p_grid_colors do
    var bounds = p_grid[color].bounds
    if bounds.hi.x < L-1 then
      var bound1 : rect2d  = { {bounds.hi.x, bounds.lo.y}, {bounds.hi.x+1, bounds.hi.y} }
      coloring.color_domain(c_bound1, color, bound1)
    end
    if bounds.hi.y < L-1 then
      var bound2 : rect2d =  { {bounds.lo.x, bounds.hi.y}, {bounds.hi.x, bounds.hi.y+1} }
      coloring.color_domain(c_bound2, color, bound2)
    end
  end

  var p_bound1 = partition(disjoint, r_grid, c_bound1, p_grid_colors)
  var p_bound2 = partition(disjoint, r_grid, c_bound2, p_grid_colors)

  -- Create partition for boundaries (For Periodic Boundaries)
  var c_edgex1 = coloring.create()
  var c_edgex2 = coloring.create()
  var c_edgey1 = coloring.create()
  var c_edgey2 = coloring.create()


  var dims = ispace(int1d, dim)
  __forbid(__parallel)
  for d in dims do
    if d == int1d {0} then
      c.printf("d == 0\n")
      var edge1 : rect2d = { {0, 0}, {0, L-1} }
      var edge2 : rect2d = { {L-1, 0}, {L-1, L-1} }
      coloring.color_domain(c_edgex1, d, edge1)
      coloring.color_domain(c_edgex2, d, edge2)
    end

    if d == int1d {1} then
      c.printf("d == 1\n")
      var edge1 : rect2d = { {0, 0}, {L-1, 0} }
      var edge2 : rect2d = { {0, L-1}, {L-1, L-1} }
      coloring.color_domain(c_edgey1, d, edge1)
      coloring.color_domain(c_edgey2, d, edge2)
    end

  end
  
  var p_xedge1 = partition(disjoint, r_grid, c_edgex1, dims)
  var p_xedge2 = partition(disjoint, r_grid, c_edgex2, dims)

  var p_edgex = partition(equal, p_xedge1[int1d {0}], p_grid_colors) --factorize(d-1)D
  var p_edgex2 = partition(equal, p_xedge2[int1d {0}], p_grid_colors)

  var p_yedge1 = partition(disjoint, r_grid, c_edgey1, dims)
  var p_yedge2 = partition(disjoint, r_grid, c_edgey2, dims)

  var p_edgey = partition(equal, p_yedge1[int1d {1}], p_grid_colors) --factorize(d-1)D
  var p_edgey2 = partition(equal, p_yedge2[int1d {1}], p_grid_colors) --factorize(d-1)D


  var t : int8 = 1
  var v : int8 = 3 
  var V : double[1][3] -- number of temps, number of variables computed
  var TS_start = c.legion_get_current_time_in_micros()  
   
  for temp = 1,t do
  var T : double = (4.0-2.0)/t*temp + 2.0
  var p : double = cmath.exp(-2/T)

  
  --Initialization
  var token : int32 = 0
  for color in p_grid_colors do
    initialize(p_grid[color],p_rng[color], p)
    token += block_task(p_grid[color])
  end
  wait_for(token)

  var iter : int32 = 0
  var maxiter : int32 = 10
 
  var mag : double = 0 

  while iter < maxiter do
    iter += 1

    for color in p_grid.colors do
      SubInitialize(p_grid[color],p_rng[color],L,p)
    end
    for color in p_grid.colors do
      Bondsx(p_grid[color], p_rng[color], p)
    end
    for color in p_grid.colors do
      Bondsy(p_grid[color], p_rng[color], p)
    end
    
    var check : int32 = 0
    for color in p_grid.colors do
      check += Local(p_grid[color])
    end
    --c.printf("Check = %d, iter = %d\n", check,iter)
    while not check == 0 do
      check = 0
      for color in p_grid.colors do
        check += Local(p_grid[color])
      end
    end 

    for color in p_bound1.colors do
      Boundsx(p_bound1[color], p_rng[color], p, L)
    end
    for color in p_edgex.colors do 
      BondEdgex(p_edgex2[color], p_edgex[color], p_rng[color], p, L)
    end
    for color in p_edgex.colors do 
      BondEdgey(p_edgey2[color], p_edgey[color], p_rng[color], p, L)
    end
    for color in p_bound2.colors do
      Boundsy(p_bound2[color], p_rng[color], p, L)
    end


    check = 1
    var checkiter = 0
    while check > 0 or checkiter > 20 do 
      c.printf('Main Check\n')
      check = 0   
      checkiter +=1 
      var checkx : int32 = 0
      var checky : int32 = 0
      var checkxe : int32 = 0
      var checkye : int32 = 0
      for color in p_bound1.colors do
        checkx += Globalx(p_bound1[color])
      end
      for color in p_bound2.colors do
        checky += Globaly(p_bound2[color])
      end

      for color in p_grid_colors do
        checkxe += Edgex(p_edgex2[color], p_edgex[color],L)
      end
      for color in p_grid_colors do
        checkye += Edgey(p_edgey2[color], p_edgey[color],L)
      end

      check += checkx
      check += checky
      check += checkxe 
      c.printf("checkx = %d, checky = %d, checkxe = %d, checkye = %d\n", checkx, checky, checkxe, checkye)
      if check > 0 then
        c.printf('Local Check\n')
	var ckloc : int32 = 1
        while ckloc > 0 do
          ckloc = 0
          for color in p_grid.colors do
            ckloc += Local(p_grid[color])
          end
        end 
      end
    end
  
    for color in p_grid_colors do
      Flip(p_grid[color])
    end

    for color in p_grid_colors do     
      mag += Mag2D(p_grid[color],L)
    end
    
  
  end --this one ends the 'while T < endtime' loop
  
  for color in p_grid_colors do
    token += block_task(p_grid[color])
  end
  wait_for(token)

  V[temp][0] = mag/t
  c.printf('T = %f, mag = %f\n', T, mag/t)


  end --this one ends the temp sweep loop
  var TS_end = c.legion_get_current_time_in_micros()
  c.printf("Done\n")
  c.printf("Total time: %.6f sec.\n", (TS_end - TS_start) * 1e-6)

  --TODO Save Grid
  --var filename = "density"
  --attach(hdf5, r_fake.U0, filename, regentlib.file_read_write)  
  --release(r_fake.U0)
  --detach(hdf5, r_fake.U0)
  --print_density(r_grid[0])

  Dump(r_grid)
end

regentlib.start(toplevel)
