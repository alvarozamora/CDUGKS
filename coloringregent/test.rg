import "regent"

local c = regentlib.c

task main()
  var N : int32[3]
  N[0] = 16
  N[1] = 1
  N[2] = 1

  
  var isp = ispace(int3d, {N[0], N[1], N[2]}) 

end

regentlib.start(main)
