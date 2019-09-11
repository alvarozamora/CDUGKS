import "regent"

local c=terralib.includec("stdio.h")

task f() 
  var r = region(ispace(int2d, {3,4}), int)

  for p in r do
    c.printf("(%d,%d)\n", p.x, p.y)
  end

end

task main()
  f()
end
regentlib.start(main)
