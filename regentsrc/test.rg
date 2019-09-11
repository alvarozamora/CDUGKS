import "regent"

local c = regentlib.c

task main()
  var a : double = 1.0

  for i = 0, 10 do
    c.printf("a = %f\n", a)
  end
end

regentlib.start(main)
