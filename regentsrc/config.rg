import "regent"

local c = regentlib.c

struct Config
{
  testproblem : int32,
  cpus  : int32,
  out : bool
}

local cstring = terralib.includec("string.h")

terra print_usage_and_abort()
  c.printf("Usage: regent edge.rg [OPTIONS]\n")
  c.printf("OPTIONS\n")
  c.printf("  -h            : Print the usage and exit.\n")
  c.printf("  -t {value}    : Set the threshold to {value}.\n")
  c.printf("  -p {value}    : Set the number of parallel tasks to {value}.\n")
  c.printf("  -o {bool}     : Output Boolean")
  c.exit(0)
end

terra Config:initialize_from_command()
  self.testproblem = -1
  self.cpus = 1 
  self.out = true

  var args = c.legion_runtime_get_input_args()
  var i = 1
  while i < args.argc do
    if cstring.strcmp(args.argv[i], "-h") == 0 then
      print_usage_and_abort()
    elseif cstring.strcmp(args.argv[i], "-t") == 0 then
      i = i + 1
      self.testproblem = c.atof(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-p") == 0 then
      i = i + 1
      self.cpus = c.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-o") == 0 then
      i = i + 1
      self.out = [bool](c.atoi(args.argv[i]))
    end
    i = i + 1
  end
end

return Config
