import "regent"

local coloring_util = {}

coloring_util.create = regentlib.c.legion_domain_point_coloring_create
coloring_util.destroy = regentlib.c.legion_domain_point_coloring_destroy
coloring_util.color_domain = regentlib.c.legion_domain_point_coloring_color_domain

return coloring_util
