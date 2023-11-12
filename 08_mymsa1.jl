# https://docs.makie.org/stable/reference/blocks/slider/index.html
# https://discourse.julialang.org/t/makie-multi-sequence-alignment-scatterplot/77932
# https://github.com/BioJulia/BioMakie.jl/blob/master/src/msa.jl
# ~/.julia/dev/BioMakie
using GLMakie
#using WGLMakie
string_length = 25
font_size = 12
cell_size = 20
x_cells = 10
y_cells = 2

axis_x = cell_size * x_cells
axis_y = cell_size * y_cells

fig_size = (2*axis_x,2*axis_y)

s1 = string.(rand('A':'Z',string_length));
s2 = string.(rand('A':'Z',string_length));

fig = Figure(fontsize = font_size, backgroundcolor = :gray95, resolution = fig_size)
ax = Axis(fig[1,1], xticks = 1:x_cells, yticks = 1:y_cells, yreversed=true, width=axix_x, hegith=axis_y)

# Scatter based
x_pos = (1:x_cells .- .5) .* cell_size
y_pos = (1:y_cells .- .5) .* cell_size

scatter!(ax, x_pos, y_pos)
# text based
text!(ax, x=1:, y=repeat([1],25), text=s1)
