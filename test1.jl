using GLMakie

# from https://discourse.julialang.org/t/makie-multi-sequence-alignment-scatterplot/77932

msa_matrix = permutedims([
    '-' 'G' 'G' '-' 'C' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' '-' '-' 'G' 'T'
    '-' 'G' 'G' 'C' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' 'G' 'T' 'G' 'G' 'T'
    'C' 'G' 'G' '-' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' 'G' '-' '-' '-' 'T'
    'C' 'G' 'G' '-' 'T' 'T' 'G' 'C' 'T' 'T' 'A' 'T' 'T' 'G' 'T' 'G' '-' '-' 'T' '-'
])

levels = Dict(
    'A' => 1,
    'C' => 2,
    'G' => 3,
    'T' => 4,
    '-' => 5,
)

colormap = tuple.([:yellow, :blue, :green, :red, :magenta], 0.5)

cellsize = (20, 20)

f = Figure(fontsize = 12, backgroundcolor = :gray95)

axis_size = size(msa_matrix) .* cellsize

ax = Axis(
    f[1, 1],
    xgridvisible = false,
    ygridvisible = false,
    xticks = 1:size(msa_matrix, 1),
    yticks = (1:size(msa_matrix, 2), ["a", "b", "c", "d"]),
    yticksvisible = false,
    yreversed = true,
    width = axis_size[1],
    height = axis_size[2],
)

hidespines!(ax)

heatmap!(ax, [levels[m] for m in msa_matrix], colormap = colormap)
text!(ax, vec(string.(msa_matrix)),
    position = vec(Point2f.(Tuple.(CartesianIndices(msa_matrix)))),
    align = (:center, :center),
    textsize = 12,
)

resize_to_layout!(f)

f
