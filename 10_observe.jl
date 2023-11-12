using GLMakie

x = 1:4
y = rand(4)
lw = Observable(2)

fig,ax = lines(x,y;linewidth=lw)

function myplot1(a)
    x = -2:.1:2
    y = a .* x.^2
    lines(x,y)
end

julia> myplot1(1)

julia> myplot1(2)

julia> a = Observable(1.0)

julia> myplot1(a) ## fails

x = Observable(-2:.1:2)
y = @lift($a .* $x.^2)

f,ax,p = lines(x,y)

julia> a[] = .1 # This works!

## So what goes into the plot must be observable.

f = Figure()
ax = Axis(f[1,1])
s = Slider(f[2,1], range = 0:.1:4, startvalue = 1)

x = Observable(-2:.1:2)
y = @lift($(s.value) .* $x.^2)

lines!(ax,x,y)

z = @lift -$y
scatter!(Axis(f[3,1]), x, z)

function myplot!(ax, x, y)
    z = @lift([z1 + z2 for  z1 in $x , z2 in $y])
    heatmap!(ax, x,x, z)
end

myplot!(Axis(f[4,1]), x,y)

function myplot1!(ax, x, s)
    z = @lift([$(s.value) * z1 * z2 for  z1 in $x , z2 in $x])
    heatmap!(ax, x,x, z)
end

myplot1!(Axis(f[5,1]), x,s)

function myplot2!(ax, x, s)
    z = @lift([sin.($(s.value) * z1)  for  z1 in $x , z2 in $x])
    heatmap!(ax,x,x, z)
#    Colorbar(f[:, end+1], hm)
end

myplot2!(Axis(f[6,1]), x,s)

function myplot3!(ax,x,s)
    y = @lift($(s.value) .* $x.^2)
    lines!(ax,x,y, color =:red)
end


myplot3!(Axis(f[7,1]), x,s)
