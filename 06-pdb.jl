using BioMakie
using GLMakie
using BioStructures

struc = retrievepdb("2vb1") |> Observable
fig = Figure()
plotstruc!(fig, struc; plottype = :ballandstick, gridposition = (1,1), atomcolors = aquacolors)
plotstruc!(fig, struc; plottype = :covalent, gridposition = (1,2))
plotstruc!(fig, struc; plottype = :spacefilling, gridposition = (1,3))
struc2 = retrievepdb("1svn") |> Observable
plotstruc!(fig, struc2; plottype = :ballandstick, gridposition = (2,1))
