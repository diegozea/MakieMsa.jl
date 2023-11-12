## https://github.com/BioJulia/BioMakie.jl
## https://biojulia.dev/BioMakie.jl/dev/msaselection/
using BioMakie
using MIToS
using MIToS.MSA, MIToS.Pfam
using GLMakie
#using WGLMakie
using Lazy
downloadpfam("pf00062")
msa1 = MIToS.MSA.read("pf00062.stockholm.gz",Stockholm)
fig = Figure(resolution = (1400,400))
plotmsa(msa1) 
plotmsa!(fig, Observable(msa1); gridposition=(2,1:3), sheetsiae=[80,20 ])
