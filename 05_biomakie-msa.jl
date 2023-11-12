## https://github.com/BioJulia/BioMakie.jl
## https://biojulia.dev/BioMakie.jl/dev/msaselection/
using BioMakie
using MIToS
using MIToS.MSA, MIToS.Pfam
using GLMakie
using Lazy
downloadpfam("pf00062")
msa1 = MIToS.MSA.read("pf00062.stockholm.gz",Stockholm)
msa2 = Observable(msa1)
plotdata = plottingdata(msa2)
fig = Figure(resolution = (1400,400))
msa = plotmsa!(fig, plotdata)
coldata = lift(plotdata[:selected]) do sel
    try
        plotdata[:matrix][][:,parse(Int,sel)]
    catch
        ["-" for i in 1:size(plotdata[:matrix][])[1]]
    end
end
allaas = [  "R", "M", "N", "E", "F",
            "I", "D", "L", "A", "Q",
            "G", "C", "W", "Y", "K",
            "P", "T", "S", "V", "H",
            "X", "-"]
sortaas = sortperm(allaas)
new_aalabels = allaas[sortaas]
hydrophobicities = [BioMakie.kideradict[new_aalabels[i]][2] for i in 1:length(new_aalabels)]
countmap1 = @lift frequencies($coldata) |> sort
aas = @lift collect(keys($countmap1))
freqs = lift(aas) do a
    collect(values(countmap1[]))
end
missingaas = @lift setdiff(allaas,$aas) |> sort
missingfreqs = @lift zeros(length($missingaas))
perm1 = @lift sortperm([$aas; $missingaas])
aafreqs = @lift ([freqs[];$missingfreqs])[$perm1]
aafreqspercent = @lift $aafreqs ./ sum($aafreqs) .* 100
new_aafreqs = @lift $aafreqspercent[sortaas]
ax = Axis(fig[1,4], xticklabelsize = 16, yticks = (0:10:100), yticklabelsize = 20,
            title = "Amino Acid Percentages",
            titlesize = 18, xticks = (1:22,new_aalabels)
)
bp = barplot!(ax, 1:22, aafreqspercent; color = hydrophobicities, strokewidth = 1,
                xtickrange=1:22, xticklabels=new_aalabels
)
ylims!(ax, (0, 100))
xlims!(ax, (0, 23))
