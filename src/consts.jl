## Clustal colors: https://www.jalview.org/help/html/colourSchemes/clustal.html
const clustal_levels = OrderedDict(
    'A' => 1,
    'B' => 9,
    'C' => 6,
    'D' => 3,
    'E' => 3,
    'F' => 1,
    'G' => 6,
    'H' => 8,
    'I' => 1,
    'J' => 9,
    'K' => 2,
    'L' => 1,
    'M' => 1,
    'N' => 4,
    'O' => 9,
    'P' => 7,
    'Q' => 4,
    'R' => 2,
    'S' => 4,
    'T' => 4,
    'U' => 9,
    'V' => 1,
    'W' => 1,
    'X' => 9,
    'Y' => 8,
    'Z' => 9,
    '-' => 9,
    '.' => 9,
    ' ' => 9,
    '*' => 9,
    'a' => 10,
    'c' => 11,
    'g' => 12,
    't' => 13,
)
const clustal_colormap = tuple.([:blue, :red, :magenta, :green, :pink, :orange, :yellow, :cyan, :white,:green, :blue, :black, :red], 0.5)


const dna_levels = OrderedDict(
    'A' => 1,
    'C' => 2,
    'G' => 3,
    'T' => 4,
)
const sanger_colormap = tuple.([:green, :blue, :black, :red], 0.5)

