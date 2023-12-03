"""
    mafft(input_file, output_file = tempname(); options = "--maxiterate 1000 --localpair", quiet = true)
    Run mafft on input_file with the options writing to output_file
    Return name of output_file
"""
function mafft(input_file, output_file = tempname(); options = "--auto", quiet = false)
  # High accuracy (for <~200 sequences x <~2,000 aa/nt):
  # % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok) L-INS-i (Probably most accurate, very slow)

  # % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out) E-INS-i (Suitable for sequences with long unalignable regions, very slow)

  # % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out) G-INS-i (Suitable for sequences of similar lengths, very slow)

    ## https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html
    ## For simple alignments, they give identical results

    if quiet
        options *= " --quiet"
    end
    _mafft = "mafft" ## TODO get from MAFFT_jll
    _options = split(options)
    cmd = `$(_mafft) $(_options)  $input_file`
    run(pipeline(cmd, stdout = output_file))
    output_file
end
