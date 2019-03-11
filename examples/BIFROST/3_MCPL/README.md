Use ssw2mcpl to generate mcpl-formatted output from the SSW data generated in step 2

Usage:

  ssw2mcpl [options] input.ssw [output.mcpl]

Converts the Monte Carlo particles in the input.ssw file (MCNP Surface
Source Write format) to MCPL format and stores in the designated output
file (defaults to "output.mcpl").

Options:

  -h, --help   : Show this usage information.
  -d, --double : Enable double-precision storage of floating point values.
  -s, --surf   : Store SSW surface IDs in the MCPL userflags.
  -n, --nogzip : Do not attempt to gzip output file.
  -c FILE      : Embed entire configuration FILE (the input deck)
                 used to produce input.ssw in the MCPL header.

- More specifically run 
ssw2mcpl -s bifrost_comblayer_EK nps1E5_with_SM_at_2m_and_sample.w 

or equivalent
