Using CINDER to estimate gamma background from activation of structural materials

The use of the Cinder'90 activation code is thouroghly described in the manual that accompanies the software.
For the specific use using the SINE2020 project, the following steps were carried out.
 - The cells constituing the experimental cave were identified: 4351, 4352 & 4351
 - Using this, and activation input file is defined: inpact
 - and run using the command: acivation inpact
 - The activation calculation results in directories containing the information on activity of the individual cells.
 - To convert into a gamma source for use in MCNPX, the following command is issued: gamma_source gamma
 - where the input file: gamma defines which irradiation step is modelled.
 - The resulting source is stored in a file called: sdef, which can then be used in a gamma MCNPX simulation
