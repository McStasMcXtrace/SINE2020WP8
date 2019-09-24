# SIMRES-MCStas integration scheme.

To enable SIMRES continuation of _any_ McStas-simulation it is sufficent to
enrich the McStas-simulation with a `MCPL_output.comp` positioned at the
switching point.
With the present releases of McStas and SIMRES the integration is one-way, but
the upcoming release of SIMRES _will_ include similar support for going the
other way around. (The software is written but awaits release)

## Example of a McStas-front end and SIMRES backend.

## Example of integrating a McStas-sample model in SIMRES.

## Example of SIMRES front-end and McStas-Mantid backend.

## Full exmaple of SIMRES-McStas-SIMRES using BEER as a vehicle

Included in the new release of SIMRES is an example that automates the connection process to take advantage of  both packages. In this example the proposed BEER instrument of ESS is simulated. Firstly the primary spectromter is simulated using a reverse Monte Carlo process. Secondly, neutron events are injected into a McStas simulation of a powder sample through and MCPL-write/read procedure. Secondly the McStas simulation lets the events interact with a the powder. Lastly, nuetron events are reinjected into SIMRES through a second MCPL-write/read pair, and the secondary spectrometer is simulated with SIMRES.  
