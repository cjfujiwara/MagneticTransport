# MagneticTransport
Author : CJ Fujiwara

This repository is a set of MATLAB code designed to calculate the magnetic transport stage for the ultracold atom experiment at the Univeristy of Toronto in the group of Joseph Thywissen.  The transport is divided into a horizontal transport and a vertical transport stage.  The two stages are requires different current traces as vertical transport respects the symmetry of a electromagnet coil while horizontal transport does not.

## References
The primary references for this calculation are :

Michael Yee. "Magnetic Trapping and Transport of Ultracold Atoms". 2009.  
David McKay. "Quantum Simulation in Strongly Correlated Optical Lattices". 2012  
Dylan Jervis. "A Fermi Gas Microcope Apparatus". 2014  
https://tiggerntatie.github.io/emagnet/offaxis/iloopoffaxis.htm  
https://www.mathworks.com/help/optim/ug/fmincon.html  

## Notes
The physics of magnetic transport is not discussed here, and this repository directs you to the relevant references.

## Important Functions
makeVerticalCoils.m                            : initialize structure arrays which describe the configuration of vertical coils  
makeHorizontalCoils.m                          : initialize structure arrays which describe the configuration of horizontal coils  
calculate_horizontal_transport.m               : calculate horizontal transport  
calculate_vertical_transport.m                 : calculate vertical transport  
fieldCoil_3D.m                                 : calculate the off-axis field of a single loop of wire  

## Horizontal Output
horizontal_output/horizontal_current.mat       : binary file of current traces for horizontal transport  
horizontal_output/horizontal_current.png       : image of current traces for horizontal transport  
horizontal_output/horizontal_field_profile.png : image of field profile for each coil pair  

## Vertical Output
vertical_output/vertical_current.mat            : binary file of current traces for vertical transport  
vertical_output/vertical_current.png            : image of current traces for vertical transport  
vertical_output/vertical_field_profile.png      : image of field profile for each vertical vertical  
