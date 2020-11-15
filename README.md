# Neutron2D_DD-SC
Solves a 2D box problem with an isotropic source using either diamond difference or step characteristics

Contents
- meshcell.m : class that defines properties and location of meshcell. Also contains step characteristics and diamond difference methods
- angle.m : function that performs angular discretization
- StepCharacteristics.m : function that produces a plot of the solution to the problem solved used step characteristics
- DiamondDifference.m : function that produces a plot of the solution to the problem solved used diamond differences
- fig/ : contains several plots of solutions. Plot format is (material type)_(no. of cells per side)_(angular discretization in x-y)_(angular in z)_(location plotted).eps


The remaining files are not important and will eventually be deleted. 
