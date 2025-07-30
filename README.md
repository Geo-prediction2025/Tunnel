Use of forward modeling program
This program is written in Fortran and is used for advanced geological prediction. It requires compiling the environment first, and then using the tetgen grid subdivision software. The anisotropy control file is in the txt file res, and the main program is in the main file.Compile and run the program on the wondow10 system.
The use of 1Tetgen
The surface measurement point node should have it, and it should be declared again in the facet
When there are multiple layers, they can be made in blocks or layers
Tetgen –pq1.2A *.poly
Tetview *.poly
The divided grid is likely to be uneven
2 one-time program usage
Need to modify input file name
Modify the source ID
You decide how to output by yourself, write it yourself
A performance does not necessarily require starting from -1
In the result process, in order to generate inversion data, the source, data, and receiver need to be deleted first. Txt consists of 3 files in total
Need to modify the size of the stack
3 Surfer draws contour maps
Click on Grid/Data to input the txt file
Click on map/new/caption to output contour map
Modify the project stack size. Right clicking on 'Project', 'Properties', linker, systerm is not a solution
If solved using a quadratic field. Define analysis in the secondary field. In F90, it is necessary to determine whether a field is in full space or half space.
