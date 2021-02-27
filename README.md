# Electrodynamics-Simultator
A simulator designed to be able to compute and visulize electric and magnetic fields,
using a mathematical expression for the charges and currents as it's input.

The following modules are required for the program to work.
1.numpy: Used for all computational work, and especially so in helpers.py
2.vispy: The plotting library of choice, the more popular option of matplotlib 
is not performace efficient for trying to visualize the behaviour of 
changing electric and magnetic fields in real time
3.pyqt5: Used for the graphical user interface and to embed the plot with a set
of controls

optional dependencies:
1.pyqt5-tools: This module contains the program 'QtDesigner' which is necessary to
edit the GUI but otherwise is not needed to simply run the program as is. 
Another option would be to convert the included .ui file into a python file and 
manually make edits to it, but the QtDesigner approach is recommended.

Issues:
The program can display stationary fields well but changing fields have not been 
implemented completely, for example the induced electric field is computed,
but the induced magnetic field is not and only a placeholder function for it is presen.
Also, the time range and the timeStep are values that need to be changed in the source code.
Ideally those would be accessible via the GUI.

The use of np.roll for the iteration in inverseLaplacian means that periodic boundary conditions
are introduced. If this is not desireable then it is recommended to keep the objects away from the
edges of the boundary

The mesh is set to be uniform in its step size on all axis, and its size needs to be 
adjusted via editing the souce code while ideally it should be available in the GUI

Since the vector plotting in vispy is based on the GPU, on some older graphics cards
or on older integrated graphics one might find that the heads of the vectors are not 
rendered properly. This is a vispy issue and to best of my knowledge,
there is no workaround.

Operation Guide:
Different types of objects:
There are two different objects one can create.
1. Or objects: These are objects which are definied by a single equation for their position.
Example: A sphere which can be defined as x^2+y^2+z^2<1 or using r<1 where r is the radial spherical co-ordinate
2. And objects: Definied by more than one equation for the position.
Example: A hemisphere definied by 'r<1 and z>0' In essence you are performaning a 3D Venn Diagram intersection
where all of the conditions need to be true at the same time

Multiple objects need to be seperated by the keyword 'or' 
Example: 'r<1 or s<1' where s is the radial cylindrical co-ordinate

About charges and currents

The charge strength at a point is a scalar so it is sufficient to define it either via a constant say '1'
or some function of the position variables such as 'z' or 'x+y+z' with the charge strength for every 
object seperated via the keyword 'or' 
Since and objects define the same object overall, one does not need to enter the charge info multiple times

Example: A valid input would be 
Position Expr: r<1 and z>0
Charge Expr: 1

Since the current at a point is a vector, one needs to define it as such and thus the input needs to be of the form,
'a,b,c' where a,b,c are either constants or some function of the coordinate variables

Example: A valid input would be
Position Expr: s<1
Current Expr: 0,0,z

The following are some settings and toggles

Normalize: If checked, it will normalize the lengths of the vectors to a fixed constant
Colour: If checked, each vector will be assigned a colour based on its magnitude

Length: The slider controls the length of the vectors via a scaling effect in both normalized and
not normalized modes
Step: Controls the step size of the generated mesh, lower meshes would be more accurate due to a higher number of vectors
but would be much more expensive computationally. 

As an example, if the step size is halfed then the number of vectors increases by a factor of 2^3 or 8
i.e for x : (1/x)^3 where x is the scaling of the step size

Notes: The program was meant to be a tool to help visualize electric and magnetic fields given a set of charges or currents
without having to compute the fields analytically and then plotting them. 

Possible Improvements: The GUI could use a reworking to also add multiple different variables to be controllable.
Such as epsilon and mu, the physical constants for charges and currents

The ability to save the results of the computation, so as to allow larger computations and then plot them or visualize 
them when needed without having to re-run the computation again

The inclusion for epsilon and mu to be variable to allow for visualizing the effect of changes in materials on the fields.
