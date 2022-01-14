Finite difference semi-vectorial wide-angle beam propagation algorithm for TE- and/or TM-Polarization of E- and/or H-fields in optical waveguide structures with arbitrary index profile. Calculation scheme is based on a semi-vectorial finite difference approach utilizing an absorbing boundary condition.

To download official releases or rate the toolbox, please visit its page on [FileExchange](https://de.mathworks.com/matlabcentral/fileexchange/105225-fd-bpm-semivec-pade11).

Theoretical Reference
=======================

I implemented this BPM-Toolbox during my time as a PhD student to characterize graded-index waveguide components manufactured by thermal ion-exchange processes. Therefore, most of the tools and instruments provided focus mainly on graded-index refractive index profiles. Many functions won't function properly or as intended if other index profiles are used. Please check if methods work correctly for other profiles. Also, the toolbox' help section at the very beginning of the main function will give you some insight into basic usage and some explanations about capabilities and restrictions.

Furthermore, this BPM-Toolbox uses a finite difference based mode solver to numerically calculate specific modes. This can be used to investigate the propagation of a specific mode the waveguide structure. Note that his mode solver is included as a git submodules. For further information on the mode solver check the according [github](https://github.com/jpvermeil/FD-Mode-Solver) page.

For further reference with respect to the numerical implementation please refer to appropriate literature. And excellent and comprehensive work is the book 'Introduction to Optical Waveguide Analysis: Solving Maxwell's Equations and the Schrödinger Equation' by K. Kawano and T. Kitoh, which served as a theoretical basis for this implementation.

Installation
==============

1. Extract the ZIP file (or clone the git repository).
2. Add the folder to your `MATLAB` path. This can be done by:
  - Adding the folder within your default path
  - Or manually by using something like `addpath(genpath('.'))` from within the folder
3. If you intend to use the FD-based mode solver, which is included as submodule, please use:
  - `git submodule init` and
  - `git submodule update`

Usage
=====

The help section at the very beginning of the main function will give you some insight into basic usage and some explanations about capabilities and restrictions. Furthermore, there are some examples for a number of different waveguide types that will show and explain basic usage.

Remarks
=======

Please bear in mind, that this is a problem oriented implementation of the according wave equations that focuses on the investigation of the phenomenological side of waveguide characterization. By no means is this toolbox viable for reliable quantification of optical properties or anything similar.

The main drawback is the lack of any other boundary condition than an absorbing boundary condition. This especially problematic for waveguide structures near the surface of the investigated area. To a certain extent, this can be avoided by creating and attaching artificial claddings as shown in one example.

The branch `TBC` provides an experimental support for the transparent boundary condition but is neither debugged nor tested properly. This is especially problematic if step-index designs are investigated, since the artificial absorber only works well for graded-index profiles.

Development
===========

While the code is likely not to be maintained or supported to any extent in the future, there are some things that I might implement:
* Proper support for the `TBC`
* Implementation of a `PML` (Perfectly Matched Layer)
* Fully-vectorial solution for the investigation of polarization coupling
* ADI-Method (Alternate Direction Implicit) for faster computation
* Option to deactivate wide-angle approach
* Python port (since I'm not sure to have access to a `MATLAB` license in the future)
