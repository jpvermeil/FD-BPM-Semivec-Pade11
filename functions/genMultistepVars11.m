function [u,v]=genMultistepVars11(dz,alpha,beta_z)
% This function generates the multistep Varialbes 'u' and 'v' for the
% Padé(1,1) finite difference BPM. For further information on the
% calculation procedure see: the book: "Introduction to optical waveguide
% analysis" by K. Kawano and T. Kitoh.
%
% SYNOPSIS
% 
% [u,v] = genMultistepVars11(dz,alpha,beta_z)
% 
% VARIABLES
% 
% dz          Step in z-direction
% alpha       Derivation paramter of Crank-Nicolson scheme
% beta_z      Propagation constant

% This function generates the multistep Varialbes 'u' and 'v' for the
% Padé(1,1) finite difference BPM. For further information on the
% calculation procedure see: the book: "Introduction to optical waveguide
% analysis" by K. Kawano and T. Kitoh. This function belongs to the program
% "Semi vectorial wide angle Padé(1,1) finite difference BPM for TE/TM E-
% and/or H-fields in 3D structures".
% Copyright (C) 2017 Jan-Philipp Roth (JanPhilipp.Roth@gmail.com)
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. You should have received a
% copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth
% Floor, Boston, MA 02110-1301  USA

u = zeros(1,1);
v = zeros(1,1);

R=2*beta_z;

C=1/R^2-1j*dz*(1-alpha)*1/R;
D=1/R^2+1j*dz*( alpha )*1/R;

u(1)= C;
v(1)= D;


end
