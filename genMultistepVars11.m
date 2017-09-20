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

u = zeros(1,1);
v = zeros(1,1);

R=2*beta_z;

C=1/R^2-1j*dz*(1-alpha)*1/R;
D=1/R^2+1j*dz*( alpha )*1/R;

u(1)= C;
v(1)= D;


end
