function [ aEastReshaped ] = padeE( n,xg,dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT)
% This function populaltes the coefficient vector that is merged into the
% system matrix A as the respective eastern diagonal of for each step in
% z-direction.
%
% SYNOPSIS
% 
% [ aEastReshaped ] = padeE(n,xg,dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT)
%
% VARIABLES
%
%   n               Refractive index profile
%   xg              Grid of x dimension. Has to match dimensions of n. 
%                   Has to be X output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   yg              Grid of y dimension. Has to match dimensions of n.
%                   Has to be Y output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   dim_y           Size of the global matrix in y-direction
%   dim_xl          Size of the local matrix in x-direction
%   dim_yl          Size of the local matrix in y-direction
%   localAdrE       Contains the local address information of eastern elements
%   POLARIZATION    String value that can either be 'TE' or 'TM'
%   FIELDCOMPONENTS Can be 'Ex', 'Ey', 'Hx' or 'Hy'

% This function populaltes the coefficient vector that is merged into the
% system matrix A as the respective eastern diagonal of for each step in
% z-direction. It belongs to the program "Semi vectorial wide angle
% Padé(1,1) finite difference BPM for TE/TM E- and/or H-fields in 3D
% structures".
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
 
    de = xg(globalAdrE + dim_y) - xg(globalAdrE);
    dw = xg(globalAdrE) - xg(globalAdrE - dim_y);
    
    % Round off error 
    de = 1e-12*round(de*1e12);
    dw = 1e-12*round(dw*1e12);

    epsilon = n.*n;
    
    if (strcmp('TM',POLARIZATION)) == 1
        
        aEast = 2./(de.*(de+dw)); 
        
    elseif (strcmp('TE',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1
        
        aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrE) ./ (epsilon(globalAdrE) + epsilon(globalAdrE + dim_y))); 
    
    elseif (strcmp('TE',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1
        
        aEast = (2./(de.*(de+dw))) .* (2*epsilon(globalAdrE + dim_y) ./ (epsilon(globalAdrE) + epsilon(globalAdrE + dim_y)));
        
    else
        
        out = 'WARNING: FD coefficients a_e could not be computed. Check POLARIZATION and FIELDOCMPONENT string.';
        disp(out);
       
    end
    
    aEastReshaped = zeros(dim_xl*dim_yl,1);
    aEastReshaped(localAdrE) = aEast;

end

