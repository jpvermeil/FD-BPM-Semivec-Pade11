function [ aSouthReshaped ] = padeS( n,yg,dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,FIELDCOMPONENT )
% This function populaltes the coefficient vector that is merged into the
% system matrix A as the respective southern diagonal of for each step in
% z-direction.
%
% SYNOPSIS
% 
% [ aSouthReshaped ] = padeS(n,xg,dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT)
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
%   localAdrS       Contains the local address information of southern elements
%   POLARIZATION    String value that can either be 'TE' or 'TM'
%   FIELDCOMPONENTS Can be 'Ex', 'Ey', 'Hx' or 'Hy'

% This function populaltes the coefficient vector that is merged into the
% system matrix A as the respective southern diagonal of for each step in
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
    
    dn = yg(globalAdrS - 1) - yg(globalAdrS);
    ds = yg(globalAdrS) - yg(globalAdrS + 1);
    
    % Round off error
    dn = 1e-12*round(dn*1e12);
    ds = 1e-12*round(ds*1e12);
    
    epsilon = n.*n; 
    
    if (strcmp('TE',POLARIZATION)) == 1
        
        a_s = (2./(ds.*(ds+dn))) .* ones(length(globalAdrS),1);
        
    elseif (strcmp('TM',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1
        
        a_s = (2./(ds.*(ds+dn))) .* (2*epsilon(globalAdrS) ./ (epsilon(globalAdrS) + epsilon(globalAdrS + 1))); 
    
    elseif (strcmp('TM',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1
        
        a_s = (2./(ds.*(ds+dn))) .* (2*epsilon(globalAdrS + 1) ./ (epsilon(globalAdrS) + epsilon(globalAdrS + 1)));
        
    else
        
        out = 'WARNING: FD coefficients a_s could not be computed. Check POLARIZATION and FIELDOCMPONENT string.';
        disp(out);
       
    end
    
    aSouthReshaped = zeros(dim_xl*dim_yl,1);
    aSouthReshaped(localAdrS) = a_s;
    
end