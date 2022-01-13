function [ PX,PY ] = checkFieldComponents( FIELDCOMPONENTS )
% This function checks which field components need to be calculated.
%
% SYNOPSIS
% 
% [PX,PY] = checkFieldComponents(FIELDCOMPONENTS)
% 
% VARIABLES
% 
% FIELDCOMPONENTS   Can be either 'Ex','Ey','Hx' oder 'Hy'

% This function checks which field components need to be calculated. It
% belongs to the program "Semi vectorial wide angle Padé(1,1) finite
% difference BPM for TE/TM E- and/or H-fields in 3D structures".
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

PX = NaN;
PY = NaN;

if (strcmp('Ex',FIELDCOMPONENTS{1}) || strcmp('Ey',FIELDCOMPONENTS{1}) || strcmp('Hx',FIELDCOMPONENTS{1}) || strcmp('Hy',FIELDCOMPONENTS{1})) ~= 1

out = ['WARNING: Invalid first FIELDCOMPONENT ''' FIELDCOMPONENTS{1} '''. Must either be ''Ex'',''Ey'',''Hx'' or ''Hy''.'];
disp(out);
return
    
end

if (strcmp('Ex',FIELDCOMPONENTS{2}) || strcmp('Ey',FIELDCOMPONENTS{2}) || strcmp('Hx',FIELDCOMPONENTS{2}) || strcmp('Hy',FIELDCOMPONENTS{2})) ~= 1
    
out = ['WARNING: Invalid second FIELDCOMPONENT ''' FIELDCOMPONENTS{2} '''. Must either be ''Ex'',''Ey'',''Hx'' or ''Hy''.'];
disp(out);
return
    
end

PX = FIELDCOMPONENTS{1};
PY = FIELDCOMPONENTS{2};

end

