%function [equal,err]=nanisequal(a,b)
%
% Compares two arrays a and b, which can contain NaNs.
% equal is a flag which is 1 if the arrays are identical.
% If equal is 0, err will contain the maximun of absolute difference
% between the arrays. 
% See isequal for other information.

% Copyright (C) 2006 Massimiliano Pirani
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License along
%  with this program; if not, write to the Free Software Foundation, Inc.,
%  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


function [equal,err]=nanisequal(a,b)

 equal=isequal(NaN2eps(a),NaN2eps(b));
 
 if equal==1
     err=0;
 else
     S=size(a);
     P=prod(S);
     err=max(abs(reshape(NaN2zero(a)-NaN2zero(b),P,1)));
 end