% function e=mse(x,y)
%
% e is the Mean Square Error between the two input arrays x and y

% Copyright (C) 2006 Massimiliano Pirani
% Copyright (C) 2021 Simone Orcioni
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

function e=mse(x,y)

if (nargin==1)
    y=0;
end

if ~(size(x)==size(y))
    y=y';
end
    
e=(x-y).^2;
e=mean(e,'all');
