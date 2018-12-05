% function swap=swapassess(N,swaptable)
% It should be used in combination with a swap table obtained with demo/swaptest
% script which generates a swaptable.mat containing a swaptable array, made
% of two columns the first of which being the vector of input lengths and the the
% second being the corresponding swap values obtained, example:
% swaptable =
% 
%          100           2
%          500           4
%         1000           5
%         5000          25
%        10000          31
%        50000          17
%       100000          10
%       500000          10
%      1000000           9
%
% swaptable can be the pathname of file containing swaptable array 
% 	 or a swaptable array ready to be used
% N is the input length for which a swap value is desired
% swap is then obtained by linear interpolation of the values in swaptable,
%     if N it's in the interval [swaptable(1,1), swaptable(end,1)], or
%     using the extremes values otherwise
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DII, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.

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

function swap=swapassess(N,swaptable)
if ischar(swaptable)
    swaptablefile=swaptable;
    F=load(swaptablefile);
    swaptable=F.swaptable;
    clear F;
end

S=size(swaptable);
if S(2)~=2, error('wrong size of swaptable array'),end

x=swaptable(:,1);
y=swaptable(:,2);

if N>=x(end)
    swap=round(y(end));
elseif N<=x(1)
    swap=round(y(1));
else
    xi=1:x(end);
    yi=interp1(x,y,xi,'linear');
    swap=round(yi(N));
end


