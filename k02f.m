% function k02=k02f(k2,A)
%
% k02 is the zero-th order component of the Wiener series corresponding to 
% the second order Wiener kernel k2
% A it's the power of input sequence
% Refer to the documentation and references provided
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

function k02=k02f(k2,A)

R=size(k2,1);
diag=zeros(R,1);
    
for tau1=1:R,diag(tau1)=k2(tau1,tau1);end
k02=-A*sum(diag);
    