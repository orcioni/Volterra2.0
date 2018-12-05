% function [k35]=k35f(k3,A)
%
% k35 is the third order component of the Wiener series corresponding to 
% the fifth order Wiener kernel k5.
% A it's the power of input sequence.
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

function k35=k35f(k5,A)

R5=size(k5,1);
diag=zeros(R5,1);
k35=zeros(R5,R5,R5);
for sgm1=1:R5
    for sgm2=1:R5;
        for sgm3=1:R5
            for tau1=1:R5,diag(tau1)=k5(tau1,tau1,sgm1,sgm2,sgm3);end
            k35(sgm1,sgm2,sgm3)=sum(diag);
        end
    end
end

k35=-10*A*k35;
