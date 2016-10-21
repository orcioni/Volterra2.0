% Copyright (C) 2006 Massimiliano Pirani
% Copyright (C) 2014 Simone Orcioni
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
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DEIT, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.
% If you are using this program for a scientific work, we encourage you to cite
% the following paper (the file cite.bib, containing the reference in bibtex
% format is also provided):
% Simone Orcioni, Massimiliano Pirani, and Claudio Turchetti. Advances in 
% Lee-Schetzen method for Volterra filter identification. Multidimensional 
% Systems and Signal Processing, 16(3):265-284, 2005.

% function [d2,d3, ...,d10]=DiagErase(h2,h3,...,h10)
%
% This function puts NaN values on diagonal points of h2 whatever the dimension
% of the array is (but less than 11), if it is provided as input alone.
% If a list of crescent order arrays/kernels is provided as argument the 
% output will be the corresponding resulting array list with diagonal points 
% marked as NaN.


function varargout=DiagErase(h2,h3,h4,h5,h6,h7,h8,h9,h10)

single=0;

switch nargin
case 1
    ord=length(size(h2));
    if ord>2,single=ord;eval(['h' num2str(ord) '=h2; clear h2']);end
otherwise
    ord=nargin+1;
end

if (ord>=2) && ((single==0)||(single==2))
    R=size(h2,1);
    for tau=1:R
        h2(tau,tau)=NaN;
    end
    varargout{1}=h2;
    clear h2;
end

if (ord>=3) && ((single==0)||(single==3))
    R=size(h3,1);
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                if (tau1==tau2)||(tau2==tau3)||(tau1==tau3)
                    h3(tau1,tau2,tau3)=NaN;
                end
            end
        end
    end
    if single==3, varargout{1}=h3;else varargout{2}=h3; clear h3;end
end

if (ord>=4) && ((single==0)||(single==4))
    R=size(h4,1);
    eqflag=0;
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                for tau4=1:R
                    ind=[tau1 tau2 tau3 tau4];
                    for i=1:4,
                        if length(find(ind(i)==ind))>1;eqflag=1;break,end
                    end
                    if eqflag==1,h4(tau1,tau2,tau3,tau4)=NaN;eqflag=0;end
                end
            end
        end
    end
    if single==4, varargout{1}=h4;else varargout{3}=h4; clear h4;end
end

if (ord>=5) && ((single==0)||(single==5))
    R=size(h5,1);
    eqflag=0;
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                for tau4=1:R
                    for tau5=1:R
                        ind=[tau1 tau2 tau3 tau4 tau5];
                        for i=1:5,if length(find(ind(i)==ind))>1;eqflag=1;break,end,end
                        if eqflag==1,h5(tau1,tau2,tau3,tau4,tau5)=NaN;eqflag=0;end
                    end
                end
            end
        end
    end
    if single==5, varargout{1}=h5;else varargout{4}=h5; clear h5;end
end

if (ord>=6) && ((single==0)||(single==6))
    R=size(h6,1);
    eqflag=0;
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                for tau4=1:R
                    for tau5=1:R
                        for tau6=1:R
                            ind=[tau1 tau2 tau3 tau4 tau5 tau6];
                            for i=1:6,if length(find(ind(i)==ind))>1;eqflag=1;break,end,end
                            if eqflag==1,h6(tau1,tau2,tau3,tau4,tau5,tau6)=NaN;eqflag=0;end
                        end
                    end
                end
            end
        end
    end
    if single==6, varargout{1}=h6;else varargout{5}=h6; clear h6;end
end

if (ord>=7) && ((single==0)||(single==7))
    R=size(h7,1);
    eqflag=0;
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                for tau4=1:R
                    for tau5=1:R
                        for tau6=1:R
                            for tau7=1:R
                                ind=[tau1 tau2 tau3 tau4 tau5 tau6 tau7];
                                for i=1:7,if length(find(ind(i)==ind))>1;eqflag=1;break,end,end
                                if eqflag==1,h7(tau1,tau2,tau3,tau4,tau5,tau6,tau7)=NaN;eqflag=0;end
                            end
                        end
                    end
                end
            end
        end
    end
    if single==7, varargout{1}=h7;else varargout{6}=h7; clear h7;end
end

if (ord>=8) && ((single==0)||(single==8))
    R=size(h8,1);
    eqflag=0;
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                for tau4=1:R
                    for tau5=1:R
                        for tau6=1:R
                            for tau7=1:R
                                for tau8=1:R
                                    ind=[tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8];
                                    for i=1:8,if length(find(ind(i)==ind))>1;eqflag=1;break,end,end
                                    if eqflag==1,h8(tau1,tau2,tau3,tau4,tau5,tau6,tau7,tau8)=NaN;eqflag=0;end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if single==8, varargout{1}=h8;else varargout{7}=h8; clear h8;end
end

if (ord>=9) && ((single==0)||(single==9))
    R=size(h9,1);
    eqflag=0;
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                for tau4=1:R
                    for tau5=1:R
                        for tau6=1:R
                            for tau7=1:R
                                for tau8=1:R
                                    for tau9=1:R
                                        ind=[tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9];
                                        for i=1:9,if length(find(ind(i)==ind))>1;eqflag=1;break,end,end
                                        if eqflag==1,h9(tau1,tau2,tau3,tau4,tau5,tau6,tau7,tau8,tau9)=NaN;eqflag=0;end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if single==9, varargout{1}=h9;else varargout{8}=h9; clear h9;end
end

if (ord>=10) && ((single==0)||(single==10))
    R=size(h10,1);
    eqflag=0;
    for tau1=1:R
        for tau2=1:R
            for tau3=1:R
                for tau4=1:R
                    for tau5=1:R
                        for tau6=1:R
                            for tau7=1:R
                                for tau8=1:R
                                    for tau9=1:R
                                        for tau10=1:R
                                            ind=[tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9 tau10];
                                            for i=1:10,if length(find(ind(i)==ind))>1;eqflag=1;break,end,end
                                            if eqflag==1,h10(tau1,tau2,tau3,tau4,tau5,tau6,tau7,tau8,tau9,tau10)=NaN;eqflag=0;end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if single==10, varargout{1}=h10;else varargout{9}=h10; clear h10;end
end
