% function [yn,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10] = 
%     BasicVoltFilt(xn, kernels, ord, delay)
%
% BasicVoltFilt is an implementation of a Volterra filter of order ord 
% (ord=10 at most), given the column input vector (xn) and a structure of 
% kernels (kernels) as follows:
%
% xn is the input vector (a column).
%
% kernels can alternatively be the name of a variable which is a structure 
% with fields h0, h1, ..., hN (N>=ord), or the name of a file (as a string)
% containing the structure itself. The N-th array hN should be a N-th 
% dimensional array with the values of a Volterra kernel.
% In this function the kernels structure will be stored in a temporary mat 
% file (namely tempor.mat) in the current working directory to let the maximum
% amount of memory available for the computations.
%
% delay (default value 0) can be used to account for a possibly causal delay 
% (formerly known or assessed) featured by the system we want to model with
% the Volterra filter. 
% If such a delay happens to be, the value of delay parameter avoids the 
% explicit computation of the zero valued points of the kernels due to 
% causality delay, indeed often reducing considerably the effective 
% computations burden. 
%
% delay can be also used as a means to trim the lower (in time lag sense) 
% part of kernels.
%
% yn will be the total ord-order Volterra filter output to the xn input 
% stimulus.
%
% Also a response for each Volterra filter series component can be obtained 
% in y0, y1, ... yN being N the order of the N-th component: 
% yn = sum(y0,...,yN).
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DII, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.
% If you are using this program for a scientific work, we encourage you to cite
% the following paper (the file cite.bib, containing the reference in bibtex
% format is also provided):
%
% Simone Orcioni. Improving the approximation ability of Volterra series 
% identified with a cross-correlation method. Nonlinear Dynamics, 2014.
%
%﻿Orcioni, S., Terenzi, A., Cecchi, S., Piazza, F., & Carini, A. (2018). 
% Identification of Volterra Models of Tube Audio Devices using 
% Multiple-Variance Method. Journal of the Audio Engineering Society, 
% 66(10), 823–838. https://doi.org/10.17743/jaes.2018.0046

% Copyright (C) 2006 Massimiliano Pirani
% Copyright (C) 2018 Simone Orcioni
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

function varargout=BasicVoltFilt(xn,kernels,ord,delay)

if ischar(kernels) 
    nomefile=kernels;
else
    nomefile='tempor.mat';
    h0=kernels.h0;save('tempor.mat','h0');
    switch ord
    case 1, h1=kernels.h1;save('tempor.mat','h0','h1');
    case 2, h1=kernels.h1;h2=kernels.h2;save('tempor.mat','h0','h1','h2');
    case 3, h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;save('tempor.mat','h0','h1','h2','h3');
    case 4, h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;h4=kernels.h4;save('tempor.mat','h0','h1','h2','h3','h4');
    case 5, h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;h4=kernels.h4;h5=kernels.h5;save('tempor.mat','h0','h1','h2','h3','h4','h5');
    case 6, h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;h4=kernels.h4;h5=kernels.h5;h6=kernels.h6;save('tempor.mat','h0','h1','h2','h3','h4','h5','h6');
    case 7, h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;h4=kernels.h4;h5=kernels.h5;h6=kernels.h6;h7=kernels.h7;save('tempor.mat','h0','h1','h2','h3','h4','h5','h6','h7');
    case 8, h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;h4=kernels.h4;h5=kernels.h5;h6=kernels.h6;h7=kernels.h7;h8=kernels.h8;save('tempor.mat','h0','h1','h2','h3','h4','h5','h6','h7','h8');
    case 9, h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;h4=kernels.h4;h5=kernels.h5;h6=kernels.h6;h7=kernels.h7;h8=kernels.h8;h9=kernels.h9;save('tempor.mat','h0','h1','h2','h3','h4','h5','h6','h7','h8','h9');
    case 10,h1=kernels.h1;h2=kernels.h2;h3=kernels.h3;h4=kernels.h4;h5=kernels.h5;h6=kernels.h6;h7=kernels.h7;h8=kernels.h8;h9=kernels.h9; h10=kernels.h10;save('tempor.mat','h0','h1','h2','h3','h4','h5','h6','h7','h8','h9','h10');
    end
   clear kernels;
   
end
    
if nargin<4
    if nargin==3
        delay=0;
    else
        error('check the correctness of function arguments');
    end
end

varargout=cell(ord+2,1);

h0=load(nomefile,'h0');h0=h0.h0;
varargout{2}=h0;

A=1;

if ord>0
    tic
    h1=load(nomefile,'h1');h1=h1.h1;
    varargout{3} = filter(h1(delay+1:end),1,xn)*A;
    clear h1;
end




if ord>1
    tic
    L=length(xn);
    h2=load(nomefile,'h2');h2=h2.h2;
    R=size(h2,1);
    R_1=R-1;
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            if J2==J1
                a=h2(J1+1,J2+1);
            else
                a=2*h2(J1+1,J2+1);
            end
            temp=temp+a*p1.*[zeros(J2,1);xn(1:end-J2)];
        end
    end
    varargout{4}=temp*A^2;
    clear h2;
end


if ord>2
    tic
    h3=load(nomefile,'h3');h3=h3.h3;
    temp=zeros(L,1);
    R=size(h3,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
             p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                temp=temp+coefs(3,6,[J1 J2 J3])*h3(J1+1,J2+1,J3+1)...
                    *p1...
                    .*p2...
                    .*[zeros(J3,1);xn(1:end-J3)];
            end
        end
    end
    varargout{5}=temp*A^3;
    clear h3;
end

if ord>3
    tic
    h4=load(nomefile,'h4');h4=h4.h4;
    temp=zeros(L,1);
    R=size(h4,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    temp=temp+coefs(4,24,[J1 J2 J3 J4])*h4(J1+1,J2+1,J3+1,J4+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*[zeros(J4,1);xn(1:end-J4)];
                 end
             end
         end
     end
        varargout{6}=temp*A^4;
        clear h4;
end
    
if ord>4
    tic
    h5=load(nomefile,'h5');h5=h5.h5;
    temp=zeros(L,1);
    R=size(h5,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                    temp=temp+coefs(5,120,[J1 J2 J3 J4 J5])*h5(J1+1,J2+1,J3+1,J4+1,J5+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*p4...
                        .*[zeros(J5,1);xn(1:end-J5)];
                    end
                 end
             end
         end
     end
        varargout{7}=temp*A^5;
        clear h5;
end    
    
if ord>5
    tic
    h6=load(nomefile,'h6');h6=h6.h6;
    temp=zeros(L,1);
    R=size(h6,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        [zeros(J5,1);xn(1:end-J5)];
                         for J6=J5:R_1
                     temp=temp+coefs(6,720,[J1 J2 J3 J4 J5 J6])*h6(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*p4...
                        .*p5...
                        .*[zeros(J6,1);xn(1:end-J6)];
                          end
                    end
                 end
             end
         end
     end
        varargout{8}=temp*A^6;
        clear h6;
end    

if ord>6
    tic
    h7=load(nomefile,'h7');h7=h7.h7;
    kernels.h7=NaN;
    temp=zeros(L,1);
    R=size(h7,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                         for J6=J5:R_1
                             p6=[zeros(J6,1);xn(1:end-J6)];
                             for J7=J6:R_1
                     temp=temp+coefs(7,5040,[J1 J2 J3 J4 J5 J6 J7])*h7(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*p4...
                        .*p5...
                        .*p6...
                        .*[zeros(J7,1);xn(1:end-J7)];
                            end
                         end
                     end
                 end
             end
         end
     end
        varargout{9}=temp^7;
        clear h7;
end 

if ord>7
    tic
    h8=load(nomefile,'h8');h8=h8.h8;
    temp=zeros(L,1);
    R=size(h8,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                         for J6=J5:R_1
                             p6=[zeros(J6,1);xn(1:end-J6)];
                             for J7=J6:R_1
                                 p7=[zeros(J7,1);xn(1:end-J7)];
                                 for J8=J7:R_1
                    temp=temp+coefs(8,40320,[J1 J2 J3 J4 J5 J6 J7 J8])*h8(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1,J8+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*p4...
                        .*p5...
                        .*p6...
                        .*p7...
                        .*[zeros(J8,1);xn(1:end-J8)];
                                 end
                            end
                         end
                     end
                 end
             end
         end
     end
        varargout{10}=temp*A^8;
        clear h8;
end 

if ord>8
    tic
    h9=load(nomefile,'h9');h9=h9.h9;
    temp=zeros(L,1);
    R=size(h9,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                         for J6=J5:R_1
                             p6=[zeros(J6,1);xn(1:end-J6)];
                             for J7=J6:R_1
                                 p7=[zeros(J7,1);xn(1:end-J7)];
                                 for J8=J7:R_1
                                     p8=[zeros(J8,1);xn(1:end-J8)];
                                     for J9=J8:R_1
                    temp=temp+coefs(9,362880,[J1 J2 J3 J4 J5 J6 J7 J8 J9])*h9(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1,J8+1,J9+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*p4...
                        .*p5...
                        .*p6...
                        .*p7...
                        .*p8...
                        .*[zeros(J9,1);xn(1:end-J9)];
                                      end
                                 end
                            end
                         end
                     end
                 end
             end
         end
     end
        varargout{11}=temp*A^9;
        clear h9;
end 

if ord>9
    tic
    h10=load(nomefile,'h10');h10=h10.h10;
    temp=zeros(L,1);
    R=size(h10,1);
    R_1=R-1;
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                         for J6=J5:R_1
                             p6=[zeros(J6,1);xn(1:end-J6)];
                             for J7=J6:R_1
                                 p7=[zeros(J7,1);xn(1:end-J7)];
                                 for J8=J7:R_1
                                     p8=[zeros(J8,1);xn(1:end-J8)];
                                     for J9=J8:R_1
                                         p9=[zeros(J9,1);xn(1:end-J9)];
                                         for J10=J9:R_1
                    temp=temp+coefs(10,3628800,[J1 J2 J3 J4 J5 J6 J7 J8 J9 J10])*h10(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1,J8+1,J9+1,J10+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*p4...
                        .*p5...
                        .*p6...
                        .*p7...
                        .*p8...
                        .*p9...
                        .*[zeros(J10,1);xn(1:end-J10)];
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
        varargout{12}=temp*A^10;
        clear h10;
end 

varargout{1}=varargout{2};
for i=3:ord+2
    varargout{1}=varargout{1}+varargout{i};
end



function c=coefs(ord,ordfact,v)
i=1;
while ~isempty(v)
    L=length(v);
    v=v(find(v(1)~=v));
    gruppi(i)=factorial(L-length(v));
    i=i+1;
end
c=ordfact/prod(gruppi);

