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

% function [ks]=symmetrize(k)
% k must be a n-th dimensional array whith N<=10.
% ks is the array with N-1 NaNs points filled symmetrically where at
% least a corresponding fundamental point in k can be found among the
% permutations of their coordianates

function [ks]=symmetrize(k)


S=size(k);
ord=length(S);
M=max(S);




switch ord
case 2
    ks=tril(k)+tril(k,-1)';


case 3
    ks=NaNmat(M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                    P=perms([x1, x2, x3 ]);
                    fill=0;
                    for i=1:6
                        val=k(P(i,1),P(i,2),P(i,3));
                        if ~isnan(val), fill=1;  break, end
                    end
                    if fill
                        for i=1:6
                            ks(P(i,1),P(i,2),P(i,3))=val;
                        end
                    end
            end
        end
    end
    
case 4
    ks=NaNmat(M,M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                for x4=x3:M
                    P=perms([x1, x2, x3, x4]);
                    fill=0;
                    for i=1:24
                    val=k(P(i,1),P(i,2),P(i,3),P(i,4));
                    if ~isnan(val), fill=1;  break, end
                    end
                    if fill
                        for i=1:24
                            ks(P(i,1),P(i,2),P(i,3),P(i,4))=val;
                        end
                    end
                end
            end
        end
    end
    
 case 5
    ks=NaNmat(M,M,M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                for x4=x3:M
                    for x5=x4:M
                    P=perms([x1, x2, x3, x4, x5]);
                    fill=0;
                    for i=1:120
                        val=k(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5));
                        if ~isnan(val), fill=1;  break, end
                    end
                        if fill
                            for i=1:120
                                ks(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5))=val;
                            end
                        end
                    end
                end
            end
        end
    end
case 6
    ks=NaNmat(M,M,M,M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                for x4=x3:M
                    for x5=x4:M
                        for x6=x5:M
                            P=perms([x1, x2, x3, x4, x5, x6]);
                            fill=0;
                            for i=1:720
                                val=k(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6));
                                if ~isnan(val), fill=1;  break, end
                            end
                            if fill
                                for i=1:720
                                    ks(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6))=val;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
case 7
    ks=NaNmat(M,M,M,M,M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                for x4=x3:M
                    for x5=x4:M
                        for x6=x5:M
                            for x7=x6:M
                                P=perms([x1, x2, x3, x4, x5, x6, x7]);
                                fill=0;
                                for i=1:5040
                                    val=k(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7));
                                    if ~isnan(val), fill=1;  break, end
                                end
                                if fill
                                    for i=1:5040
                                        ks(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7))=val;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
case 8
    ks=NaNmat(M,M,M,M,M,M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                for x4=x3:M
                    for x5=x4:M
                        for x6=x5:M
                            for x7=x6:M
                                for x8=x7:M
                                    P=perms([x1, x2, x3, x4, x5, x6, x7, x8]);
                                    fill=0;
                                    for i=1:40320
                                        val=k(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7),P(i,8));
                                        if ~isnan(val), fill=1;  break, end
                                    end
                                    if fill
                                        for i=1:40320
                                            ks(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7),P(i,8))=val;
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
case 9
    ks=NaNmat(M,M,M,M,M,M,M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                for x4=x3:M
                    for x5=x4:M
                        for x6=x5:M
                            for x7=x6:M
                                for x8=x7:M
                                    for x9=x8:M
                                        P=perms([x1, x2, x3, x4, x5, x6, x7, x8, x9]);
                                        fill=0;
                                        for i=1:362880
                                            val=k(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7),P(i,8),P(i,9));
                                            if ~isnan(val), fill=1;  break, end
                                        end
                                        if fill
                                            for i=1:362880
                                                ks(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7),P(i,8),P(i,9))=val;
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
    end  
case 10
    ks=NaNmat(M,M,M,M,M,M,M,M,M,M);
    for x1=1:M
        for x2=x1:M
            for x3=x2:M
                for x4=x3:M
                    for x5=x4:M
                        for x6=x5:M
                            for x7=x6:M
                                for x8=x7:M
                                    for x9=x8:M
                                        for x10=x9:M
                                            P=perms([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]);
                                            fill=0;
                                            for i=1:3628800
                                                val=k(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7),P(i,8),P(i,9),P(i,10));
                                                if ~isnan(val), fill=1;  break, end
                                            end
                                            if fill
                                                for i=1:3628800
                                                    ks(P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6),P(i,7),P(i,8),P(i,9),P(i,10))=val;
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
        end
    end  
    
end  
