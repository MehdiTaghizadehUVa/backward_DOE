function [y] = wingweight(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WING WEIGHT FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT AND INPUT:
%
% y  = wing weight
% xx = [Sw, Wfw, A, LamCaps, q, lam, tc, Nz, Wdg, Wp]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sw      = xx(1) * 50 + 150;
Wfw     = xx(2) * 80 + 220;
A       = xx(3) * 4 + 6;
LamCaps = (xx(4) * 20 - 10) * (pi/180);
q       = xx(5) * 29 + 16;
lam     = xx(6) * 0.5 + 0.5;
tc      = xx(7) * 0.1 + 0.08;
Nz      = xx(8) * 3.5 + 2.5;
Wdg     = xx(9) * 800 + 1700;
Wp      = xx(10) * 0.055 + 0.025;

fact1 = 0.036 * Sw^0.758 * Wfw^0.0035;
fact2 = (A / ((cos(LamCaps))^2))^0.6;
fact3 = q^0.006 * lam^0.04;
fact4 = (100*tc / cos(LamCaps))^(-0.3);
fact5 = (Nz*Wdg)^0.49;

term1 = Sw * Wp;

y = fact1*fact2*fact3*fact4*fact5 + term1;

end