function sys_square = square_system_inputs_pp(sys_orig,ZERO)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Purpose:
%   Square-up a non-square linear MIMO system where the 
%   system is tall, i.e. number of outputs is more than
%   the number of inputs.
% 
%   Assumes D = 0 of appropriate dimensions
%
%   Reference:
%   Misra, P., "A Computational Algorithm for Squaring-Up
%   Part I: Zero Input-Output Matrix", Proceedings of 
%   1992 CDC
%
% Input:
%   sys_orig - original unsquared system with sizes:
%       A       - nxn system state matrix
%       B       - nxm system input matrix
%       C       - pxn system output matrix (p>m)
%       D       - assumed to be zero of dim(p,m)
%       ZERO    - place the transmission zero here
%
% Output:
%   sys_square - squared-up system object
%
% Created:
%   05-31-12, Ross Gadient
%
% Modified:
% K. A. Wise Modified to use poleplacement instead of LQR
% to place transmission zero at ZERO
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eta_c = 0.;

% System matrices
A = sys_orig.a;
B = sys_orig.b;
C = sys_orig.c;

% Extract nominal system sizes
np = size(A,1);
mp = size(B,2);
pp = size(C,1);

if isequal(mp,pp)
    disp('System is already square!');
    sys_square = sys_orig;
    return
end
if mp>pp
    error('More inputs than outputs, need to square OUTPUTS');
end

% Enforce controllability
CC = ctrb(A,B);
svd_CC = svd(CC);
very_small_svd_CC = find(svd_CC <= 1e-06);
if ~isempty(very_small_svd_CC)
    disp('System NOT Controllable');
end;
% Enfore full column rank of B
if rank(B) < mp
    disp('B not full column rank');
end
% Enforce rank of CB
if rank(C*B) < mp
    disp('CB rank deficient');
end

% Do SVD decomposition to transform system
[~,~,V] = svd(C);

% Orthogonal state coorinate transformation
T = V';
Abar = T*A*inv(T);
Bbar = T*B;
Cbar = C*inv(T);
% Verify C matrix partition
very_small_C = find(Cbar(pp+1:end,:) >= 1e-06);
if ~isempty(very_small_C)
    disp('C matrix partition may be incorrect');
end

% Transformed system matrices partitioned
A11 = Abar(1:pp,1:pp);
A12 = Abar(1:pp,pp+1:end);
A21 = Abar(pp+1:end,1:pp);
A22 = Abar(pp+1:end,pp+1:end);
B11 = Bbar(1:pp,:);
B21 = Bbar(pp+1:end,:);
C1  = Cbar(:,1:pp);

% Form B12 so that it lies in nullspace of B11.  This ensures full rank of
% matrix B1, which must be inverted later
rank_B11 = rank(B11);
[UU,~,~] = svd(B11);
B12 = UU(end-(size(B11,1)-rank_B11-1):end,:)';
B1 = [B11,B12];
% Verify B1 is full rank
if rank(B1) < pp
    disp('B1 not full rank');
end

% Form matrices need to do pole placement for B22
B2_wiggle = [B21,zeros(np-pp,pp-mp)];
A22_wiggle = A22-B2_wiggle*inv(B1)*A12;
pseudo_B = A12'*(inv(B1)');

% Use LQR to place poles
% Mod
A22_wiggle_mod = A22_wiggle + eta_c*eye(np-pp);
%B22 = lqr(A22_wiggle_mod',pseudo_B(:,end-(pp-mp-1):end),Q,R)';
B22 = place(A22_wiggle_mod',pseudo_B(:,end-(pp-mp-1):end),ZERO)';

% Formed transformed squared-up input matrix
squared_B = [B11,B12;B21,B22];

% Change back to original system coordinates B matrix
squared_B_original = inv(T)*squared_B;

% Normalize squared-up columns to prevent high gain
for ii=mp+1:pp
    squared_B_original(:,ii) = squared_B_original(:,ii)./norm(squared_B_original(:,ii));
end

% Return square-up system
sys_square = ss(A,squared_B_original,C,0*C*squared_B_original);

return
