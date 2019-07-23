function [x, alpha, q] = AlphaRQSweep(H, S, F, iV, x, q)

%[x,alpha,q] = AlphaRQSweep(H,S,F,iV,x,q) is the one-sweep 
%alpha-eigenvalue Rayleigh Quotient Fixed Point. The function 
%performs one sweep of the system.

%Input: H - Leakage/Transport Matrix
%       S - Scattering Matrix
%       F - Fission Matrix
%       iV - Inverse Velocity Matrix
%       x - Previous Angular Flux Vector Iterate
%       q - Previous Scattering and Fission Source Vector

%Output: x - New Angular Flux Iterate
%        alpha - New Alpha-Eigenvalue Iterate
%        q - New Scattering and Fission Source Vector

    %Check if source is nonzero
    if ( sum(q) == 0 )
        alpha = 0;
    else
        %Alpha-eigenvalue RQ update
        alpha = (x'*(S+F)*x - x'*q)/(x'*iV*x);
    end
    
    %Set new source
    q = (-alpha*iV + S + F)*x;
    
    %Transport Sweep
    x = H\q;
    
return
