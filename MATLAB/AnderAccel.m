function [alpha,x,residual] = AnderAccel(H,S,F,iV,maxits,tol,mmax,fpiters,beta)

%[alpha,x,residual] = AnderAccel(H,S,F,iV,maxits,tol,mmax,fpiters,beta) is
%the Anderson acceleration implementation for the alpha-eigenvalue Rayleigh
%Quotient Fixed Point method.

%Input: H - Leakage/Transport Matrix
%       S - Scattering Matrix
%       F - Fission Matrix
%       iV - Inverse Velocity Matrix
%       maxits - Maximum Number of Iterations
%       tol - L2 Norm Residual Tolerance
%       mmax - Maximum Number of Residual Vectors used by Anderson
%              acceleration
%       fpiters - Initial fixed-point iteration function evaluations
%       beta - Relaxation parameter

%Output: alpha - Converged Alpha-Eigenvalue
%        x - Converged Angular Flux Vector
%        residual - Residual vector

    %Initialize initial angular flux guess and source
    x = ones(length(H),1);
    q = zeros(length(H),1);

    fold = 0;

    maa = 0;
    G = [];

    %Initial fixed point iterations
    for i = 1:fpiters
    [x,alpha,q] = AlphaRQSweep(H,S,F,iV,x,q);
    end

    %Start Anderson Acceleration
    for k = 0:maxits
    
        %Fixed Point Iteration Evaluation
        xold = x;
        [gcur, alpha, q] = AlphaFcn(H,S,F,iV,x,q);
        fcur = gcur - x;
    
        %Form residual and function matrices
        if ( k > 0 )
            dF = fcur - fold;
            dG = gcur - gold;
            if ( maa < mmax )
                G = [G,dG];
            else
                G = [G(:,2:maa),dG];
            end
            maa = maa + 1;
        end
        fold = fcur;
        gold = gcur;
    
        %Gram-Schmidt Orthogonalization
        if ( maa == 0 )
            x = gcur;
        else
            if ( maa == 1 )
                Q(:,1) = dF/norm(dF); R(1,1) = norm(dF);
            else
                if ( maa > mmax )
                    [Q,R] = qrdelete(Q,R,1); %Delete matrix column
                    maa = maa - 1;
                end
                for i = 1:maa - 1
                    R(i,maa) = Q(:,i)'*dF;
                    dF = dF - R(i,maa)*Q(:,i); %QR decomposition
                end
                Q(:,maa) = dF/norm(dF); R(maa,maa) = norm(dF);
            end
            gamma = R\(Q'*fcur); %Solve for gamma coefficients
            x = gcur - G*gamma; %Set next iterate
            x = x - (1-beta)*(fcur - Q*R*gamma); %Relaxation Coefficient
        end
    
        residual(k+1) = norm(x-xold)/norm(x);
    
        fprintf('Iter: %i alpha = %e residual = %e \n',...
            k+1, alpha, residual(k+1));
    
        %Termination Criterion
        if ( residual(k+1) < tol )
            break
        end
    
        %Break if NaN
        if (isnan(alpha) == 1 )
            break
        end
    
    end
        
return