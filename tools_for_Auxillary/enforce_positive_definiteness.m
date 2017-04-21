function Mout = enforce_positive_definiteness(M)
% Force the square matrix M to be positive definite.
%
% Author: Jason E. Hill, Ph.D.
% Date: 31 AUG 2016

nR = size(M,1);
if nR ~= size(M,2)
   error('Input matrix must be square!');
end

[~,p] = chol(M);
firstFlag = true;
while p > 0 
    % need to force M to be positive semidefinite
    [V,D] = eig(M);
    evals = diag(D);
    negEvals = (evals<0);
    [~,firstPosEval] = max(1./diag(D));
    for evali = 1:nR
        if negEvals(evali)
            D(evali,evali) = D(firstPosEval,firstPosEval);
            V(:,evali) = V(:,evali);
        end
    end
    M = V*D*V';
    if firstFlag
        M(M> 1) = 1.0;
        M(M<-1) = -1.0;
    end
    [~,p] = chol(M);
    firstFlag = false;
end 

Mout = M;

end

