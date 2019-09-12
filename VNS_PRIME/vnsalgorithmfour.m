function [ best_cdes, objvalue ] = vnsalgorithmfour( parents, objvalue, n, d, m, sset, nset, possibJhc, vecbigm )
% VNS Algorithm for improving the J4-characteristics of a concatenated design.
% Version for parent designs with a number of basic factors that is a prime
%
% INPUTS:
% m         Number of factors to keep fixed.
% n         Run size of the parent designs. 
% d         Number of parent designs.
% parents   The n-by-m-by-d array with the two-level parent designs with 
%           the fixed factor columns.
% objvalue  Starting objective value. 
% nset      Number of 4-factor sets with J4-characteristics equal to N. 
% sset      The nset-by-4 matrix of 4-factor sets.
% possibJhc Possible values of the J4-characteristics.
% vecbigm   The 1-by-nelF4 vector with the weight for each component of the
%           F4 vector, where nelF4 is the number of elements in the 
%           F4 vector of a strength-3 design.
%
% OUTPUTS:
% best_cdes The (n*d)-by-m improved design. 
% objval    The best objective value found.
%==========================================================================

ii = 1; 
while ii <= (d-1) 
    ii = ii + 1;
    factors = randperm(m); % Random Order of the elements in the neighborhood.
    for f = factors    
        parents(:,f,ii) = -1*parents(:,f,ii); % Fold-over factor f.--------
        cdes = concatenate(parents, n, m, d); % Construct concatenated design.
        % Evaluate design.-------------------------------------------------
        resvalue = objective(cdes, sset, nset, possibJhc, vecbigm);
        
        if resvalue < objvalue 
               objvalue = resvalue; 
               best_cdes = cdes; 
               ii = 1; % Go back to N1.----------------
               break 
        end
        parents(:,f,ii) = -1*parents(:,f,ii); % Re-store parents.----------
    end 
end 
end

