function [t_Neumann] = Neumannboundary(t, p, t_Neumann, axes, logicalcond, 'x')
% DESCRIPTION:
% This function defines the nodes where the distribuited forces will be
% applied.
for e = 1 : size(t,1)
        elnodes = t(e, :);
        
        I = p(1, elnodes) == 0 | p(1, elnodes) == 1000 ;%This should apply the force on the last x value?
        if( sum(I) == 2)
            t_Neumann = [t_Neumann; elnodes(I)];
        end
        
end
end