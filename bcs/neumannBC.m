%this class implements a Neumann boundary condition
%it contains the sideSet (set of edges) where the conditions must be
%applied, as well as the values that a variable must assume at these edges.
classdef neumannBC < handle

    properties
        mesh
        sideSet %true or false style
        valuesSet %include all nodes, if no condition, set to 0
    end
    
    methods
        %this function receives a residual and a jacobian after the
        %elements contributions are assemble. It then modifies it to
        %account for the boundary condition
        function [Res, Jac] = apply(Res, Jac)
            
            error("Neumann BC not yet implemented");
        
        end
        
    end


end