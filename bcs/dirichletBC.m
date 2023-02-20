%this class implements a Dirichlet boundary condition
%it contains the nodeSet (set of nodes) where the conditions must be
%applied, as well as the values that a variable must assume at these nodes.
classdef dirichletBC < handle

    properties
        mesh
        nodeSet %true or false style
        valuesSet %include all nodes, if no condition, set to 0
    end
    
    methods
        %this function receives a residual and a jacobian after the
        %elements contributions are assemble. It then modifies it to
        %account for the boundary condition
        function [Res, Jac] = apply(Res, Jac)
            
            % number of degrees of freedom
            ndof = msh.numnod;
            
            % use this sparse multiplication procedure to use matlab's
            % built-in sparse matrix operations
            % loop over nodes to set essential boundary conditions
            BC_vector_plus = ones(1,ndof);
            BC_vector_minus = zeros(1,ndof);
            for n=1:ndof
                if (obj.nodeSet(n) == 1)
                  BC_vector_plus(n) = 0.0;
                  BC_vector_minus(n) = 1.0;
                  Res(n) = obj.valuesSet(n);
                end 
            end
            BC_editor_mult = sparse(1:ndof,1:ndof,BC_vector_plus);
            BC_editor_sum = sparse(1:ndof,1:ndof,BC_vector_minus);
            Jac = BC_editor_mult*Jac + BC_editor_sum;
            
        
        end
        
    end


end