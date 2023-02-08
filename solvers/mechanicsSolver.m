%handle classes in MATLAB are accesed by reference 
%check docs for value (standard) vs handle class in matlab
%https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html
classdef mechanicsSolver < handle 
    
    properties
        mesh
        constitutiveModel
        boundaryConditions
        nonlinearTolerance
        displacements
    end
    
    methods
        %constructor of damage solver
        function obj = mechanicsSolver(msh, mat, deg_fun, tol)
            obj.mesh = msh;
            obj.material = mat;
            obj.degradation_function = deg_fun;
            obj.nonlinear_tolerance = tol;
            %initialize damage
            obj.displacements = zeros(2*obj.mesh.numnod,1); % number of nodes
            %initialize fixed_dofs, fixed_values - maybe used dof numbers
            %instead of 0s and 1s
            obj.fixed_dofs = [];%zeros(obj.mesh.numnod,1);
            obj.fixed_values = [];%zeros(obj.mesh.numnod,1);
        end
        
        %parallel assembly routine, that computes the residual and jacobian
        %at a given step
        %TODO: use parfor
        function [R, J] = assemble(obj)
            %set quadrature rule
            gauss = [-3^(-0.5), 3^(-0.5)];
            ndof = 2*obj.mesh.numnod;
            % init list of position and entries for jacobian
            row_indices_J = [];
            column_indices_J = [];
            values_J = [];
            % init list of position and entries for residual
            indices_R = [];
            values_R = [];
            %begin element loop
            %creates stress and stiffness fields at quadrature points -
            %stress(e,i,j) - element e, quadrature point i,j
            for e=1:obj.mesh.numele
                
              %This term computes the residual and jacobian contribution of
              %(stress, grad_test) in weak form 
              [stressDivR, stressDivJ] = getElemStressTerm(obj.mesh,gauss,e,obj.constitutiveModel,obj.displacements);
              %This computes the contibution of the body force to the
              %residual -(b, test) in weak form
              bodyForceR = getElemBodyForceResTerm(obj.mesh,gauss,e);              

              % assemble in J and R
              % this loop fill in the I,J,K notation for sparse matrix
              % assembly
              for i=1:4 %loop over elem nodes
                rbk = obj.mesh.ien(e,i+1);
                indices_R = [indices_R; rbk];
                values_R = [values_R; stressDivR(i)+bodyForceR(i)];  
                for j=1:4 %inner loop over elem nodes
                  cbk = obj.mesh.ien(e,j+1);     
                  row_indices_J = [row_indices_J; rbk];
                  column_indices_J = [column_indices_J ; cbk];
                  values_J = [values_J; stressDivJ(i,j)];            
                end
              end
            end % end element loop
            %assemble sparse residual and jacobian
            %TODO: is this sparse residual actually efficient?
            R = sparse(indices_R, ones(size(indices_R)), values_R, ndof, 1);
            J = sparse(row_indices_J, column_indices_J, values_J, ndof, ndof);
        end
               
        function obj = nonlinearSolve(obj)
            %solve system
            %for loop Newton method
            maxIter = 50; %nonlinear iterations
            tol = obj.nonlinearTolerance;
            %set initial guess - current value of displacements;
            u = obj.displacements;
            convergence_flag = 0;

            for iter = 1:maxIter

                % assemble residual and jacobian
                [R, J] = obj.assemble(u);
                % apply all boundary conditions
                for bc = obj.boundaryConditions
                    [R, J] = bc.apply(R,J,u);
                end
                % simple L2 norm of the residual
                % TODO: consider other norms
                residualNorm = norm(R);
                fprintf('Iteration %d | Residual norm: %.2e\n', iter, residualNorm)
                if resiualdNorm < tol %absolute tolerance criteria
                    convergence_flag = iter;
                    break
                end
                % compute damage increment
                delta_u = J\R;
                u = u - delta_u;

            end
            if convergence_flag == 0
                error('Newton loop did not converge\n')
            end
            fprintf('Newton loop converged in %d iterations\n', convergence_flag)
            obj.displacements = d;
            
        end
        
    end
    
end