%handle classes in MATLAB are accesed by reference 
%check docs for value (standard) vs handle class in matlab
%https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html
classdef damageSolver < handle 
    
    properties
        mesh
        constitutiveModel
        boundaryConditions
        irreversibilityEnforcement
        nonlinearTolerance
        damage
    end
    
    methods
        %constructor of damage solver
        function obj = damageSolver(msh, constitutive, irrev, tol, d)
            obj.mesh = msh;
            obj.constitutiveModel = constitutive;
            obj.irreversibilityEnforcement = irrev;
            obj.nonlinearTolerance = tol;
            %initialize damage
            obj.damage = d; % number of nodes
        end
        
        %parallel assembly routine, that computes the residual and jacobian
        %at a given step
        function [R, J] = assemble(obj,d)
            %set quadrature rule
            gauss = [-3^(-0.5), 3^(-0.5)];
            ndof = obj.mesh.numnod;
            % init list of position and entries for jacobian
            row_indices_J = [];
            column_indices_J = [];
            values_J = [];
            % init list of position and entries for residual
            indices_R = [];
            values_R = [];
            %begin element loop
            for e=1:obj.mesh.numele
                
              %The residual is separated into many contributions, that are called separetely
              r_d = getElemLocalDissipationResTerm(obj.mesh,obj.material,gauss,obj.local_dissipation,d,e);
              r_grad = getElemGradientResTerm(obj.mesh,obj.material,gauss,obj.local_dissipation,d,e);
              r_driv = getElemDrivingForceResTerm(obj.mesh,gauss,obj.degradation_function,obj.active_energy,d,e);
              
              %The jacobian is separated into many contributions, that will 
              %be called separately
              j_d = getElemLocalDissipationJacTerm(obj.mesh, obj.material, gauss, obj.local_dissipation, d, e);
              j_grad = getElemGradientJacTerm(obj.mesh, obj.material, gauss, obj.local_dissipation, e);
              j_driv = getElemDrivingForceJacTerm(obj.mesh, gauss, obj.degradation_function, obj.active_energy, d, e);
    
              % assemble in J and R
              for i=1:4
                rbk = obj.mesh.ien(e,i+1);
                indices_R = [indices_R; rbk];
                values_R = [values_R; r_d(i)+r_grad(i)+r_driv(i)];  
                for j=1:4
                  cbk = obj.mesh.ien(e,j+1);     
                  row_indices_J = [row_indices_J; rbk];
                  column_indices_J = [column_indices_J ; cbk];
                  values_J = [values_J; j_d(i,j) + j_grad(i,j) + j_driv(i,j)];            
                end
              end
            end % end element loop
            R = sparse(indices_R, ones(size(indices_R)), values_R, ndof, 1);
            J = sparse(row_indices_J, column_indices_J, values_J, ndof, ndof);
        end
        
        function obj = nonlinearSolve(obj)
            %solve system
            %for loop Newton method
            maxIter = 50;
            tol = obj.nonlinearTolerance;
            %set initial guess;
            d = obj.damage;
            convergence_flag = 0;

            for iter = 1:maxIter

                % assemble residual and jacobian
                [R, J] = obj.assemble(d);
                for bc = obj.boundaryConditions
                    [R, J] = bc.apply(R, J);
                end
                fprintf('Iteration %d | Residual norm: %.2e\n', iter, norm(R))
                if norm(R) < tol %absolute tolerance criteria
                    convergence_flag = iter;
                    break
                end
                % compute damage increment
                delta_d = J\R;
                d = d - delta_d;

            end
            if convergence_flag == 0
                error('Newton loop did not converge\n')
            end
            fprintf('Newton loop converged in %d iterations\n', convergence_flag)
            obj.damage = d;
            
        end
        
    end
    
end