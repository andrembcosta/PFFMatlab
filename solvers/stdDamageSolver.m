%handle classes in MATLAB are accesed by reference 
%check docs for value (standard) vs handle class in matlab
%https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html
classdef stdDamageSolver < handle 
    
    properties
        mesh
        material
        active_energy
        irreversibility_enforcement
        local_dissipation
        degradation_function
        fixed_dofs
        fixed_values
        nonlinear_tolerance
        damage
    end
    
    methods
        %constructor of damage solver
        function obj = stdDamageSolver(msh, mat, irrev, local_dis, deg_fun, tol, dr_force)
            obj.mesh = msh;
            obj.material = mat;
            obj.irreversibility_enforcement = irrev;
            obj.local_dissipation = local_dis;
            obj.degradation_function = deg_fun;
            obj.nonlinear_tolerance = tol;
            obj.active_energy = dr_force;
            %initialize damage
            obj.damage = zeros(obj.mesh.numnod,1); % number of nodes
            %initialize fixed_dofs, fixed_values - maybe used dof numbers
            %instead of 0s and 1s
            obj.fixed_dofs = [];%zeros(obj.mesh.numnod,1);
            obj.fixed_values = [];%zeros(obj.mesh.numnod,1);
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
        
        function [R, J] = applyBCs(obj, R, J)
            ndof = obj.mesh.numnod;
            bc_vector_ones = ones(1,ndof); %this will be a diagonal matrix with 1s almost everywhere, and 0s where a dof is fixed
            bc_vector_zeros = zeros(1,ndof); %this is be 0s everywhere and 1 where a dof is fixed
            for n = obj.fixed_dofs
                bc_vector_ones(n) = 0;
                bc_vector_zeros(n) = 1;
            end
%             for n=1:ndof
%                 if (ifix(n) == 1)
%             %       for i=1:ndof
%             %         F(i) = F(i) - K(i,n);
%             %       end
%                   BC_vector_plus(n) = 0.0;
%                   BC_vector_minus(n) = 1.0;
%                   F(n) = 1.0;
%                 end 
%             end
            bc_remove_fixed = sparse(1:ndof,1:ndof,bc_vector_ones);
            bc_add_fixed = sparse(1:ndof,1:ndof,bc_vector_zeros);
            J = bc_remove_fixed*J + bc_add_fixed;
            R = bc_remove_fixed*R;
        end
        
        function obj = nonlinearSolve(obj)
            %solve system
            %for loop Newton method
            max_nonlinear_iterations = 50;
            tol = obj.nonlinear_tolerance;
            %set initial guess;
            d = obj.damage;
            d(obj.fixed_dofs) = obj.fixed_values;
            convergence_flag = 0;

            for iter = 1:max_nonlinear_iterations

                % assemble residual and jacobian
                [R, J] = obj.assemble(d);
                fprintf('Iteration %d | Residual norm: %.2e\n', iter, norm(R))
                if norm(R) < tol %absolute tolerance criteria
                    convergence_flag = iter;
                    break
                end
                % apply boundary conditions
                [R, J] = obj.applyBCs(R, J);
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