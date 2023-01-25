%handle classes in MATLAB are accesed by reference 
%check docs for value (standard) vs handle class in matlab
%https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html
classdef mechanicsSolver < handle 
    
    properties
        mesh
        material
        degradation_function
        fixed_dofs
        fixed_values
        nonlinear_tolerance
        displacements
        stress_calculator
        stiffness_calculator
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
        function [R, J] = assemble(obj,u)
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
            stress = obj.stress_calculator(obj.mesh, obj.material, obj.displacements, gauss);
            Kmat = obj.stiffness_calculator(obj.mesh, obj.material, obj.displacements, gauss);
            for e=1:obj.mesh.numele
                
              %The residual is separated into many contributions, that are called separetely
              [r_stress_div, j_stress_div] = getElemStressTerm(obj.mesh,gauss,stress,Kmat,e);
              r_body_force = getElemBodyForceResTerm(obj.mesh,obj.material,gauss,obj.local_dissipation,u,e);
              r_surface_force = getElemSurfaceForceResTerm(obj.mesh,gauss,obj.degradation_function,obj.active_energy,u,e);
              

              % assemble in J and R
              for i=1:4
                rbk = obj.mesh.ien(e,i+1);
                indices_R = [indices_R; rbk];
                values_R = [values_R; r_stress_div(i)+r_body_force(i)+r_surface_force(i)];  
                for j=1:4
                  cbk = obj.mesh.ien(e,j+1);     
                  row_indices_J = [row_indices_J; rbk];
                  column_indices_J = [column_indices_J ; cbk];
                  values_J = [values_J; j_stress_div(i,j)];            
                end
              end
            end % end element loop
            R = sparse(indices_R, ones(size(indices_R)), values_R, ndof, 1);
            J = sparse(row_indices_J, column_indices_J, values_J, ndof, ndof);
        end
        
        function [R, J] = applyBCs(obj, R, J)
            ndof = 2*obj.mesh.numnod;
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
            u = obj.displacements;
            u(obj.fixed_dofs) = obj.fixed_values;
            convergence_flag = 0;

            for iter = 1:max_nonlinear_iterations

                % assemble residual and jacobian
                [R, J] = obj.assemble(u);
                fprintf('Iteration %d | Residual norm: %.2e\n', iter, norm(R))
                if norm(R) < tol %absolute tolerance criteria
                    convergence_flag = iter;
                    break
                end
                % apply boundary conditions
                [R, J] = obj.applyBCs(R, J);
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