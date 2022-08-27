classdef stdDamageSolver
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
        function obj = stdDamageSolver(msh, mat, decomp, irrev, local_dis, deg_fun, tol)
            mesh = msh;
            material = mat;
            irrevesibility_enforcement = irrev;
            local_dissipation = local_dis;
            degradation_function = deg_fun;
            nonlinear_tolerance = tol;
            %initialize damage
            d = zeros(mesh.numnod,1); % number of nodes
            %initialize fixed_dofs, fixed_values
            fixed_dofs = zeros(mesh.numnod,1);
            fixed_values = zeros(mesh.numnod,1);
        end
        
        %parallel assembly routine, that computes the residual and jacobian
        %at a given step
        function [R, J] = assemble()
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
            parfor e=1:obj.msh.numele
                
              %The residual is separated into many contributions, that are called separetely
              r_d = getElemLocalDissipationRes(obj.mesh,obj.material,gauss,obj.local_dissipation,e);
              r_grad = getElemGradientResTerm(obj.mesh,obj.material,gauss,obj.local_dissipation,d,e);
              r_driv = getElemDrivingForceResTerm(obj.mesh,gauss,obj.degradation_function,obj.active_energy,e,d);
              
              %The jacobian is separated into many contributions, that will 
              %be called separately
              j_d = getElemLocalDissipationJacTerm(obj.mesh, obj.material, obj.local_dissipation, e);
              j_grad = getElemGradientJacTerm(obj.mesh, obj.material, e);
              j_driv = getElemDrivingForceJacTerm(obj.mesh, obj.material, obj.degradation_function, e);
    
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
            R = sparse(indices_R, values_R, ndof, 1);
            J = sparse(row_indices_J, column_indices_J, values_J, ndof, ndof);
        end
        
        function applyBCs()
            %set fixed_dofs
            %set fixed_values
        end
        function nonlinearSolve()
            %solve system
        end
    end
    
end