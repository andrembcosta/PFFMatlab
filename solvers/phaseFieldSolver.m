classdef phaseFieldSolver < handle
   
    properties
       FEProblem;
       damageSolver;
       mechanicsSolver;
       degradation_function;
       decomposition;
       AM_tolerance;
       AM_max_iter;
    end
    
    methods
        function obj = phaseFieldSolver(FEProblem)
           %create damageSolver
           
           %create mechanicsSolver
           
           
        end
        
        function alternateMinimization(obj)
            
            tol = 1e15;
            num_iterations = 0;
            while tol > obj.AM_tolerance
                obj.mechanicsSolver.solve();
                obj.damageSolver.solve();
                num_iterations = num_iterations+1;
                tol = L2norm();
            end
            
        
        end
        
        
    end
    
end