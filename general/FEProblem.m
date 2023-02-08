%handle classes in MATLAB are accesed by reference 
%check docs for value (standard) vs handle class in matlab
%https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html
classdef FEProblem < handle
   
    properties
        mesh
        physics
        materials
        solvers
        boundary_conditions
        initial_conditions
        output_path
    end
    
    methods
        function obj = FEProblem()
        end
        
        function init(inputStruct)
            obj.mesh = setupMesh(inputStruct);
            obj.physics = inputStruct.physics;
            obj.materials = inputStruct.materials;
            obj.solvers = getSolversFromPhysics(obj, inputStruct);
            obj.boundary_conditions = setupBCs(inputStruct);
            obj.initial_conditions = setupICs(inputStruct);
            obj.output_path = inputStruct.output_path;
        end
        
        function mesh = setupMesh(inputStruct)
            name = convertStringsToChars(inputFileStruct.mesh);
            if (convertCharsToStrings(name(end-3:end))==".msh")
                mesh = gmshParser(name);
            else %mesh is not from msh file
                mesh = generateMesh(inputStruct.mesh_subdivisions, inputStruct.mesh_size, 0);
            end
        end
        
        function solvers = getSolversFromPhysics(obj, inputStruct)
            solvers=[];
            for ph = obj.physics
                if ph=="Phase-field"
                    pfSolver = PhaseFieldSolver();
                    damageSolver = stdDamageSolver();
                    solvers = [solvers, pfSolver, damageSolver];
                    %catch pf options
                    pfSolver.degradation_function = inputStruct.pf_degradation_function; 
                    pfSolver.strain_decomposition = inputStruct.pf_strain_decomposition;
                    %catch damage options
                    damageSolver.local_dissipation=inputStruct.damage_local_dissipation;
                    damageSolver.irreversibility=inputStruct.damage_irreversibility;
                end
            end
            for ph = obj.physics
                if ph=="Mechanics"
                    mechSolver = mechanicsSolver();
                end
            end
            for ph = obj.physics
                if ph=="Poromechanics"
                    mechSolver = PoroMechanicsSolver();
                end
            end
            solvers = [solvers, mechSolver];
        end
        
        function parseInput(filename)
            fprintf('Parsing input...\n\n');
            obj.mesh = gmshParser('./mesh/tension0.msh', dim , nnod);
        end
   
        
        
    end
    
end