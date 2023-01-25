%test damage solver

%create damage solver object
addpath('../mesh/')
addpath('../functions/')
addpath('./DamageKernels/')
mesh = generateNbyNMesh(1,1,0);
%mesh.ien = [1 2 1 3 4];
mat = materials;
irrev = 'none';
local_dissipation = quadratic;
degradation_function = quadratic_degradation;
newton_tolerance = 1e-6;
%active_energy = [1 2 3 4]; %quadrature point quantity
active_energy = [1 1 1 1]; %quadrature point quantity
dSolver = stdDamageSolver(mesh, mat, irrev, local_dissipation, degradation_function, newton_tolerance, active_energy);
dSolver.nonlinearSolve();
disp("damage: ")
disp(dSolver.damage)
disp("validation: ")
disp(active_energy(1)/(1+mat.Gc/(local_dissipation.c0*mat.l)))