clear all;clc;
%this is an example input file
mesh = "my_mesh_file.msh";
%get variables, if ends with msh, call gmsh parser, if "internal_generator"
%check for mesh_type cracked or uncracked, mesh_size and mesh_subdivisions

%physics
physics = ["mechanics_only","phase-field","poromechanics","porous_flow","fracture_flow"]; 
%this doesnt require many options, can build defauls constructors of
%solvers

%solvers
pf_degradation_function = "lorentz"; %or "quadratic"
pf_strain_decomposition = "spectral";
damage_local_dissipation = "linear";
damage_irreversibility = "none";
damage_nonlinear_tol = 1e-9;
mechanics_nonlinear_tol = 1e-10;

%materials
materials = ["shale", "water"];

%boundary_conditions
bc_surface_names = ["left", "right"];
bc_node_nums = [[1,2,3,4],[56,59,60,75]];
bc_type = ["Neumann","Dirichlet"];
bc_variable = ["u","d","p"];
bc_value = [[1,1]; 2, 3]; %only constants for now

%initial_conditions
%initial_crack = [[x1 y1],[x2 y2],[x3 y3]]; %several initial cracks

%output file
output = "output/problemName/problemName_out";

%run program
%MatFrac;