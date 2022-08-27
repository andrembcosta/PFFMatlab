classdef materials
  
  properties
    lambda; % bulk modulus [N/mm2]
    mu;     % shear modulus [N/mm2]
    Gc;     % toughness [N/mm]
    l;      % length scale
    k;      % stability parameter
    alpha;  % biot coefficients
    intrinsic_permeability;
    fluid_viscosity;
    M; 
  end
  
  methods   
    function obj = materials
      %obj.mu     = 8.077e4;%4
      %obj.lambda = 1.212e5;%5
%       obj.mu     = 8.3333e9;
%       obj.lambda = 5.5556e9 - 2*obj.mu/3;
      obj.mu     = 80.77e3;
      obj.lambda = 121.15e3;
      %obj.Gc     = 27;
      obj.Gc     = 2.7;
      %obj.l      = 0.05;
      obj.l      = 0.015;
      obj.k      = 0.0;
      obj.alpha  = 0; %biot coefficient
      obj.intrinsic_permeability = 0;
      obj.fluid_viscosity = 0;
      obj.M = 0;
    end
  end
  
end