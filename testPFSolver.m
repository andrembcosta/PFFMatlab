%script to test phase-field solver

msh = generateBlock();
Gc = 1;
ell = 1;
split = "none";
degradation = quadratic;
localDis = quadratic;
FEProblem.mesh = msh;
FEProblem.materialModel = damageSolidModel(FEProblem.getMesh(), Gc, ell, split, degradation, localDis);

obj = phaseFieldSolver(FEProblem);
obj.alternateMinimization();