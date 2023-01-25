%MatFrac is a matlab reservoir simulator 
%authors: Andre Costa, Matteo Cusini, Rudy Geelen
fprintf('Thanks for using MatFrac!');
%build big struct with all the data defined on input file
inputFileCell=whos;
inputFileStruct=cell2struct({inputFileCell.name}.',{inputFileCell.name});
inputFileStruct=eval(structvars(inputFileStruct,0).');
%allocate FEProblem
FEProblem(); 
%initialize FEProblem with input data
FEProblem.init(inputFileStruct);
fprintf('Simulation ended.');
