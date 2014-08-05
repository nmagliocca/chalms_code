function Sfarmmap=load_farmmap(farmmapfname)
myVars={'AGLAYER'};
Sfarmmap=load(farmmapfname,myVars{:});
end