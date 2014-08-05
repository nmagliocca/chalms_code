function [Sdist2cbd]=load_DIST2CBD_east(distfname)
myVars={'dist2cbd'};
Sdist2cbd=load(distfname,myVars{:});
end