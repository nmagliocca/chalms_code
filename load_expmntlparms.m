%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Experimental Parameter File   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [am0,am_slope,ampref_max,ampref_min,maxPflood,highrisk,stormfreq,...
    Cdam,Cmit,miteff,AVGFARMRETURN,STDFARMRETURN,coastvalue,midvalue,...
    inlandvalue,milecost,milestraveled]=load_expmntlparms(EXPTRUNS,NPARMS)
% Coastal Amenity
am0=reshape(repmat(linspace(200000,800000,EXPTRUNS)',1,EXPTRUNS^(NPARMS-1))',...
    1,EXPTRUNS^NPARMS);
% am0=500000*ones(1,EXPTRUNS^NPARMS);
am_slope=reshape(repmat(linspace(0.025,0.125,EXPTRUNS),EXPTRUNS,EXPTRUNS),...
    1,EXPTRUNS^NPARMS);
% am_slope=0.1*ones(1,EXPTRUNS^NPARMS);

% Consumer preferences
ampref_max=0.9*ones(1,EXPTRUNS^NPARMS);
ampref_min=0.1*ones(1,EXPTRUNS^NPARMS);

% Storm Impacts and Mitigation
maxPflood=0.7*ones(1,EXPTRUNS^NPARMS);
highrisk=30*ones(1,EXPTRUNS^NPARMS);
stormfreq=0.1*ones(1,EXPTRUNS^NPARMS);
% stormthresh=15*ones(1,EXPTRUNS^NPARMS);
Cdam=0.5*ones(1,EXPTRUNS^NPARMS);
Cmit=3000*ones(1,EXPTRUNS^NPARMS);
miteff=1*ones(1,EXPTRUNS^NPARMS);
    
% Initial land value
AVGFARMRETURN=2486.3*ones(1,EXPTRUNS^NPARMS);
STDFARMRETURN=0.10*AVGFARMRETURN.*ones(1,EXPTRUNS^NPARMS);
lvset=mat2cell(repmat([3; 2; 1]*linspace(0.5,1.5,EXPTRUNS),1,EXPTRUNS^(NPARMS-1)),...
    ones(3,1),EXPTRUNS^NPARMS);
% lvset=[3; 2; 1]*ones(1,EXPTRUNS^NPARMS);
% coastvalue=3*ones(1,EXPTRUNS^NPARMS);
% midvalue=2*ones(1,EXPTRUNS^NPARMS);
% inlandvalue=1*ones(1,EXPTRUNS^NPARMS);
coastvalue=lvset{1};
midvalue=lvset{2};
inlandvalue=ones(1,EXPTRUNS^NPARMS);

% Travel costs
milecost=1.30*ones(1,EXPTRUNS^NPARMS);
milestraveled=500*ones(1,EXPTRUNS^NPARMS);

end