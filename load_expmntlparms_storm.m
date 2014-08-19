%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Experimental Parameter File   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [am0,am_slope,ampref_max,ampref_min,maxPflood,highrisk,stormfreq,maxdam,...
    Cmit,miteff,AVGFARMRETURN,STDFARMRETURN,coastvalue,midvalue,...
    inlandvalue,milecost,milestraveled]=load_expmntlparms_storm(EXPTRUNS)
% Coastal Amenity
am0=linspace(200000,800000,EXPTRUNS);
% am0=500000*ones(1,EXPTRUNS);
am_slope=linspace(0.025,0.125,EXPTRUNS);
% am_slope=0.1*ones(1,EXPTRUNS);

% Consumer preferences
ampref_max=0.9*ones(1,EXPTRUNS);
ampref_min=0.1*ones(1,EXPTRUNS);

% Storm Impacts and Mitigation
maxPflood=0.7*ones(1,EXPTRUNS);
highrisk=30*ones(1,EXPTRUNS);
maxdam=ones(1,EXPTRUNS);
stormfreq=ones(1,EXPTRUNS);
% stormthresh=15*ones(1,EXPTRUNS);
% Cdam=0.5*ones(1,EXPTRUNS);
Cmit=3000*ones(1,EXPTRUNS);
miteff=1*ones(1,EXPTRUNS);
    
% Initial land value
AVGFARMRETURN=2486.3*ones(1,EXPTRUNS);
STDFARMRETURN=0.10*AVGFARMRETURN.*ones(1,EXPTRUNS);
lvset=mat2cell([3; 2; 1]*linspace(0.5,1.5,EXPTRUNS),ones(3,1),EXPTRUNS);
% lvset=[3; 2; 1]*ones(1,EXPTRUNS);
% coastvalue=3*ones(1,EXPTRUNS);
% midvalue=2*ones(1,EXPTRUNS);
% inlandvalue=1*ones(1,EXPTRUNS);
coastvalue=lvset{1};
midvalue=lvset{2};
inlandvalue=ones(1,EXPTRUNS);

% Travel costs
milecost=1.30*ones(1,EXPTRUNS);
milestraveled=500*ones(1,EXPTRUNS);

end