%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Experimental Parameter File   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coastal Amenity
am0=500000*ones(1,EXPTRUNS);
am_slope=linspace(0.025,0.125,EXPTRUNS);

% Consumer preferences
ampref_max=0.9*ones(1,EXPTRUNS);
ampref_min=0.1*ones(1,EXPTRUNS);

% Storm Impacts and Mitigation
maxPflood=0.7*ones(1,EXPTRUNS);
highrisk=30*ones(1,EXPTRUNS);
stormfreq=0.1*ones(1,EXPTRUNS);
% stormthresh=15*ones(1,EXPTRUNS);
Cdam=0.5*ones(1,EXPTRUNS);
Cmit=3000*ones(1,EXPTRUNS);
miteff=1*ones(1,EXPTRUNS);
    
% Initial land value
AVGFARMRETURN=2486.3*ones(1,EXPTRUNS);
STDFARMRETURN=0.10*AVGFARMRETURN.*ones(1,EXPTRUNS);
coastvalue=3*ones(1,EXPTRUNS);
midvalue=2*ones(1,EXPTRUNS);
inlandvalue=1*ones(1,EXPTRUNS);

% Travel costs
milecost=1.30*ones(1,EXPTRUNS);
milestraveled=500*ones(1,EXPTRUNS);