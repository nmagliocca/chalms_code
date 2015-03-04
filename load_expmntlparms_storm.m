%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Experimental Parameter File   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [am0,am_slope,ampref_max,ampref_min,maxPflood,highrisk,stormfreq,maxdam,...
    Cmit,miteff,AVGFARMRETURN,STDFARMRETURN,coastvalue,midvalue,...
    inlandvalue,milecost,milestraveled]=load_expmntlparms_storm(EXPTRUNS)
% Coastal Amenity
% am0=linspace(200000,800000,EXPTRUNS);
am0=500000*ones(1,EXPTRUNS);        %baseline
am_slope_parms=[0.025 0.05 0.075 0.1 0.125 0.15 0.175];
am_slope=repmat(am_slope_parms,1,4);
% % am_slope=repmat(reshape(repmat(am_slope_parms,MRUNS,1),MRUNS*...
% %     length(am_slope_parms),1),4,1);
% am_slope=0.1*ones(1,EXPTRUNS);      %baseline

% Consumer preferences
ampref_max=0.9*ones(1,EXPTRUNS);    %baseline
ampref_min=0.1*ones(1,EXPTRUNS);    %baseline

% Storm Impacts and Mitigation
maxPflood=0.7*ones(1,EXPTRUNS);     %baseline
highrisk=30*ones(1,EXPTRUNS);       %baseline
maxdam=ones(1,EXPTRUNS);            %baseline
% stormfreq=ones(1,EXPTRUNS);         %baseline
stormfreq_parms=[1 2 3 4];
stormfreq=reshape(repmat(stormfreq_parms,length(am_slope_parms),1),...
    length(am_slope_parms)*length(stormfreq_parms),1)'; 
% % stormfreq=reshape(repmat(stormfreq_parms,MRUNS*length(am_slope_parms),1),...
% %     MRUNS*length(am_slope_parms)*length(stormfreq_parms),1);         

% stormthresh=15*ones(1,EXPTRUNS);
% Cdam=0.5*ones(1,EXPTRUNS);
Cmit=3000*ones(1,EXPTRUNS);         %baseline
miteff=1*ones(1,EXPTRUNS);          %baseline
    
% Initial land value
AVGFARMRETURN=2486.3*ones(1,EXPTRUNS);
STDFARMRETURN=0.10*AVGFARMRETURN.*ones(1,EXPTRUNS);
% lvset=mat2cell([3; 2; 1]*linspace(0.5,1.5,EXPTRUNS),ones(3,1),EXPTRUNS);    
% coastvalue=lvset{1};
% midvalue=lvset{2};
% inlandvalue=ones(1,EXPTRUNS);
coastvalue=3*ones(1,EXPTRUNS);
midvalue=2*ones(1,EXPTRUNS);
inlandvalue=1*ones(1,EXPTRUNS);


% Travel costs
milecost=1.30*ones(1,EXPTRUNS);     %baseline
milestraveled=500*ones(1,EXPTRUNS); %baseline

end

