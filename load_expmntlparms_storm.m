%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Experimental Parameter File   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [am0,am_slope,ampref_max,ampref_min,maxPflood,highrisk,stormfreq,maxdam,...
    Cmit,miteff,AVGFARMRETURN,STDFARMRETURN,coastvalue,midvalue,...
    inlandvalue,milecost,milestraveled]=load_expmntlparms_storm(EXPTRUNS)
% Coastal Amenity
% am0=linspace(200000,800000,EXPTRUNS);
am0=500000*ones(1,EXPTRUNS);        %baseline
% am_slope_parms=[0.05 0.075 0.1 0.125 0.15];
% am_slope=repmat(am_slope_parms,1,1);
% am_slope=repmat(reshape(repmat(am_slope_parms,MRUNS,1),MRUNS*...
%     length(am_slope_parms),1),4,1);
am_slope=0.1*ones(1,EXPTRUNS);      %baseline

% Consumer preferences
ampref_max=0.9*ones(1,EXPTRUNS);    %baseline
ampref_min=0.1*ones(1,EXPTRUNS);    %baseline

% Storm Impacts and Mitigation
maxPflood=0.7*ones(1,EXPTRUNS);     %baseline
highrisk=30*ones(1,EXPTRUNS);       %baseline
maxdam=ones(1,EXPTRUNS);            %baseline
stormfreq_parms=ones(1,EXPTRUNS);         %baseline
% stormfreq_parms=[1 2 3 4];
% stormfreq=reshape(repmat(stormfreq_parms,length(am_slope_parms),1),...
%     length(am_slope_parms)*length(stormfreq_parms),1)';         
stormfreq=stormfreq_parms;

% stormthresh=15*ones(1,EXPTRUNS);
% Cdam=0.5*ones(1,EXPTRUNS);
Cmit=3000*ones(1,EXPTRUNS);         %baseline
miteff=1*ones(1,EXPTRUNS);          %baseline
    
% Initial land value
AVGFARMRETURN=2486.3*ones(1,EXPTRUNS);
STDFARMRETURN=0.10*AVGFARMRETURN.*ones(1,EXPTRUNS);
coastvalue=2*ones(1,EXPTRUNS);      %baseline
midvalue=1.5*ones(1,EXPTRUNS);        %baseline
inlandvalue=1*ones(1,EXPTRUNS);     %baseline
% coastvalue=[1 2 3 4 5];
% coastvalue=reshape(repmat(coastvalue,length(am_slope_parms),1),...
%     length(am_slope_parms)*length(coastvalue),1)';
% midvalue=[1 1.5 2 2.5 3];
% midvalue=reshape(repmat(midvalue,length(am_slope_parms),1),...
%     length(am_slope_parms)*length(midvalue),1)';
% inlandvalue=[1 1 1 1 1];
% inlandvalue=reshape(repmat(inlandvalue,length(am_slope_parms),1),...
%     length(am_slope_parms)*length(inlandvalue),1)';


% Travel costs
milecost=1.30*ones(1,EXPTRUNS);     %baseline
milestraveled=500*ones(1,EXPTRUNS); %baseline

end

