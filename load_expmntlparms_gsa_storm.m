%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Experimental Parameter File   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [am0,am_slope,ampref_max,ampref_min,maxPflood,highrisk,stormfreq,maxdam,...
    Cmit,miteff,AVGFARMRETURN,STDFARMRETURN,coastvalue,midvalue,...
    inlandvalue,milecost,milestraveled,dscnt]=load_expmntlparms_gsa_storm(erun,MRUNS)

 cd X:\model_results\CHALMS_coast_gsa_baseline
 load gsaparms.mat
% Coastal Amenity
% am0=linspace(200000,800000,EXPTRUNS);
% am0=500000*ones(1,EXPTRUNS);        %baseline
am0=X{erun,3};
% am_slope_parms=[0.025 0.05 0.075 0.1 0.125];
% am_slope=repmat(am_slope_parms,1,1);
% am_slope=repmat(reshape(repmat(am_slope_parms,MRUNS,1),MRUNS*...
%     length(am_slope_parms),1),4,1);
% am_slope=0.025*ones(1,EXPTRUNS);      %baseline
am_slope=X{erun,4};
% Consumer preferences
% ampref_max=0.9*ones(1,EXPTRUNS);    %baseline
ampref_max=X{erun,5};
ampref_min=0.1*ones(1,MRUNS);    %baseline

% Storm Impacts and Mitigation
maxPflood=0.7*ones(1,MRUNS);     %baseline
highrisk=30*ones(1,MRUNS);       %baseline
% maxdam=ones(1,MRUNS);            %baseline
maxdam=X{erun,2};            %baseline
% stormfreq=ones(1,EXPTRUNS);         %baseline
stormfreq=X{erun,1};         %baseline
% stormfreq_parms=[1 2 3 4];
% stormfreq=reshape(repmat(stormfreq_parms,length(am_slope_parms),1),...
%     length(am_slope_parms)*length(stormfreq_parms),1)';         

% stormthresh=15*ones(1,EXPTRUNS);
% Cdam=0.5*ones(1,EXPTRUNS);
Cmit=3000*ones(1,MRUNS);         %baseline
miteff=1*ones(1,MRUNS);          %baseline
    
% Initial land value
AVGFARMRETURN=2486.3*ones(1,MRUNS);
STDFARMRETURN=0.10*AVGFARMRETURN.*ones(1,MRUNS);
coastvalue=3*ones(1,MRUNS);      %baseline
midvalue=2*ones(1,MRUNS);        %baseline
inlandvalue=1*ones(1,MRUNS);     %baseline
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
% milecost=1.30*ones(1,EXPTRUNS);     %baseline
milecost=X{erun,6}; 
milestraveled=500*ones(1,MRUNS); %baseline
dscnt=X{erun,7};
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\base-chalms-code

end
