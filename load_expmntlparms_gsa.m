%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Experimental Parameter File   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [am0,am_slope,ampref_max,ampref_min,maxPflood,highrisk,stormfreq,maxdam,...
    Cmit,miteff,AVGFARMRETURN,STDFARMRETURN,coastvalue,midvalue,...
    inlandvalue,milecost,milestraveled,discountrate]=load_expmntlparms_gsa(EXPTRUNS,MRUNS)

load 'C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files\gsaparms.mat'

% Coastal Amenity
am0=cat(1,X{:,3});
% am0=500000*ones(1,EXPTRUNS);        %baseline
am_slope=cat(1,X{:,4});
% am_slope=0.1*ones(1,EXPTRUNS);      %baseline

discountrate=cat(1,X{:,7});
% discountrate=0.05*ones(1,EXPTRUNS*MRUNS);   %baseline

% Consumer preferences
ampref_max=cat(1,X{:,5});
% ampref_max=0.9*ones(1,EXPTRUNS);    %baseline
ampref_min=0.1*ones(1,EXPTRUNS*MRUNS);    %baseline

% Storm Impacts and Mitigation
maxPflood=0.7*ones(1,EXPTRUNS*MRUNS);     %baseline
highrisk=30*ones(1,EXPTRUNS*MRUNS);       %baseline
maxdam=cat(1,X{:,2}); 
% maxdam=ones(1,EXPTRUNS);            %baseline
% stormfreq=ones(1,EXPTRUNS);         %baseline
stormfreq=cat(1,X{:,1});         

% stormthresh=15*ones(1,EXPTRUNS);
% Cdam=0.5*ones(1,EXPTRUNS);
Cmit=3000*ones(1,EXPTRUNS*MRUNS);         %baseline
miteff=1*ones(1,EXPTRUNS*MRUNS);          %baseline
    
% Initial land value
AVGFARMRETURN=2486.3*ones(1,EXPTRUNS*MRUNS);
STDFARMRETURN=0.10*AVGFARMRETURN.*ones(1,EXPTRUNS*MRUNS);
% lvset=mat2cell([3; 2; 1]*linspace(0.5,1.5,EXPTRUNS),ones(3,1),EXPTRUNS);    
% lvset=[3; 2; 1]*ones(1,EXPTRUNS);`%baseline
coastvalue=3*ones(1,EXPTRUNS*MRUNS);
midvalue=2*ones(1,EXPTRUNS*MRUNS);
inlandvalue=1*ones(1,EXPTRUNS*MRUNS);
% coastvalue=lvset{1};
% midvalue=lvset{2};
% inlandvalue=ones(1,EXPTRUNS);

% Travel costs
milecost=cat(1,X{:,6});
% milecost=1.30*ones(1,EXPTRUNS);     %baseline
milestraveled=500*ones(1,EXPTRUNS*MRUNS); %baseline

end

