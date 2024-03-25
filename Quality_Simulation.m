% Clear workspace and command window
clear all; close('all'); clc;

% Load toolkit paths
start_toolkit;

% Load the network
inpname = 'skeletonized_intodbp_pilot_prepared.inp';
G = epanet(inpname);

% Load MSX file with reactions
msxname = 'THMs_HAAs.msx';
G.loadMSXFile(msxname);

%% Sensor locations
link_description={'WTP outlet pipe','Before Second Chlorination pipe','Zone2 outlet pipe'};
link_virtual_id = {'N24','N31','p-1748'};
link_index = G.getLinkIndex(link_virtual_id);

sensor_description={'DMA point','DMA Entrance' ,'Zone2 Tank','WTP reservoir'};
sensor_id = {'dist412','dist1268','T_Zone2','WTP'};
sensor_index = G.getNodeIndex(sensor_id);

%% Get Network data
demand_pattern = G.getPattern;
roughness_coeff = G.getLinkRoughnessCoeff;

%% Injection locations
Zone2Inlet={'N31'};
Zone2Chlorination={'414'};
Zone2Level={'T_Zone2'};
Zone2TO4={'232'};
WTP={'WTP'};

IndexTank=G.getNodeIndex(Zone2Level);
IndexChlorination=G.getNodeIndex(Zone2Chlorination);

IndexInlet=G.getLinkIndex(Zone2Inlet);

IndexZ4=G.getLinkIndex(Zone2TO4);
IndexWTP=G.getNodeIndex(WTP);

node_id = G.getNodeNameID;
resInd = G.getNodeReservoirIndex;

%% Simulation Setup
t_d = 7;% day
t_sim = t_d*24*60*60;
G.setTimeSimulationDuration(t_sim);

%% Load Sensor Data
load DMAcenterorganics1622102023
load DMAnodeorganics1622102023
load DWTPorganics1622102023
load Reservoirorganics1622102023


%% Load Measurements Data
load reservoir_chlorine.mat

%% Chlorine and TOC initial conditions

%Adjust initial conditions on reservoir
InletChlorine=DWTPorganics1622102023{1:t_d*12*24,10};
InletTOC=DWTPorganics1622102023{1:t_d*12*24,3};
MeanInletChlorine=mean(InletChlorine);
MeanInletTOC=mean(InletTOC);

%Adjust initial conditions on tank
InletChlorineTank=Reservoirorganics1622102023{1:t_d*12*24,10};
InletTOCTank=Reservoirorganics1622102023{1:t_d*12*24,3};
MeanInletChlorineTank=mean(InletChlorineTank);
MeanInletTOCTank=mean(InletTOCTank);
InitialChlorine=InletChlorineTank(1);
InitialTOC=InletTOCTank(1);

speciesIndex=[1;2];
values = G.getMSXNodeInitqualValue;
values{IndexTank}(speciesIndex)=[InitialChlorine;InitialTOC];
G.setMSXNodeInitqualValue(values) 
values = G.getMSXNodeInitqualValue;

%Add chlorination
G.setMSXSources(node_id(IndexWTP), 'CL2' , 'Concen' ,MeanInletChlorine, 'CL2PAT' ); % Specify Chlorine injection source
G.setMSXSources(node_id(IndexWTP), 'TOC' , 'Concen' ,MeanInletTOC, 'TOCPAT' );

G.setMSXSources(node_id(IndexTank), 'CL2' , 'Mass' ,200.24, 'CL2PATT' ); % Specify Chlorine injection source
% G.setMSXSources(node_id(IndexTank), 'TOC' , 'Concen' ,MeanInletTOCTank, 'TOCPATT' );

%% Set Patterns

CL2_pat = InletChlorine/MeanInletChlorine;
TOC_pat = InletTOC/MeanInletTOC;
G.setMSXPattern('CL2PAT', CL2_pat); % Set pattern for Chlorine injection
G.setMSXPattern('TOCPAT', TOC_pat); % Set pattern for Chlorine injection

% CL2_patt = InletChlorineTank/MeanInletChlorineTank;
% TOC_patt = InletTOCTank/MeanInletTOCTank;
% G.setMSXPattern('CL2PATT', CL2_patt); % Set pattern for Chlorine injection
% G.setMSXPattern('TOCPATT', TOC_patt); % Set pattern for Chlorine injection

%%
% Set msx time step
msx_timestep = 300; % seconds

% Get species names
species = G.getMSXSpeciesNameID;

% Get quality for specific species type (nodes and links).
MSX_comp=cell(length(species),1); % initialize a cell array to store the results;
for i = 1: length (species) 
    MSX_comp{i} = G.getMSXComputedQualitySpecie(species{i});
end

% Allspecies=G.getMSXComputedQualityNode;

% Save results
% G.saveInputFile('savenet.inp');

%% Plots
% Sensor measurements
SensorDateTime=DWTPorganics1622102023{1:t_d*12*24,1}; 

MeasuredTOC(:,1)=DMAcenterorganics1622102023{1:t_d*12*24,3};
MeasuredTOC(:,2)=DMAnodeorganics1622102023{1:t_d*12*24,3};
MeasuredTOC(:,3)=Reservoirorganics1622102023{1:t_d*12*24,3};
MeasuredTOC(:,4)=DWTPorganics1622102023{1:t_d*12*24,3};

MeasuredChlorine(:,1)=DMAcenterorganics1622102023{1:t_d*12*24,10};
MeasuredChlorine(:,2)=DMAnodeorganics1622102023{1:t_d*12*24,10};
MeasuredChlorine(:,3)=Reservoirorganics1622102023{1:t_d*12*24,10};
MeasuredChlorine(:,4)=DWTPorganics1622102023{1:t_d*12*24,10};

%Sampling measurements
% StartTime=SensorDateTime(1); EndTime=SensorDateTime(end);
% 
% Tinlet=WBLZone2reservoirfreechlorineresidualsampling{:,1};
% Toutlet=WBLZone2reservoirfreechlorineresidualsampling{:,6};
% tf1 =isbetween(Tinlet,StartTime,EndTime);
% tf2=isbetween(Toutlet,StartTime,EndTime);
% tfpos1=find(tf1==1);
% tfpos2=find(tf2==1);

InletResDateTime=WBLZone2reservoirfreechlorineresidualsampling{144,6};
InletResChlorine=mean(WBLZone2reservoirfreechlorineresidualsampling{72:73,4});

OutletResDateTime=WBLZone2reservoirfreechlorineresidualsampling{144,6};
OutletResChlorine=WBLZone2reservoirfreechlorineresidualsampling{144,9};
%%
figure
k=1;
for i=length(sensor_index):-1:1

subplot(4,1,k)
plot(SensorDateTime,MeasuredTOC(:,i))
hold on
plot(SensorDateTime, MSX_comp{2}.NodeQuality(1:length(SensorDateTime),sensor_index(i)))
ylabel('mg/L') 
xlabel('Datetime');
title(sensor_description{i})
k=k+1;
end
legend({'Sensor Measured','Simulated'})
sgtitle('TOC')


figure
k=1;
for i=length(sensor_index):-1:1

subplot(4,1,k)
plot(SensorDateTime,MeasuredChlorine(:,i))
hold on
plot(SensorDateTime, MSX_comp{1}.NodeQuality(1:length(SensorDateTime),sensor_index(i)))
ylabel('mg/L') 
% xlabel('Datetime');
title(sensor_description{i})
k=k+1;
end
legend({'Sensor Measured','Simulated'})
sgtitle('Chlorine')

figure
plot(SensorDateTime, MSX_comp{1}.LinkQuality(1:length(SensorDateTime),link_index(1)))
hold on
plot(SensorDateTime,MeasuredChlorine(:,4))
ylabel('mg/L') 
ylim([0.15 0.5])
title(link_description{1})
legend({'Simulated','Sensor Measured'})

figure
plot(SensorDateTime, MSX_comp{1}.LinkQuality(1:length(SensorDateTime),link_index(2)))
hold on
scatter(InletResDateTime,InletResChlorine,100,'k','x')
ylabel('mg/L') 
ylim([0.15 0.5])
title(link_description{2})
legend({'Simulated','Sampling'})

figure
plot(SensorDateTime, MSX_comp{1}.LinkQuality(1:length(SensorDateTime),link_index(3)))
hold on
plot(SensorDateTime,MeasuredChlorine(:,3))
hold on
scatter(OutletResDateTime,OutletResChlorine,100,'k','x')
ylabel('mg/L') 
ylim([0.05 0.5])
title(link_description{3})
legend({'Simulated','Sensor Measured','Sampling'})


% subplot(4,1,3)
% plot(SensorDateTime, MSX_comp{3}.NodeQuality(1:length(SensorDateTime),sensor_index))
% ylabel('ug/L') 
% xlabel('Datetime');
% legend({'Sensor Simulated'})
% title('THMs')


%% Map plot
MaxCL=max(MSX_comp{1}.NodeQuality);
MaxTOC=max(MSX_comp{2}.NodeQuality);
MaxTHM=max(MSX_comp{3}.NodeQuality);
MaxHAAs=max(MSX_comp{4}.NodeQuality);
coords = G.getNodeCoordinates;

%CL
G.plot('legendposition', 'best')
hold on
for i=1:length(MaxCL)
    scatter(coords{1}(i),coords{2}(i),35, MaxCL(i), 'filled')
end
title('Max Chlorine concentration')
c = colorbar('southoutside');
c.Label.String = ['Chlorine (', G.MSXSpeciesUnits{1}, ')'];

%TOC
G.plot('legendposition', 'best')
hold on
for i=1:length(MaxTOC)
    scatter(coords{1}(i),coords{2}(i),35, MaxTOC(i), 'filled')
end
title('Max TOC concentration')
c = colorbar('southoutside');
c.Label.String = ['TOC (', G.MSXSpeciesUnits{2}, ')'];

%THMS
G.plot('legendposition', 'best')
hold on
for i=1:length(MaxTHM)
    scatter(coords{1}(i),coords{2}(i),35, MaxTHM(i), 'filled')
end
title('Max THMs concentration')
c = colorbar('southoutside');
c.Label.String = ['THMs (', G.MSXSpeciesUnits{3}, ')'];

%HAAs
G.plot('legendposition', 'best')
hold on
for i=1:length(MaxHAAs)
    scatter(coords{1}(i),coords{2}(i),35, MaxHAAs(i), 'filled')
end
title('Max HAAs concentration')
c = colorbar('southoutside');
c.Label.String = ['HAAs (', G.MSXSpeciesUnits{3}, ')'];

 %Unload libraries
G.unloadMSX;
G.unload;
