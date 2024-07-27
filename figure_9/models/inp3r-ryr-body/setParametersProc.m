function out = setParametersProc(proc_outfile,varargin)
%%%% setParametersProc.m %%%%

% Script to set, compute, and store parameter sets for
%   calcium waves in a cellular process model
% User adjusts user-settable parameters in script
% When run, script prompts for unique output filename for .mat file

%% Parse input parameters
p = inputParser;
%kbCf: Ca buffering forward rate. 
addOptional(p,'kbCf',1.9);% Multiplies [Ca] & [CalB]; units s^-1.uM^-1 ****
%rmIf: rate constant per unit area, [Ca]-dependent production of InP3
addOptional(p,'rmIf',1.4E-15);% units umol.s^-1.um^-2.uM^-1
addOptional(p,'rIg',70);% gating function delay rate constant. Units s^-1
addOptional(p,'RrIb',10);
parse(p,varargin{:});
kbCf=p.Results.kbCf;
rmIf=p.Results.rmIf;
rIg=p.Results.rIg;
RrI=p.Results.RrIb;%Always set RrI=RrIb

%filename = input('Specify name of output file: ','s');
filename = 'parameters_proc';
%% Reference level for ER calcium (fixed)

Ce0 = 400; % units uM

%% User-settable (at present) model parameter values

% MORPHOLOGICAL PARAMETERS

% volume ratio, cytosol/ER 
volR = 0.1;

% process radius
%rd = 0.5; % units um
rd = 1.5; % units um

% process length
%ld = 100; % units um ****
% process lengths for body simulations
ld = 10; % units um ****
% minimal test process
%ld = 2.5; % units um; applies when dx = 0.5um

% Ca TRANSPORT/REACTION PARAMETERS

%%% Universal calcium parameters %%%

% Calcium membrane leakage rates per unit area
% plasma membrane (assumed constant)
rlx = 5E-16; % units umol.s^-1.um^-2
% ER membrane
rle = 5E-16/Ce0; % units umol.s^-1.um^-2.uM^-1

% Ca diffusion reduction parameters
fCac = 0.3; % fraction to reduce cyto. DCa due to intracellular crowding
fCae = 0.6; % fraction to reduce ER DCa due to intracellular crowding

% Ca buffering forward rate
%kbCf = 1.9; % multiplies [Ca] & [CalB]; units s^-1.uM^-1 ****
% Ca buffering release rate
kbCb = 19; % multiplies [CaCalB]; units s^-1 ****

%%% Calcium pump parameters %%%

% Calcium NCX pump density
RpN = 14; % units um^-2

% Calcium PMCA pump density
RpP = 500; % units um^-2

% Calcium SERCA pump density
RpS =  4000;  % units um^-2

%%% InP3 receptor and related parameters %%%

% InP3R density
%RrI = 10; % units um^-2

% InP3 diffusion reduction parameter
fInP3 = 0.5; % fraction to reduce DInP3 due to intracellular crowding

% gating function delay rate constant
%rIg = 70; % units s^-1

% rate constant per unit area, [Ca]-dependent production of InP3
%rmIf = 1.4E-15; % units umol.s^-1.um^-2.uM^-1
% rate constant, degradation of InP3
kmIb = 2.5; % units s^-1 (multiplies [InP3]) ****

% constant, inhibition of InP3 production by activated PKC
kiPI = 0.0943; % multiplies [PKCv]; units uM^-1 ****

% available PKC concentration
PKC0 = 1.0; % units uM ****

% rate for activation of PKC
kiPf = 0.6; % multiplies [Ca] & [PKC]; units s^-1.uM^-1 ****
% rate for deactivation of PKCv
kiPb = 0.5; % multiplies [PKCv]; units s^-1 ****

%%% Ryanodine receptor parameters %%%

% RyR receptor density
RrR = 0.4; % units um^-2

%save([filename,'_assigned']);

%% Fixed (at present) model parameter values

%%% Universal calcium parameters %%%

% free Ca diffusion coefficient in aqueous solution
DCa = 240; % units um^2.s^-1

% available Calbindin concentration 
CalB0 = 40; % units = uM ****

%%% Calcium pump parameters %%%

% NCX pump constants
% specific reference current, individual pump
IpN = 2.5E-15; % units umol.s^-1
% saturation concentration
kpNC = 1.8; % units uM ****

% PMCA pump constants
% specific reference current, individual pump
IpP = 1.7E-17; % units umol.s^-1
% saturation concentration
kpPC = 0.24; % units uM
kpPC2pw = kpPC^2; % square of saturation concentration; units uM^2 ****

% SERCA pump constants
% specific reference current parameter, individual pump
IpS =  7.0E-18*Ce0;  % units umol.s^-1.uM
% saturation concentration
kpSC = 0.10; % units uM ****

%%% InP3 receptor and related parameters %%%

% free InP3 diffusion coefficient in aqueous solution
DInP3 = 300; % units um^2.s^-1

% InP3R channel reference current parameter
IrI = 6.5E-13/Ce0; % units umol.s^-1.uM^-1
%IrI = 0; % units umol.s^-1.uM^-1

%%% Ryanodine receptor parameters %%%

% RyR channel reference current 
%IrR = 3.0E-12/Ce0; % units umol.s^-1.uM^-1
IrR = 0; % units umol.s^-1.uM^-1

%% Computation-related parameters

% MORPHOLOGICAL PARAMETERS

% ER radius
re = sqrt(volR/(volR+1))*rd; % units um

% GRID PARAMETERS

% center-to-center axial spacing of grid compartments (um)
dx = 0.5; % SHOULD GO EVENLY INTO ld! ****
% number of compartments in axial direction
ng = ld/dx; % ****

% spacing for 3-compartment radial grid
drp = (rd-re)/3;
% put bounding radii into array
rxp = re;
for k=1:3
    rxp = [rxp, rxp(end)+drp];
end

% membrane areas
ae = 2*pi()*re*dx; % ER membrane area
ac = 2*pi()*rd*dx; % plasma membrane area

% cross-sectional areas of radial elements
for k=1:3
    axcp(k) = pi()*(rxp(k+1)^2 - rxp(k)^2);
end

% total compartment volumes
vole = pi()*re*re*dx; % ER compartment volume
volc = pi()*rd*rd*dx - vole; % cytosol compartment volume

% radial compartment (annulus) volumes
for k=1:3
    volr(k) = axcp(k)*dx;
end
vi = volr(1);  % volume of innermost radial compartment
vo = volr(3);  % volume of outermost radial compartment

% volume ratio, inner radial compartment to ER
volRie = vi/vole; % ****
% volume ratio, inner radial compartment to entire cytoplasm
volRic = vi/volc; % ****
% volume ratio, outer radial compartment to entire cytoplasm
volRoc = vo/volc; % ****

% Constants for axial diffusion computations (DERIVED WITH MOLE)
kD = 2; % order of Laplacian computation ****
LP = lap(kD, ng, dx); % matrix to compute Laplacian ****
% get boundary condition coefficients
LB = robinBC(kD, ng, dx, 0, 1); % kth-order Neumann BC constants
% row matrix to compute right side grad BC, with BV coeff. separated out
JP = LB(end,:);
JPb = full(JP(end));
JP(end) = 0;
% vector to compute boundary value for 0 slope, right side
BR = -JP./JPb;
% vector to compute boundary value for 0 slope, left side
BL = fliplr(BR);

% Constants for (2nd order) radial diffusion computations
  % compose matrix f
R0 = [  -2*rxp(2)/(rxp(1)+rxp(2)), 2*rxp(2)/(rxp(1)+rxp(2)),      0     ];
R0 = [ R0;  2*rxp(2)/(rxp(2)+rxp(3)),   -2,    2*rxp(3)/(rxp(2)+rxp(3)) ];
R0 = [ R0;  0  ,  2*rxp(3)/(rxp(3)+rxp(4)),   -2*rxp(3)/(rxp(3)+rxp(4)) ];
RD = R0./drp^2; % matrix to compute radial diffusion ****

%% Derived parameters

kvol = 1E15; % volumetric conversion factor, um^3/l

%%% Universal calcium parameters %%%

% Calcium molar membrane leakage rates
% rate to cytosol from extracellular space
klx = rlx * ac/vo * kvol; % units uM.s^-1 ****
% rate to cytosol from ER
kle = rle * ae/vi * kvol; % units s^-1 ****

% Ca diffusion coefficients
% in cytosol
DCac = DCa*fCac; % units um^2.s^-1 ****
% in ER
DCae = DCa*fCae; % units um^2.s^-1 ****

%%% Calcium pump parameters %%%

% Calcium NCX pump constant
NpN = round(RpN*ac);
kpN = NpN*IpN /vo * kvol; % units uM.s^-1 ****

% Calcium PMCA pump constant
kpP = RpP*IpP * ac/vo * kvol; % units uM.s^-1 ****

% Calcium SERCA pump constants
kpS =  RpS*IpS * ae/vi * kvol;  % units uM.s^-1 ****

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3 = DInP3*fInP3;  % units um^2.s^-1 ****

% Ca molar flow rate through InP3R channels
NrI = round(RrI*ae);
krI = (NrI*IrI) /vi * kvol; % units s^-1. ****

% rate constant, [Ca]-dependent production of InP3
kmIf = rmIf * ac/vo * kvol; % multiplies [Ca]; units s^-1 ****

%%% Ryanodine receptor parameters %%%

% Set RyR molar flow rate parameters according to density
NrR = RrR*ae;
MrR = round(1/NrR); % period: compartments per RyR (for NrR < 0.5) ****
NrR = round(NrR);   % frequency: RyRs per compartment (for NrR > 0.5) ****
%Drflag = ( NrR>=1 ); % indicator of multiple RyRs per compartment ****
Drflag = ( NrR>1 ); % indicator of multiple RyRs per compartment ****
if Drflag % if there is more than 1 RyR per compartment
    krR = (NrR*IrR) /vi * kvol; % units uM.s^-1 ****
else % if there is less than 1 RyR per compartment
    % krR for a single receptor
    krR = IrR /vi * kvol; % units uM.s^-1 ****
    % assign an RyR to one out of every MrR compartments 
    IRyR = zeros(ng,1);
    Mr0 = max(floor(MrR/2),1)+1;
    for k = Mr0:MrR:ng
        IRyR(k,1) = 1; % indicator array for compartments with RyRs ****
    end
    
end

%%% coefficients for contribution of input flux to boundary values %%%
%   Note: vol/area for cytosolic annuli = dx
BCc = dx/(JPb*DCac);
BInP3 = dx/(JPb*DInP3);

%% Write output

clear k LB R0 re ae ac vole volc vi vo kD LB R0 kvol varargin

save(proc_outfile);
