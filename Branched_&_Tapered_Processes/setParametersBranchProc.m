%%%% setParametersMultiSegProc.m %%%%

% Script to set, compute, and store parameter sets for
%   calcium waves in a in a multi-segment (possibly branching)
%   cellular process model
% Distal (smaller diameter) segments to the left,
%   proximal (larger diameter) segments to the right.
% User adjusts user-settable parameters in script
% When run, script prompts for unique output filename for .mat file

clear;

filename = input('Specify name of output file: ','s');

%% Reference level for ER calcium (fixed)

Ce0 = 400; % units uM

%% User-settable (at present) model parameter values

% MORPHOLOGICAL PARAMETERS

% volume ratio, cytosol/ER 
volR = 0.1;

% segment radii: must have 1 column per segment in model
% for check runs
%rd = 0.1667.*[1,1]; % pre-junction size
%rd = 0.1667*[sqrt(2),sqrt(2)]; % post-junction size
% for examining branches
rd = 0.1667.*[1,1,sqrt(2)]; % units um

% ns = number of segments
ns = length(rd);

% segment lengths must have 1 column per segment in model
% for check runs
%ld = [32,32];
% for examining branches
ld = [32,32,32]; % units um

% cnx matrix must have 1 column per junction in model
%   junction = interconnection of 2 or 3 segments 
%   entries in ea column = indices of proximal, distal, distal segments
%   0 in row 3 allows simple junction of 2 segments only
% for check runs
%cnx = [2;1;0];
% for examining branches
cnx = [3;1;2];

% nx = number of junctions
[k,nx] = size(cnx);

% identify dendrite ends for sealed boundary conditions (non-junctions)

sbcl = true(1,ns);
sbcr = true(1,ns);

for k=1:nx
    
    if cnx(1,k)
        sbcl(cnx(1,k)) = false;
    end
    if cnx(2,k)
        sbcr(cnx(2,k)) = false;
    end
    if cnx(3,k)
        sbcr(cnx(3,k)) = false;
    end
        
end

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
kbCf = 1.9; % multiplies [Ca] & [CalB]; units s^-1.uM^-1 ****
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
RrI = 10; % units um^-2

% InP3 diffusion reduction parameter
fInP3 = 0.5; % fraction to reduce DInP3 due to intracellular crowding

% gating function delay rate constant
rIg = 70; % units s^-1

% rate constant per unit area, [Ca]-dependent production of InP3
rmIf = 1.4E-15; % units umol.s^-1.um^-2.uM^-1
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

% cross-sectional area
ax = pi().*rd.*rd; % ****

% ER radius
re = sqrt(volR/(volR+1)).*rd; % units um

% Numerical/run control parameters

% GRID PARAMETERS

% center-to-center axial spacing of grid compartments (um): universal
%dx = 0.5; % SHOULD GO EVENLY INTO ld!
dx = 0.4; % SHOULD GO EVENLY INTO ld!
% number of compartments in axial direction
ng = ld./dx; % ****

% spacing for 3-compartment radial grid
dr = (rd-re)./3; 
r2 = re+dr; % outer radii of innermost radial compartments
r3 = r2+dr; % inner radii of outermost radial compartments

% membrane areas
ae = 2*pi().*re.*dx; % ER membrane area
ac = 2*pi().*rd.*dx; % plasma membrane area

% total compartment volumes
vole = pi().*re.*re.*dx; % ER compartment volumes
volc = pi().*rd.*rd.*dx - vole; % cytosol compartment volumes

% inner & outer radial compartment volumes
vi = pi().*(r2.*r2-re.*re).*dx; % volumes of innermost radial compartments
vo = pi().*(rd.*rd-r3.*r3).*dx; % volumes of outermost radial compartments

% volume ratios, inner radial compartments to ER
volRie = vi./vole; % ****
% volume ratios, inner radial compartments to entire cytoplasm
volRic = vi./volc; % ****
% volume ratios, outer radial compartments to entire cytoplasm
volRoc = vo./volc; % ****

% Constants for axial diffusion computations (DERIVED WITH MOLE)

% allocate cell arrays
LP = cell(1,ns);
JPL = cell(1,ns);
JPR = cell(1,ns);
BL = cell(1,ns);
BR = cell(1,ns);
RD = cell(1,ns);

kD = 2; % order of Laplacian computation

% get boundary value coefficient (common for all segments)
LB = robinBC(kD, 2*kD+1, dx, 0, 1); % Neumann BC constants
JPb = full(LB(1,1));

% get Laplacian and boundary condition-related matrices/vectors
for k=1:ns
    
    LP{1,k} = lap(kD, ng(k), dx); % matrix to compute Laplacian ****
     
    % get boundary condition coefficients
    LB = robinBC(kD, ng(k), dx, 0, 1); % Neumann BC constants
    % row matrix to compute left side grad BC, w/ BV coeff. separated out
    JPL{1,k} = [ 0, LB(1,2:end) ];
    % row matrix to compute right side grad BC, w/ BV coeff. separated out
    JPR{1,k} = [ LB(end,1:end-1), 0 ];
    % row matrix to compute boundary value for 0 slope, right side
    BL{1,k} = -JPL{1,k}./JPb;
    % row matrix to compute boundary value for 0 slope, left side
    BR{1,k} = -JPR{1,k}./JPb;
end

% modify Laplacian:
%   account for smaller cross-section of upstream segments at junctions
for k=1:nx
    LP1 = LP{1,cnx(1,k)}(2,1);
    if cnx(3,k) % rescale leading coefficient for branch
        LP1 = (ax(cnx(2,k))+ax(cnx(3,k)))/ax(cnx(1,k)) * LP1;
    else % rescale leading coefficient for simple junction
        LP1 = ax(cnx(2,k))/ax(cnx(1,k)) * LP1;
    end
    LP{1,cnx(1,k)}(2,1) = LP1; % replace leading coefficient
    % modify 2nd coefficient to maintain conservation of mass
    LP{1,cnx(1,k)}(2,2) = -( LP1 + LP{1,cnx(1,k)}(2,3) );    
end


% Constants for radial diffusion computations
RD = cell(1,ns);
for k=1:ns
    % compose matrix
    R0 = [    -2*r2(k)/(re(k)+r2(k)),  2*r2(k)/(re(k)+r2(k)),     0      ];
    R0 = [ R0; 2*r2(k)/(r2(k)+r3(k)),    -2,       2*r3(k)/(r2(k)+r3(k)) ];
    R0 = [ R0;     0  ,    2*r3(k)/(r3(k)+rd(k)), -2*r3(k)/(r3(k)+rd(k)) ];
    
    RD{1,k} = R0./dr(k)^2; % matrix to compute radial diffusion ****;
    
end

%% Derived parameters

kvol = 1E15; % volumetric conversion factor, um^3/l

%%% Universal calcium parameters %%%

% Calcium molar membrane leakage rates
% rate to cytosol from extracellular space
klx = rlx .* ac./vo .* kvol; % units uM.s^-1 ****
% rate to cytosol from ER
kle = rle .* ae./vi .* kvol; % units s^-1 ****

% Ca diffusion coefficients
% in cytosol
DCac = DCa*fCac; % units um^2.s^-1 ****
% in ER
DCae = DCa*fCae; % units um^2.s^-1 ****

%%% Calcium pump parameters %%%

% Calcium NCX pump constant
% NpN = round(RpN.*ac);
% kpN = NpN.*IpN ./vo .* kvol; % units uM.s^-1 ****
% following: non-discretized NCX
kpN = RpN.*IpN .* ac./vo .* kvol; % units uM.s^-1 ****

% Calcium PMCA pump constant
kpP = RpP.*IpP .* ac./vo .* kvol; % units uM.s^-1 ****

% Calcium SERCA pump constants
kpS =  RpS.*IpS .* ae./vi .* kvol;  % units uM.s^-1 ****

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3 = DInP3*fInP3;  % units um^2.s^-1 ****

% Ca molar flow rate through InP3R channels
% NrI = round(RrI.*ae);
% krI = (NrI.*IrI) ./vi .* kvol; % units s^-1. ****
% following: non-discretized InP3Rs
krI = RrI.*IrI * ae./vi .* kvol; % units s^-1. ****

% rate constant, [Ca]-dependent production of InP3
kmIf = rmIf .* ac./vo .* kvol; % multiplies [Ca]; units s^-1 ****

%%% Ryanodine receptor parameters %%%

% Set RyR molar flow rate parameters according to density
NrR = RrR.*ae;
MrR = round(1./NrR); % period: compartments per RyR (for NrR < 0.5) ****
NrR = round(NrR);   % frequency: RyRs per compartment (for NrR > 0.5) ****
Drflag = ( NrR>=1 ); % indicator of multiple RyRs per compartment ****
IRyR = cell(1,ns);
for k=1:ns
    if Drflag(k) % if there is more than 1 RyR per compartment
        krR(k) = (NrR(k)*IrR) /vi(k) * kvol; % units uM.s^-1 ****
    else % if there is less than 1 RyR per compartment
        % krR for a single receptor
        krR(k) = IrR /vi(k) * kvol; % units uM.s^-1 ****
        % assign an RyR to one out of every MrR compartments
        IRk = zeros(ng(k),1);
        Mr0 = max(floor(MrR(k)/2),1)+1;
        for kk=Mr0:MrR(k):ng(k)
            IRk(kk,1) = 1; % indicator array for compartments w/ RyRs ****
        end
        IRyR{1,k} = IRk;
    end
end

%%% coefficients for contribution of input flux to boundary values %%%
%   Note: vol/area for cytosolic annuli = dx
BCc = dx/(JPb*DCac);
BInP3 = dx/(JPb*DInP3);

%% Write output

%clear k LB R0
clear k re dr r2 r3 ae ac ve vc vi vo LB R0 kvol NpN Mr0 IrK

save(filename);
