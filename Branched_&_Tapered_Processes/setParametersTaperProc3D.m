%%%% setParametersTaperProc3D.m %%%%

% Script to set, compute, and store parameter sets for
%   calcium waves in a in a multi-segment (possibly branching)
%   cellular process model w/ 3D diffusion
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

% seed for random circumferential placement of RyRs
nseed = 8;

% volume ratio, cytosol/ER 
volR = 0.10;

% distal radius
rd0 = 0.15; % units um
% proximal radius
rdf = 1.5; % units um

% taper ratio
tp = 0.05;
% nominal tapered length
ld = (rdf-rd0)/tp; % units um

% length for initial & terminal untapered regions
ld0 = 15; % units um

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
IrR = 3.0E-12/Ce0; % units umol.s^-1.uM^-1
%IrR = 0; % units umol.s^-1.uM^-1

% scaling factor for open probability (1>=Rscale>=0)
Rscale = 1.0;
%Rscale = IrI/IrR;

%% Computation-related parameters

% GRID PARAMETERS

% center-to-center axial spacing of grid compartments (um): universal
dx = 0.5;
% number of tapered compartments in axial direction
ng = round(ld/dx); % ****
% number of untapered axial compartments in terminal segments
ng0 = round(ld0/dx); % ****
% total number of axial compartments
ngt = ng + 2*ng0; % ****
% index into first element of taper
kt1 = ng0 + 1; % ****
% index into last element of taper
kt2 = ng0 + ng; % ****

% set length to accommodate ng tapered compartments exactly
ld = ng * dx; % ****
% set terminal lengths to accomodate ng0 tapered compartments exactly
ld0 = ng0 * dx; % ****
% total process length
ldt = ld + 2*ld0; % ****

% # vol elements in circumferential direction
nc = 6; 
thc = 2*pi()/nc; % angle subtended by vol elements (units radians)

% (GRID-DEPENDENT MORPHOLOGICAL PARAMETERS)

% mean radii of segments
rd = ( rd0 + tp*(0:ng-1)*dx )';
rd = [ (rd0-dx*tp)*ones(ng0,1); rd; (rd(end)+dx*tp)*ones(ng0,1) ];

% squares of radii
rdpow2 = rd.*rd;

% mean ER radius relative to mean process radius
rel = sqrt(volR/(volR+1));

% relative spacing for 3-compartment radial grid
dr = (1-rel)./3;
r2 = rel+dr; % relative outer radius of innermost radial compartment
r3 = r2+dr; % relative inner radius of outermost radial compartments

% annular dr^-2 for radial elements
drm2pw = ( dr * rd ).^-2;

% relative radii for midlines of radial elements
rm = [ rel+dr/2; r2+dr/2; r3+dr/2 ];

% membrane areas
ac = 2*pi()*rd*dx; % plasma membrane areas
ae = rel * ac; % ER membrane areas

% inner & outer radial compartment volumes
vi = pi()*(r2^2-rel^2)*rdpow2*dx; % volumes of innermost radial compartments
vo = pi().*(1-r3^2)*rdpow2*dx; % volumes of outermost radial compartments

% volume ratios, inner radial compartments to ER
volRie = (r2^2-rel^2)/rel^2; % ****
% volume ratio, inner volume elements to ER
volRce = volRie/nc; % ****
% volume ratios, inner radial compartments to entire cytoplasm
volRic = (r2^2-rel^2)/(1-rel^2); % ****
% volume ratios, outer radial compartments to entire cytoplasm
volRoc = (1-r3^2)/(1-rel^2); % ****

% Constants for axial diffusion computations (DERIVED WITH MOLE)
kD = 2; % order of Laplacian computation ****
% get boundary condition coefficients
LB = robinBC(kD, ngt, dx, 0, 1); % kth-order Neumann BC constants
% row matrix to compute right side grad BC, with BV coeff. separated out
JP = LB(end,:);
JPb = full(JP(end));
JP(end) = 0;
% vector to compute boundary value for 0 slope, right side
BR = -JP./JPb;
% vector to compute boundary value for 0 slope, left side
BL = fliplr(BR);
LP = lap(kD, ngt, dx); % matrix to compute Laplacian
% modify Laplacian to account for taper
for k=kt1+1:kt2+1
    LP(k, k-1) = ( 1 - tp*dx/(2*rd(k)) ) * LP(k, k-1);
    LP(k, k+1) = ( 1 + tp*dx/(2*rd(k)) ) * LP(k, k+1);
end

% Constants for circumferential diffusion computations (DERIVED WITH MOLE) 
L1=lap(kD,2*kD+1,1); % minimal-sized sparse Laplacian matrix, order kD
L1 = full(L1);
L1 = L1(kD+2,:); % middle row of L1
L1 = L1(3:end-2); % vector of (2*kD-1) non-zero Laplacian coefficients
kk = 1:2*kD-1;
L0 = zeros(nc,nc);
for k=1:nc
    L0(k, mod((kk+k-kD-1),nc)+1) = L1;
end
% scaling by (circumferential 'length')^-2 (reuse L1; becomes output)
L1 = cell(3,ngt); % cell array to hold scaled Laplacian coefficients
for kk=1:3
    dcm2pw = (thc * rm(kk))^-2;
    for k=1:ngt
        L1{kk,k} = drm2pw(k).*dcm2pw.*L0;
    end
end

% Constants for radial diffusion computations
% compose matrix
R0 = [    -2*r2/(rel+r2),  2*r2/(rel+r2),     0       ];
R0 = [ R0; 2*r2/(r2+r3),      -2,        2*r3/(r2+r3) ];
R0 = [ R0;     0  ,        2*r3/(r3+1),  -2*r3/(r3+1) ];
% scaling by (radial length)^-2
R1 = cell(ngt);
for k=1:ngt
    R1{k} = R0.*drm2pw(k);
end

%% Derived parameters

kvol = 1E15; % volumetric conversion factor, um^3/l

%%% Universal calcium parameters %%%
    % Note all parameters derived with annular volumes and areas, except krR

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
NrR = RrR*ae; % RyRs per annulus (non-integer)
MrR = round(1./NrR); % period: annuli per RyR (for NrR < 1) ****

rng(nseed) % initiate random number generator
% assign an RyR to one of every MrR axial locations; alternating circumf.
%   NOTE: it is assumed that MrR is always >= 1 
IRyR = false(ngt,1); % indicator vector for axial RyR placement
KRyR = false(ngt,nc); % indicator array for axial & circumf. RyR placement
k0 = 5; % first RyR placement
IRyR(k0) = true; % set indicator, first RyR
ncr = ceil( nc * rand(1) ); % choose axial placement randomly
KRyR(k0,ncr) = true;
kl = k0; % indicator of last RyR placement
% place remaining RyRs
for k=k0:ngt
    if (k-kl) >= MrR(k) % place another RyR when local MrR is exceeded
        kl = k;
        IRyR(k) = true;
        ncr = ceil( nc * rand(1) ); % choose axial placement randomly
        KRyR(k,ncr) = true;
    end
end
nRyR = sum(IRyR); % total number of RyRs

% assemble axial indices evoked by KRyR, in order
L=[];
for k=1:nc
    L=[L,[1:ngt]'];
end
LRyR = L(KRyR);

% compute krR for each compartment with RyRs
krR = Rscale * IrR * nc./vi(LRyR) * kvol; % units uM.s^-1 ****;

%%% coefficients for contribution of input flux to boundary values %%%
%   Note: vol/area for cytosolic annuli = dx
BCc = dx/(JPb*DCac);
BInP3 = dx/(JPb*DInP3);

%% Write output

clear k LB re dr r2 r3 ae ac vi vo kD LB L0 dc kvol k0 kl L

save(filename);
