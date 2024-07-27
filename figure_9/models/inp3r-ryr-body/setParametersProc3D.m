function out = setParametersProc3D(proc_outfile,varargin)
%%%% setParametersProc3D.m %%%%

% Script to set, compute, and store parameter sets for
%   calcium waves in a cellular process model
% User adjusts user-settable parameters in script
% When run, script prompts for unique output filename for .mat file

%% Reference level for ER calcium (fixed)
Ce0 = 400; % units uM

%% Parse input parameters
p = inputParser;
addOptional(p,'rd',0.5);% process radius
addOptional(p,'IrR',3.0E-12/Ce0);% Ryr current. Used to turn RyR on/off
addOptional(p,'volR',0.1);% volume ratio, cytosol/ER 
addOptional(p,'RrI',10);% InP3R density, units um^-2
addOptional(p,'RrR',1.5);% RyR receptor density, units um^-2
addOptional(p,'rlx',5E-16);%plasma membrane Ca leakage, units umol.s^-1.um^-2
addOptional(p,'Ce',200);% Initial ER Ca, uM

parse(p,varargin{:});
rd=p.Results.rd;%
IrR=p.Results.IrR;%units umol.s^-1.uM^-1
volR = p.Results.volR;
RrI = p.Results.RrI;
RrR = p.Results.RrR;
rlx = p.Results.rlx;
Ce = p.Results.Ce;

%filename = input('Specify name of output file: ','s');
filename = proc_outfile%'parameters_proc';



%% User-settable (at present) model parameter values

% Small calcium concentrations used to ensure convergence to equilibrium
Cc = 0.01; % cytosolic Ca; units uM
%Ce = 200; % ER Ca; units uM

% MORPHOLOGICAL PARAMETERS

% volume ratio, cytosol/ER 
%volR = 0.1;

% process radius
%rd = 0.5; % units um
%rd = 1.5; % value for attaching process; units um

% process length for dendrite sims
ld = 100; % units um ****
%ld = 25; % units um ****
% process lengths for body simulations
%ld = 10; % units um ****
% minimal test process
%ld = 2.5; % units um; applies when dx = 0.5um

% seed for random circumferential placement of RyRs
nseed = 2;

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
%RrI = 10; % units um^-2

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
%RrR = 0.4; % units um^-2
%RrR = 1.5; % units um^-2

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
%IrR = 0; % units umol.s^-1.uM^-1

% scaling factor for open probability (1>=Rscale>=0)
Rscale = 1.0;
%Rscale = IrI/IrR;

%% Computation-related parameters

% MORPHOLOGICAL PARAMETERS

% ER radius
re = sqrt(volR/(volR+1))*rd; % units um

% GRID PARAMETERS

% center-to-center axial spacing of grid compartments (um)
dx = 0.5; % SHOULD GO EVENLY INTO ld! ****
%dx = 0.25; % SHOULD GO EVENLY INTO ld! ****
% number of compartments in axial direction
ng = ld/dx; % ****

% spacing for 3-compartment radial grid
drp = (rd-re)/3;
% put bounding radii into array
rxp = re;
for k=1:3
    rxp = [rxp, rxp(end)+drp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # vol elements in circumferential direction
nc = 6; 
thc = 2*pi()/nc; % angle subtended by vol elements (units radians)

% membrane areas of annular elements
ae = 2*pi()*re*dx; % ER membrane area
ac = 2*pi()*rd*dx; % plasma membrane area

% cross-sectional areas of radial (annular) elements
for k=1:3
    axcp(k) = pi()*(rxp(k+1)^2 - rxp(k)^2);
end

% total compartment volumes
vole = pi()*re*re*dx; % ER compartment volume
volc = pi()*rd*rd*dx - vole; % cytosol annulus volume

% radial compartment (annulus) volumes
for k=1:3
    volr(k) = axcp(k)*dx;
end
vi = volr(1);  % volume of innermost radial annulus
vo = volr(3);  % volume of outermost radial annulus

% volume ratio, inner radial annulus to ER
volRie = vi/vole; % ****
% volume ratio, inner volume elements to ER
volRce = volRie/nc; % ****
% volume ratio, inner radial annulus to entire cytoplasm
volRic = vi/volc; % ****
% volume ratio, outer radial annulus to entire cytoplasm
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% scaling by (circumferential 'length')^-2
LC = cell(3,1); % cell array to hold scaled Laplacian coefficients
for k=1:3
    dc = thc * (rxp(k) + 0.5*drp);
    LC{k} = L0./dc^2;
end

% Constants for (2nd order) radial diffusion computations
  % compose matrix f
R0 = [  -2*rxp(2)/(rxp(1)+rxp(2)), 2*rxp(2)/(rxp(1)+rxp(2)),      0     ];
R0 = [ R0;  2*rxp(2)/(rxp(2)+rxp(3)),   -2,    2*rxp(3)/(rxp(2)+rxp(3)) ];
R0 = [ R0;  0  ,  2*rxp(3)/(rxp(3)+rxp(4)),   -2*rxp(3)/(rxp(3)+rxp(4)) ];
RD = R0./drp^2; % matrix to compute radial diffusion ****

%% Derived parameters

kvol = 1E15; % volumetric conversion factor, um^3/l

%%% Universal calcium parameters %%%
    % Note all parameters based on annular volumes and areas, except krR

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
NrR = RrR*ae; % RyRs per annulus
MrR = round(1/NrR); % period: annuli per RyR (for NrR < 1) ****
NrR = round(NrR);   % frequency: RyRs per annulus (for NrR > 1) ****
% Drflag = ( NrR>=1 ); % indicator of multiple RyRs per compartment ****
Drflag = ( NrR>1 ); % indicator of multiple RyRs per compartment ****
if Drflag % if there is more than 1 RyR per compartment
    disp('More than one RyR per axial compartment!');
else % if there is less than 1 RyR per compartment
    rng(nseed) % initiate random number generator
    % krR for a single receptor: note scaling for single circmf. cmprt.
    krR = Rscale * IrR * nc/vi * kvol; % units uM.s^-1 ****
    % assign an RyR to one of every MrR compartments; alternating circumf. 
    IRyR = false(ng,1); % indicator vector for axial RyR placement
    KRyR = false(ng,nc); % indicator array for axial+circumf. RyR placement
    Mr0 = max(floor(MrR/2),1)+1; %  axial index of first RyR
    for k = Mr0:MrR:ng
        IRyR(k,1) = true; % axial RyR location ****
        ncr = ceil( nc * rand(1) ); % choose axial placement randomly
        KRyR(k,ncr) = true;  % axial+circumferential RyR location ****
    end
    nRyR = sum(IRyR); % total number of RyRs
    % assemble axial indices evoked by KRyR, in order
    L=[];
    for k=1:ng
        L=[L,[1:ng]'];
    end
    LRyR = L(KRyR);
end

MrR0 = 5;   % MrR for baseline parameter set

%%% coefficients for contribution of input flux to boundary values %%%
%   Note: vol/area for cytosolic annuli = dx
BCc = dx/(JPb*DCac);
BInP3 = dx/(JPb*DInP3);

%% Write output

clear k LB R0 re ae ac volc vi vo kD LB L1 L0 dc R0 L

save(filename);
