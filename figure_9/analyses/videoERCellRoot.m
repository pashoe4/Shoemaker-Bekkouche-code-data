%%%% videoERCellRoot.m %%%%
% plots ER calcium from output of a BodyCaDynamics simulation
% contains proximal portions of input processes
% creates a video file

clear;

%% load data for plotting
filename = input('Specify name of output file for plotting: ','s');
load(filename);

%% create the video writer
v = VideoWriter('video_ERCellBodyRoot.avi');
v.FrameRate = 10;

% open the video writer
open(v);

% define figure
figure(23);
ax = gca();

COLOR = colormap;

%% plot data

kdr = 5; % further decimation ratio for plotting

iend = length(Cebt); % # time steps

% plot parameters
dph = 2*pi()/np; % azimuthal angle between processes
ns = ng+ngr; % total number of grid elements / process
ds = ld/ns; % length of grid element
es = 0.8*ds; % distance to 'protrude' process into cell (for visualization)
rb0 = rb-es; % radius at which to start depicting process (visualization)
nex = round(0.5*ns); % fraction of processes to exclude from display

dth = 2*pi()/nc; % azimuthal angle subtended by body element
nos = ncp/2; % offset to align body with processes (visualization)

lwb = 10; % line width, body elements (for visualization)
lwr = 41 + 41*0.6/(ngr*1.5)*(1:ngr); % line widths, root elements (viz)
lwp = [ 41*ones(1,ng), lwr ]; % line width, process elts. incl. root (viz)

Cemax = 400; % maximum expected calcium value (for colormap scaling)

% loop over time (i.e., frames)
for it=1:iend
    
    if ~mod(it,kdr) % decimate data to be plotted
        
        % plot body calcium
        %%%   OPTION 1: Average ER calcium over z-dimension (thickness)
            Ceb = zeros(nc,nr);
            for kt = 1:nt
                Ceb = Ceb + Cebt{it}(:,:,kt);
            end
            Ceb = Ceb/nt;
        %%%   OPTION 2: Display ER calcium at a single z-value
%             kt = 5;
%             Ceb = Cebt{it}(:,:,kt);
        % compute indices into colormap
        iCeb = round( 255*Ceb/Cemax ) + 1;
        iCeb(iCeb>256) = 256; % just in case Cemax is exceeded
        % plot
        for kr=1:nr
            r = rc(kr);
            for kc = 1:nc
                th = [ (kc-nos-1.05)*dth, (kc-nos)*dth ];
                plot( r*sin(th), r*cos(th),...
                    'Color', COLOR(iCeb(kc,kr),:),...
                    'LineWidth', lwb );
                axis square
                axis([-13 13 -13 13]);
                hold on;
            end
        end
        
        % plot calcium in processes
        for kp=1:np % loop over processes
            % load ER Ca for this process
            iCep = Cept{it,kp}(2:ns+1);
            % compute indices into colormap 
            iCep = round( 255*iCep/Cemax ) + 1;
            iCep(iCep>256) = 256;  % just in case Cemax is exceeded
            ph = (kp-1)*dph;
            % plot
            for ks=nex+1:ns
                sk = ns-ks;
                r = [ rb0+sk*ds, rb0+(sk+1)*ds ];
                plot( r*sin(ph), r*cos(ph),...
                    'Color', COLOR(iCep(ks),:),...
                    'LineWidth', lwp(ks) );
                axis square
                axis([-13 13 -13 13]);
                hold on;
            end
        end
        
        hold off
        
        % For making video
        frame = getframe(ax);
        frame.cdata = imresize(frame.cdata,[600, 600]);
        writeVideo(v,frame);
        
    end
    
end

%% close the writer object

close(v);