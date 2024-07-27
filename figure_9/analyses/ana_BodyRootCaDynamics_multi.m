function out=ana_BodyRootCaDynamics_multi(results_name_base, var_str_list, intervals,fVar,tuned,proc_nr,default_i,var_str_leg)
fig1=figure(1);
var_str_list_lab = var_str_list
if ~isempty(tuned)
    for i = tuned(1:end-1)
        var_str_list_lab{i} = [var_str_list_lab{i}, '_{tuned}'];
        %var_str_leg{i} = [var_str_leg{i}, '_{tuned}'];
    end
end
min_peak_val = 0.04;
delay_mat = zeros(length(var_str_list),length(intervals{1}));
linestyle_list = {'-','-','-','-','-','-',':','--'};
marker_list = {'o','x','*','.','diamond','>','pentagram',"square"};
for i=1:length(var_str_list)
    var_str=var_str_list{i};
    %var_str_lab=var_str_list_lab{i};
    var_str_lab=var_str_leg{i};
    interval=intervals{i};
    home=pwd;
    cd results
    if sum(i==tuned)>0
        results_name=['exp_BodyRootCaDynamics_',var_str,'_tuned_',results_name_base];
    else
        results_name=['exp_BodyRootCaDynamics_',var_str,'_',results_name_base];
    end
    cd(results_name)
    input_file_proc='parameters_proc';
    input_file_body='parameters_body';
    delays=[];
    j=1;
    for var_value=interval
        input_file_body_var=[input_file_body,'_', var_str, num2str(var_value),'.mat'];
        output=['OUT_',input_file_body_var];
        load(output);
        [nrows, ncolumns, nzdim] = size(Ccpt{1});
        Ccpt_vec=[];Ccpt_vec_start=[];
        for ii=1:length(Ccpt)
            %[Ccpt_vec Ccpt{time,process_index}(process_segment_index,process_radial_depth)]
            %process_index, 6 processes
            %process_segment_index 1-26 but only use 2-25, (1 and 26 are
            %needed for MOLE to compute Laplacian
            %process_radial_depth, 3 levels
            Ccpt_vec=[Ccpt_vec Ccpt{ii,proc_nr}(22,2)];%26 
            Ccpt_vec_start=[Ccpt_vec_start Ccpt{ii,1}(22,2)];%26 
        end

        %[max_val max_val_i]=max(Ccpt_vec);
        [max_vals,max_vals_i] = findpeaks(Ccpt_vec);
        select = max_vals>min_peak_val;
        max_vals = max_vals(select);
        max_vals_i = max_vals_i(select);
        
        max_val = nan;
        max_val_i = nan;
        if length(max_vals)>1
            max_val = max_vals(1);
            max_val_i = max_vals_i(1);
        elseif length(max_vals)==1
            max_val = max_vals;
            max_val_i = max_vals_i;
        end

        [max_val_start max_val_i_start]=max(Ccpt_vec_start);
        
        

        if max_val_start>min_peak_val%Don't consider start wave unless alot above base level
            max_val_start_time=data_time(max_val_i_start);
        else
            max_val_start_time=nan;
        end
        
        if max_val>min_peak_val%Don't consider wave unless alot above base level
            max_val_time=data_time(max_val_i);
            delta_time = max_val_time-max_val_start_time;
            delays=[delays, delta_time];
        else
            delta_time=nan;
            delays=[delays, delta_time];
        end

        delay_mat(i,j) = delta_time;
        j=j+1;
    end
    %total length of roots (outer to cell body) = 2um
    %radius = r = 6um
    if proc_nr == 4
        trav_dist = 2*6*pi/2 + 2*2;% 2*r*pi/2 + 2*2um
    elseif proc_nr == 3
        trav_dist = 2*6*pi/3 + 2*2;% 2*r*pi/3 + 2*2um
    elseif proc_nr == 2
        trav_dist = 2*6*pi/6 + 2*2;% 2*r*pi/6 + 2*2um
    end
    velocity = trav_dist./delays;
    var_value_mid=interval(default_i);%interval(ceil(end/2));
    %plot(fVar,velocity,'marker',marker_list{i},'LineStyle',linestyle_list{i},'linewidth',1,'DisplayName',[var_str_lab,'=',num2str(var_value_mid)]);hold on
    if ~strcmp(var_str_lab,'R_{rI}')
        plot(fVar,velocity,'marker',marker_list{i},'LineStyle',linestyle_list{i},'linewidth',1,'DisplayName',[var_str_lab]);hold on
    end
    grid on
    cd(home)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig1
xlabel('Parameter Scaling (relative to baseline)')
ylabel('Ca^{2+} wave velocity (µm/s)')
%ylim([0 0.8]);
ylim([0 160]);
if proc_nr == 4
    %legend('Location','best');
    %leg = legend('Location','northwest');
    leg = legend('Position',[0.2 0.65 0.15 0.0869])
    
    %h = findobj('type', 'axes');  % Find all sets of axes
    %set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
    %set( leg ,'edgecolor','none' )%'fontsize' , 20 , 'location' , 'east', 'color' , 'g' , 'box' , 'on','edgecolor','g'
end

%Save figure1
resolution=300;
output_size = [(1800) (1200)];
set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
print(fig1, '-dpng', ['delays',num2str(proc_nr),'.png'], ['-r' num2str(resolution)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig2
fig2 = figure(2);
[M, N]=size(delay_mat)
x = repmat(1:N,M,1); % generate x-coordinates
y = repmat(1:M,N,1)'; % generate y-coordinates
% Generate Labels
delay_mat2 = round(delay_mat,2)
t = num2cell(delay_mat2); % extact values into cells
%t = cellfun(@(t)num2str, t, 'UniformOutput', false,'precision',2); % convert to string
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% Draw Image and Label Pixels
colormap parula%gray%jet%hot %pink, parula
imagesc(delay_mat)
h=gca; h.YAxis.TickLength = [0 0];
text(x(:), y(:), t, 'HorizontalAlignment', 'Center','FontWeight','bold','Color','white')
text(x(:), y(:), t, 'HorizontalAlignment', 'Center','FontWeight','normal','Color','black')
yticklabels([])
if proc_nr==3
    xwidth=500;
elseif proc_nr==4
    yticklabels(var_str_list_lab)
    xwidth=600;
elseif proc_nr==2
    hcb=colorbar;
    xwidth=650;
end
xwidth=xwidth+500;
xticks([1,2,3,4])
xticklabels(fVar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%caxis([0.48, 1.42]);
%caxis([0.13, 0.75]);
caxis([0.045, 0.765]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hcb.Title.String = "Delay (s)";
%hcb.Label.String = "asd"
%xlabel('Factor change from baseline (Legend: baseline [1] value)')

%Save figure2
resolution=300;
output_size = [(xwidth) (1200)];%[(1800) (1200)];
set(fig2,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
print(fig2, '-dpng', ['delay_mat',num2str(proc_nr),'.png'], ['-r' num2str(resolution)]);

