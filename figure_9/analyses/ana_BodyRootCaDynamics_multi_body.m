function out=ana_BodyRootCaDynamics_multi_body(results_name_base, var_str_list, intervals,fVar,tuned,elem_nr,default_i,var_str_leg)
fig1=figure(1);
var_str_list_lab = var_str_list
if ~isempty(tuned)
    for i = tuned(1:end-1)
        var_str_list_lab{i} = [var_str_list_lab{i}, '_{tuned}'];
        %var_str_leg{i} = [var_str_leg{i}, '_{tuned}'];
    end
end
min_peak_val = 0.05;
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
        results_name=['exp_BodyRootCaDynamics_',var_str,'_tuned',results_name_base];
    else
        results_name=['exp_BodyRootCaDynamics_',var_str,results_name_base];
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
        [nrows, ncolumns, nzdim] = size(Ccbt{1});
        Ccbt_vec=[];Ccbt_vec_start=[];
        for ii=1:length(Ccbt)
            %time (ii)
            %circumferential element index (elem_nr,next to process 2,3,4)
            % process index to element index: 1(1), 2(9), 3(17), 4(25)
            %radial segment index (6, outermost)
            %axial segment index (5,middle)
            Ccbt_vec=[Ccbt_vec Ccbt{ii}(elem_nr,6,5)];%5 is middle
            Ccbt_vec_start=[Ccbt_vec_start Ccbt{ii}(1,6,5)];
        end

        %[max_val max_val_i]=max(Ccbt_vec);
        [max_vals,max_vals_i] = findpeaks(Ccbt_vec);
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

        [max_val_start max_val_i_start]=max(Ccbt_vec_start);
        
        

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
    if elem_nr == 25
        trav_dist = 2*6*pi/2;% 2*r*pi/2
    elseif elem_nr == 17
        trav_dist = 2*6*pi/3;% 2*r*pi/3
    elseif elem_nr == 9
        trav_dist = 2*6*pi/6;% 2*r*pi/6
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
%ylim([0 160]);
ylim([0 300]);
if elem_nr == 25
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
print(fig1, '-dpng', ['delays_body',num2str(elem_nr),'.png'], ['-r' num2str(resolution)]);

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
if elem_nr==17
    xwidth=500;
elseif elem_nr==25
    yticklabels(var_str_list_lab)
    xwidth=600;
elseif elem_nr==9
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
print(fig2, '-dpng', ['delay_mat_body',num2str(elem_nr),'.png'], ['-r' num2str(resolution)]);

