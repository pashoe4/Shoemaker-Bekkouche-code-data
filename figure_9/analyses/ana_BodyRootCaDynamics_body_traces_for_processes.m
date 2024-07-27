function out=ana_BodyRootCaDynamics_body_traces_for_processes(results_name_base, var_str, val, fVar, tuned, elem_nrs,var_str_leg)
fig1=figure(1);
var_str_list_lab = var_str
if ~isempty(tuned)
    var_str_lab = [var_str, '_{tuned}'];
end

cd results
if tuned ~= []
        results_name=['exp_BodyRootCaDynamics_',var_str,'_tuned_',results_name_base];
    else
        results_name=['exp_BodyRootCaDynamics_',var_str,'_',results_name_base];
    end
cd(results_name)

legend_labs={};
for elem_nr=elem_nrs
    var_str_lab=var_str_list_lab;
    home=pwd;
    
    

    input_file_proc='parameters_proc';
    input_file_body='parameters_body';
    val=val;
    input_file_body_var=[input_file_body,'_', var_str, num2str(val),'.mat'];
    output=['OUT_',input_file_body_var];
    load(output);
    [nrows, ncolumns, nzdim] = size(Ccbt{1});
    Ccbt_vec=[];Ccbt_vec_start=[];
    for ii=1:length(Ccbt)
        Ccbt_vec=[Ccbt_vec Ccbt{ii}(elem_nr,6,5)];%
        Ccbt_vec_start=[Ccbt_vec_start Ccbt{ii}(1,6,5)];% 
    end
    
    if elem_nr == 1
        plot(data_time, Ccbt_vec,'LineStyle','--','LineWidth',1,'Color','black');hold on
    else
        plot(data_time, Ccbt_vec,'LineWidth',1);hold on
    end
%     elseif proc_nr == 4
%         plot(data_time, Ccbt_vec,'LineWidth',1,'Color',[0.8500, 0.3250, 0.0980]);hold on
%     elseif proc_nr == 2
%         plot(data_time, Ccbt_vec,'LineWidth',1,'Color',[0.4660, 0.6740, 0.1880]);hold on
%     elseif proc_nr == 3
%         plot(data_time, Ccbt_vec,'LineWidth',1,'Color',[0.6350, 0.0780, 0.1840]);hold on
%     end
    meas=[1,2,3,4];
    elem_nr_i = meas(elem_nr==[1,9,17,25]);
    %legend_labs{end+1}=['Process ',num2str(elem_nr)];
    legend_labs{end+1}=['Measurement ',num2str(elem_nr_i)];
    %val_mid=interval(ceil(end/2));
    %plot(fVar,delays,'marker','o','linewidth',1,'DisplayName',[var_str_lab,'=',num2str(val_mid)]);hold on
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig1
xlabel('Time (s)')
ylabel('Distal root [Ca^{2+}] (μM)')
legend(legend_labs,'Location','best')
ylim([0, 3.7])
xlim([0, 1.5])
grid on

%Save figure1
resolution=300;
output_size = [(1800) (1200)];%[(1800) (1200)];
set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
print(fig1, '-dpng', ['Ccbt_processes_' var_str_lab '.png'], ['-r' num2str(resolution)]);

% resolution=300;
% output_size = [(1800) (700)];%[(1800) (1200)];
% set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
% print(fig1, '-dpng', ['Ccbt_processes_' var_str_lab '.png'], ['-r' num2str(resolution)]);

cd(home)

disp("end")

