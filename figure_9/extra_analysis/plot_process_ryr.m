clear all
close all

%ParamDepend_Graphs.xlsx
run("Speed_RyR.m")
run("Ampl_RyR.m")
run("Speed_noRyR.m")
run("Ampl_noRyR.m")

data={speed_ryr,ampl_ryr,speed_noryr,ampl_noryr};
data_lab = {'Speed_RyR','Ampl_RyR','Speed_noRyR','Ampl_noRyR'};
ylabs = {'Ca^{2+} wave velocity (µm/s)','Ca^{2+} wave amplitude (µM)','Ca^{2+} wave velocity (µm/s)','Ca^{2+} wave amplitude (µM)'};
ylims = {[90,250],[7.9,17],[0,120],[0.9,5.3]};
xlims = {[0,3],[0,3],[0,3.5],[0,3.5]};
legpos = {[0.68 0.76 0.15 0.0869],[0.68 0.76 0.15 0.0869],[0.68 0.3 0.15 0.0869],[0.68 0.3 0.15 0.0869]};
linestyle_list = {'-','-','-','-','-','-',':','--','-','-.','-.','-'};
marker_list = {'o','x','*','.','diamond','>','pentagram',"square",'<','|','pentagram','.'};
Ce0=400;
%rlx CeCtl, rle          fCac,fCae,fInP3, Rp_PM, RpS,  kbCf, rmIf,   kiPI,  rIg, Ce0 
baselines = [5.0E-16,   1.25E-18, 0.3, 0.6, 0.5,   500,   4000, 1.9,  1.4E-15,0.0943,70 , Ce0];
legend_labs = {'r_{lx,CeCtl}', 'r_{le}','D_{Cac}','D_{Cae}','D_{InP3}', 'R_{pP}', 'R_{pS}',  'k_{bCf}', 'r_{mIf}',   'k_{iPI}',  'r_{Ig}', 'Ce_{0}'};


for m = [1,2,3,4]
    fig1 = figure(m);
    cur_data = data{m};
    cur_lab = data_lab{m};
    for i = 1:2:24
        ii=(i+1)/2
        leg_label = legend_labs{ii};%speed_ryr.Properties.VariableNames{i};
        leg_label_file = cur_data.Properties.VariableNames{i};
        if strcmp(leg_label_file, 'Ce0')
            %plot(cur_data{:,i}, cur_data{:,i+1},'marker',marker_list{ii},'LineStyle',linestyle_list{(i+1)/2},'linewidth',2, 'DisplayName', [leg_label,'=',num2str(baselines(ii))],'Color',[128 128 128]/255);hold on
            plot(cur_data{:,i}, cur_data{:,i+1},'marker',marker_list{ii},'LineStyle',linestyle_list{(i+1)/2},'linewidth',2, 'DisplayName', [leg_label],'MarkerSize',12,'Color',[128 128 128]/255);hold on
        elseif ~strcmp(leg_label, 'kiPI')
            plot(cur_data{:,i}, cur_data{:,i+1},'marker',marker_list{ii},'LineStyle',linestyle_list{(i+1)/2},'linewidth',1, 'DisplayName', [leg_label]);hold on
        end
        
    end
    ylabel(ylabs{m});
    xlabel('Parameter Scaling (relative to baseline)');
    ylim(ylims{m});
    xlim(xlims{m});
    
    if m==1
        %legend('Location','best','NumColumns',2);
        %legend('Location','northeast','NumColumns',2);
        %legend('Position',[0.68 0.76 0.15 0.0869],'NumColumns',2)
        legend('Position',legpos{m},'NumColumns',2);
    end    
    
    grid on

    %Save figure1
    resolution=300;
    output_size = [(1800) (1200)];
    set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
    print(fig1, '-dpng', [cur_lab,'.png'], ['-r' num2str(resolution)]);
end