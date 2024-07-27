function out=ana_BodyRootCaDynamics(results_name, var_str, interval, kerval_path,var_str_leg)
home=pwd;
set_tuned_vals=true;
var_str_lab=var_str;
if ~strcmp(kerval_path, '')
    var_str_lab = [var_str,'_{tuned}'];
end

results_path = fullfile('results',results_name)
if not(isfolder(results_path))
    mkdir(results_path)
end
cd(results_path)

input_file_proc='parameters_proc';
input_file_body='parameters_body';
fig1=figure(1);
if strcmp(var_str,'')
    %Don't use this function perhaps
    input_file_body_var=[input_file_body,'.mat'];
    output=['OUT_',input_file_body_var];
    load(output);
    [nrows, ncolumns, nzdim] = size(Ccpt{1});
    Ccpt_vec=[];
    for i=1:length(Ccpt)
        %Ccpt_vec=[Ccpt_vec Ccpt{i}(round(nc/2),round(nr/2))];
        %Ccpt_vec=[Ccpt_vec Ccpt{i}(1,nr)];
        Ccpt_vec=[Ccpt_vec Ccpt{ii,proc_nr}(22,2)];%26
    end
    plot(data_time,Ccpt_vec);hold on
else
    for var_value=interval
        input_file_body_var=[input_file_body,'_', var_str, num2str(var_value),'.mat'];
        output=['OUT_',input_file_body_var];
        load(output);
        [nrows, ncolumns, nzdim] = size(Ccpt{1});
        Ccpt_vec=[];
        for i=1:length(Ccpt)
            %Ccpt_vec=[Ccpt_vec Ccpt{i}(round(nc/2),round(nr/2))];
            %Ccpt_vec=[Ccpt_vec Ccpt{i}(round(nrows/2),round(ncolumns/2),round(nzdim/2))];
            %Ccpt_vec=[Ccpt_vec Ccpt{i}(1,1,round(nzdim/2))];
            %Ccpt_vec=[Ccpt_vec Ccpt{i,4}(22,1)];%previosly Ccpt{i,4}(26,1)
            proc_nr=round(np/2+1);%Always the one on the opposite side of stimulus at proc_nr=1
            %[Ccpt_vec Ccpt{time,process_index}(process_segment_index,process_radial_depth)]
            Ccpt_vec=[Ccpt_vec Ccpt{i,proc_nr}(22,2)]; 
        end
        plot(data_time,Ccpt_vec,'DisplayName',[var_str_leg,'=',num2str(var_value)],'LineWidth',1);hold on
        grid on
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ylim([0, 3.5])%if all are plotted
if ~strcmp(var_str, 'np')
    ylim([0, 1.8])%For RrIb
    xlim([0, 1.5])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time (s)')
ylabel('Distal root [Ca^{2+}] (μM)')
legend('Location','best')

%Save figure
resolution=300;
output_size = [(1800) (1200)];%[(1800) (1200)];
set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
print(fig1, '-dpng', ['Ccpt_vs_' var_str_lab '.png'], ['-r' num2str(resolution)]);

% resolution=300;
% output_size = [(1800) (700)];%[(1800) (1200)];
% set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
% print(fig1, '-dpng', ['Ccpt_vs_' var_str_lab '.png'], ['-r' num2str(resolution)]);

cd(home)