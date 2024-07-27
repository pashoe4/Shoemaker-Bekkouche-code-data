function out=ana_BodyRootCaDynamics_body(results_name, var_str, interval, kerval_path,var_str_leg)
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
if strcmp(var_str,'')%
    %Don't use this if-statement (else is ok), should be updated or removed
    input_file_body_var=[input_file_body,'.mat'];
    output=['OUT_',input_file_body_var];
    load(output);
    [nrows, ncolumns, nzdim] = size(Ccbt{1});
    Ccbt_vec=[];
    for i=1:length(Ccbt)
        Ccbt_vec=[Ccbt_vec Ccbt{i}(elem_nr,6,5)];
    end
    plot(data_time,Ccbt_vec);hold on
else
    for var_value=interval
        input_file_body_var=[input_file_body,'_', var_str, num2str(var_value),'.mat'];
        output=['OUT_',input_file_body_var];
        load(output);
        [nrows, ncolumns, nzdim] = size(Ccbt{1});
        Ccbt_vec=[];
        for i=1:length(Ccbt)
            elem_nr=25;%round(np/2+1);
            Ccbt_vec=[Ccbt_vec Ccbt{i}(elem_nr,6,5)]; 
        end
        plot(data_time,Ccbt_vec,'DisplayName',[var_str_leg,'=',num2str(var_value),'μm^{-2}'],'LineWidth',1);hold on
        %plot(data_time,Ccbt_vec,'DisplayName',[var_str_leg],'LineWidth',1);hold on
        grid on
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ylim([0, 3.5])%if all are plotted
% if ~strcmp(var_str, 'np')
%     ylim([0, 1.8])%For RrIb
%     xlim([0, 1.5])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time (s)')
ylabel('Distal root [Ca^{2+}] (μM)')
legend('Location','best')

%Save figure
resolution=300;
output_size = [(1800) (1200)];%[(1800) (1200)];
set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
print(fig1, '-dpng', ['Ccbt_vs_' var_str_lab '.png'], ['-r' num2str(resolution)]);

% resolution=300;
% output_size = [(1800) (700)];%[(1800) (1200)];
% set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
% print(fig1, '-dpng', ['Ccbt_vs_' var_str_lab '.png'], ['-r' num2str(resolution)]);

cd(home)