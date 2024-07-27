function out=ana_BodyCaDynamics(results_name, var_str, interval)
home=pwd;
cd results
if not(isfolder(results_name))
    mkdir(results_name)
end
cd(results_name)

input_file_proc='parameters_proc';
input_file_body='parameters_body';
fig1=figure(1);
for var_value=interval
    input_file_body_var=[input_file_body,'_', var_str, num2str(var_value),'.mat'];
    output=['OUT_',input_file_body_var];
    load(output);
    Ccbt_vec=[];
    for i=1:length(Ccbt)
        %Ccbt_vec=[Ccbt_vec Ccbt{i}(round(nc/2),round(nr/2))];
        Ccbt_vec=[Ccbt_vec Ccbt{i}(1,nr)];
    end
    plot(data_time,Ccbt_vec,'DisplayName',[var_str,'=',num2str(var_value)]);hold on
end
xlabel('Time (s)')
ylabel('Cell body Ca2+ conc. (μM)')
legend('Location','best')

%Save figure
resolution=300;
output_size = [(1800) (1200)];
set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
print(fig1, '-dpng', ['Ccbt_vs_' var_str '.png'], ['-r' num2str(resolution)]);

cd(home)