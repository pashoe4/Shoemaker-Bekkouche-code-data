function out=exp_BodyCaDynamics(network_name, var_str, interval)
home=pwd;
cd results
if not(isfolder(network_name))
    mkdir(network_name)
end
cd(network_name)

input_file_proc='parameters_proc';
input_file_body='parameters_body';
rateFunctionsSCS
setParametersProc
for var_value=interval
    input_file_body_var=[input_file_body,'_', var_str, num2str(var_value),'.mat'];
    setParametersBody(input_file_proc,input_file_body_var,var_str,var_value);
    BodyCaDynamics(input_file_body_var);
end
cd(home)
