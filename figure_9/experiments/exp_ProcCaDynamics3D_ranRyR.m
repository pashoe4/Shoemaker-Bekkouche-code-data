function out=exp_ProcCaDynamics3D_ranRyR(network_name, var_str, interval,var_location,useRyr)

home=pwd;
cd results
if not(isfolder(network_name))
    mkdir(network_name)
end
cd(network_name)

input_file_proc=['parameters_proc_useRyr',num2str(useRyr)];
rateFunctionsSCS

if strcmp(var_str,'')
    input_file_proc_var=[input_file_proc,'.mat'];
    setParametersProc3D(input_file_proc,'useRyr')
    ProcCaDynamics3D_ranRyR(input_file_body_var);
else
    for var_value=interval
        input_file_proc_var=[input_file_proc,'_', var_str, num2str(var_value),'.mat'];
        if useRyr
            setParametersProc3D(input_file_proc_var, var_str, var_value);
        else
            setParametersProc3D(input_file_proc_var, var_str, var_value,'IrR',0);
        end

        ProcCaDynamics3D_ranRyR(input_file_proc_var);
    end
end
cd(home)

%Order to run original scripts
% rateFunctionsSCS
% setParametersProc
% ProcCaDynamics