function out=exp_BodyRootCaDynamics(network_name, var_str, interval,var_location,kerval_path)

home=pwd;
cd results
if not(isfolder(network_name))
    mkdir(network_name)
end
cd(network_name)

varargin='';
ker_val_map=[];
if strcmp(kerval_path, '')
    set_tuned_vals=false;
else
    set_tuned_vals=true;
    load(kerval_path)
end

input_file_proc='parameters_proc';
input_file_body='parameters_body';
rateFunctionsSCS

if strcmp(var_str,'')
    input_file_body_var=[input_file_body,'.mat'];
    setParametersProc(input_file_proc)
    setParametersBodyRoot(input_file_proc,input_file_body_var);
    BodyRootCaDynamics(input_file_body_var, false);
else
    for var_value=interval
        input_file_body_var=[input_file_body,'_', var_str, num2str(var_value),'.mat'];
        if strcmp(var_location,'setParametersProc')
            setParametersProc(input_file_proc, var_str, var_value)
            setParametersBodyRoot(input_file_proc,input_file_body_var);
        elseif strcmp(var_location,'setParametersBodyRoot')
            setParametersProc(input_file_proc)
            setParametersBodyRoot(input_file_proc,input_file_body_var, var_str, var_value);
        else
            disp('Variable function undefined!')
        end
        if set_tuned_vals
            select = ker_val_map(:,1) == var_value;
            ker_value = ker_val_map(select,2);
            BodyRootCaDynamics(input_file_body_var, 'ker',ker_value, 'tend_sim',1.5);
        else
            BodyRootCaDynamics(input_file_body_var, 'tend_sim',1.5);
        end
    end
end
cd(home)

% rateFunctionsSCS
% setParametersProc
% setParametersBodyRoot
% BodyRootCaDynamics 