function out=exp_BodyRootCaDynamics_tuning(network_name, var_str, interval,var_location,min_start,max_start)

home=pwd;

if not(isfolder('results'))
    mkdir('results')
end
cd('results')
if not(isfolder(network_name))
    mkdir(network_name)
end
cd(network_name)

input_file_proc='parameters_proc';
input_file_body='parameters_body';
rateFunctionsSCS
ker_vals=[];
ce_vals=[];
precision = 1;
target_ce = 415;
target_ce_upper = target_ce + precision;
target_ce_lower = target_ce - precision;
for var_value = interval
    min=min_start;
    max=max_start;
    ker_val = min+(max-min)/2;
    for i=1:10%10 attempts
        disp('Attempt nr ' +string(i))
        input_file_body_var=[input_file_body,'_', var_str, num2str(var_value),'.mat'];
        if strcmp(var_location,'setParametersProc')
            setParametersProc(input_file_proc, var_str, var_value)
            setParametersBodyRoot(input_file_proc,input_file_body_var);
        elseif strcmp(var_location,'setParametersBodyRoot')
            setParametersProc(input_file_proc, var_str, var_value)
            setParametersBodyRoot(input_file_proc,input_file_body_var, var_str, var_value);
        else
            disp('Variable function undefined!')
        end
        cur_ce = BodyRootCaDynamics(input_file_body_var, 'tuning_only', true, 'ker', ker_val);
        disp(['ce = ', num2str(cur_ce)])
        if cur_ce > target_ce_upper %Half-interval search
            max = ker_val;
            ker_val=min+(ker_val-min)/2;
        elseif cur_ce < target_ce_lower
            min = ker_val;
            ker_val=ker_val+(max-ker_val)/2;
        else
            disp('Converged ker_val = ' + string(ker_val) + ', with ce value= ' +string(cur_ce) + ', for ' + var_str + '=' + string(var_value));
            converged=true;
            break
        end
        disp(['Trying ker_val = ', num2str(ker_val)])
    end
    ker_vals = [ker_vals ker_val];
    ce_vals = [ce_vals cur_ce];
end
disp(['ker_vals = ',num2str(ker_vals),' for ',var_str,' = ',num2str(interval)])
%ker_vals = -2.8516           0      3.2812 for RrIb = 5  10  15
ker_val_map = [interval',ker_vals',ce_vals'];
save('ker_val_map.mat','ker_val_map')

%disp('Made tuning attempts without success')
disp('End of tuning')
cd(home)
