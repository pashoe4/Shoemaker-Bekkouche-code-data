function out=ana_ProcCaDynamics3D_ranRyR_multi(results_name_base, var_str, interval)
fig1=figure(1);
useRyRList=[0 1];
min_peak_val = 0.04;
delay_mat = zeros(2,length(interval));%with & without ryr
linestyle_list = {'-','-','-','-','-','-',':','--'};
marker_list = {'o','x','*','.','diamond','>','pentagram',"square"};
i=1;
for useRyR=useRyRList
    useRyRLab = num2str(useRyR);
    home=pwd;
    cd results
    cd(results_name_base)
    delays=[];
    j=1;
    for var_value=interval
        input_file_proc=['parameters_proc_useRyr',useRyRLab,'_', var_str, num2str(var_value)];
        output=['OUT_',input_file_proc];
        load(output);
        [nrows, ncolumns, nzdim] = size(Ccpt{1});
        Ccpt_vec=[];Ccpt_vec_start=[];
        for time_i=1:length(Ccpt)
            %One process 
            %process_segment_index (length) 1-202 (select close to first and last 3 to 200)
            %needed for MOLE to compute Laplacian
            %3 process_radial_depth levels (select nr2, middle)
            %6 elements in circumferential direction (select nr3, middle)
            
            Ccpt_vec=[Ccpt_vec Ccpt{time_i}(200,2,3)];
            Ccpt_vec_start=[Ccpt_vec_start Ccpt{time_i}(3,2,3)];
        end
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[max_val max_val_i]=max(Ccpt_vec);
        [max_vals,max_vals_i] = findpeaks(Ccpt_vec);

        %Don't consider start wave unless alot above base level
        select = max_vals>min_peak_val;
        max_vals = max_vals(select);
        max_vals_i = max_vals_i(select);
        
        %Make sure only the first is choosen if there are more than 1
        max_val = nan;
        max_val_i = nan;
        if length(max_vals)>1
            max_val = max_vals(1);
            max_val_i = max_vals_i(1);
        elseif length(max_vals)==1
            max_val = max_vals;
            max_val_i = max_vals_i;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[max_val_start max_val_i_start]=max(Ccpt_vec_start);
        [max_vals_start,max_vals_i_start] = findpeaks(Ccpt_vec_start);
        
        %Don't consider start wave unless alot above base level
        select = max_vals_start>min_peak_val;
        max_vals_start = max_vals_start(select);
        max_vals_i_start = max_vals_i_start(select);
        
        %Make sure only the first is choosen if there are more than 1
        max_val_start = nan;
        max_val_i_start = nan;
        if length(max_vals_start)>1
            max_val_start = max_vals_start(1);
            max_val_i_start = max_vals_i_start(1);
        elseif length(max_vals_start)==1
            max_val_start = max_vals_start;
            max_val_i_start = max_vals_i_start;
        end
        
        %A bit unnecessary code since it is already done, but I'll leave it
        %for now
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
    var_value_mid=interval(ceil(end/2));
    plot(interval,delays,'marker',marker_list{i},'LineStyle',linestyle_list{i},'linewidth',1,'DisplayName',['useRyr=',useRyRLab]);hold on
    cd(home)
    i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig1
xlabel(var_str)
ylabel('Wave travel time (s)')
%ylim([0 0.8]);

legend('Location','best');

%Save figure1
resolution=300;
output_size = [(1800) (1200)];
set(fig1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
print(fig1, '-dpng', ['results/',results_name_base,'/','delays_ryr_',var_str,'.png'], ['-r' num2str(resolution)]);



