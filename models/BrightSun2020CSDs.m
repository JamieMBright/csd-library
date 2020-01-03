function [csd, csd_tricomponent, ghics_optimised,dnics_optimised,difcs_optimised] = BrightSun2020CSDc(zen,ghi,ghics,dif,difcs,LST,plot_figure)
% The Bright-Sun clear-sky detection methodology, 2020.
% Cloudless sky CSDc
%
% ## References ##
%
% M. Alia-Martinez, J. Antonanzas, R. Urraca, F. J. Martinez-De-Pison, and
% F. Antonanzas-Torres, ‘Benchmark of algorithms for solar clear- sky
% detection’, J. Renew. Sustain. Energy, vol. 8, no. 3, 2016.
%
% Ellis, 2018, PVLIB detect_ghics at commit
% f4e4ad3bbe335fa49bf300ce53784d69d719ca98
% https://github.com/pvlib/pvlib-python/pull/596/commits
%
% Reno, M. J. and Hansen, C. W. 2016. Identification of periods of clear
% sky irradiance in time series of GHI measurements. Renewable Energy. 90,
% 520-531.
%
% Shen, Yu; Wei, Haikun; Zhu, Tingting; Zhao, Xin; Zhang, Kanjian. 2018. A
% Data-driven Clear Sky Model for Direct Normal Irradiance. IOP Conf.
% Series: Journal of Physics: Conf. Series. 1072.
%
% Inman, Rich H; Edson, James G and Coimbra, Carlos F M. 2015. Impact of
% local broadband turbidity estimation on forecasting of clear sky direct
% normal irradiance Solar Energy. 117, 125-138
%
% Gueymard, C.A; Bright, J.M.; Lingfors, D.; Habte, A. and Sengupta, M.
% 2019. A posteriori clear-sky identification methods in solar irradiance
% time series: Review and preliminary validation using sky imagers.
% Renewable and Sustainable Energy Reviews. In review.
%
% Larraneta, M; Reno, M J; Lillo-bravo, I; Silva-p, M A. 2017. Identifying
% periods of clear sky direct normal irradiance. Renewable Energy. 113,
% 756-763.
%
% Quesada-Ruiz, Samuel Linares-Rodriguez, Alvaro Ruiz-arias, José A
% Pozo-Vázquez, David. 2015. ScienceDirect An advanced ANN-based method to
% estimate hourly solar radiation from multi-spectral MSG imagery. Solar
% Energy. 115, 494-504.
%
% Zhang, Wenqi; Florita, Anthony R; Hodge, Bri-mathias; Mather, Barry.
% 2018. Modeling and Simulation of High Frequency Solar Irradiance.
% Preprint to Journal of Photovoltaics.
%
% Zhang, Wenqi Kleiber, William Florita, Anthony R Hodge, Bri-mathias
% Mather, Barry. 2018. A stochastic downscaling approach for generating
% high-frequency solar irradiance scenarios. Solar Energy. 176, 370-379.
%
%
% ------------------------------------------------------------------------
%                              INPUTS
% ------------------------------------------------------------------------
% ASSUMPTION: all data is in 1-min resolution.
% zen = Zenith angle in degrees. Corresponding to all inputs
%
% ghi = Global horizontal irradiance, column vector. Not necessary
%       to be continuous so long as the ghics perfectly corresponds.
%
% ghics = Clear-sky global horizontal irradiance, column vector.
%
% dif = Diffuse horizontal irradiance, column vector. Not necessary
%       to be continuous so long as the ghics perfectly corresponds.
%
% difcs = Clear-sky diffuse horizontal irradiance, column vector.
%
%
% LST = the local solar time, datevector. This can be calculated as
%       LST = datevec(datenum(UTC) + (lon*15)/24);
%
% plot_figure = if this variable is defined by any class, then a figure is
%               plotted illustrating the outcome of the CSD method, e.g. 1;
% ------------------------------------------------------------------------
%                              OUTPUTS
% ------------------------------------------------------------------------
% csd = clear-sky detection for cloudless sky. a flag of 1 means that
% clouds were suspected whereas a flag of 0 means that the hour is clear.
%
% csd_tricomponent = clear-sky detection from the reno tri component
% optimisation, this returns many many more suspected clear-sky periods
% without the durational filters etc.
%
% dnics/difcs_optimised - the resulting optimised clear-sky curves. Note
% that closure is not maintained or guaranteed with this approach, GHIcs
% should be recalculated to maintian closure.

%% safety checks
if ~iscolumn(zen)
    zen = zen';
    if ~iscolumn(zen)
        error('var zen must be a column vector')
    end
end
if ~iscolumn(ghi)
    ghi = ghi';
    if ~iscolumn(ghi)
        error('var ghi must be a column vector')
    end
end
if ~iscolumn(ghics)
    ghics = ghics';
    if ~iscolumn(ghics)
        error('var ghi must be a column vector')
    end
end
if ~iscolumn(dif)
    dif = dif';
    if ~iscolumn(dif)
        error('var dif must be a column vector')
    end
end
if ~iscolumn(difcs)
    difcs = difcs';
    if ~iscolumn(difcs)
        error('var difcscs must be a column vector')
    end
end
if length(unique([length(ghi),length(ghics),length(zen),length(dif),length(difcs),length(datenum(LST))]))~=1
    error('vars must be equal in length')
end

% turn off hankel diagonal conflict warning. This happens if the last
% element in ghi/dni/dif is NaN.
warning('off','MATLAB:hankel:AntiDiagonalConflict')

%% calculate DNI
% we could also request DNI as an input.

% make sure dif is not greater than ghi
dif(dif>ghi)=ghi(dif>ghi);
difcs(difcs>ghics)=ghics(difcs>ghics);

% calculate dni and dnics from closure
dnics = (ghics-difcs)./cosd(zen);
dni = (ghi-dif)./cosd(zen);

%% Clean the data for NaN values
% find all the indices where there are missing data of ghi, dni or dif.
idxs_nan = (isnan(ghi) | isnan(dni) | isnan(dif));
% find the corresponding index of where the data is not nan.
not_nan_indxs = (~isnan(ghi) & ~isnan(dni) & ~isnan(dif));
% find the length of the required ouput for sanity.
n_output = length(ghi);

% Now delete all cases of missing data. We cannot perform on this.
ghi(idxs_nan)=[];
dni(idxs_nan)=[];
dif(idxs_nan)=[];
ghics(idxs_nan)=[];
dnics(idxs_nan)=[];
difcs(idxs_nan)=[];
zen(idxs_nan)=[];
LST(idxs_nan,:)=[];
% what remains is perfectly continuous time series of all the variables.
% The tests do not matter should time skip because all tests are
% comparative to the clear-sky irradiance and they too also skip. 

%% Parameterisation
% The methodology applies the Reno 2016 methodology to each sub-component
% of irradiance. There have been many proposed parametrisations in
% literature---all contained here and commented out if not used.
% This methodology mixes and matches the different parametrisations and
% also proposes new limits to certain tests. 

% Reno 2016 Parameterisations for Optimisation
reno_ghi_mean_diff=75;
reno_ghi_max_diff=75;
reno_ghi_lower_line_length=-5;
reno_ghi_upper_line_length=10;
reno_ghi_var_diff=0.1;% old reno 0.005;
reno_ghi_slope_dev=8;

% optimisation lower irradiance limit threshold (Wm-2)
opt_thres = 30; 

% Bright 2019 parametrisations for modified-Reno
% GHI
BrightSunCriteria.ghi_c1_big_zenith = 0.5;
BrightSunCriteria.ghi_c1_mean_diff=0.125;
BrightSunCriteria.ghi_c1_small_zenith = 0.25;

BrightSunCriteria.ghi_c2_big_zenith = 0.5;
BrightSunCriteria.ghi_c2_max_diff=0.125;
BrightSunCriteria.ghi_c2_small_zenith = 0.25;

BrightSunCriteria.ghi_c3_lower_line_length=-0.5;
BrightSunCriteria.ghi_c3_small_zenith = -7;

BrightSunCriteria.ghi_c4_var_diff=0.4; 

BrightSunCriteria.ghi_c5_slope_dev=15; 
BrightSunCriteria.ghi_c5_small_zenith = BrightSunCriteria.ghi_c5_slope_dev*3;

%DIF
BrightSunCriteria.dif_c1_big_zenith = 0.5;
BrightSunCriteria.dif_c1_mean_diff=0.5;
BrightSunCriteria.dif_c1_small_zenith = 0.25;

BrightSunCriteria.dif_c2_big_zenith = 0.5;
BrightSunCriteria.dif_c2_max_diff=0.5;
BrightSunCriteria.dif_c2_small_zenith = 0.25;

BrightSunCriteria.dif_c3_lower_line_length=-1.7;
BrightSunCriteria.dif_c3_small_zenith = -6;

BrightSunCriteria.dif_c4_var_diff=0.2;

BrightSunCriteria.dif_c5_slope_dev=8;
BrightSunCriteria.dif_c5_small_zenith = BrightSunCriteria.dif_c5_slope_dev*3;


% set parametrisation for the turnaround in the zenith tuned parameters
% At this zenith parameter, the criteria 1,2,4,5 for Bright-Sun will be
% exactly as set.
zenith_turn_around = 30;

% the window length is as standard across all the Reno variants at 10 mins
window_length = 10;
% the max iterations comes from Ellis 2018 whereby optimisations are
% acheived up until 20 attempts or convergence is reached.
max_iterations=20;

% assert a minimum and maximum allowable modification to the clear-sky
% irradiance curve, where 1 = no change. 
upper_alpha_limit = 1.5;
lower_alpha_limit = 0.7;


% Duration filter strictness
% First durational window. Higher time but lower strictness.
% The window is the duration size of the window in mins. 
window_1st_filter = 1;
% Tolerance is the permissible number of CSD periods within the window. 
tolerance_1st_filter = 1;

% Second durational window. Lower time but higher strictness.
window_2nd_filter = 1;
tolerance_2nd_filter = 1;

%% Initial CSD guess using Reno for Optimisation Only.
%%%
%%% GHI
%%%

% produce Hankel matrices for GHI and GHIcs as defined by window_length. a
% Hankel matrix is essentially a staggered set of time series where each
% column is the same data offset by 1 time step, this means that each row
% has a "window" length time series, and each column has the entire time
% series offset by the column number.
ghi_window = hankel(ghi,[ghi(end),NaN(1,window_length-1)]);
ghics_window = hankel(ghics,[ghics(end),NaN(1,window_length-1)]);

% Calculate the statistics for use in the criteria. This approach is lifted
% from Ellis 2018, Inman 2015 and Reno 2012/2016.
% calculate measurement statistics for ghi.
meas_mean = nanmean(ghi_window,2); % DIM=2, so the mean is taken on a row-wise basis
% [M,NDX] = NANMAX(A,[],DIM) operates along the dimension DIM.
meas_max = nanmax(ghi_window,[],2); % the [] is just matlab specific operations.
% diff reduces the matrix size by 1 due to delta between. The resultant is
% increaasaed by a row of NaNs representing that there is no t+1 for the
% final time step.
meas_diff = [diff(ghi_window);NaN(1,window_length)];
meas_slope = [diff(ghi_window);NaN(1,window_length)];
% Y = STD(X,W,DIM), we want DIM=2.
% Pass in 0 for W to use the default normalization by N-1, or 1 to use N.
% use the omitnan option to ignore NaN values, this automatically updates N
meas_slope_nstd = std(meas_slope,0,2,'omitnan') ./ meas_mean;
meas_line_length = nansum(sqrt(meas_diff .* meas_diff + 1^2),2);

% Repeat the above for GHIcs statistics
clear_mean = nanmean(ghics_window,2);
clear_max = nanmax(ghics_window,[],2);
clear_diff = [diff(ghics_window);NaN(1,window_length)];
clear_slope = [diff(ghics_window);NaN(1,window_length)];
clear_line_length = nansum(sqrt(1^2 + clear_diff.^2),2);

% the final statistic is the line difference that requires both GHIcs and
% GHI statistics combined.
line_diff = meas_line_length - clear_line_length;

% initialise the criteria matrices where 1 = cloud. The default assumption
% is that cloud is present, satisfaction of all criteria results in an
% assertion of 0 for that time step.
c1 = ones(size(ghi));
c2 = ones(size(ghi));
c3 = ones(size(ghi));
c4 = ones(size(ghi));
c5 = ones(size(ghi));
c6 = ones(size(ghi));

% perform criteria. We use the Ellis parameterisation for GHI.
% these criteria are very well documented in all Reno, Ellis and Inman
% papers. Particularly in the figure 5 of Inman 2015! 

% if clear>meas, this is ok
c1(abs(meas_mean - clear_mean) < reno_ghi_mean_diff)=0;
c2(abs(meas_max - clear_max) < reno_ghi_max_diff)=0;
c3((line_diff > reno_ghi_lower_line_length) & (line_diff < reno_ghi_upper_line_length))=0;
c4(meas_slope_nstd < reno_ghi_var_diff)=0;
c5(nanmax(abs(meas_slope -  clear_slope),[],2) < reno_ghi_slope_dev)=0;
% this 6th criteria only exists in the Ellis 2018 coded version in PVlib.
% It is a sensibility check for NaN results, though probably redundant due
% to the prior screening.
c6((clear_mean ~= 0) & ~isnan(clear_mean))=0;

% should any criteria be 1(cloud), then the time-step is deemed cloudy
csd_initial = c1 + c2 + c3 + c4 + c5 + c6;
% csd is now the sum of failed tests and so should all 6 criteria have
% failed, csd_1=6. A value of 0 is the only case where clear sky is
% detected. We are unconcerned with how many tests failed, so long as 1 did
% fail, it can then be represented as a binary and so we set all failures
% to be 1.
csd_ghiz = zeros(size(ghi));
csd_ghiz(csd_initial>0)=1;
% initialise the first CSD guess.
% zeros indicate that the period is clear. 1s indicate cloud
% if all three components do not corroborate, then we assume the period was
% cloudy
csd_initial = csd_ghiz;

if exist('plot_figure','var')
    % % plot an example of the first guess CSD.
    lst_plot = datetime(datevec(LST));
    figure('name','Initial CSD from Reno 2016 example','color','w');
    hold on
    plot(lst_plot,ghi)
    plot(lst_plot,dni)
    plot(lst_plot,dif)
    ghiCSD = ghi; ghiCSD(csd_initial==1)=NaN;
    plot(lst_plot,ghiCSD,'linewidth',2,'color','k')
    plot(lst_plot,ghics,'k:')
    dniCSD = dni; dniCSD(csd_initial==1)=NaN;
    difCSD = dif; difCSD(csd_initial==1)=NaN;
    plot(lst_plot,dniCSD,'linewidth',2,'color','k')
    plot(lst_plot,difCSD,'linewidth',2,'color','k')
    plot(lst_plot,dnics,'k:')
    plot(lst_plot,difcs,'k:')
    ylim([1,1300])
    hold off
    legend({'GHI','DNI','DIF','CSD','Clear-sky irradiance'},'fontname','times')
    ylabel(gca,'Irradiance [Wm$^{-2}$]','Interpreter','latex')
    % set(gca,'xtick',[]);
    xlabel('Time','Interpreter','latex')
    set(gca,'FontName','times')
end

%% Find all the unique days
% find the unique days so that optimisation can be done on a daily basis.
dtn = floor(datenum(LST));
unique_days = unique(dtn);

% initial clear-sky optimisation is that the clear-sky irradiance that is
% input is already perfect (e.g. alpha=1), with each iteration and
% optimisation, we will adjust 1 based on those periods identified as clear
% by the modified Zhang. 
alpha_ghi = ones(size(unique_days));
alpha_dif = ones(size(unique_days));
alpha_dni = ones(size(unique_days));

%% Optimisation of clear-sky irradiance
% This is a similar principal as employed by Alia-Martinez 2016 and Ellis
% 2018.

if exist('plot_figure','var')
    old_ghics = ghics;
    old_dnics = dnics;
    old_difcs = difcs;
end

% First, we look through every single day identified in the input data
for d = 1:length(unique_days)
    % find the indices within the data that correspond to this day
    idxs = find(dtn == unique_days(d));
    
    % isolate the variables to only this day
    ghid = ghi(idxs);
    ghicsd = ghics(idxs);
    difd = dif(idxs);
    difcsd = difcs(idxs);
    dnid = dni(idxs);
    dnicsd = dnics(idxs);
    csdd = csd_initial(idxs);
    
    % apply the optimisation lower limit threshold
    csdd(difd<opt_thres | ghid<opt_thres)=1;
    
    % isolate clear measured ghi and corresponding clear-sky ghi from within this day.
    test_ghi = ghid(csdd==0);
    test_ghics = ghicsd(csdd==0);
    test_dif = difd(csdd==0);
    test_difcs = difcsd(csdd==0);
    test_dni = dnid(csdd==0);
    test_dnics = dnicsd(csdd==0);
    
    % define the rmse functions for optimisation strategy.
    rmse_ghi = @(a) sqrt(nanmean((test_ghi - a*test_ghics).^2));
    rmse_dif = @(a) sqrt(nanmean((test_dif - a*test_difcs).^2));
    rmse_dni = @(a) sqrt(nanmean((test_dni - a*test_dnics).^2));
    
    % if there were at least 60 clear periods to optimise with, then we can
    % proceed. Fewer sites may lead to really strange optimisation, which
    % is undesireable. 60 seems a solid amount for trustworthy
    % optimsiations.
    if sum(csdd==0)>60
        % Begining of optimisation
        
        % initialise the "current_alpha", which is always 1.
        current_alpha = alpha_ghi(d);
        % set the previous_alpha as NaN so that it cannot pass the "if"
        % statement featured in the following "while" loop.
        previous_alpha = NaN;
        % initialise the iteration count
        iter = 0;
        
        % enter the while loop. This means that the process will continue
        % until the condition is satisfied. The condition is satisfied if
        % the optimisation has converged within 0.00001 or if 20 iterations
        % has occured.
        while (iter<max_iterations && round(current_alpha*10000) ~= round(previous_alpha*10000))
            % run a Multidimensional unconstrained nonlinear minimization
            % (Nelder-Mead) search of the function=rmse minimising through alpha.
            previous_alpha = current_alpha;
            current_alpha = fminsearch(rmse_ghi,current_alpha);
            
            % update the iteration count
            iter = iter + 1;
        end
        alpha_ghi(d) = current_alpha;
        
        % repeat process for DIF
        current_alpha = alpha_dif(d);
        previous_alpha = NaN;
        iter = 0;
        while (iter<max_iterations && round(current_alpha*10000) ~= round(previous_alpha*10000))
            previous_alpha = current_alpha;
            current_alpha = fminsearch(rmse_dif,current_alpha);
            iter = iter + 1;
        end
        alpha_dif(d) = current_alpha;
        
        % repeat process for DNI
        current_alpha = alpha_dni(d);
        previous_alpha = NaN;
        iter = 0;
        while (iter<max_iterations && round(current_alpha*10000) ~= round(previous_alpha*10000))
            previous_alpha = current_alpha;
            current_alpha = fminsearch(rmse_dni,current_alpha);
            iter = iter + 1;
        end
        alpha_dni(d) = current_alpha;
    end
    
    % Occasionally, the Reno fails at identifying clear periods (such as in
    % extreme latitudes in the antarctic), as such, we assert some upper
    % limits to the possible alpha correction.

    %ghi
    temp = alpha_ghi(d);
    temp(temp>upper_alpha_limit) = upper_alpha_limit;
    temp(temp<lower_alpha_limit) = lower_alpha_limit;
    alpha_ghi(d)=temp;
    %dni
    temp = alpha_dni(d);
    temp(temp>upper_alpha_limit) = upper_alpha_limit;
    temp(temp<lower_alpha_limit) = lower_alpha_limit;
    alpha_dni(d)=temp;
    %dif
    temp = alpha_dif(d);
    temp(temp>upper_alpha_limit) = upper_alpha_limit;
    temp(temp<lower_alpha_limit) = lower_alpha_limit;
    alpha_dif(d)=temp;
        
    % Apply the clear-sky correction factors to the clear-sky estimates for
    % this day. Should the estimate already be ideal, alpha=1 and no change
    % will occur.
    ghics(idxs) = ghics(idxs).*alpha_ghi(d);
    difcs(idxs) = difcs(idxs).*alpha_dif(d);
    dnics(idxs) = dnics(idxs).*alpha_dni(d);
end

% plot example of the optimisation working
if exist('plot_figure','var')
    lst_plot = datetime(datevec(LST));
    figure('Name','Example of clear-sky optimisation','color','w')
    plot(lst_plot,ghi)
    hold on
    plot(lst_plot,dni)
    plot(lst_plot,dif)
    plot(lst_plot,ghiCSD,'linewidth',2,'color','k')
    plot(lst_plot,old_ghics,'b:','linewidth',2)
    plot(lst_plot,ghics,'r:','linewidth',2)
    
    plot(lst_plot,dniCSD,'linewidth',2,'color','k')
    plot(lst_plot,difCSD,'linewidth',2,'color','k')   
    plot(lst_plot,old_dnics,'b:','linewidth',2)
    plot(lst_plot,old_difcs,'b:','linewidth',2)
    plot(lst_plot,dnics,'r:','linewidth',2)
    plot(lst_plot,difcs,'r:','linewidth',2)
    hold off
    legend({'GHI','DIF','DNI','Reno CSD','Clear-sky original','Clear-sky optimised'},'FontName','times')
    set(gca,'FontName','Times')
    xlabel('Time','FontName','Times')
    ylabel('Irradiance [Wm$^{-2}$]','FontName','Times','Interpreter','latex')
end

%% perform tri-component analysis.
% tri component analysis is a fancy way of saying we perform criteria on
% all irradiance components (global, direct and diffuse). Each subcomponent
% undergoes its own CSD methodology and only if all component CSDs
% corroborate that it is infact cloud free does the tri-component CSD time
% series register as clear. 


%%%
%%% GHI
%%%

% produce Hankel matrices for GHI and GHIcs as defined by window_length. a
% Hankel matrix is essentially a staggered set of time series where each
% column is the same data offset by 1 time step, this means that each row
% has a "window" length time series, and each column has the entire time
% series offset by the column number.
ghi_window = hankel(ghi,[ghi(end),NaN(1,window_length-1)]);
ghics_window = hankel(ghics,[ghics(end),NaN(1,window_length-1)]);

% Calculate the statistics for use in the criteria. This approach is lifted
% from Ellis 2018, Inman 2015 and Reno 2012/2016.
% calculate measurement statistics for ghi.
meas_mean = nanmean(ghi_window,2); % DIM=2, so the mean is taken on a row-wise basis
% [M,NDX] = NANMAX(A,[],DIM) operates along the dimension DIM.
meas_max = nanmax(ghi_window,[],2); % the [] is just matlab specific operations.
% diff reduces the matrix size by 1 due to delta between. The resultant is
% increaasaed by a row of NaNs representing that there is no t+1 for the
% final time step.
meas_diff = [diff(ghi_window);NaN(1,window_length)];
meas_slope = [diff(ghi_window);NaN(1,window_length)];
% Y = STD(X,W,DIM), we want DIM=2.
% Pass in 0 for W to use the default normalization by N-1, or 1 to use N.
% use the omitnan option to ignore NaN values, this automatically updates N
meas_slope_nstd = std(meas_slope,0,2,'omitnan') ./ meas_mean;
meas_line_length = nansum(sqrt(meas_diff .* meas_diff + 1^2),2);

% Repeat the above for GHIcs statistics
clear_mean = nanmean(ghics_window,2);
clear_max = nanmax(ghics_window,[],2);
clear_diff = [diff(ghics_window);NaN(1,window_length)];
clear_slope = [diff(ghics_window);NaN(1,window_length)];
clear_line_length = nansum(sqrt(1^2 + clear_diff.^2),2);

% the final statistic is the line difference that requires both GHIcs and
% GHI statistics combined.
line_diff = (meas_line_length - clear_line_length)./clear_line_length;

% initialise the criteria matrices where 1 = cloud. The default assumption
% is that cloud is present, satisfaction of all criteria results in an
% assertion of 0 for that time step.
c1 = ones(size(ghi));
c2 = ones(size(ghi));
c3 = ones(size(ghi));
c4 = ones(size(ghi));
c5 = ones(size(ghi));
c6 = ones(size(ghi));

% We also introduce the zenith flexibility introduced by Larraneta 2017.
% We observe in many occasions events where high zenith angles behave
% smoothly as if clear sky, however, are often lower than the clear sky
% curve (particularly if not a very good clearsky model. For that reason,
% we must apply some flexibility with zenith angle.
z = (20:0.01:90)';
% produce a linearly spaced correction factor as proposed by Larraneta
% 2017. Note that they used three fixed bins whereas  we smooth this
% through interpolation

c1_lim = flip([linspace(BrightSunCriteria.ghi_c1_big_zenith,BrightSunCriteria.ghi_c1_mean_diff,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.ghi_c1_mean_diff,BrightSunCriteria.ghi_c1_small_zenith,length(z)-find(z==zenith_turn_around))']);
% find the indices where the zenith angles correspond.
inds = knnsearch(z,zen);
% find the kc_limimts for those associated time steps
c1_lim = c1_lim(inds);

% same principle for low and high zenith.
c2_lim = flip([linspace(BrightSunCriteria.ghi_c2_big_zenith,BrightSunCriteria.ghi_c2_max_diff,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.ghi_c2_max_diff,BrightSunCriteria.ghi_c2_small_zenith,length(z)-find(z==zenith_turn_around))']);
c2_lim = c2_lim(inds);

% apply a relxation for very low zenith to c3 and 5.
c3_lim_lower = flip([linspace(BrightSunCriteria.ghi_c3_lower_line_length,BrightSunCriteria.ghi_c3_lower_line_length,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.ghi_c3_lower_line_length,BrightSunCriteria.ghi_c3_small_zenith,length(z)-find(z==zenith_turn_around))']);
c3_lim_lower = c3_lim_lower(inds);
c3_lim_upper = abs(c3_lim_lower);

c5_lim_lower = flip([linspace(BrightSunCriteria.ghi_c5_slope_dev,BrightSunCriteria.ghi_c5_slope_dev,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.ghi_c5_slope_dev,BrightSunCriteria.ghi_c5_small_zenith,length(z)-find(z==zenith_turn_around))']);
c5_lim_lower = c5_lim_lower(inds);

% perform criteria. We use the Ellis parameterisation for GHI.
% these criteria are very well documented in all Reno, Ellis and Inman
% papers. Particularly in the figure 5 of Inman 2015! 

% if clear>meas, this is ok
c1(abs(meas_mean - clear_mean)./clear_mean < c1_lim)=0;
c2(abs(meas_max - clear_max)./clear_max < c2_lim)=0;
c3((line_diff > c3_lim_lower) & (line_diff < c3_lim_upper))=0;
c4(meas_slope_nstd < BrightSunCriteria.ghi_c4_var_diff)=0;
c5(nanmax(abs(meas_slope -  clear_slope),[],2) < c5_lim_lower)=0;
% this 6th criteria only exists in the Ellis 2018 coded version in PVlib.
% It is a sensibility check for NaN results, though probably redundant due
% to the prior screening.
c6((clear_mean ~= 0) & ~isnan(clear_mean))=0;

% should any criteria be 1(cloud), then the time-step is deemed cloudy
csd_1 = c1 + c2 + c3 + c4 + c5 + c6;
% csd is now the sum of failed tests and so should all 6 criteria have
% failed, csd_1=6. A value of 0 is the only case where clear sky is
% detected. We are unconcerned with how many tests failed, so long as 1 did
% fail, it can then be represented as a binary and so we set all failures
% to be 1.
csd_ghi = zeros(size(ghi));
csd_ghi(csd_1>0)=1;

% Keep diagnostics for demonstration plots
criteria_ghi = [c1,c2,c3,c4,c5,c6];

% The same process is repeated for DIF
%%%
%%% DIF
%%%
dif_window = hankel(dif,[dif(end),NaN(1,window_length-1)]);
difcs_window = hankel(difcs,[difcs(end),NaN(1,window_length-1)]);
meas_mean = nanmean(dif_window,2); % DIM=2.
meas_max = nanmax(dif_window,[],2);
meas_diff = [diff(dif_window);NaN(1,window_length)];
meas_slope = [diff(dif_window);NaN(1,window_length)];
meas_slope_nstd = std(meas_slope,0,2,'omitnan') ./ meas_mean;
meas_line_length = nansum(sqrt(meas_diff .* meas_diff + 1),2);

clear_mean = nanmean(difcs_window,2);
clear_max = nanmax(difcs_window,[],2);
clear_diff = [diff(difcs_window);NaN(1,window_length)];
clear_slope = [diff(difcs_window);NaN(1,window_length)];
clear_line_length = nansum(sqrt(1^2 + clear_diff.^2),2);

line_diff = (meas_line_length - clear_line_length)./clear_line_length;

c1 = ones(size(dif));
c2 = ones(size(dif));
c3 = ones(size(dif));
c4 = ones(size(dif));
c5 = ones(size(dif));
c6 = ones(size(dif));

z = (20:0.01:90)';
c1_lim = flip([linspace(BrightSunCriteria.dif_c1_big_zenith,BrightSunCriteria.dif_c1_mean_diff,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.dif_c1_mean_diff,BrightSunCriteria.dif_c1_small_zenith,length(z)-find(z==zenith_turn_around))']);
inds = knnsearch(z,zen);
c1_lim = c1_lim(inds);
c2_lim = flip([linspace(BrightSunCriteria.dif_c2_big_zenith,BrightSunCriteria.dif_c2_max_diff,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.dif_c2_max_diff,BrightSunCriteria.dif_c2_small_zenith,length(z)-find(z==zenith_turn_around))']);
c2_lim = c2_lim(inds);
c3_lim_lower = flip([linspace(BrightSunCriteria.dif_c3_lower_line_length,BrightSunCriteria.dif_c3_lower_line_length,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.dif_c3_lower_line_length,BrightSunCriteria.dif_c3_small_zenith,length(z)-find(z==zenith_turn_around))']);
c3_lim_lower = c3_lim_lower(inds);
c3_lim_upper = abs(c3_lim_lower);
c5_lim_lower = flip([linspace(BrightSunCriteria.dif_c5_slope_dev,BrightSunCriteria.dif_c5_slope_dev,find(z==zenith_turn_around))';...
    linspace(BrightSunCriteria.dif_c5_slope_dev,BrightSunCriteria.dif_c5_small_zenith,length(z)-find(z==zenith_turn_around))']);
c5_lim_lower = c5_lim_lower(inds);


% note that the BrightSunCriteria.dif parameters are now used.
c1(abs(meas_mean - clear_mean)./clear_mean < c1_lim)=0;
c2(abs(meas_max - clear_max)./clear_max < c2_lim)=0;
c3((line_diff > c3_lim_lower) & (line_diff < c3_lim_upper))=0;
c4(meas_slope_nstd < BrightSunCriteria.dif_c4_var_diff)=0;
c5(nanmax(abs(meas_slope -  clear_slope),[],2) < c5_lim_lower)=0;
c6((clear_mean ~= 0) & ~isnan(clear_mean))=0;

csd_2 = c1 + c2 + c3 + c4 + c5 + c6;
csd_dif = zeros(size(ghi));
csd_dif(csd_2>0)=1;

% Keep diagnostics for demonstration plots
criteria_dif = [c1,c2,c3,c4,c5,c6];

%%%
%%% DNI
%%%
% Apply the Quesada-Ruiz 2015 methodology, though due to optimisation of
% the clear-sky curves, we can afford to be more conservative and apply
% 0.95 instead. 
% We also introduce the zenith flexibility introduced by Larraneta 2017.
% create zenith bins from 30:85 at a 0.01 interval size
z = (30:0.01:90)';
% produce a linearly spaced correction factor as proposed by Larraneta
% 2017. Note that they used three fixed bins, we smooth this through
% interpolation
kc_lims = flip(linspace(0.5,0.9,length(z)))';
% find the indices where the zenith angles correspond.
inds = knnsearch(z,zen);
% find the kc_limimts for those associated time steps
kc_lim = kc_lims(inds);

% find the clear sky index for DNI/beam (kcb)
kcb = dni./dnics;
% pre allocate the DNI CSD time series with assumption of cloudy=1 as
% initial guess.
csd_dni=ones(size(dni));
% should the clear-sky beam index be below the limt, we assert that this
% DNI period is clear---with a 0.
csd_dni(kcb>kc_lim)=0;

% keep diagnostics for dni
criteria_dni = csd_dni;

%% Overall criteria
% pre allocate assuming all cloudy
csd_overall = ones(size(ghi));
% if all componentats are clear, then clear.
csd_overall(csd_ghi+csd_dif+csd_dni==0)=0;

if exist('plot_figure','var')
    % % Exploration plot of the tri-component analysis. Zoom in to have a look.
    f = figure('name','Demonstrate that the tri-componenet is working','color','w');
    f.Position(3:4) = [1000,250];
    hold on
    plot(ghi)
    plot(dni)
    plot(dif)
    ghiCSD = ghi; ghiCSD(csd_ghi==1)=NaN;
    dniCSD = dni; dniCSD(csd_dni==1)=NaN;
    difCSD = dif; difCSD(csd_dif==1)=NaN;
    allCSD = ghi; allCSD(csd_overall==1)=NaN;
    plot(ghiCSD,'linewidth',4,'color','k')
    plot(ghics,'k:')
    plo t(allCSD,'LineWidth',2,'Color','g')
    plot(dnics,'k:')
    plot(difcs,'k:')
    plot(dniCSD,'linewidth',2,'color','k')
    plot(difCSD,'linewidth',2,'color','k')
    hold off
    leg = legend('GHI','DNI','DIF','CSD','CS','combined CSD');
    ylabel(gca,'Irradiance [Wm$^{-2}$]','FontName','Times','Interpreter','latex')
    set(gca,'xtick',[]);
    ylim([0,1300])
    set(gca,'FontName','Times')
    xlabel('Time','FontName','Times')
    set(gca,'Position',[0.0564000000000000,0.0944000000000000,0.922800000000000,0.873800000000002])
end

%% Make a plot of the reno & DNI diagnostics.
if exist('plot_figure','var')
    major_line_width = 20;
    line_sep = (major_line_width - 3)/5;
    cols = [199,233,180;...
        127,205,187;...
        65,182,196;...
        44,127,184;...
        37,52,148]./255;
    % % Exploration plot of the tri-component analysis. Zoom in to have a look.
    f = figure('name','explore the nature of Reno criteria','color','w');
    f.Position(3:4) = [1000,250];
    hold on
    plot(ghi)
    plot(dni)
    plot(dif)
    plot(ghics,'k:')
    plot(dnics,'k:')
    plot(difcs,'k:')
    for c = 1:5
        eval(['c',num2str(c),'=ghi;']);
        crit = criteria_ghi(:,c);
        eval(['c',num2str(c),'(crit==1)=NaN;']);
        eval(['plot(c',num2str(c),',''linewidth'',major_line_width-line_sep*',num2str(c),',''color'',cols(',num2str(c),',:));']);
    end
    csdghi=ghi;csdghi(csd_ghi==1)=NaN;plot(csdghi,'m');
    csdghi=ghi;csdghi(csd_overall==1)=NaN;plot(csdghi,'r','linewidth',2);
    
    csddni=dni;csddni(csd_dni==1)=NaN;plot(csddni,'m');
    csddni=dni;csddni(csd_overall==1)=NaN;plot(csddni,'r','linewidth',2);
   
    for c = 1:5
        eval(['c',num2str(c),'=dif;']);
        crit = criteria_dif(:,c);
        eval(['c',num2str(c),'(crit==1)=NaN;']);
        eval(['plot(c',num2str(c),',''linewidth'',major_line_width-line_sep*',num2str(c),',''color'',cols(',num2str(c),',:));']);
    end
    csddif=dif;csddif(csd_dif==1)=NaN;plot(csddif,'m');
    csddif=dif;csddif(csd_overall==1)=NaN;plot(csddif,'r','linewidth',2);
    hold off
    leg = legend('GHI','DNI','DIF','GHIcs','DNIcs','DIFcs','c1','c2','c3','c4','c5','CSD_{g/b/d}','tri-CSD');
    ylabel(gca,'Irradiance [Wm$^{-2}$]','FontName','Times','Interpreter','latex')
    set(gca,'xtick',[]);
    ylim([0,1300])
    set(gca,'FontName','Times')
    xlabel('Time','FontName','Times')
    set(gca,'Position',[0.0564000000000000,0.0944000000000000,0.922800000000000,0.873800000000002])
end


%% Duration criteria
% we build an hour durational filter looking ahead and behind for 45
% minutes. should there not have been a continuous CSD for an hour, then we
% reject all those instances.

% safety check on windows
if window_1st_filter<=1
    if tolerance_1st_filter==window_1st_filter
        tolerance_1st_filter = 2;
    end
    window_1st_filter=2;
end
if window_2nd_filter<=1
    if tolerance_2nd_filter==window_2nd_filter
        tolerance_2nd_filter = 2;
    end
    window_2nd_filter=2;
end

% First Duration Filter
csdh_1st_duration = nansum(hankel(csd_overall,[csd_overall(end),NaN(1,window_1st_filter-1)]),2);
csdh_1st_duration = [NaN(window_1st_filter/2,1);csdh_1st_duration(1:end-window_1st_filter/2)];
csd_1st_duration=zeros(size(ghi));
csd_1st_duration(csdh_1st_duration>tolerance_1st_filter)=1;
% this test makes morning and evening during sunlight hours impossible,
% therefore we relax the CSD at low zenith. However, in polar climates,
% some whole weeks can be >80.
A = (1:length(zen))';
B = find(round(zen)==85);
[~, distance_to_nearest_sunrise_set] = knnsearch(B,A);
csd_1st_duration(distance_to_nearest_sunrise_set<window_1st_filter)=0;

% Second Duration Filter
csdh_2nd_duration = nansum(hankel(csd_overall,[csd_overall(end),NaN(1,window_2nd_filter-1)]),2);
csdh_2nd_duration = [NaN(window_2nd_filter/2,1);csdh_2nd_duration(1:end-window_2nd_filter/2)];
csd_2nd_duration=zeros(size(ghi));
csd_2nd_duration(csdh_2nd_duration>tolerance_2nd_filter)=1;

% As sun-down is considered cloudy (which it should not be), the second
% filter unfavourably elimiates these periods. For that reason, the scond
% filter during sun down proximity (defined by the
% distance_to_nearest_sunset principle) is overridden by a less strict
% duration filter.
% Proximity to Sundown Duration Filter
window_3rd_duration = 10;
tolerance_3rd_duration = 2;
csdh_3rd_duration= nansum(hankel(csd_overall,[csd_overall(end),NaN(1,window_3rd_duration-1)]),2);
csdh_3rd_duration = [NaN(window_3rd_duration/2,1);csdh_3rd_duration(1:end-window_3rd_duration/2)];
csd_3rd_duration=zeros(size(ghi));
csd_3rd_duration(csdh_3rd_duration>tolerance_3rd_duration)=1;

% override the 2nd duration filter at peak proximity to sunrise defined as
% within the first filter window mins of sunrise and where the 10min filter
% found it permissable.
csd_2nd_duration(csd_3rd_duration==0 & distance_to_nearest_sunrise_set<window_1st_filter) = 0;

% CSD corroboration from the tri-component and two duration filters
csd_nan_indexed = ones(size(ghi));
csd_nan_indexed (csd_overall+csd_1st_duration+csd_2nd_duration==0)=0;

%% outputs that are indexed.
csd = ones(n_output,1);
csd(not_nan_indxs)=csd_nan_indexed;

% output the GHI/DNI/DIF cs optimised
ghics_optimised = NaN(n_output,1);
dnics_optimised = NaN(n_output,1);
difcs_optimised = NaN(n_output,1);

ghics_optimised(not_nan_indxs) = ghics;
dnics_optimised(not_nan_indxs) = dnics;
difcs_optimised(not_nan_indxs) = difcs;

% output the CSD from tri component;
csd_tricomponent = NaN(n_output,1);
csd_tricomponent(not_nan_indxs) = csd_overall;


%% plots
if exist('plot_figure','var')
    % % Exploration plot of the Duration Filters. Zoom in to look
    f = figure('name','durational window example','color','w');
    f.Position(3:4) = [1000,250];
    csdtri = ghi;csdtri(csd_overall==1)=NaN;
    csd60w = ghi;csd60w(csd_1st_duration==1)=NaN;
    csd10w = ghi;csd10w(csd_2nd_duration==1)=NaN;
    csdall = ghi;csdall(csd_nan_indexed==1)=NaN;
    hold on
    plot(ghi)
    plot(csdtri,'color','k','LineWidth',10)
    plot(csd60w,'color',[215,48,31]./255,'LineWidth',7.5)
    plot(csd10w,'color',[254,196,79]./255,'LineWidth',5)
    plot(csdall,'color',[120,198,121]./255,'LineWidth',2.5)
    hold off
    leg=legend({'GHI','Tri-Component CSD','Passed 90-min Filter','Passed 30-min Filter','Bright-Sun CSD'},'FontName','Times');
    leg.Position = [0.914866927054168,0.723066661478678,0.0762000008392335,0.239600005187988];
    set(gca,'FontName','Times')
    ylim([0 1300])
    set(gca,'XTickLabel',[]);
    xlabel('Time','FontName','Times')
    ylabel('Irradiance [Wm$^{-2}$]','FontName','Times','Interpreter','latex')
    set(gca,'Position',[0.0564000000000000,0.0944000000000000,0.922800000000000,0.873800000000002])
end

%% plot an example csd from the input data.
if exist('plot_figure','var')
    %% COMPLEX FIGURE WITH ALL CSD STEPS AND CLEAR_SKY CURVES
    colours = [215,48,39;...
        254,224,139;...
        145,207,96;...
        120,198,121]./255;
    
    figure('name','Bright-Sun CSD 2019 example','color','w');
    hold on
    plot(ghi)
    plot(dni)
    plot(dif)
    ghiCSD = ghi; ghiCSD(csd_ghi==1)=NaN;
    dniCSD = dni; dniCSD(csd_dni==1)=NaN;
    difCSD = dif; difCSD(csd_dif==1)=NaN;
    plot(ghics,'k:')
    plot(ghiCSD,'linewidth',12.5,'color','k')
    ghiCSD = ghi; ghiCSD(csd_overall==1)=NaN;
    plot(ghiCSD,'LineWidth',10,'color',colours(1,:))
    ghiCSD = ghi; ghiCSD(csd_1st_duration==1)=NaN;
    plot(ghiCSD,'LineWidth',7.5,'color',colours(2,:))    
    ghiCSD = ghi; ghiCSD(csd_2nd_duration==1)=NaN;
    plot(ghiCSD,'LineWidth',5,'color',colours(3,:))
    ghiCSD = ghi; ghiCSD(csd_nan_indexed==1)=NaN;
    plot(ghiCSD,'LineWidth',2,'color','m')
    
    plot(dnics,'k:')
    plot(difcs,'k:')
    plot(dniCSD,'linewidth',2,'color','k')
    plot(difCSD,'linewidth',2,'color','k')
    hold off
    if length(ghi)>20000
        xlim([1,20000])
    end
    legend({'GHI','DNI','DIF','Clear-sky irrad.','Tri-component CSD','Corroborating CSD','90min window','30min window CSD','Final CSD'},'fontname','times')
    ylabel(gca,'Irradiance [Wm$^{-2}$]','Interpreter','latex')
    set(gca,'xtick',[]);
    xlabel('Time','Interpreter','latex')
    set(gca,'fontname','times')

    colours = [215,48,39;254,224,139;145,207,96]./255;
    figure('name','Bright-Sun CSD 2019 example','color','w');
    hold on
    plot(ghi)
    plot(ghics,'k:')
    ghiCSD = ghi; ghiCSD(csd(not_nan_indxs)==1)=NaN;
    plot(ghiCSD,'LineWidth',4,'color',colours(3,:))
    hold off
    legend('GHI','GHIcs','CSD')
    ylabel(gca,'Irradiance [Wm$^{-2}$]','Interpreter','latex')
    set(gca,'xtick',[]);
    xlabel('Time')

end
end