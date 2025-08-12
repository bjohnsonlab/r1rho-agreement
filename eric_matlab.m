% MATLAB code for generating Figure 3

% Written by Eric G. Keeler
% July 2022
% ekeeler@nysbc.org; aem5@columbia.edu
% Cite as: E.G.Keeler, A.E.McDermott, ChemRev 2022
%% Figure 3(a) Plots
clear
close all
% Load the colormap for contour plots
load('../Fig2_3_4_Input_files/colorscheme.mat')
% Load the Experimental DMS data taken from Quinn and McDermott, JMR 2012
load('../Fig2_3_4_Input_files/DMS_exp_val.mat')

%Set the temperatures run in the data
DMS.T=[37,57,67,77]+273;

% Fit variables for 37°C point
wI=99.8e6*2*pi;
MAS=10e3;
CSA=-37e-6;
S=0.3295; %Initial S2 guess
T=1e-6; %Initial T guess

% Fit the DMS data at 37°C using the equations from Kurbanov et al., JCP 2011
tic
    try [fitresult,gof,jacobian]=Kurbanov_fit_v3(DMS_irr_power_37.*1e3.*2.*pi,DMS_rate_37-(1.8),CSA,MAS,wI,T,S);
            c_values=coeffvalues(fitresult);
            S2=c_values(1);
            tau=c_values(2);
            c_int=confint(fitresult,0.95);
            ci_tau = c_int(:,2);
            ci_S2 = c_int(:,1);
            rsquared = gof.rsquare;
            rmse = gof.rmse;

    catch ME
            disp(ME);
            fprintf('Curve failed to fit: %s\n', ME.message);
            S2=NaN;
            tau=NaN;
            ci_tau=NaN;
            ci_S2=NaN;
            rsquared = NaN;
            rmse = NaN;
    end

% Fix the S2 value to 0.3295 and fit again
S = 0.3295; %Fixed S2 value
    try [fitresult,gof,jacobian_fixed]=Kurbanov_fit_v3_fixedS(DMS_irr_power_37.*1e3.*2.*pi,DMS_rate_37-(1.8),CSA,MAS,wI,T,S);
            c_values=coeffvalues(fitresult);
            tau_fixed=c_values(1);
            c_int=confint(fitresult,0.95);
            ci_tau_fixed = c_int(:,1);
            rmse_fixed = gof.rmse;

    catch ME
            fprintf('Curve failed to fit', ME.message);
            tau_fixed=NaN;
            ci_tau_fixed=NaN;
            rmse_fixed = NaN;
    end
toc
clear c_int c_values fitresult gof ME q

% Save the results for the 37°C fit (variable S2)
    fprintf("37°C, fixed S2:\ntau = %s\n", tau);
    DMS.t7.tau=tau;
    DMS.t7.k=1/(2*tau);
    DMS.t7.S2=S2;
    DMS.t7.hopang=acosd(sqrt(4.*S2-1)./sqrt(3));
    DMS.t7.tau_err=DMS.t7.tau-ci_tau;
    DMS.t7.S2_err=DMS.t7.S2-ci_tau;
    DMS.t7.tau_gof_r2=rsquared;
    DMS.t7.rmse=rmse;
    DMS.t7.jacobian=jacobian;
% Save the results for the 37°C fit (fixed S2)
    DMS.t7.tau_fixed=tau_fixed;
    DMS.t7.k_fixed=1/(2*tau_fixed);
    DMS.t7.hopang_fixed=109;
    DMS.t7.tau_err_fixed=DMS.t7.tau_fixed-ci_tau_fixed;
    DMS.t7.rmse_fixed=rmse_fixed;
    DMS.t7.jacobian_fixed=jacobian_fixed;

clear tau k S2 tau_err S2_err tau_gof rmse jacobian jacobian_fixed


% Fit variables for 57°C point
wI=99.8e6*2*pi;
MAS=10e3;
CSA=-37e-6;
S=0.3295; %Initial S2 guess
T=1e-5; %Initial T guess

% Fit the DMS data at 57°C using the equations from Kurbanov et al., JCP 2011
tic
     try [fitresult,gof,jacobian]=Kurbanov_fit_v3(DMS_irr_power_57.*1e3.*2.*pi,DMS_rate_57-(1.8),CSA,MAS,wI,T,S);
            c_values=coeffvalues(fitresult);
            S2=c_values(1);
            tau=c_values(2);
            c_int=confint(fitresult,0.95);
            ci_tau = c_int(:,2);
            ci_S2 = c_int(:,1);
            rsquared = gof.rsquare;
            rmse = gof.rmse;

    catch ME
            fprintf('Curve failed to fit', ME.message);
            S2=NaN;
            tau=NaN;
            ci_tau=NaN;
            ci_S2=NaN;
            rsquared = NaN;
            rmse = NaN;
     end

% Fix the S2 value to 0.3295 and fit again
S=0.3295; %Fixed S2 value
    try [fitresult,gof,jacobian_fixed]=Kurbanov_fit_v3_fixedS(DMS_irr_power_57.*1e3.*2.*pi,DMS_rate_57-(1.8),CSA,MAS,wI,T,S);
            c_values=coeffvalues(fitresult);
            tau_fixed=c_values(1);
            c_int=confint(fitresult,0.95);
            ci_tau_fixed = c_int(:,1);
            rmse_fixed = gof.rmse;

    catch ME
            fprintf('Curve failed to fit', ME.message);
            tau_fixed=NaN;
            ci_tau_fixed=NaN;
            rmse_fixed = NaN;
    end
toc
clear c_int c_values fitresult gof ME q

% Save the results for the 57°C fit (variable S2)
    DMS.f7.tau=tau;
    DMS.f7.k=1/(2*tau);
    DMS.f7.S2=S2;
    DMS.f7.hopang=acosd(sqrt(4.*S2-1)./sqrt(3));
    DMS.f7.tau_err=DMS.f7.tau-ci_tau;
    DMS.f7.S2_err=DMS.f7.S2-ci_tau;
    DMS.f7.tau_gof_r2=rsquared;
    DMS.f7.rmse=rmse;
    DMS.f7.jacobian=jacobian;
% Save the results for the 57°C fit (fixed S2)
    DMS.f7.tau_fixed=tau_fixed;
    DMS.f7.k_fixed=1/(2*tau_fixed);
    DMS.f7.hopang_fixed=109;
    DMS.f7.tau_err_fixed=DMS.f7.tau_fixed-ci_tau_fixed;
    DMS.f7.rmse_fixed=rmse_fixed;
    DMS.f7.jacobian_fixed=jacobian_fixed;

clear tau k S2 tau_err S2_err tau_gof rmse jacobian jacobian_fixed


% Fit variables for 67°C point
wI=99.8e6*2*pi;
MAS=10e3;
CSA=-37e-6;
S=0.3295; %Initial S2 guess
T=1e-5; %Initial T guess

% Fit the DMS data at 67°C using the equations from Kurbanov et al., JCP 2011
tic
     try [fitresult,gof,jacobian]=Kurbanov_fit_v3(DMS_irr_power_67.*1e3.*2.*pi,DMS_rate_67-(1.8),CSA,MAS,wI,T,S);
            c_values=coeffvalues(fitresult);
            S2=c_values(1);
            tau=c_values(2);
            c_int=confint(fitresult,0.95);
            ci_tau = c_int(:,2);
            ci_S2 = c_int(:,1);
            rsquared = gof.rsquare;
            rmse = gof.rmse;

    catch ME
            fprintf('Curve failed to fit', ME.message);
            S2=NaN;
            tau=NaN;
            ci_tau=NaN;
            ci_S2=NaN;
            rsquared = NaN;
            rmse = NaN;
     end

% Fix the S2 value to 0.3295 and fit again
S=0.3295; %Fixed S2 value
    try [fitresult,gof,jacobian_fixed]=Kurbanov_fit_v3_fixedS(DMS_irr_power_67.*1e3.*2.*pi,DMS_rate_67-(1.8),CSA,MAS,wI,T,S);            c_values=coeffvalues(fitresult);
            tau_fixed=c_values(1);
            c_int=confint(fitresult,0.95);
            ci_tau_fixed = c_int(:,1);
            rmse_fixed = gof.rmse;

    catch ME
            fprintf('Curve failed to fit', ME.message);
            tau_fixed=NaN;
            ci_tau_fixed=NaN;
            rmse_fixed = NaN;
    end
toc
clear c_int c_values fitresult gof ME q

% Save the results for the 67°C fit (variable S2)
    DMS.si7.tau=tau;
    DMS.si7.k=1/(2*tau);
    DMS.si7.S2=S2;
    DMS.si7.hopang=acosd(sqrt(4.*S2-1)./sqrt(3));
    DMS.si7.tau_err=DMS.si7.tau-ci_tau;
    DMS.si7.S2_err=DMS.si7.S2-ci_tau;
    DMS.si7.tau_gof_r2=rsquared;
    DMS.si7.rmse=rmse;
    DMS.si7.jacobian=jacobian;
% Save the results for the 67°C fit (fixed S2)
    DMS.si7.tau_fixed=tau_fixed;
    DMS.si7.k_fixed=1/(2*tau_fixed);
    DMS.si7.hopang_fixed=109;
    DMS.si7.tau_err_fixed=DMS.si7.tau_fixed-ci_tau_fixed;
    DMS.si7.rmse_fixed=rmse_fixed;
    DMS.si7.jacobian_fixed=jacobian_fixed;

clear tau k S2 tau_err S2_err tau_gof rmse jacobian jacobian_fixed


% Fit variables for 77°C point
wI=99.8e6*2*pi;
MAS=10e3;
CSA=-37e-6;
S=0.3295; %Initial S2 guess
T=6e-5; %Initial T guess

% Fit the DMS data at 77°C using the equations from Kurbanov et al., JCP 2011
tic
     try [fitresult,gof,jacobian]=Kurbanov_fit_v3(DMS_irr_power_77.*1e3.*2.*pi,DMS_rate_77-(1.8),CSA,MAS,wI,T,S);
            c_values=coeffvalues(fitresult);
            S2=c_values(1);
            tau=c_values(2);
            c_int=confint(fitresult,0.95);
            ci_tau = c_int(:,2);
            ci_S2 = c_int(:,1);
            rsquared = gof.rsquare;
            rmse = gof.rmse;

    catch ME
            fprintf('Curve failed to fit', ME.message);
            S2=NaN;
            tau=NaN;
            ci_tau=NaN;
            ci_S2=NaN;
            rsquared = NaN;
            rmse = NaN;
     end

% Fix the S2 value to 0.3295 and fit again
S=0.3295; %Fixed S2 value
    try [fitresult,gof,jacobian_fixed]=Kurbanov_fit_v3_fixedS(DMS_irr_power_77.*1e3.*2.*pi,DMS_rate_77-(1.8),CSA,MAS,wI,T,S);            c_values=coeffvalues(fitresult);
            tau_fixed=c_values(1);
            c_int=confint(fitresult,0.95);
            ci_tau_fixed = c_int(:,1);
            rmse_fixed = gof.rmse;

    catch ME
            fprintf('Curve failed to fit', ME.message);
            tau_fixed=NaN;
            ci_tau_fixed=NaN;
            rmse_fixed = NaN;
    end
toc
clear c_int c_values gof ME q

% Save the results for the 77°C fit (variable S2)
    DMS.se7.tau=tau;
    DMS.se7.k=1/(2*tau);
    DMS.se7.S2=S2;
    DMS.se7.hopang=acosd(sqrt(4.*S2-1)./sqrt(3));
    DMS.se7.tau_err=DMS.se7.tau-ci_tau;
    DMS.se7.S2_err=DMS.se7.S2-ci_tau;
    DMS.se7.tau_gof_r2=rsquared;
    DMS.se7.rmse=rmse;
    DMS.se7.jacobian=jacobian;
% Save the results for the 77°C fit (fixed S2)
    DMS.se7.tau_fixed=tau_fixed;
    DMS.se7.k_fixed=1/(2*tau_fixed);
    DMS.se7.hopang_fixed=109;
    DMS.se7.tau_err_fixed=DMS.se7.tau_fixed-ci_tau_fixed;
    DMS.se7.rmse_fixed=rmse_fixed;
    DMS.se7.jacobian_fixed=jacobian_fixed;

clear tau k S2 tau_err S2_err tau_gof rmse jacobian jacobian_fixed

% Plot the results of the DMS fits

% Set the variables for the fixed S2 plot
MAS=10e3;
CSA=-37e-6;
S = 0.3295;
T = DMS.t7.tau_fixed;
X=linspace(25,42,35)*1e3*2*pi;
wI=6.283185307179586e+08;

% Generate the R1rho curves for the plots
R1rho_plot_37_fixed = (3/4).*(CSA.*wI)^2.*((2/5)*...
    (1-S).*(T./(1+((wI).*T).^2)))+(((1/6)*...
    (CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))));
T = DMS.f7.tau_fixed;
R1rho_plot_57_fixed = (3/4).*(CSA.*wI)^2.*((2/5)*...
    (1-S).*(T./(1+((wI).*T).^2)))+(((1/6)*...
    (CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))));
T = DMS.si7.tau_fixed;
R1rho_plot_67_fixed = (3/4).*(CSA.*wI)^2.*((2/5)*...
    (1-S).*(T./(1+((wI).*T).^2)))+(((1/6)*...
    (CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))));
T = DMS.se7.tau_fixed;
R1rho_plot_77_fixed = (3/4).*(CSA.*wI)^2.*((2/5)*...
    (1-S).*(T./(1+((wI).*T).^2)))+(((1/6)*...
    (CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))));

tau37 = DMS.t7.tau_fixed;
tau57 = DMS.f7.tau_fixed;
tau67 = DMS.si7.tau_fixed;
tau77 = DMS.se7.tau_fixed;
disp(R1rho_plot_77_fixed);
fprintf("tau37: %s", tau37);
fprintf("tau57: %s", tau57);
fprintf("tau67: %s", tau67);
fprintf("tau77: %s", tau77);
save("r1rho_curves.mat", "X", "tau37", "R1rho_plot_37_fixed", "tau57", "R1rho_plot_57_fixed", "tau67", "R1rho_plot_67_fixed", "tau77", "R1rho_plot_77_fixed");
X = linspace(1.49, 49.02, 98) * 1e3 * 2 *pi;
disp(X);
T = DMS.se7.tau_fixed;
disp(T); disp(log10(T));
r1rho = (3/4).*(CSA.*wI)^2.*((2/5)*...
    (1-S).*(T./(1+((wI).*T).^2)))+(((1/6)*...
    (CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))));
disp(mat2str(r1rho));
return;

% Plot the R1rho curves for fixed S2 results
f2=figure(1);
axes2 = axes('Parent',f2);
plot(DMS_irr_power_37,DMS_rate_37,'co')
hold on
plot(linspace(25,42,35),R1rho_plot_37_fixed,'c--')
hold on
plot(DMS_irr_power_57,DMS_rate_57,'go')
hold on
plot(linspace(25,42,35),R1rho_plot_57_fixed,'g--')
hold on
plot(DMS_irr_power_67,DMS_rate_67,'ro')
hold on
plot(linspace(25,42,35),R1rho_plot_67_fixed,'r--')
hold on
plot(DMS_irr_power_77,DMS_rate_77,'bo')
hold on
plot(linspace(25,42,35),R1rho_plot_77_fixed,'b--')
hold on
legend('37°C, exp','37°C MF fit','57°C, exp','57°C MF fit','67°C, exp','67°C MF fit','77°C, exp','77°C MF fit')
title('Model Free with S^2 = 0.3295 (\beta = 109°)')
xlabel('\nu_1 (kHz)')
ylabel('R_1_\rho (s^-^1)')
xlim([24 42])
set(axes2,'BoxStyle','full','FontSize',14,'Layer','top',...
    'TickDir','out','XMinorTick','on','YMinorTick','on');

%% Fit to a single S2 with varying tau_c using all temperature data at once
f3=figure(3);
% Set the field strengths for each set of temperature points (some points
% are repeated to make the length of each dataset equal to 7)
X = [25;28;32;35;38;25;28;25.6;31.3;35.7;41.7;25.6;31.3;35.7;25.6;28.8;31.3...
    ;33.3;35.7;38.5;41.7;25.6;31.3;35.7;38.5;41.7;25.6;31.3];
% Set the rate constant for each set of DMS data from Quinn and McDermott,
% JMR 2012 (some points are repeated to make the length of each dataset
% equal to 7)
Y = [11.8;5.9;5.1;3.5;3.15;11.8;5.9;35;19.7;13.9;9.3;35;19.7;13.9;75.1;50.7...
    ;39.2;28.9;24.6;23.2;18.2;104.5;58.7;40.3;32.6;30.8;104.5;58.7];
% Correct the rate constant for the R_2^0 value
Y = Y-1.8;
% Put the field strengths in the correct values
X=X*2*pi*1e3;
% Set a variable that defines the length of each data set to 7
dsid = repelem((1:4)',7);
% Run the fit and plot the results
h=gscatter(X,Y,dsid,'cgrb','oooo');
x = X;
wI=6.283185307179586e+08;
CSA= -37e-6;
f0 = @(T,S,x) (3/4)*(CSA*wI)^2.*((2/5)*(1-S).*...
    (T(dsid)./(1+((wI).*T(dsid)).^2)))+(((1/6)*...
    (CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*(T(dsid)./...
    (1+((x-2.*10e3*2*pi).*T(dsid)).^2)))+((2/5)*(1-S).*(T(dsid)./...
    (1+((x-10e3*2*pi).*T(dsid)).^2)))+((2/5)*(1-S).*(T(dsid)./...
    (1+((x+10e3*2*pi).*T(dsid)).^2)))+((1/2).*(2/5)*(1-S).*(T(dsid)./...
    (1.+((x+2*10e3*2*pi).*T(dsid)).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*(T(dsid)./...
    (1+((wI).*T(dsid)).^2)))));
f1 = @(params,x) f0(params(1:4),params(5),x);
start = [1e-4;1e-4;1e-4;1e-4;0.3295];
lb=[0;0;0;0;0];
ub=[Inf;Inf;Inf;Inf;1];
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
b = lsqcurvefit(f1,start,x,Y,lb,ub,options);
f2= @(T,S,x) (3/4)*(CSA*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))+(((1/6)*(CSA*wI)^2.*...
    (((1/2).*(2/5)*(1-S).*(T./(1+((x-2.*10e3*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((x-10e3*2*pi).*T).^2)))+((2/5)*(1-S).*(T./(1+((x+10e3*2*pi).*T).^2)))+...
    ((1/2).*(2/5)*(1-S).*(T./(1.+((x+2*10e3*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*(T./(1+((wI).*T).^2)))));

DMS.allfit.k=1./(2.*b(1:4));
DMS.allfit.S=b(5);

xgrid = linspace(24,42,35)*2*pi*1e3;
for j=1:length(h)
    hh = line(xgrid,f2(b(j),b(5),xgrid));
    hh.Color = h(j).Color;
end
xlim([24 42]*2*pi*1e3)
ylim([0 200])
%

%% Figure 3(d)

% Angle between the sigma_33 component of the CA CSA tensors in °
angle=109;
% Determine the fixed S2 for that angle.
fixed_S2=1/4*(3*cos(angle*pi/180)^2+1);
% Generate the Plot for ln(k) vs 1/T for Figure 3(d)
f4=figure('Name','ln(k) vs 1/T Plot');
ax4=axes('Parent',f4);
% Plot the Numerical Simulations from Quinn and McDermott
plot(1./DMS.T,log([300,1600,4500,7500]),'ko')
hold on
% Plot the Global S2 fit
plot(1./DMS.T,log([DMS.allfit.k]),'ks')
hold on
% Plot the fixed S2 fit
plot(1./DMS.T,log(1./(2.*[DMS.t7.tau_fixed,DMS.f7.tau_fixed,DMS.si7.tau_fixed,DMS.se7.tau_fixed])),'rs')
hold on
% Plot the Spinach Simulation results from Keeler and McDermott, ChemRev 2022
plot(1./DMS.T,log([315, 1406, 4446, 7498]),'cd')
hold on
globalstr = sprintf('Global S^2 = %.4f ',DMS.allfit.S);
fixedstr = sprintf('Fixed S^2 = %.4f ',fixed_S2);
legend('Exp',globalstr,'Fixed S^2=0.3295','Numerical Simulation','Location','southwest')
xlabel('1/T (K^-^1)');
ylabel('ln(k) (s^-^1)');
set(ax4,'BoxStyle','full','FontSize',14,'Layer','top',...
    'TickDir','out','XMinorTick','on','YMinorTick','on','yscale','lin','xscale','lin');
ylim([0 12]);
xlim([2.65e-3 3.75e-3]);
clear f4 ax4

%% Calculate and plot the RMSE for a set of S2 and tau_c data points for 77°C

% Set variables for the plots and calculations
MAS=10e3;
CSA=-37e-6;
s = logspace(-2.5,0,50);
s = 1-s;
s(end) = 0;
t = logspace(-8,-1,100);
wI=6.283185307179586e+08;
% Save variables into a data space to use in plots
Results.experiment=DMS_rate_77-1.8;
Results.experiment_2=DMS_rate_77-1.8;
y=DMS_irr_power_77*1e3*2*pi;
Results.R1rho_sim=zeros(length(y),length(s),length(t));
% Generate the R1rho curves from Kurbanov et al. JCP 2011 equations
for l=1:length(y)
    for i=1:length(s)
        for k=1:length(t)
            Results.R1rho_sim(l,i,k) = (3/4).*(CSA.*wI)^2.*((2/5)*...
                (1-s(i)).*(t(k)./(1+((wI).*t(k)).^2)))+(((1/6)*...
                (CSA*wI)^2.*(((1/2).*(2/5)*(1-s(i)).*(t(k)./...
                (1+((y(l)-2.*MAS*2*pi).*t(k)).^2)))+((2/5)*(1-s(i)).*...
                (t(k)./(1+((y(l)-MAS*2*pi).*t(k)).^2)))+((2/5)*(1-s(i)).*...
                (t(k)./(1+((y(l)+MAS*2*pi).*t(k)).^2)))+((1/2).*(2/5)*...
                (1-s(i)).*(t(k)./(1.+((y(l)+2*MAS*2*pi).*t(k)).^2)))))...
                -0.5*((3/4)*(CSA*wI)^2.*((2/5)*(1-s(i)).*(t(k)./...
                (1+((wI).*t(k)).^2)))));
        end
    end
end
% Generate the residuals of the Kurbanov curve in comparison to the
% experimental results on DMS
Results.residuals=zeros(length(y),length(s),length(t));
for l=1:length(y)
    for i=1:length(s)
        for k=1:length(t)
                Results.residuals(l,i,k)=Results.experiment_2(l)-Results.R1rho_sim(l,i,k);
        end
    end
end
Results.residuals(isnan(Results.residuals))=0;
% Use the jacobian of the DMS fit from above to adjust the residuals
    jacob=full(DMS.se7.jacobian.Jacobian);
    [Q,~]=qr(jacob,0);
    h=min(0.9990,sum(Q.*Q,2));
    adj_fac=1./sqrt(1-h);
    Results.adjust.residuals=Results.residuals.*adj_fac;
% Bisquare weight RMSE for 77°C by adjusting the residuals due to
% leverage using MF fit from above
    BW=zeros(length(y),length(s),length(t));
    R=zeros(length(y),length(s),length(t));
    Results.BW.residuals=zeros(length(y),length(s),length(t));
    Results.BW.RMSE=zeros(length(s),length(t));
    RS = sort(abs(Results.adjust.residuals));
    MAR = median(RS(2:end,:,:));
    SIGMA = MAR./0.6745;
    tune=4.685;
    tSIGMA=tune.*SIGMA;

for l = 1:length(y)
    for i = 1:length(s)
        for k = 1:length(t)
        %Calculate the Bisquare weighting function
        if abs(Results.residuals(l,i,k)) > tSIGMA(:,i,k)
            BW(l,i,k) = 0;
        else
            R(l,i,k) = Results.residuals(l,i,k)/tSIGMA(:,i,k);
            BW(l,i,k) = (abs(R(l,i,k))<1).*(1-(R(l,i,k)).^2).^2;
        end
        Results.BW.residuals(l,i,k) = Results.residuals(l,i,k) .* BW (l,i,k);
        end
    end
end
clear i l k
% Generate the RMSE and the Bisquare Weighted RMSE
for i=1:length(s)
    for k=1:length(t)
        Results.RMSE(i,k)=sqrt(mean(Results.residuals(:,i,k).^2));
        Results.BW.RMSE(i,k)=sqrt(mean(Results.BW.residuals(:,i,k).^2));
    end
end

% Plot the Bisquare Weighted RMSE normalized by the average rate constant
% of the R1rho dispersion curve of the 77°C data - FIGURE 3(b, left)
f1=figure('Name', 'Contour Plot of 77°C, BW v2');
axes1 = axes('Parent',f1);
contourf(t,1-s,log10((squeeze(Results.BW.RMSE(:,:)))./squeeze(mean(Results.R1rho_sim))),25,'LineStyle','none')
hold on
% Minimum of fit with only 77°C data
scatter(DMS.se7.tau,1-DMS.se7.S2,'ko','filled')
hold on
% Minimum of fit with fixed S2 = 0.3295
scatter(DMS.se7.tau_fixed,1-0.3295,'ro','filled')
hold on
colormap(flipud(colorscheme));
set(axes1,'BoxStyle','full','FontSize',14,'Layer','top',...
'TickDir','out','XMinorTick','on','YMinorTick','on','yscale','lin','xscale','log');
c=colorbar('peer',axes1);
c.TickDirection='out';
xlabel('\tau_c (s)')
ylabel('1-S^2')
zlabel('RMSE of MF fit')
xlim([2e-8 1e-2])
caxis([-1.5 3.8])
title('Contour Plot of RMSE for 77°C, BW')

clear f1 axes1

%% RMSE calculations for Numerical Simulation Grid

% Load the Experimental DMS data taken from Quinn and McDermott, JMR 2012
    load('Fig2_3_4_Input_files/DMS_exp_val.mat')
    DMS_rate_77=DMS_rate_77-1.8;
    clear DMS_rate_37 DMS_rate_57 DMS_rate_67
% Load the numerical simulations of varying S2 and tau_c data
    load('Fig2_3_4_Input_files/Fig2_data.mat')
% Extract the numerical simulations points for each temperature point of
% experimental DMS data so that they align
    Results.NumSim.irr_power_77=Results.NumSim.irr_power([2,5,9,11,12],1);
    Results.NumSim.curve_77=Results.NumSim.curve([2,5,9,11,12],:,:);

% Determine the Residuals and RMSE for the numerical simulation data set
% with respect to the experimental data for 77°C data.
Results.NumSim.se7.residuals=zeros(length(DMS_irr_power_77),length(Results.NumSim.T),length(Results.NumSim.S));
for l=1:length(DMS_irr_power_77)
            for h=1:length(Results.NumSim.T)
                for q=1:length(Results.NumSim.S)
                    Results.NumSim.se7.residuals(l,h,q)=DMS_rate_77(l)-Results.NumSim.curve_77(l,h,q);
                    Results.NumSim.se7.residuals(1,k)=0.67*Results.NumSim.se7.residuals(1,k);
                end
            end
end
%
Results.NumSim.se7.residuals(isnan(Results.NumSim.se7.residuals))=0;

for i=1:length(DMS_irr_power_77)
        for h=1:length(Results.NumSim.T)
            for q=1:length(Results.NumSim.S)
                Results.NumSim.se7.RMSE(h,q)=sqrt(mean(Results.NumSim.se7.residuals(:,h,q).^2));
            end
        end
end

%Calculate Bisquare Weighted residuals and RMSE
[Results.NumSim.se7.BW.residuals]=bisquare_residuals_3D(DMS.se7.jacobian.Jacobian,Results.NumSim.se7.residuals,DMS_irr_power_77,Results.NumSim.T,Results.NumSim.S);

for k=1:length(Results.NumSim.T)
    for q=1:length(Results.NumSim.S)
        Results.NumSim.se7.BW.RMSE(k,q)=sqrt(mean(Results.NumSim.se7.BW.residuals(:,k,q).^2));
    end
end

% Plot the Bisquare Weighted RMSE normalized by the average rate constant
% of the R1rho dispersion curve of the 77°C data - FIGURE 3(b, right)
f1=figure('Name', 'Contour Plot of 77°C - BW -v2');
axes1 = axes('Parent',f1);
A=Results.NumSim.se7.BW.RMSE';
B=squeeze(mean(Results.NumSim.curve_77))';
contourf(Results.NumSim.T,1-Results.NumSim.S,log10(A./B),25,'LineStyle','none')
hold on
% Minimum of fit with fixed S2 = 0.3295
scatter(DMS.se7.tau_fixed,1-0.3295,'ro','filled')
hold on
colormap(flipud(colorscheme));
set(axes1,'BoxStyle','full','FontSize',14,'Layer','top',...
    'TickDir','out','XMinorTick','on','YMinorTick','on','yscale','lin','xscale','log');
c=colorbar('peer',axes1);
c.TickDirection='out';
xlabel('\tau_c (s)')
ylabel('1-S^2')
zlabel('RMSE of MF fit')
ylim([min(1-Results.NumSim.S) 1]);
xlim([2e-8 1e-2])
caxis([-1.5 3.8])
title('Contour Plot of Numerical Simulation RMSE for 77°C - BW')

%%
% Plot of 77°C with data focused around the minimum
load('Fig2_3_4_Input_files/Fig3_data.mat');
disp(wI);
[X,S,T]=meshgrid(Results.NumSim.FixedS.full.irr_power,Results.NumSim.FixedS.full.S,Results.NumSim.FixedS.full.T);
Results.NumSim.FixedS.R1rho_Kurbanov = (3/4).*(CSA.*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))+(((1/6)*(CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*(T./(1+((wI).*T).^2)))));
Results.NumSim.FixedS.R1rho_Kurbanov = squeeze(Results.NumSim.FixedS.R1rho_Kurbanov);
Results.NumSim.FixedS.full.se7.residuals=zeros(length(Results.NumSim.FixedS.full.irr_power),length(Results.NumSim.FixedS.full.T));
Results.NumSim.FixedS.full.se7.RMSE=zeros(length(Results.NumSim.FixedS.full.T),1);
for l=1:length(Results.NumSim.FixedS.full.irr_power)
        for k=1:length(Results.NumSim.FixedS.full.T)
                Results.NumSim.FixedS.full.se7.residuals(l,k)=(DMS_rate_77(l)-1.8)-Results.NumSim.FixedS.full.curve(l,k);
        end
end
Results.NumSim.FixedS.full.se7.residuals(isnan(Results.NumSim.FixedS.full.se7.residuals))=0;

%Bisquare Weighting of the residuals
[Results.NumSim.FixedS.full.se7.BW.residuals]=bisquare_residuals(...
    DMS.se7.jacobian_fixed.Jacobian,Results.NumSim.FixedS.full.se7.residuals,...
    Results.NumSim.FixedS.full.irr_power,Results.NumSim.FixedS.full.T);

    for k=1:length(Results.NumSim.FixedS.full.T)
        Results.NumSim.FixedS.full.se7.BW.RMSE(k)=sqrt(mean(Results.NumSim.FixedS.full.se7.BW.residuals(:,k).^2));
        Results.NumSim.FixedS.full.se7.RMSE(k)=sqrt(mean(Results.NumSim.FixedS.full.se7.residuals(:,k).^2));
    end

% Same for narrow data set (to define the curve better)
[X,S,T]=meshgrid(Results.NumSim.FixedS.narrow.irr_power,Results.NumSim.FixedS.narrow.S,Results.NumSim.FixedS.narrow.T);
Results.NumSim.FixedS.R1rho_Kurbanov = (3/4).*(CSA.*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))+(((1/6)*(CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*...
    (CSA*wI)^2.*((2/5)*(1-S).*(T./(1+((wI).*T).^2)))));
Results.NumSim.FixedS.R1rho_Kurbanov = squeeze(Results.NumSim.FixedS.R1rho_Kurbanov);
Results.NumSim.FixedS.narrow.se7.residuals=zeros(length(Results.NumSim.FixedS.narrow.irr_power),length(Results.NumSim.FixedS.narrow.T));
Results.NumSim.FixedS.narrow.se7.RMSE=zeros(length(Results.NumSim.FixedS.narrow.T),1);
for l=1:length(Results.NumSim.FixedS.narrow.irr_power)
        for k=1:length(Results.NumSim.FixedS.narrow.T)
                Results.NumSim.FixedS.narrow.se7.residuals(l,k)=(DMS_rate_77(l)-1.8)-Results.NumSim.FixedS.narrow.curve(l,k);
        end
end
Results.NumSim.FixedS.narrow.se7.residuals(isnan(Results.NumSim.FixedS.narrow.se7.residuals))=0;

% Bisquare Weighting of the residuals
[Results.NumSim.FixedS.narrow.se7.BW.residuals]=bisquare_residuals(...
    DMS.se7.jacobian_fixed.Jacobian,Results.NumSim.FixedS.narrow.se7.residuals,...
    Results.NumSim.FixedS.narrow.irr_power,Results.NumSim.FixedS.narrow.T);

    for k=1:length(Results.NumSim.FixedS.narrow.T)
        Results.NumSim.FixedS.narrow.se7.BW.RMSE(k)=sqrt(mean(Results.NumSim.FixedS.narrow.se7.BW.residuals(:,k).^2));
        Results.NumSim.FixedS.narrow.se7.RMSE(k)=sqrt(mean(Results.NumSim.FixedS.narrow.se7.residuals(:,k).^2));
    end

% Generate the R1rho residuals for the MF using the Kurbanov et al., JCP
% 2011 equations
single_s = 0.3295;
single_t = logspace(-8,-1,250);
x=DMS_irr_power_77*1e3*2*pi;
[X,S,T]=meshgrid(x,single_s,single_t);
wI=6.283185307179586e+08;
Results.singleS.R1rho_Kurbanov = (3/4).*(CSA.*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))+(((1/6)*(CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((X-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((X+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((X+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*(CSA*wI)^2.*...
    ((2/5)*(1-S).*(T./(1+((wI).*T).^2)))));
Results.singleS.R1rho_Kurbanov = squeeze(Results.singleS.R1rho_Kurbanov);
Results.singleS.residuals=zeros(length(y),length(single_s),length(single_t));
Results.singleS.RMSE=zeros(length(single_t),1);
for l=1:length(y)
for k=1:length(single_t)
Results.singleS.residuals(l,k)=(DMS_rate_77(l)-1.8)-Results.singleS.R1rho_Kurbanov(l,k);
end
end
Results.singleS.residuals(isnan(Results.singleS.residuals))=0;
for k=1:length(single_t)
Results.singleS.RMSE(k)=sqrt(mean(Results.singleS.residuals(:,k).^2));
end
% Determine the Jacobian for the Bisquare Weighting of the residuals
jacob=full(DMS.se7.jacobian.Jacobian);
[Q,~]=qr(jacob,0);
h=min(0.9990,sum(Q.*Q,2));
adj_fac=1./sqrt(1-h);
Results.singleS.adjust.residuals=Results.singleS.residuals.*adj_fac;
BW=zeros(length(y),length(single_s),length(single_t));
R=zeros(length(y),length(single_s),length(single_t));
Results.singleS.BW.residuals=zeros(length(y),length(single_s),length(single_t));
Results.singleS.BW.RMSE=zeros(length(single_s),length(single_t));
RS = sort(abs(Results.singleS.adjust.residuals));
MAR = median(RS(2:end,:,:));
SIGMA = MAR./0.6745;
tune=4.685;
tSIGMA=tune.*SIGMA;
% Bisquare Weighting of the residuals
for l = 1:length(y)
    for i=1:length(single_s)
        for k = 1:length(single_t)
        %Calculate the Bisquare weighting function
        if abs(Results.singleS.residuals(l,i,k)) > tSIGMA(:,k)
            BW(l,i,k) = 0;
        else
            R(l,i,k) = Results.singleS.residuals(l,i,k)/tSIGMA(:,i,k);
            BW(l,i,k) = (abs(R(l,i,k))<1).*(1-(R(l,i,k)).^2).^2;
        end
        Results.singleS.BW.residuals(l,i,k) = Results.singleS.residuals(l,i,k) .* BW (l,i,k);
        end
    end
end
clear i l k
    for k=1:length(single_t)
        for i=1:length(single_s)
        Results.singleS.RMSE(i,k)=sqrt(mean(Results.singleS.residuals(:,i,k).^2));
        Results.singleS.BW.RMSE(i,k)=sqrt(mean(Results.singleS.BW.residuals(:,i,k).^2));
        end
    end
    Results.singleS.BW.RMSE=squeeze(Results.singleS.BW.RMSE);

% Figure 3(c)
f1=figure();
axes1=axes('Parent',f1);
semilogx(single_t,log10(Results.singleS.BW.RMSE(:)./squeeze(mean(Results.singleS.R1rho_Kurbanov)')),'r-')
hold on
semilogx(Results.NumSim.FixedS.full.T,...
    log10(Results.NumSim.FixedS.full.se7.BW.RMSE(:)./...
    squeeze(mean(Results.NumSim.FixedS.full.curve))'),'k-')
hold on
semilogx(Results.NumSim.FixedS.narrow.T,...
    log10(Results.NumSim.FixedS.narrow.se7.BW.RMSE(:)./...
    squeeze(mean(Results.NumSim.FixedS.narrow.curve))'),'k-')
set(axes1,'BoxStyle','full','FontSize',14,'Layer','top',...
    'TickDir','out','XMinorTick','on','YMinorTick','on','yscale','lin','xscale','log');
xlabel('\tau_c (s)')
ylabel('RMSE of MF fit')
ylim([-1.5 2.2]);
xlim([2e-8 1e-2])
clear f1 axes1

%% Functions for the above code
% Function for the fit of the Kurbanov et al., JCP 2011 equations
function [fitresult, gof, Jacobian] = Kurbanov_fit_v3(omega1, fid,CSA,MAS,wI,T,S)
%% Mono-exponential fit
% Prepare the data
[xData, yData] = prepareCurveData( omega1, fid );
% Set up fittype and options.
ft = fittype( @(S,T,CSA,wI,MAS,x) (3/4)*(CSA*wI)^2.*((2/5)*(1-S).*...
    (T./(1+((wI).*T).^2)))+(((1/6)*(CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
    (T./(1+((x-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((x-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
    (T./(1+((x+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
    (T./(1.+((x+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*(CSA*wI)^2.*...
    ((2/5)*(1-S).*(T./(1+((wI).*T).^2))))), 'problem',{'CSA','wI','MAS'}, 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Algorithm = 'Trust-Region';
opts.Robust = 'Bisquare';
opts.StartPoint = [S*(1-(rand(1))*0.025), T+(T*((-1+(1+1)*rand(1))*0.08))];
opts.Lower = [0.01,0];
opts.Upper = [1.0,Inf];
opts.TolX = 10^-8;
opts.DiffMaxChange = 0.01;
opts.MaxFunEvals = 750;
opts.MaxIter = 500;
% Fit model to data.
[fitresult, gof, Jacobian] = fit( xData, yData, ft,'problem',{CSA, wI, MAS}, opts);
end
% Function for the fit of the Kurbanov et al., JCP 2011 equations with
% fixed S2
function [fitresult, gof, Jacobian] = Kurbanov_fit_v3_fixedS(omega1, fid,CSA,MAS,wI,T,S)
%% Mono-exponential fit
    % Prepare the data
    [xData, yData] = prepareCurveData( omega1, fid );

    % Set up fittype and options.
    ft = fittype( @(T,CSA,S,wI,MAS,x) (3/4)*(CSA*wI)^2.*((2/5)*(1-S).*...
        (T./(1+((wI).*T).^2)))+(((1/6)*(CSA*wI)^2.*(((1/2).*(2/5)*(1-S).*...
        (T./(1+((x-2.*MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
        (T./(1+((x-MAS*2*pi).*T).^2)))+((2/5)*(1-S).*...
        (T./(1+((x+MAS*2*pi).*T).^2)))+((1/2).*(2/5)*(1-S).*...
        (T./(1.+((x+2*MAS*2*pi).*T).^2)))))-0.5*((3/4)*(CSA*wI)^2.*...
        ((2/5)*(1-S).*(T./(1+((wI).*T).^2))))), 'problem',{'CSA','S','wI','MAS'}, 'independent', 'x', 'dependent', 'y');
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Algorithm = 'Trust-Region';
    % opts.Robust = 'LAR';
    opts.StartPoint = T+(T*((-1+(1+1)*rand(1))*0.08));
    opts.Lower = 0;
    opts.Upper = Inf;
    opts.TolX = 10^-8;
    opts.DiffMaxChange = 0.01;
    opts.MaxFunEvals = 750;
    opts.MaxIter = 500;

    % Fit model to data.
    [fitresult, gof, Jacobian] = fit( xData, yData, ft,'problem',{CSA, S, wI, MAS}, opts);
end
% Function for Bisquare Weighted residuals using Jacobian from MF fits
function [res_BW] = bisquare_residuals(jacobian,residuals,y,x)
    jacobian=full(jacobian);
    [Q,~]=qr(jacobian,0);
    h=min(0.9999,sum(Q.*Q,2));
    adj_fac=1./sqrt(1-h);
    res_adj=residuals.*adj_fac;

BW=zeros(length(y),length(x));
R=zeros(length(y),length(x));
res_BW=zeros(length(y),length(x));


RS = sort(abs(res_adj));
MAR = median(RS(2:end,:));
SIGMA = MAR./0.6745;
tune=4.685;
tSIGMA=tune.*SIGMA;
        for k = 1:length(x)
            for l = 1:length(y)
                if abs(residuals(l,k)) > tSIGMA(k)
                    BW(l,k) = 0;
                else
                    R(l,k) = residuals(l,k)/tSIGMA(k);
                    BW(l,k) = (abs(R(l,k))<1).*(1-(R(l,k)).^2).^2;
                end
                    res_BW(l,k) = residuals(l,k) .* BW (l,k);
            end
        end
end
% Function for Bisquare Weighted residuals using Jacobian from MF fits
function [res_BW] = bisquare_residuals_3D(jacobian,residuals,y,x,z)
    jacobian=full(jacobian);
    [Q,~]=qr(jacobian,0);
    h=min(0.9999,sum(Q.*Q,2));
    adj_fac=1./sqrt(1-h);
    res_adj=residuals.*adj_fac;

BW=zeros(length(y),length(x),length(z));
R=zeros(length(y),length(x),length(z));
res_BW=zeros(length(y),length(x),length(z));


RS = sort(abs(res_adj));
MAR = median(RS(2:end,:,:));
SIGMA = MAR./0.6745;
tune=4.685;
tSIGMA=tune.*SIGMA;
    for q = 1:length(z)
        for k = 1:length(x)
            for l = 1:length(y)
                if abs(residuals(l,k,q)) > tSIGMA(:,k,q)
                    BW(l,k,q) = 0;
                else
                    R(l,k,q) = residuals(l,k,q)/tSIGMA(:,k,q);
                    BW(l,k,q) = (abs(R(l,k,q))<1).*(1-(R(l,k,q)).^2).^2;
                end
                    res_BW(l,k,q) = residuals(l,k,q) .* BW (l,k,q);
            end
        end
    end
end
