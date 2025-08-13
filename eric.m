j% MATLAB code for generating Figure 3

% Written by Eric G. Keeler
% July 2022
% ekeeler@nysbc.org; aem5@columbia.edu
% Cite as: E.G.Keeler, A.E.McDermott, ChemRev 2022

% August 2025
% Adapted by Simon G. Hulse for checking consistency with RING NMR
% Load the Experimental DMS data taken from Quinn and McDermott, JMR 2012
load('DMS_exp_val.mat')

wI=99.8e6*2*pi;
wR=10e3;
delta=-37e-6;
ssq=0.3295;

powers = {DMS_irr_power_37 DMS_irr_power_57 DMS_irr_power_67 DMS_irr_power_77};
rates = {DMS_rate_37 DMS_rate_57 DMS_rate_67 DMS_rate_77};
inital_taus = {1e-6 1e-5 1e-5 6e-5};
w1 = 2 * pi * 1e3 * linspace(1.49, 49.02, 98);

taus = zeros(4);
curves = cell(1, 4);

for i = 1:4
    power = 2 * pi * 1e3 * powers{i};
    rate = rates{i} - 1.8;
    initial_tau = inital_taus{i};

    [fitresult, ~, ~] = Kurbanov_fit_v3_fixedS(power, rate, delta, wR, wI, initial_tau, ssq);
    coeffs = coeffvalues(fitresult);
    tau = coeffs(1);
    taus(i) = tau;

    curve = (3/4).*(delta.*wI)^2.*((2/5)*...
        (1-ssq).*(tau./(1+((wI).*tau).^2)))+(((1/6)*...
        (delta*wI)^2.*(((1/2).*(2/5)*(1-ssq).*...
        (tau./(1+((w1-2.*wR*2*pi).*tau).^2)))+((2/5)*(1-ssq).*...
        (tau./(1+((w1-wR*2*pi).*tau).^2)))+((2/5)*(1-ssq).*...
        (tau./(1+((w1+wR*2*pi).*tau).^2)))+((1/2).*(2/5)*(1-ssq).*...
        (tau./(1.+((w1+2*wR*2*pi).*tau).^2)))))-0.5*((3/4)*...
        (delta*wI)^2.*((2/5)*(1-ssq).*...
        (tau./(1+((wI).*tau).^2)))));
    curves{i} = curve;
end

tau37 = taus(1); tau57 = taus(2); tau67 = taus(3); tau77 = taus(4);
curve37 = curves{1}; curve57 = curves{2}; curve67 = curves{3}; curve77 = curves{4};
save("data/eric.mat", "tau37", "tau57", "tau67", "tau77", "w1", "curve37", "curve57", "curve67", "curve77");

function [fitresult, gof, Jacobian] = Kurbanov_fit_v3_fixedS(omega1,fid,CSA,MAS,wI,T,S)
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
