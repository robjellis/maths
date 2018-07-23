function [res] = irace_eval

% evaluate the error induced by sampling rates
%
% version = 2012.10.19

means = input('\n Enter the target means; "99" = 250:250:1000: ');
if means == 99
   means = 250:250:1000;
end

nmeans = numel(means);

cvs = input(' Enter the target CVs; "99" = [1 2.5 5 10]: ');
if cvs == 99
   cvs = [1 2.5 5 10];
end

ncvs = numel(cvs);

alphas = input(' Enter the alpha values for 1/f noise (e.g., [-2 -1 0]): ');
nalphas = numel(alphas);

ser = input(' Enter the ITI series length: ');
Hz = input(' Enter the sampling rate in Hz (e.g., 60): ');
iter = input(' Enter the number of iterations (e.g., 5000): ');


obsCV = nan(3,nmeans);
obsDFA = nan(3,nmeans);

for a = 1:nalphas
    for c = 1:ncvs
        for m = 1:nmeans
            [output] = sd_vs_rms(iter,ser,means(m),cvs(c),alphas(a),Hz);
            cur.CV = cvs(c);
            cur.alpha = alphas(a);
            obsCV(1,m) = output.cvsd_err_minus;
            obsCV(2,m) = output.cvsd_err_med;
            obsCV(3,m) = output.cvsd_err_plus;
            
            obsDFA(1,m) = output.dfa_err_minus;
            obsDFA(2,m) = output.dfa_err_med;
            obsDFA(3,m) = output.dfa_err_plus;
            
  
        end
        cur
        obsCV
        %obsDFA  % requires a seq. length of > 65
    end

end
        
