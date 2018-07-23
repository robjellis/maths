function output = mc_pairedsamp(data,iter,nsamp,seed,plot_figs)

% paired-sample tests (paired T, Wilcoxon) using N/2 random subsamples
%
% updated = 2014.04.26

if iter == 0
    do_sub = 0;
    iter = 1;
elseif iter > 0
    do_sub = 1;
end

% get rid of NaN cases
is_nan = isnan(data);
all_nan = sum(is_nan,2);
data = data(all_nan==0,:);

% new: change Inf to the highest observed value
   ex_vec = data(:);
       maxd = max(ex_vec(ex_vec<Inf));
      isinf = data == Inf;
data(isinf) = maxd;

% size of D?
nsub = size(data,1);
ncol = size(data,2);


% summary plot
   if plot_figs == 1
       figure(201)
       for m = 1:ncol
           d = data(:,m);
           plot_this = 0;
           fout = ecdf_manual(d,plot_this);

           if m == 1
                plot(fout.x,fout.f,'r','LineWidth',2)
           else
                plot(fout.x,fout.f,'b')
           end
           hold on
              
       end
       hold off
       
       xlabel('Data values')
       ylabel('Cumulative Probability')
   end


%% subsample analysis
% have to do this for each PAIR independently

%Wval  = nan(iter,ncol); 
Wsig  = nan(iter,ncol); % 1 if significant, 0 if not significant
Tval  = nan(iter,ncol); 
Tsig  = nan(iter,ncol); % 1 if significant, 0 if not significant

%progressbar(0) % initialize

% take the smaller value just in case we want to explore smaller subsamples    
n = min(nsamp,floor(nsub/2));

rho = nan(1,ncol);

progressbar(0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% important: where do we want to take a new sample of subjects? Inside or
% outside the SF loop (columns)?

% set the random stream value; only have to do this once
s = RandStream('mt19937ar','Seed',seed);

for i = 1:iter
    %progressbar([],i/iter)
    progressbar(i/iter)
        
    if do_sub == 0
        ind = 1:nsub; i1 = i1(:);
    elseif do_sub == 1
        % permute the full set
        rp = randperm(s,nsub);        
        ind = rp(1:n); % for CO
    end
    
    for m = 2:ncol
    %progressbar(m/ncol,[])
        
    % Spearman correlation for the full set of values
    if i == 1
        rho(m) = corr(data(:,1), data(:,m),'type','Spearman');      
    end

    % ************************************
    % get the subsample
    p1 = data(ind,1);
    p2 = data(ind,m);


    % paired T test: ttest(B,A) will be positive if mean B > mean A
    [h xxx xxx stats] = ttest(p2,1); 
    %Tsub(i,m) = stats.tstat; % keep the sign so that we can tell about direction of inflation of Ex relative to Co!
    Tval(i,m) = abs(stats.tstat); % abs to make things simpler
    Tsig(i,m) = h;


    % Wilcoxon test
    wout = signrank_rje(p2,p1); % two sample test, so order doesn't matter
    Wsig(i,m) = wout.sig;   

    end % column loop
    
end % iter loop

% now get the observed FPR as a PROPORTION
Tfpr = sum(Tsig) / iter;
Wfpr = sum(Wsig) / iter;

    %% output

    output.n      = n;
    output.iter   = iter;
    output.Tsig = Tsig;
    output.Wsig = Wsig;
    output.Tfpr = Tfpr;
    output.Wfpr = Wfpr;
    output.rho = rho;    
    
end % do_sub 

