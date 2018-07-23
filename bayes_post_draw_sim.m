function out = bayes_post_draw_sim(a0,b0,s1,n1,s2,n2,simN,do_plot)

%
% compute the probability Pr(item2 > item1) using simulation of beta
% distribution values
%
% where a0 = hyperparameter alpha
%		b0 = hyperparameter beta
%		s_i = success rate (0 to 1) for item i
%		n_i = total number of observances for item i
%
% RJE | 8 Feb 2018

if nargin < 7
	simN = 10000;
end

if nargin < 8
	do_plot = 1;
end

tic;

% define hits (# successes) and misses (#failures)
h1 = s1*n1;
m1 = n1 - h1;

h2 = s2*n2;
m2 = n2 - h2;

% find the posterior params for item 1 and item 2
a1 = a0 + h1;
b1 = b0 + m1;

a2 = a0 + h2;
b2 = b0 + m2;

% get the posterior expected value
ev1 = a1 / (a1 + b1);
ev2 = a2 / (a2 + b2);

% pull values from beta distribution
v1 = betarnd(a1,b1,simN,1);
v2 = betarnd(a2,b2,simN,1);

% compute rate of v1 > v2
d = v2 - v1;

Pr = mean(d>0); % taking average is same as taking sum and dividing by simN

% alternative - nope, it's not intuitive
% d_sort = sort(v2) - sort(v1);
% d_std  = std(d_sort);
% 
% figure(38)
% plot(sort(v1),'b')
% hold on
% plot(sort(v2),'r')
% hold off
% 
% figure(39)
% ecdf(v1)
% hold on
% ecdf(v2)
% hold off
% 
% figure(40)
% ecdf(d_sort)


%% plots (optional)
if do_plot == 1
	% compute PDFs and CDFs for plotting purposes
	x = 0:.001:1;
	pdf0 = betapdf(x,a0,b0);
	
	pdf1 = betapdf(x,a1,b1);
	cdf1 = betacdf(x,a1,b1);
	
	pdf2 = betapdf(x,a2,b2);
	cdf2 = betacdf(x,a2,b2);
	
	figure(10)
	subplot(1,2,1)
	ecdf(v1)
	hold on
	plot(x,cdf1)
	hold off

	subplot(1,2,2)
	ecdf(v2)
	hold on
	plot(x,cdf2)
	hold off
	
	figure(12)

	thr = .001; % to make the plots nicer, don't plot PDF when it's near 0

	plot(x(pdf0>thr),pdf0(pdf0>thr),'g')
	hold on
	plot(x(pdf1>thr),pdf1(pdf1>thr),'b','LineWidth',1.2)
	plot(x(pdf2>thr),pdf2(pdf2>thr),'r')
	hold off
	xlabel('Success rate')
	ylabel('Probability density')
	title('PDF')
	
	% rescale x-axis for beter visibility
	minx = min([min(x(pdf0>thr)) min(x(pdf1>thr)) min(x(pdf2>thr))]);
	maxx = max([max(x(pdf0>thr)) max(x(pdf1>thr)) max(x(pdf2>thr))]);
	xlim([minx maxx]);
	
	% observed differences between the two distributions
	figure(13)
	ecdf(d)
	xlabel('Item 2 - Item 1 simulated scores')

end

tocc = toc;

out.Pr = Pr;
out.ev1 = ev1;
out.ev2 = ev2;
out.ev_diff = ev2 - ev1;
out.d_avg = mean(d);
out.process_sec = tocc;
%out.d_std = d_std;


