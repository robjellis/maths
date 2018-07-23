function tm = trimean(x)

% take the trimean of x

q1 = prctile(x,25);
q2 = prctile(x,50);
q3 = prctile(x,75);

tm = (q1 + 2*q2 + q3) / 4;