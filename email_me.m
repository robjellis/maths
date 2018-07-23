function email_me(proc_name,recip,starttime,attach)

% function email_me(proc_name,recip,starttime,attach)
%
% send an email from the msg.matlab account to alert the user that a
% process ("proc_name") has finished
%
% proc_name = a text_string, e.g., 'data_simulation.m'
% recip     = the recipient's email address, e.g., 'robjellis@gmail.com'
% start     = a starting time (determined by the user); e.g., datestr(clock) is
%             the current time
% attach    = a file name, e.g., 'output.mat'
%
% code from: http://obasic.net/how-to-send-e-mail-from-matlab
% modified by rob ellis, august 2011

setpref('Internet', 'E_mail', 'msg.matlab@gmail.com');
setpref('Internet', 'SMTP_Username', 'msg.matlab@gmail.com');
setpref('Internet', 'SMTP_Password', 'matlab123');
setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port', '465');

endtime = datestr(clock); % the current time
% Send the email
subj   = 'MATLAB Process is Complete';
msg    = ['Your MATLAB process <', proc_name, '> is complete! // start time: ', starttime, ' // end time: ', endtime];


sendmail(recip, subj, msg, attach);

fprintf(['\n Email sent to: ', recip '\n']);