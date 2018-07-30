function er_senderror(to,ME,verbose)
%ER_SENDERROR(TO,ME,VERBOSE) sends error messages defined in MEexpection
% object to the mail address (specified in TO) along with the function that
% caused the error (as attachment). 
%
% - TO is either a string specifying a single address, or a cell array of
%   addresses.
% - ME (optional) is a matlab MException object.
% - VERBOSE (optional): 1/0. If verbose == 1 all stacked files will be sent
%   by mail, otherwise only the one that caused the error. 
%
% If ME is not provided, the mail will be sent with unspecified error.
%__________________________________________________________________________
%
%   Version 1.0
%   Copyright (C) 2018   danielemascali@gmail.com

%------------CONFIGURE SMTP SERVER and SENDER MAIL-------------------------
mail_service = 'matlab.communications@gmail.com';
password = 'Qwert90.';
setpref('Internet','E_mail',mail_service);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail_service);
setpref('Internet','SMTP_Password',password);
% for unix system:
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
%--------------------------------------------------------------------------

if nargin == 0
    error('Not enough input arguments. At least one mail address is required.');
end

if nargin < 3
    verbose = 0;
end

ST = dbstack('-completenames');
callingFunction = ST(end);

subject = ['<',callingFunction.name,'> FAILURE'];


if nargin == 1 || isempty(ME)
    msg_str = 'has just failed with non specified error.';
    error_str = '';
else
    msg_str = sprintf('has just failed with the following error:\n');
    error_str = getReport(ME,'extended','hyperlinks','off');
end

if ispc
    username = getenv('username');
    computername = getenv('computername');
elseif isunix
    username = getenv('USER');
    [~, computername] = system('hostname');
    computername = strtrim(computername);
end

%varius info
message_final{1} = sprintf(['Dear %s,\n\n',...
                           'matlab function: %s\n',...
                           'full path: %s\n',...
                           'on machine: %s\n',...
                           '%s'],username,callingFunction.name,callingFunction.file,computername,msg_str);
%error if available                       
message_final{2} = error_str;

%add signature
message_final{3} = sprintf(['____________________________\n',...
                    'Automatic Mail Service created by\n',...
                    'Daniele Mascali, PhD\n',...
                    '"Enrico Fermi" Centre\n',...
                    'MARBILab\n',...
                    'Email: daniele.mascali@gmail.com']);

                
if nargin == 1 || isempty(ME)    
    sendmail(to,subject,message_final)
else
    if verbose > 0
        for l = 1:length(ME.stack)
            attachments{l} = ME.stack(l).file;
        end
        sendmail(to,subject,message_final,attachments)
    else
        sendmail(to,subject,message_final,ME.stack(1).file)
    end
end

return
end