function progress(n,first,last,tint)
%This is a periodic update function created to provide occasional updates
%on the progress of an iterative process.  Inputs are "n" - the iteration
%number, "first" - the first iteration value, "last" - the last iteration
%value, and "tint" - the time interval in seconds upon which the program is
%to give an update on its progress.  The function outputs this display 
%every 'tint' seconds: "...Processing: 'n' on 'first'-'last': 'percentage'%"
%
%If this function is to be placed into a loop of unknown value "last", input
%"last" as infinity(Inf) and program will display: 
%"...Processing: 'n' on 'first'-unknown".
%
%Written by: Martin Kearney-Fischer
%       3/12/2007

global time1

if isempty(time1)||(n==first)
    time1 = clock;
end
if abs(etime(time1,clock))>tint
    time1 = clock;
    if last ~= Inf
        if n-1==first
            percent = (n-first)/(last-first)*100;
        else
            percent = (n-1-first)/(last-first)*100;
        end
        disp(['...Processing: ' num2str(n) ' on ' num2str(first) '-' num2str(last) ':  ' num2str(percent,'%.2f') '%'])
    else
        disp(['...Processing: ' num2str(n) ' on ' num2str(first) '- unkown'])
    end
end