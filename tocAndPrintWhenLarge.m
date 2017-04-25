function tocTime = tocAndPrintWhenLarge(minSeconds2print,processDescriptor,tocTime)
% Similar to the toc function of matlab, but displays the passed time only
% when it is larger than minSeconds2print, in seconds.
% Input:
% minSeconds2print - The number of seconds above which the elapsed time
% will be printed.
% processDescriptor - A description that will be printed in case the time
% is printed.
% givenTocTime - If given, this will be considered as elapsed time and not
% the actual toc time.
% 
% Output:
% tocTime - The elapsed time, in seconds.
if nargin<3
    tocTime = toc;
end
if nargin<2
    processDescriptor = [];
else
    processDescriptor = [': ' processDescriptor];
end
% if nargin>2 && nargout>0
%     warning('The output toc time is the one given by the user')
% end

if tocTime>minSeconds2print
    if tocTime<60
        tocString = sprintf('%.2f [SS]',tocTime);
    elseif tocTime<3600
        tocString = sprintf('%d:%02.f [MM:SS]',floor(tocTime/60),round(rem(tocTime,60)));
    else
        tocString = sprintf('%d:%02.f:%02.f [HH:MM:SS]',floor(tocTime/3600),floor(rem(tocTime/60,60)),round(rem(rem(tocTime,3600),60)));
    end
    display(sprintf('%s passed%s',tocString,processDescriptor))
end