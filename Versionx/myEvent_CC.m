function [value,isterminal,direction] = myEvent_CC(~,~,~,~,~,~)

%   Define the timeout in seconds
TimeOut = 10;

%   The solver runs until this VALUE is negative (does not change the sign)
value = toc-TimeOut;

%   The function should terminate the execution, so
isterminal = 1;

%   The direction does not matter
direction = 0;

end

