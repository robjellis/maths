function plot_index(x,y,fignum)

% return the index associated with a specific point (x,y) in vector [X,Y]
%
% RJE + John CAI | 2013.03.19


fig = figure(fignum);

plot(x,y,'.'); % rje function

datacursormode on;
h = datacursormode(fig);
set(h,'UpdateFcn',@myupdatefcn);

function [indtxt] = myupdatefcn(obj,event_obj)
    pos = get(event_obj,'Position');

    index = get(event_obj,'DataIndex');
    
    indtxt = {['Ind: ' num2str(index)];...
              ['  X: ' num2str(pos(1))];...
              ['  Y: ' num2str(pos(2))]};
end

%%%%%%%%%%%% Extension Point %%%%%%%%%%%%%%%%%%%%%%
%         DO SOMETHING HERE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


