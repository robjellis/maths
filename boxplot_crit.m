function boxplot_crit(data,critu)

% boxplot with critical values as determined by the user
% will plot each column in DATA as a separate box
       

       ncol = size(data,2);
       xrange = [0 ncol+1];
       
       figure(90)
       
       if nargin == 2
           critu = [critu critu];
           critl = -1*critu;
           plot(xrange,critu,'g','LineWidth',2)

           hold on
           plot(xrange,critl,'g','LineWidth',2)
           
       end
       
       boxplot(data,'notch','off','whisker',1.5,'symbol','r') % whisker = 1.5 is standard value for matlab
           
       hold off