function [coefs] = DRPlookup(n,m,h)
%Looks up coefficients for 4th Order, 1st Derivative DRP scheme
    a{33} = [-0.02084314277031176 ...
                0.166705904414580469 ...
                -0.77088238051822552 ...
                0 ...
                0.77088238051822552 ...
                -0.166705904414580469 ...
                0.02084314277031176]';
            
    a{42} = [0.02636943100 ...
                -0.16613853300 ...
                0.518484526d0 ...
                -1.27327473700 ...
                0.47476091400 ...
                0.46884035700 ...
                -0.049041958d0]';
        
    a{51} = [-0.048230454 ...
                0.281814650 ...
                -0.768949766 ...
                1.388928322 ...
                -2.147776050 ...
                1.084875676 ...
                0.209337622]';
            
    a{60} = [0.203876371 ...
                -1.128328861 ...
                2.833498741 ...
                -4.461567104 ...
                5.108851915 ...
                -4.748611401 ...
                2.192280339]';
    
    a{24} = -flipud(a{42});
    a{15} = -flipud(a{51});
    a{06} = -flipud(a{60});
            
    coefs = a{str2double(strcat([num2str(n),num2str(m)]))}/h;
end