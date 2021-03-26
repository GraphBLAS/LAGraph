
% get the latest results
clear
% system ('./mget > res.m') ;
res
% whos
nhosts = length (Results) ;

method = {
    'Dot:auto '
    'Dot2:auto'
    } ;

nalgo = 2 ;

% find the fastest run time for each of the algorithms
tbest = inf (nalgo,1) ;
for host = 1:nhosts
    t = Results {host}  ;
    tmin = min (t, [ ], 2) ;
    tbest = min (tbest, tmin) ;
end

% find the fastest run time for each of the algorithms
for host = 1:nhosts
    hostname = Host {host} ;
    hostid = hostname (5:13) ;
    t = Results {host}  ;
    tmin = min (t, [ ], 2) ;
    tmax = max (t, [ ], 2) ;
    tmean = mean (t, 2) ;
    tmedian = median (t, 2) ;
    fprintf ('host: %s ', hostid) ;
    fprintf ('min: %8.3f %8.3f ', tmin ./ tbest) ;
    fprintf ('mean: %8.3f %8.3f ', tmean ./ tbest) ;
    fprintf ('med: %8.3f %8.3f ', tmedian ./ tbest) ;
    fprintf ('max %8.3f %8.3f\n', tmax ./ tbest) ;
end

fprintf ('\nbest run times:\n') ;
for m = 1:nalgo
    fprintf ('%s : %g\n', method {m}, tbest (m)) ;
end
fprintf ('\n\n') ;

