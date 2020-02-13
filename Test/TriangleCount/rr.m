
% get the latest results
clear
% system ('./mget > res.m') ;
res
% whos
nhosts = length (Results) ;

method = {
    'Dot:- '
    'Dot:0 '
    'Dot:+ '
    'Dot2:-'
    'Dot2:0'
    'Dot2:+' } ;

% find the fastest run time for each of the 6 algorithms
tbest = inf (1,6) ;
for host = 1:nhosts
    hostname = Host {host} ;
    R = Results {host} ;
    nresults = length (R) ;
    % fprintf ('%s: %d\n', hostname, nresults) ;
    for k = 1:nresults
        t = Results {host}{k} ;
        tbest = min (tbest, t) ;
    end
end

fprintf ('\nbest run times:\n') ;
for m = 1:6
    fprintf ('%s : %g\n', method {m}, tbest (m)) ;
end
fprintf ('\n\n') ;

% convert results to relative times
quality = inf (nhosts,1) ;
minq = inf (nhosts,1) ;
maxq = zeros (nhosts,1) ;

for host = 1:nhosts
    hostname = Host {host} ;
    R = Results {host} ;
    nresults = length (R) ;
%   fprintf ('%s: %d\n', hostname, nresults) ;
    tmax = zeros (1,6) ;
    tmin = inf   (1,6) ;
    ttot = zeros (1,6) ;
    for k = 1:nresults
        t = Results {host}{k} ;
        t = (t ./ tbest - 1) * 100 ;
        Results {host}{k} = t ;
%       fprintf ('      ') ;
%       for m = 1:6
%           if (t (m) == 0)
%               fprintf ('  best ') ;
%           else
%               fprintf (' %6.0f', t (m)) ;
%           end
%       end
%       fprintf ('\n') ;
        ttot = ttot + t ;
        tmin = min (tmin, t) ;
        tmax = max (tmax, t) ;
    end
    if (nresults > 0)
        ttot = ttot / nresults ;
%       fprintf ('  avg:') ;
%       for m = 1:6
%           fprintf (' %6.1f', ttot (m)) ;
%       end
%       fprintf ('\n') ;
        quality (host) = sum (ttot) / 6;
        minq (host) = min (tmin) ;
        maxq (host) = max (tmax) ;
    end
%   fprintf ('\n') ;
end

% convert results to relative times
[~, hosts] = sort (quality, 'ascend') ;

for host = hosts' % 1:nhosts
    hostname = Host {host} ;
    R = Results {host} ;
    nresults = length (R) ;
    fprintf ('%s: %d\n', hostname, nresults) ;
    tmax = zeros (1,6) ;
    tmin = inf   (1,6) ;
    ttot = zeros (1,6) ;
    for k = 1:nresults
        t = Results {host}{k} ;
%       t = (t ./ tbest - 1) * 100 ;
%       Results {host}{k} = t ;
        fprintf ('      ') ;
        for m = 1:6
            if (t (m) == 0)
                fprintf ('  best ') ;
            else
                fprintf (' %6.0f', t (m)) ;
            end
        end
        fprintf ('\n') ;
        ttot = ttot + t ;
        tmin = min (tmin, t) ;
        tmax = max (tmax, t) ;
    end
    if (nresults > 0)
        ttot = ttot / nresults ;
        fprintf ('  min:') ;
        for m = 1:6
            fprintf (' %6.1f', tmin (m)) ;
        end
        fprintf (' :   overall: %6.1f\n', min (tmin)) ;

        fprintf ('  avg:') ;
        for m = 1:6
            fprintf (' %6.1f', ttot (m)) ;
        end
        fprintf (' :   overall: %6.1f\n', sum (ttot)/6) ;

        fprintf ('  max:') ;
        for m = 1:6
            fprintf (' %6.1f', tmax (m)) ;
        end
        fprintf (' :   overall: %6.1f\n', max (tmax)) ;

    end
    fprintf ('\n') ;
end


nhosts = length (find (isfinite (quality)))
q = quality (hosts) ;
qmin = minq (hosts) ;
qmax = maxq (hosts) ;
q = q (1:nhosts) ;
qmin = qmin (1:nhosts) ;
qmax = qmax (1:nhosts) ;

figure (1)
errorbar (1:nhosts, q, qmin-q, qmax-q)
xlabel ('host (sorted by quality lo to hi)') ;
ylabel ('run time / best run time (lower is better)') ;


