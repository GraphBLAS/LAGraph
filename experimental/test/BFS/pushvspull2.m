function pushvspull (name, allpush_results, allpull_results, pushpull_results, n, nz, fig)
% pushvspull ('kron', allpush_results, allpull_results, n, nz, 1) ;

d = nz / n  ;
ntrials = length (allpush_results) ;
format long g
figure (fig)
clf (fig)

fprintf ('--------------- %s:\n', name) ;

t_pull_tot = 0 ;
t_push_tot = 0 ;
t_auto_tot = 0 ;
t_best_tot = 0 ;
t_pushpull_tot = 0 ;

for k = 1:ntrials

    % fprintf ('%d:\n', k) ;

    push_result = allpush_results {k}  ;
    pull_result = allpull_results {k}  ;
    pushpull_result = pushpull_results {k}  ;

    edges_unexplored = nz ;
    last_nq = 0 ;
    do_push = true ;

    nq = push_result (:,2) ;
    nvisited = push_result (:,3) ;
    edges_in_frontier = pushpull_result (:,4) ;
    t_pull = pull_result (:,5) ;
    t_push = push_result (:,5) ;
    t_pushpull = pushpull_result (:,5) ;

    reltime = t_push ./ t_pull ;

    tt = [t_pull t_push]  ;
    t_best = min (tt, [ ], 2) ;
    nlevels = size (push_result, 1) ;
    t_auto = nan (nlevels, 1) ;
    edges_unexploreds = nan (nlevels, 1) ;
    growings = nan (nlevels, 1) ;

    % good for kron, urand, twitter:
    alpha = 12 ;
    beta1 = 8 ;
    beta2 = 100 ;

    % good for orkut, web.
    % seems OK for the other matrices too
    alpha = 8 ;
    beta1 = 8 ;
    beta2 = 500 ;     % better for web
    % beta2 = 100 ;

    npull = 0 ;
    npush = 0 ;

    any_pull = false ;
    for level = 1:nlevels
        % select push or pull
        this_nq = nq (level) ;

        growing = this_nq > last_nq ;
        shrinking = this_nq < last_nq ;
        growings (level) = growing ;

        if (do_push)
            if (any_pull)
                switch_to_pull = (growing && this_nq > n / beta1) ;
            elseif (edges_unexplored < n)
                % turn off the pull
                fprintf ('turn off pull\n') ;
            else
                edges_unexplored = edges_unexplored - edges_in_frontier (level);
                edges_unexploreds (level) = edges_unexplored ;
                big_frontier =  ...
                    alpha * edges_in_frontier (level) > edges_unexplored  ;
                switch_to_pull = (growing && big_frontier) ;
            end
%           big_frontier = ...
%               ((edges_in_frontier (level)) ./ edges_unexplored) > alpha ;
            if (switch_to_pull && growing)
                do_push = false ;
            end
        else
            if ((this_nq < n / beta2) && shrinking)
                do_push = true ;
            end
        end

        if (do_push)
            npush = npush + 1 ;
            t_auto (level) = t_push (level);
        else
            any_pull = true ;
            npull = npull + 1 ;
            t_auto (level) = t_pull (level);
        end
        bad_choice = (t_auto (level) > t_best (level)) ;

if (0)
        if (bad_choice)
            fprintf ('   level %2d: push %10.4f pull %10.3f ', ...
                level, t_push (level), t_pull (level)) ;
            if (do_push)
                fprintf ('      ');
            else
                fprintf (' pull ');
            end
            fprintf (' %15.4f ', edges_unexplored / edges_in_frontier (level)) ;
            fprintf (' %15.4f ', 1/(this_nq / n)) ;
            if (bad_choice)
                fprintf (' (oops %d)', growing) ;
            end
            fprintf ('\n') ;
        end
end

        last_nq = this_nq ;
    end

if (0)
    fprintf ('trial %2d : ', k) ;
    fprintf ('npush %4d npull %4d ', npush, npull) ;
    fprintf ('push %10.4f ', sum (t_push)) ;
    fprintf ('pull %10.4f ', sum (t_pull)) ;
    fprintf ('auto %10.4f ', sum (t_auto)) ;
    fprintf ('actual push/pull %10.4f ', sum (t_pushpull)) ;
    fprintf ('(%10.3f) ', sum (t_auto) / sum (t_pushpull)) ;
    fprintf (': best %10.4f ', sum (t_best)) ;
    fprintf ('(%10.3f) ', sum (t_auto) / sum (t_push)) ;
    fprintf ('(%10.3f)\n', sum (t_auto) / sum (t_best)) ;
end

    t_pull_tot = t_pull_tot + sum (t_pull) ;
    t_push_tot = t_push_tot + sum (t_push) ;
    t_auto_tot = t_auto_tot + sum (t_auto) ;
    t_best_tot = t_best_tot + sum (t_best) ;
    t_pushpull_tot = t_pushpull_tot + sum (t_pushpull) ;

    relalpha =   edges_unexploreds ./ (edges_in_frontier) ;
    relnq = nq / n ;

    % relalpha
    % reltime
    % growings

    subplot (1,3,1)
    loglog ( relalpha (growings ~= 0), reltime (growings ~= 0), 'o', ...
             [1e-4 200], [1 1], 'g-', ...
             [alpha alpha], [1e-3 10], 'g-') ;
    ylabel ('pushtime / pulltime') ;
    xlabel ('growing: (edges unexplored - edges in frontier)/(edges in frontier) ') ;
    hold on
    title (name) ;

    subplot (1,3,2)
    loglog ( relnq (growings ~= 0), reltime (growings ~= 0), 'o', ...
             [1e-6 1], [1 1], 'g-', ...
             [1/beta1 1/beta1], [1e-3 10], 'g-') ;
    ylabel ('pushtime / pulltime') ;
    xlabel ('growing: nq/n') ;
    hold on
    title (name) ;

    subplot (1,3,3)
    loglog ( relnq (growings == 0), reltime (growings == 0), 'o', ...
             [1e-6 1], [1 1], 'g-', ...
             [(1/beta2) (1/beta2)], [1e-3 10], 'g-') ;
    ylabel ('pushtime / pulltime') ;
    xlabel ('shrinking: nq/n') ;
    title (name) ;
    hold on
    drawnow

%{
    subplot (2,2,2)
    loglog ( pushwork, t_push, 'o') ;
    ylabel ('push time') ;
    xlabel ('push work') ;
    hold on

    subplot (2,2,3)
    loglog ( pullwork, t_pull, 'o') ;
    ylabel ('pull time') ;
    xlabel ('pull work') ;
    hold on
%}

%   subplot (2,2,2)
%   semilogy ( nvisited/n, reltime, 'o') ;
%   ylabel ('nvisited/n') ;
%   ylabel ('pushtime / pulltime') ;

%   subplot (2,2,3)
%   x = (nvisited -nq) / n; % ./(nvisited-nq) ;
%   semilogy ( x, t_push, 'o', x, t_pull, 'rx') ;

%   subplot (1,2,2)
%   loglog ( ...
%       relwork (:), t_push (:), 'o', ...
%       relwork (:), t_pull (:), 'rx', ...
%            [1 1], [1e-3 10], 'g-') ;

end
hold off

    t_pull_tot = t_pull_tot / ntrials ;
    t_push_tot = t_push_tot / ntrials ;
    t_auto_tot = t_auto_tot / ntrials ;
    t_best_tot = t_best_tot / ntrials ;
    t_pushpull_tot = t_pushpull_tot / ntrials ;

% fprintf ('\n-------\n') ;
fprintf ('all trials: allpush: %10.3f allpull: %10.3f auto: %10.3f pushpull %10.3f best: %10.3f (%10.3f)\n', ...
    t_push_tot, t_pull_tot, t_auto_tot, t_pushpull_tot, t_best_tot, t_pushpull_tot / t_best_tot) ;

