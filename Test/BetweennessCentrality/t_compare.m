
clear
t_nomask ;
t_allmask ;
whos

figure (1)
clear figure

% phase1:  C<!M>=A*B
% phase2:  C<M>=A*B

ss = 1 ;

for id = 1:5
    matrix = Matrix {id} ;
    fprintf ('\nmatrix: %s\n', matrix) ;


    for phase = 1:2

        tno   = NoMask_Result {id}{phase} ;
        tno_tot   = tno (:,1) ;     % total flops
        tno_mwork = tno (:,2) ;
        tno_axb   = tno (:,3) ;
        tno_time  = tno (:,4) ;
        k = size (tno,1) ;

        % mask disabled so mwork must be zero
        assert (sum (tno_mwork) == 0) ;

        tall = Mask_Result {id}{phase} ;
        tall_tot   = tall (:,1) ;
        tall_mwork = tall (:,2) ;
        tall_axb   = tall (:,3) ;
        tall_time  = tall (:,4) ;
        assert (k == size (tall ,1)) ;

        tall_mwork = max (tall_mwork,1) ;

        % x axis: (tall_axb ./ tall_mwork)

        % y axis: log10 (tno_time ./ tall_time) ;

        subplot (5,2, ss) ;
        loglog (tall_axb ./ tall_mwork, tno_time ./ tall_time, 'o', ...
                [1e-8 1e6], [1 1], 'r-', ...
                [1 1], [1e-4 1e1], 'k-' ) ;
        ss = ss + 1 ;

        title (sprintf ('%s:%d', matrix, phase)) ;

        tbest = min (tno_time, tall_time) ;

        if (phase == 1)
            tol = 0.01 ;
        else
            tol = 0.01 ;
        end

        % heuristic:
        tguess = zeros (k,1) ;
        for i = 1:k
            mw  = tall_mwork (i) ;
            axb = tall_axb (i) ;
            % if ( axb / mw > 0.01)
            if (axb > mw * tol)
                % use the mask
                tguess (i) = tall_time (i) ;
            else
                % do not use the mask
                tguess (i) = tno_time (i) ;
            end
%           if (id < 5)
%               fprintf ('       %4d: axb/mw %7.2e   tno %7.2e tmask %7.2e tno/tmask %7.2g : guess/best %7.2g\n', ...
%               i, axb/mw, tno_time(i), tall_time(i), tno_time(i)/tall_time(i), tguess(i)/tbest(i)) ;
%           end
        end

        fprintf ('\n  phase: %d\n', phase) ;
        fprintf ('  tno   time: %10.4g\n', sum (tno_time)) ;
        fprintf ('  tall  time: %10.4g\n', sum (tall_time)) ;
        fprintf ('  best  time: %10.4g\n', sum (tbest)) ;
        fprintf ('  guess time: %10.4g  : %10.4g\n', sum (tguess), sum(tguess)/sum(tbest)) ;

    end

end
    
