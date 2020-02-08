
for node = 1:234
    filename = sprintf ('myjob_TC_urand_n%03d', node) ;
    f = fopen (filename, 'w') ;
    fprintf (f, '#PBS -l nodes=s001-n%03d\n', node) ;
    fprintf (f, '#PBS -l walltime=24:00:00\n') ;
    fprintf (f, 'cd $PBS_O_WORKDIR\n') ;
    fprintf (f, 'cat /proc/cpuinfo | sort | uniq -c\n') ;
    fprintf (f, 'cat /proc/meminfo\n') ;
    fprintf (f, 'numastat -v\n') ;
    fprintf (f, 'numactl -l\n') ;
    for k = 1:8
        fprintf (f, './build/gap_tc_test ~/GAP/GAP-urand/GAP-urand.grb\n') ;
    end
    fclose (f) ;
end

