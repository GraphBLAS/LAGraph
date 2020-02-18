
for node = 3 % 1:234
    filename = sprintf ('myjob_TC_urand_n%03db', node) ;
    f = fopen (filename, 'w') ;
    fprintf (f, '#PBS -l nodes=s001-n%03d\n', node) ;
    fprintf (f, '#PBS -l walltime=24:00:00\n') ;
    fprintf (f, 'cd $PBS_O_WORKDIR\n') ;
    fprintf (f, 'cat /proc/cpuinfo | sort | uniq -c\n') ;
    fprintf (f, 'cat /proc/meminfo\n') ;
    fprintf (f, 'numastat -v\n') ;
    fprintf (f, 'numactl -s\n') ;
    fprintf (f, 'export OMP_PLACES=cores\n') ;
    fprintf (f, 'export OMP_PROC_BIND=spread\n') ;
    fprintf (f, 'export OMP_DISPLAY_AFFINITY=true\n') ;
    fprintf (f, 'export OMP_AFFINITY_FORMAT="Thrd Lev=%%3L, thrd num=%%5n, thrd aff=%%15A"\n') ;
    for k = 1:8
        fprintf (f, './build/gap_tc_test ~/GAP/GAP-urand/GAP-urand.grb\n') ;
    end
    fclose (f) ;
end

f = fopen ('runb_001_to_100', 'w') ;
for node = 1:100
    filename = sprintf ('myjob_TC_urand_n%03db', node) ;
    fprintf (f, 'qsub %s\n', filename) ;
end

f = fopen ('runb_101_to_200', 'w') ;
for node = 101:200
    filename = sprintf ('myjob_TC_urand_n%03db', node) ;
    fprintf (f, 'qsub %s\n', filename) ;
end

f = fopen ('runb_201_to_234', 'w') ;
for node = 201:234
    filename = sprintf ('myjob_TC_urand_n%03db', node) ;
    fprintf (f, 'qsub %s\n', filename) ;
end

