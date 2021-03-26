
BEGIN {
printf "host = host + 1 ;\n" 
printf "Host {host} = hostname ;\n\n"
i = 1 ;
j = 0 ;
}

/Message/ {
    is_urand = match ($0, "urand")
    m = $0
    method = $4
    if (match ($0, "auto")) {
        sort = "auto"
    }
    else if (match ($0, "none")) {
        sort = "none"
    }
    else if (match ($0, "ascend")) {
        sort = "+"
    }
    else if (match ($0, "descend")) {
        sort = "-"
    }
}

/time:/ {
    t = $6
    nthreads = $4
    if (is_urand) {
        if (match (method, "Dot:")) {
            j = j + 1 ;
            i = 1 ;
        }
        if (nthreads == 64) {
            printf "    t(%d,%2d) = %10.3f ; %%  %s %s\n", i, j, t, method, sort
            i = i + 1 ;
        }
    }
}

END {
    printf ("Results {host} = t ; clear t\n\n") ;
}
