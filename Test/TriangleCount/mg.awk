
BEGIN {
printf "host = host + 1 ;\n" 
printf "Host {host} = hostname ;\n\n"
i = 1 ;
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
    if (is_urand) {
        if (match (method, "Dot:")) {
            if (match (sort, "\-")) {
                i = 1 ;
            }
        }
        printf "    t(%d) = %10.3f ; %%  %s %s\n", i, t, method, sort
        i = i + 1 ;

        if (match (method, "Dot2")) {
            if (match (sort, "\+")) {
                printf ("\nk = k + 1 ;\n") ;
                printf ("Results {host}{k} = t ;\n\n") ;
            }
        }
    }
}
