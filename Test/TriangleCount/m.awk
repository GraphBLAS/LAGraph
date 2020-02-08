

/Message/ {
    is_urand = match ($0, "urand")
    m = $0
    method = $4
    if (match ($0, "none")) {
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
        printf "%10.3f   %%  %s %s\n", t, method, sort

        if (match (method, "Dot2")) {
            if (match (sort, "\+")) {
                printf ("\n") ;
            }
        }
    }
}
