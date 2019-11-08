//
// Created by Scott Kolodziej on 11/7/19.
//

#include "sssp_test.h"
#include <stdio.h>

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A);                  \
    GrB_free (&Abool);              \
    GrB_free (&path_lengths);       \
}

int main (int argc, char **argv)
{
    FILE *results_file;
    FILE *check_file;
    char* results_filename = argv[1];
    char* check_filename = argv[2];

    results_file = fopen(results_filename, "r");
    if (results_file == NULL)
    {
        fprintf(stderr, "Could not open results file %s", results_filename);
        return EXIT_FAILURE;
    }

    check_file = fopen(check_filename, "r");
    if (check_file == NULL)
    {
        fprintf(stderr, "Could not open check file %s", check_filename);
        fclose(results_file);
        return EXIT_FAILURE;
    }

    bool tests_pass = true;
    double x_result;
    double x_check;

    while (fscanf(results_file, "%lf\n", &x_result) != EOF &&
           fscanf(check_file, "%lf\n", &x_check) != EOF)
    {
        tests_pass &= (x_result == x_check);
    }

    fclose(results_file);
    fclose(check_file);

    fprintf (stderr,
             "------------------------------------------------------------\n\n");
    fprintf (stderr, "sssp_test: ");
    if (tests_pass)
    {
        fprintf (stderr, "all tests passed\n");
        printf("all tests passed\n");
    }
    else
    {
        fprintf (stderr, "TEST FAILURE\n");
        printf("TEST FAILURE\n");
    }
    fprintf (stderr,
             "------------------------------------------------------------\n\n");
}