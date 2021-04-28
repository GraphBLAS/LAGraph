//------------------------------------------------------------------------------
// LAGraphs.h: include file for LAGraph experimental code
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#ifndef LAGRAPHX_H
#define LAGRAPHX_H


//****************************************************************************
// ascii header prepended to all *.grb files
#define LAGRAPH_BIN_HEADER 512

// LAGraph_BinRead: read a matrix from a binary file
int LAGraph_BinRead         // returns 0 if successful, -1 if failure
(
    GrB_Matrix *A,          // matrix to read from the file
    GrB_Type   *A_type,     // type of the scalar stored in A
    FILE       *f,          // file to read it from, already open
    char       *msg
) ;

#endif
