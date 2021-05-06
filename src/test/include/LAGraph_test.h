//-----------------------------------------------------------------------------
// LAGraph/src/test/include/LAGraph_test.h:  defintions for testing LAGraph
//-----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#ifndef LAGRAPH_TEST_H
#define LAGRAPH_TEST_H

#include <LAGraph.h>
#include <acutest.h>
#include <graph_zachary_karate.h>

// The tests are compiled with -DLGDIR=/home/me/LAGraph, or whatever, where
// /home/me/LAGraph is the top-level LAGraph source directory.
#define LG_XSTR(x) LG_STR(x)
#define LG_STR(x) #x
#define LG_SOURCE_DIR LG_XSTR (LGDIR)

#define LG_DATA_DIR LG_SOURCE_DIR "/data/"

#define OK(method) TEST_CHECK (method == 0)

#endif

