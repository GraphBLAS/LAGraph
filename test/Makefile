#-------------------------------------------------------------------------------
# LAGraph/test/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause
#
# See additional acknowledgments in the LICENSE file,
# or contact permission@sei.cmu.edu for the full terms.

#-------------------------------------------------------------------------------

JOBS ?= 1

# build all the tests
default:
	( cd BreadthFirstSearch    ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd BetweennessCentrality ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd PageRank              ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd TriangleCount         ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd ConnectedComponents   ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd SingleSourceShortestPath ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )

clean: distclean

purge: distclean

# remove all files not in the distribution.
distclean:
	( cd BreadthFirstSearch    ; $(MAKE) distclean )
	( cd BetweennessCentrality ; $(MAKE) distclean )
	( cd PageRank              ; $(MAKE) distclean )
	( cd TriangleCount         ; $(MAKE) distclean )
	( cd ConnectedComponents   ; $(MAKE) distclean )
	( cd SingleSourceShortestPath ; $(MAKE) distclean )

