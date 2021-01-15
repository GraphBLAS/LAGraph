#-------------------------------------------------------------------------------
# LAGraph/Test2/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
# SPDX-License-Identifier: BSD-2-Clause

#-------------------------------------------------------------------------------

JOBS ?= 1

# build all the tests
default:
	( cd BFS                   ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd BetweennessCentrality ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd PageRank              ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )

clean: distclean

purge: distclean

# remove all files not in the distribution.
distclean:
	( cd BFS                   ; $(MAKE) distclean )
	( cd BetweennessCentrality ; $(MAKE) distclean )
	( cd PageRank              ; $(MAKE) distclean )

