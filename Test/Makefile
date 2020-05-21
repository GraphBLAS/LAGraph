#-------------------------------------------------------------------------------
# LAGraph/Test/Makefile
#-------------------------------------------------------------------------------

# LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
# http://graphblas.org  See LAGraph/LICENSE for license.

#-------------------------------------------------------------------------------

# simple Makefile for LAGraph, relies on cmake to do the actual build.  Use
# the CMAKE_OPTIONS argument to this Makefile to pass options to cmake.

# Install GraphBLAS before trying to compile LAGraph.

JOBS ?= 1

# build all the tests
default:
	( cd Complex               ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd MXM                   ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd LCC                   ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd AllKTruss             ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd KTruss                ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd TriangleCount         ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd PageRank              ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd PageRank2             ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd PageRank3             ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd BFS                   ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd BetweennessCentrality ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd CC                    ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd CDLP                  ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd MatrixMarket          ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd BuildMatrix           ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd BinRead               ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd DenseRelabel          ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd SSSP                  ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )
	( cd MatrixExtractKeepDimensions ; $(MAKE) --jobs=$(JOBS) CC=$(CC) CXX=$(CC) )

clean: distclean

purge: distclean

# remove all files not in the distribution
distclean:
	( cd MXM                   ; $(MAKE) distclean )
	( cd LCC                   ; $(MAKE) distclean )
	( cd AllKTruss             ; $(MAKE) distclean )
	( cd KTruss                ; $(MAKE) distclean )
	( cd TriangleCount         ; $(MAKE) distclean )
	( cd PageRank              ; $(MAKE) distclean )
	( cd PageRank2             ; $(MAKE) distclean )
	( cd PageRank3             ; $(MAKE) distclean )
	( cd BFS                   ; $(MAKE) distclean )
	( cd BetweennessCentrality ; $(MAKE) distclean )
	( cd CC                    ; $(MAKE) distclean )
	( cd CDLP                  ; $(MAKE) distclean )
	( cd MatrixMarket          ; $(MAKE) distclean )
	( cd BuildMatrix           ; $(MAKE) distclean )
	( cd BinRead               ; $(MAKE) distclean )
	( cd DenseRelabel          ; $(MAKE) distclean )
	( cd SSSP                  ; $(MAKE) distclean )
	( cd MatrixExtractKeepDimensions ; $(MAKE) distclean )
