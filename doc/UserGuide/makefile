pdf: headers
	pdflatex LAGraph_UserGuide.tex
	pdflatex LAGraph_UserGuide.tex
	pdflatex LAGraph_UserGuide.tex

headers:
	python func_extract.py ../../include/LAGraph.h  > LAGraph_h.tex
	python func_extract.py ../../include/LAGraphX.h > LAGraphX_h.tex

clean:
	rm -f *.toc *.log *.aux *.out *.lof *.lot

cleanall: clean
	rm -f LAGraph_UserGuide.pdf
