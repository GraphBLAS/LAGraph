all:
	gcc -o main.o -lgraphblas -llagraph main.c
	gcc -o test.o -lgraphblas -llagraph test.c
	echo "Running algorithms on test graphs..."; 	    \
	for number in 1000 2003 2500; do 		    \
		./test.o "test$$number.mtx" $$number;    \
	done

run:
	./main.o $(iFile) $(oFile) $(size)

clean:
	rm -rf *.o

