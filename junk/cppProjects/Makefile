lookup: testLUT.o lookupTables.o
	g++  testLUT.o lookupTables.o -o lookup
	
lookupTables.o: lookupTables.C lookupTables.H
	g++ -c lookupTables.C
	
testLUT.o: lookupTables.H
	g++ -c testLUT.C

clean:
	rm lookup *.o
