compile:
	clang++ -Wall -O3 sim.cpp -o sim
run:
	./sim | tee res.csv
	R 
full: compile run
