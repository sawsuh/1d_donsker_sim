compile:
	clang++ -Wall -O3 sim.cpp -o sim
run:
	./sim | tee res.csv
full: compile run