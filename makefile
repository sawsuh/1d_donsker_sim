compile:
	clang++ -Wall -g -O3 sim.cpp -o sim
run:
	time ./sim
full: compile run
