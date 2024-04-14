compile:
	clang++ -Wall -O3 sim.cpp -o sim
run:
	time (./sim)
full: compile run
