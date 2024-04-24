check: compile_debug run_check
time: compile run_check
compile:
	clang++ -fopenmp -std=c++20 -Wall -g -rdynamic -pedantic -O3 sim.cpp -o sim
compile_debug:
	clang++ -fopenmp -std=c++20 -Wall -g -rdynamic -pedantic -O3 -D _DEBUG sim.cpp -o sim
run_check:
	(time ./sim 2> ./log.txt && Rscript ./analyse.R) || cat log.txt
run:
	./sim
