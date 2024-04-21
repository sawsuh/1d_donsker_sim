check: compile_debug run_check
time: compile run_check
compile:
	clang++ -fopenmp -std=c++20 -Wall -g -O3 sim.cpp -o sim
compile_debug:
	clang++ -fopenmp -std=c++20 -Wall -g -O3 -D _DEBUG sim.cpp -o sim
run_check:
	export OMP_NUM_THREADS=16
	time ./sim 2> ./log.txt
	Rscript ./analyse.R
run:
	export OMP_NUM_THREADS=16
	./sim
