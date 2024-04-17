time_check: compile_debug run_check
time: compile run_check
check: compile_debug run
compile:
	clang++ -Wall -g -O3 sim.cpp -o sim
compile_debug:
	clang++ -Wall -g -O3 -D _DEBUG sim.cpp -o sim
run:
	./sim
run_check:
	time ./sim
	Rscript ./analyse.R
