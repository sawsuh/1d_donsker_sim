check: compile_debug run_check
time: compile run_check
compile:
	clang++ -std=c++20 -Wall -g -O3 sim.cpp -o sim
compile_debug:
	clang++  -std=c++20 -Wall -g -O3 -D _DEBUG sim.cpp -o sim
run_check:
	(time ./sim 2> ./log.txt && Rscript ./analyse.R) || cat log.txt
run:
	./sim
