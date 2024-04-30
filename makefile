SHELL := /bin/bash
check: compile_debug run
time: compile run
compile:
	clang++ -fopenmp -std=c++20 -Wall -g -rdynamic -pedantic -O3 sim.cpp -o sim
compile_debug:
	clang++ -fopenmp -std=c++20 -Wall -g -rdynamic -pedantic -O3 -D _DEBUG sim.cpp -o sim
run:
	(time ./sim 2> ./log.txt && Rscript ./analyse.R) || cat log.txt
figs:
	for i in {1..5}; \
	do \
		echo "working on $${i}" && \
		clang++ -fopenmp -std=c++20 -Wall -g -rdynamic -pedantic -O3 -D "_FIG$$i" -D _DEBUG sim.cpp -o sim && \
		./sim && \
		Rscript ./analyse.R --vanilla "_FIG$$i" && \
		mv Rplots.pdf "fig_$$i.pdf"; \
	done
