compile:
	clang++ -Wall -g -O3 sim.cpp -o sim
compile_debug:
	clang++ -Wall -g -O3 -D _DEBUG sim.cpp -o sim
run:
	./sim
time:
	time ./sim
full: compile_debug run compile time
