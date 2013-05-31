all: ga
ga: ga.hh ga.cc
	g++ -o ga ga.cc
debug: ga.hh ga.cc
	g++ -g -o ga-debug ga.cc
run: ga
	./ga
clear: ga
	rm ga
