all:a.out

a.out: force.cc
	g++ -O3 force.cc

clean:
	rm -f a.out
