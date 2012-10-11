CXX = g++
CXXFLAGS = -g -Wall -Wextra -Werror
OBJECTS = map.o point.o list.o
TESTS = point_test.o map_test.o list_test.o

all: tsp

tsp: main.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) main.cpp -o tsp

test_runner: $(OBJECTS) $(TESTS) gtest-all.o
	$(CXX) $(OBJECTS) $(TESTS) gtest-all.o gtest_main.cc -lpthread \
	-o test_runner

check: test_runner
	./test_runner


test: tsp
	./tsp < ../utils/cities.dat > ans1.tmp
	sort ans1.tmp | uniq -u | wc -l
	./tsp < ../utils/random_cities.dat > ans2.tmp
	sort ans2.tmp | uniq -u | wc -l
	./tsp < ../utils/big_random_cities.dat > ans3.tmp
	sort ans3.tmp | uniq -u | wc -l
        

submit: tsp
	python kattis_submit.py -f --problem=oldkattis:tsp \
		   main.cpp  \
		   point.h   \
		   point.cpp \
		   map.h	 \
		   map.cpp   \
		   list.h    \
		   list.cpp


.o :
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf *.o tsp test_runner *.tmp
