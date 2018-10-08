HTS_INC=-Ihtslib/htslib
example: main.o htslib/libhts.a
	g++ -o $@ $^ -lz -lm -lbz2 -llzma -lpthread

%.o: %.cpp
	g++ -g -std=c++11 -o $@ -c $< -I bamtools $(HTS_INC)

htslib/htslib.a:
	cd htslib && make

