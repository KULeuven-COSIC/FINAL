CCX = g++
CCXFLAGS = -O3 -funroll-loops -march=native -std=c++11 -pthread -I. -I./include
DEPS = -lntl -lgmp -lfftw3 -lm

all: clean test
	 
clean:
	$(RM) test test.o lwehe.o ntruhe.o fft.o sampler.o keygen.o

test: include/params.h ntruhe.o lwehe.o keygen.o fft.o sampler.o
	$(CCX) $(CCXFLAGS) -o test test.cpp ntruhe.o lwehe.o keygen.o fft.o sampler.o $(DEPS)

ntruhe.o: include/ntruhe.h keygen.o sampler.o lwehe.o src/ntruhe.cpp
	$(CCX) $(CCXFLAGS) -c src/ntruhe.cpp

lwehe.o: include/lwehe.h keygen.o sampler.o src/lwehe.cpp
	$(CCX) $(CCXFLAGS) -c src/lwehe.cpp

keygen.o: include/keygen.h sampler.o fft.o src/keygen.cpp
	$(CCX) $(CCXFLAGS) -c src/keygen.cpp

fft.o: include/fft.h
	$(CCX) $(CCXFLAGS) -c src/fft.cpp

sampler.o: include/sampler.h include/params.h src/sampler.cpp
	$(CCX) $(CCXFLAGS) -c src/sampler.cpp