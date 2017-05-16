CC=g++ -std=c++11 -pthread
OPT=-O2
DBG=-g
.PHONY: default all clean debug
default: all
all: NB.run
NB.run: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp
	$(CC) $(OPT) Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run -lboost_filesystem -lboost_system -lboost_program_options
	chmod +x NB.run
debug: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp
	$(CC) $(DBG) Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run -lboost_filesystem -lboost_system -lboost_program_options
	chmod +x NB.run
clean:
	@rm NB.run
