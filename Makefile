CXX=g++ -std=c++11 -pthread
OPT=-O3
CPPFLAGS=-I$(BOOSTINCLUDEDIR)
LDFLAGS=-L$(BOOSTLIBDIR) -lboost_filesystem -lboost_system -lboost_program_options
DBG=-g

.PHONY: default all clean debug

default: all

all: NB.run

proteus: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp
	$(CXX) $(OPT) $(CPPFLAGS) Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run $(LDFLAGS)
	chmod +x NB.run

debug: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp
	$(CXX) $(OPT) $(DBG) $(CPPFLAGS) Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run $(LDFLAGS)
	chmod +x NB.run

NB.run: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp
	$(CXX) $(OPT) $(CPPFLAGS) Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run $(LDFLAGS)
	chmod +x NB.run
clean:
	@rm -f NB.run
