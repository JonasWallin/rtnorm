
PRO = rtnorm

CXX = g++
CXXFLAG = -O3 -fopenmp -std=c++11 -W -Wall -g
LIB = -lm  -lgslcblas -lgomp
LIBRARY_DIRS = -L/opt/local/lib/
SRC1 = src/rtnorm.cpp 
OBJ1 = tmp/rtnorm.o

SRC2 = src/main.cpp
OBJ2 = tmp/main.o

SRC3 = src/qnorm.cpp
OBJ3 = tmp/qnorm.o

HDRS = src/rtnorm.hpp

all : $(PRO)

$(PRO): $(OBJ1) $(OBJ2) $(OBJ3)
	$(CXX) -o $(PRO) $(OBJ1) $(OBJ2) $(LIBRARY_DIRS) $(LIB)

$(OBJ1) : $(SRC1) $(HDRS)
	$(CXX) -o $(OBJ1) -c $(SRC1) $(CXXFLAG)

$(OBJ2) : $(SRC2) $(HDRS)
	$(CXX) -o $(OBJ2) -c $(SRC2) $(CXXFLAG)

$(OBJ3) : $(SRC3) $(HDRS)
	$(CXX) -o $(OBJ3) -c $(SRC3) $(CXXFLAG)

clean : 
	rm -rf $(PRO) ./tmp/*.o
