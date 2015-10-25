CXX = g++-mp-4.8  
CXXFLAG_CPP = -O3 -fopenmp -std=c++11 -W -Wall -g
 
OPENMP = -fopenmp
MATLAB_ROOT_BIN = /Applications/MATLAB_R2013a.app/bin/
INCL     =  -I/opt/local/include 
LIBS 	 = -L/usr/local/lib/ -L/opt/local/lib/
CXXFLAGS = -Isrc/ -std=c++11 $(OPENMP)
LDFLAGS = -std=c++11

MEX = $(MATLAB_ROOT_BIN)mex -v CXXFLAGS="\$$CXXFLAGS $(CXXFLAGS)" LDFLAGS="\$$LDFLAGS $(OPENMP)"    
INCL += -I../  


PRO = rtnorm

all: $(PRO)


#rtnorm_mex: rtnorm_mex.cpp 
#	$(MEX) $^  $(INCL) $(LIBS) -o $@


SRC1 = src/rtnorm.cpp 
OBJ1 = tmp/rtnorm.o

SRC2 = src/qnorm.cpp
OBJ2 = tmp/qnorm.o

SRC3 = rtnorm_mex.cpp

HDRS = src/rtnorm.hpp


$(PRO): $(SRC3)  $(OBJ1)
	$(MEX) $^  $(INCL) $(LIBS) -o $(PRO)

$(OBJ1) : $(SRC1) $(HDRS)
	$(CXX) -o $(OBJ1) -c $(SRC1) $(CXXFLAG_CPP)

$(OBJ2) : $(SRC2) $(HDRS)
	$(CXX) -o $(OBJ2) -c $(SRC2) $(CXXFLAG_CPP)

.PHONY: clean
clean:
	rm -f *.*~ *.o  tmp/* *.mex* 

