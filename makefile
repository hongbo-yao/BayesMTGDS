CXX					:=  mpicxx
CXXFLAGS				:=  -std=c++11 -O3 -fopenmp

INCLUDES				:= $(wildcard src/*.h)
SRCS					:= $(wildcard src/*.cpp) 
OBJS					:= $(patsubst %.cpp, %.o, $(SRCS))

BayesMTGDS_INCLUDE	:= -I src/

%.o : %.cpp 
	@echo "BayesMTGDS is compiling "$<"..."
	@$(CXX) $(CXXFLAGS) $(BayesMTGDS_INCLUDE) -c $< -o $@

BayesMTGDS: $(OBJS) $(INCLUDES)
	@$(CXX) $(CXXFLAGS) -o BayesMTGDS $(OBJS)
	
clean:	
	@rm -rf *.o $(OBJS) BayesMTGDS
