# Makefile for program: Filament Dynamcis
# Benjamin Dalton: 18-09-2016
# need to  export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"
 
export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"

EX_NAME = memoryextraction
CPP_FLAGS = -c -std=c++14 -O1

OBJS = main.o histogram_constructor.o correlation.o figure_plot.o

$(EX_NAME) : $(OBJS)
	@ export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"; \
	g++ -o $(EX_NAME) $(OBJS) -I/usr/include/python3.9 -I /home/physikqy/.conda/envs/clib/lib/python3.9/site-packages/numpy/core/include -L /home/physikqy/.conda/envs/clib/lib -lpython3.9 -lpthread -lutil -ldl -Xlinker -export-dynamic -lpython3.9 -lfftw3 -lfftw3_threads -lm   


main.o : main.cpp
	@ export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"; \
	g++ $(CPP_FLAGS) main.cpp


histogram_constructor.o : histogram_constructor.cpp potential.h
	g++ $(CPP_FLAGS) histogram_constructor.cpp 

correlation.o : correlation.cpp correlation.h
	g++ $(CPP_FLAGS) correlation.cpp -lfftw3 -lfftw3_threads -lm   
	
figure_plot.o: figure_plot.cpp figure_plot.h 
	@ export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"; \
	g++ $(CPP_FLAGS) figure_plot.cpp -I/usr/include/python3.9 -I /home/physikqy/.conda/envs/clib/lib/python3.9/site-packages/numpy/core/include -L /home/physikqy/.conda/envs/clib/lib -lpython3.9 -lpthread -lutil -ldl -Xlinker -export-dynamic -lpython3.9 
clean:
	rm -f core $(EX_NAME) $(OBJS)
	
