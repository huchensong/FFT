# Makefile for fft project

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -pthread -g -lm -fopenmp
src_dir  = ./src
test_dir = ./test
script_dir = ./script
# Source files
SOURCES = $(src_dir)/fft.cpp $(src_dir)/main.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Target executable
TARGET = fft_app

# Rule to build the executable
$(TARGET): $(OBJECTS) $(src_dir)/fft.h
	$(CXX) $(CXXFLAGS) $^ -o $(src_dir)/$@ 
	./$(src_dir)/$(TARGET)

# Rule to build object files
%.o: %.cpp $(src_dir)/fft.h
	$(CXX) $(CXXFLAGS)  -c $< -o $@
# $@表示目标
# Rule to clean the project
clean:
	rm -f $(TARGET) $(OBJECTS)
complex_test:
	$(CXX) $(CXXFLAGS) $(test_dir)/complex_test.cpp -o $(test_dir)/complex_test
generate_w: $(script_dir)/generate_w_table.cpp $(src_dir)/fft.o
	$(CXX) -c $(script_dir)/generate_w_table.cpp -o generate_w_table.o
	 $(CXX) $(CXXFLAGS) $(src_dir)/fft.o generate_w_table.o -o generate_w
	./generate_w
