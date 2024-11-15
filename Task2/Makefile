# Compiler settings
CC = icc  # Default compiler
CFLAGS = $(OPTFLAGS) -c -g          # Flags to compile without linking
LDFLAGS = -lgsl -lm -lopenblas  # Link with OpenBLAS library
EXEC = spmv           # Output executable name
OPTFLAGS = -O0  # Default optimization level

# Source files
SRC = timer.c spmv.c my_dense.c my_sparse.c
OBJ = $(SRC:.c=.o)    # Object files corresponding to source files

# Default target to build and run the project
all: $(EXEC) run

# Rule to link the object files and create the executable
$(EXEC): $(OBJ)
	echo "Linking object files to create the executable..."
	$(CC) $(OBJ) $(LDFLAGS) -o $(EXEC)

# Rule to compile .c files into .o files
%.o: %.c
	echo "Compiling $<..."
	$(CC) $(CFLAGS) $<

# Run the executable
run: $(EXEC)
	echo "Running the program..."
	./$(EXEC)

# Clean up object files and the executable
clean:
	echo "Cleaning up..."
	rm -f $(OBJ) $(EXEC)

#make #to run
#make clean #to clean