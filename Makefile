# Compiler settings
CC = gcc              # Compiler
CFLAGS = -c -g           # Flags to compile without linking
LDFLAGS = -lopenblas  # Link with OpenBLAS library
EXEC = spmv           # Output executable name

# Source files
SRC = timer.c spmv.c my_dense.c my_sparse.c
OBJ = $(SRC:.c=.o)    # Object files corresponding to source files

# Default target to build and run the project
all: $(EXEC) run

# Rule to link the object files and create the executable
$(EXEC): $(OBJ)
	@echo "Linking object files to create the executable..."
	$(CC) $(OBJ) $(LDFLAGS) -o $(EXEC)

# Rule to compile .c files into .o files
%.o: %.c
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) $<

# Run the executable
run: $(EXEC)
	@echo "Running the program..."
	./$(EXEC)

# Clean up object files and the executable
clean:
	@echo "Cleaning up..."
	rm -f $(OBJ) $(EXEC)


#make : to compile and run
#make clean : to remove the compile and executable file