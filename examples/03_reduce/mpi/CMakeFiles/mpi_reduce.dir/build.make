# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.4

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/projects/x86-64-haswell/cmake/3.4.3/bin/cmake

# The command to remove a file.
RM = /home/projects/x86-64-haswell/cmake/3.4.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/evalen/public/Fenix

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/evalen/public/Fenix

# Include any dependencies generated for this target.
include examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/depend.make

# Include the progress variables for this target.
include examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/progress.make

# Include the compile flags for this target's objects.
include examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/flags.make

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o: examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/flags.make
examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o: examples/03_reduce/mpi/mpi_reduce.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/evalen/public/Fenix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o"
	cd /home/evalen/public/Fenix/examples/03_reduce/mpi && /home/evalen/devtools/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o   -c /home/evalen/public/Fenix/examples/03_reduce/mpi/mpi_reduce.c

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mpi_reduce.dir/mpi_reduce.c.i"
	cd /home/evalen/public/Fenix/examples/03_reduce/mpi && /home/evalen/devtools/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/evalen/public/Fenix/examples/03_reduce/mpi/mpi_reduce.c > CMakeFiles/mpi_reduce.dir/mpi_reduce.c.i

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mpi_reduce.dir/mpi_reduce.c.s"
	cd /home/evalen/public/Fenix/examples/03_reduce/mpi && /home/evalen/devtools/bin/mpicc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/evalen/public/Fenix/examples/03_reduce/mpi/mpi_reduce.c -o CMakeFiles/mpi_reduce.dir/mpi_reduce.c.s

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.requires:

.PHONY : examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.requires

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.provides: examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.requires
	$(MAKE) -f examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/build.make examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.provides.build
.PHONY : examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.provides

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.provides.build: examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o


# Object files for target mpi_reduce
mpi_reduce_OBJECTS = \
"CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o"

# External object files for target mpi_reduce
mpi_reduce_EXTERNAL_OBJECTS =

examples/03_reduce/mpi/mpi_reduce: examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o
examples/03_reduce/mpi/mpi_reduce: examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/build.make
examples/03_reduce/mpi/mpi_reduce: examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/evalen/public/Fenix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable mpi_reduce"
	cd /home/evalen/public/Fenix/examples/03_reduce/mpi && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpi_reduce.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/build: examples/03_reduce/mpi/mpi_reduce

.PHONY : examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/build

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/requires: examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/mpi_reduce.c.o.requires

.PHONY : examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/requires

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/clean:
	cd /home/evalen/public/Fenix/examples/03_reduce/mpi && $(CMAKE_COMMAND) -P CMakeFiles/mpi_reduce.dir/cmake_clean.cmake
.PHONY : examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/clean

examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/depend:
	cd /home/evalen/public/Fenix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/evalen/public/Fenix /home/evalen/public/Fenix/examples/03_reduce/mpi /home/evalen/public/Fenix /home/evalen/public/Fenix/examples/03_reduce/mpi /home/evalen/public/Fenix/examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/03_reduce/mpi/CMakeFiles/mpi_reduce.dir/depend

