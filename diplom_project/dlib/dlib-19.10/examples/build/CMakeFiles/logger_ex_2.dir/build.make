# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.11.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.11.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/artem/Documents/диплом/dlib/dlib-19.10/examples

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/build

# Include any dependencies generated for this target.
include CMakeFiles/logger_ex_2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/logger_ex_2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/logger_ex_2.dir/flags.make

CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.o: CMakeFiles/logger_ex_2.dir/flags.make
CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.o: ../logger_ex_2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/artem/Documents/диплом/dlib/dlib-19.10/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.o -c /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/logger_ex_2.cpp

CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/logger_ex_2.cpp > CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.i

CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/logger_ex_2.cpp -o CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.s

# Object files for target logger_ex_2
logger_ex_2_OBJECTS = \
"CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.o"

# External object files for target logger_ex_2
logger_ex_2_EXTERNAL_OBJECTS =

logger_ex_2: CMakeFiles/logger_ex_2.dir/logger_ex_2.cpp.o
logger_ex_2: CMakeFiles/logger_ex_2.dir/build.make
logger_ex_2: dlib_build/libdlib.a
logger_ex_2: /usr/local/lib/libpng.dylib
logger_ex_2: /usr/lib/libz.dylib
logger_ex_2: /usr/lib/libcblas.dylib
logger_ex_2: /usr/lib/liblapack.dylib
logger_ex_2: /usr/lib/libsqlite3.dylib
logger_ex_2: CMakeFiles/logger_ex_2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/artem/Documents/диплом/dlib/dlib-19.10/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable logger_ex_2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/logger_ex_2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/logger_ex_2.dir/build: logger_ex_2

.PHONY : CMakeFiles/logger_ex_2.dir/build

CMakeFiles/logger_ex_2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/logger_ex_2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/logger_ex_2.dir/clean

CMakeFiles/logger_ex_2.dir/depend:
	cd /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/artem/Documents/диплом/dlib/dlib-19.10/examples /Users/artem/Documents/диплом/dlib/dlib-19.10/examples /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/build /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/build /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/build/CMakeFiles/logger_ex_2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/logger_ex_2.dir/depend

