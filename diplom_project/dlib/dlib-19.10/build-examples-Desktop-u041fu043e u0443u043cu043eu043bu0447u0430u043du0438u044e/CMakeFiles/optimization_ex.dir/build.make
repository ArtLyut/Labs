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
CMAKE_BINARY_DIR = "/Users/artem/Documents/диплом/dlib/dlib-19.10/build-examples-Desktop-u041fu043e u0443u043cu043eu043bu0447u0430u043du0438u044e"

# Include any dependencies generated for this target.
include CMakeFiles/optimization_ex.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/optimization_ex.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/optimization_ex.dir/flags.make

CMakeFiles/optimization_ex.dir/optimization_ex.cpp.o: CMakeFiles/optimization_ex.dir/flags.make
CMakeFiles/optimization_ex.dir/optimization_ex.cpp.o: /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/optimization_ex.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/artem/Documents/диплом/dlib/dlib-19.10/build-examples-Desktop-u041fu043e u0443u043cu043eu043bu0447u0430u043du0438u044e/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/optimization_ex.dir/optimization_ex.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optimization_ex.dir/optimization_ex.cpp.o -c /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/optimization_ex.cpp

CMakeFiles/optimization_ex.dir/optimization_ex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization_ex.dir/optimization_ex.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/optimization_ex.cpp > CMakeFiles/optimization_ex.dir/optimization_ex.cpp.i

CMakeFiles/optimization_ex.dir/optimization_ex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization_ex.dir/optimization_ex.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/artem/Documents/диплом/dlib/dlib-19.10/examples/optimization_ex.cpp -o CMakeFiles/optimization_ex.dir/optimization_ex.cpp.s

# Object files for target optimization_ex
optimization_ex_OBJECTS = \
"CMakeFiles/optimization_ex.dir/optimization_ex.cpp.o"

# External object files for target optimization_ex
optimization_ex_EXTERNAL_OBJECTS =

optimization_ex: CMakeFiles/optimization_ex.dir/optimization_ex.cpp.o
optimization_ex: CMakeFiles/optimization_ex.dir/build.make
optimization_ex: dlib_build/libdlib.a
optimization_ex: /usr/local/lib/libpng.dylib
optimization_ex: /usr/lib/libz.dylib
optimization_ex: /usr/lib/libcblas.dylib
optimization_ex: /usr/lib/liblapack.dylib
optimization_ex: /usr/lib/libsqlite3.dylib
optimization_ex: CMakeFiles/optimization_ex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/artem/Documents/диплом/dlib/dlib-19.10/build-examples-Desktop-u041fu043e u0443u043cu043eu043bu0447u0430u043du0438u044e/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable optimization_ex"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/optimization_ex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/optimization_ex.dir/build: optimization_ex

.PHONY : CMakeFiles/optimization_ex.dir/build

CMakeFiles/optimization_ex.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/optimization_ex.dir/cmake_clean.cmake
.PHONY : CMakeFiles/optimization_ex.dir/clean

CMakeFiles/optimization_ex.dir/depend:
	cd "/Users/artem/Documents/диплом/dlib/dlib-19.10/build-examples-Desktop-u041fu043e u0443u043cu043eu043bu0447u0430u043du0438u044e" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/artem/Documents/диплом/dlib/dlib-19.10/examples /Users/artem/Documents/диплом/dlib/dlib-19.10/examples "/Users/artem/Documents/диплом/dlib/dlib-19.10/build-examples-Desktop-u041fu043e u0443u043cu043eu043bu0447u0430u043du0438u044e" "/Users/artem/Documents/диплом/dlib/dlib-19.10/build-examples-Desktop-u041fu043e u0443u043cu043eu043bu0447u0430u043du0438u044e" "/Users/artem/Documents/диплом/dlib/dlib-19.10/build-examples-Desktop-u041fu043e u0443u043cu043eu043bu0447u0430u043du0438u044e/CMakeFiles/optimization_ex.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/optimization_ex.dir/depend

