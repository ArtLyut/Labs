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
CMAKE_SOURCE_DIR = /Users/artem/Documents/диплом/dlib/dlib-19.10/dlib/cmake_utils/test_for_neon

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/artem/Documents/диплом/dlib/dlib-19.10/build/temp.macosx-10.13-x86_64-3.6/dlib_build/neon_test_build

# Include any dependencies generated for this target.
include CMakeFiles/neon_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/neon_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/neon_test.dir/flags.make

CMakeFiles/neon_test.dir/neon_test.cpp.o: CMakeFiles/neon_test.dir/flags.make
CMakeFiles/neon_test.dir/neon_test.cpp.o: /Users/artem/Documents/диплом/dlib/dlib-19.10/dlib/cmake_utils/test_for_neon/neon_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/Users/artem/Documents/диплом/dlib/dlib-19.10/build/temp.macosx-10.13-x86_64-3.6/dlib_build/neon_test_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/neon_test.dir/neon_test.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/neon_test.dir/neon_test.cpp.o -c /Users/artem/Documents/диплом/dlib/dlib-19.10/dlib/cmake_utils/test_for_neon/neon_test.cpp

CMakeFiles/neon_test.dir/neon_test.cpp.i: cmake_force
	@echo "Preprocessing CXX source to CMakeFiles/neon_test.dir/neon_test.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/artem/Documents/диплом/dlib/dlib-19.10/dlib/cmake_utils/test_for_neon/neon_test.cpp > CMakeFiles/neon_test.dir/neon_test.cpp.i

CMakeFiles/neon_test.dir/neon_test.cpp.s: cmake_force
	@echo "Compiling CXX source to assembly CMakeFiles/neon_test.dir/neon_test.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/artem/Documents/диплом/dlib/dlib-19.10/dlib/cmake_utils/test_for_neon/neon_test.cpp -o CMakeFiles/neon_test.dir/neon_test.cpp.s

# Object files for target neon_test
neon_test_OBJECTS = \
"CMakeFiles/neon_test.dir/neon_test.cpp.o"

# External object files for target neon_test
neon_test_EXTERNAL_OBJECTS =

libneon_test.a: CMakeFiles/neon_test.dir/neon_test.cpp.o
libneon_test.a: CMakeFiles/neon_test.dir/build.make
libneon_test.a: CMakeFiles/neon_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/Users/artem/Documents/диплом/dlib/dlib-19.10/build/temp.macosx-10.13-x86_64-3.6/dlib_build/neon_test_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libneon_test.a"
	$(CMAKE_COMMAND) -P CMakeFiles/neon_test.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/neon_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/neon_test.dir/build: libneon_test.a

.PHONY : CMakeFiles/neon_test.dir/build

CMakeFiles/neon_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/neon_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/neon_test.dir/clean

CMakeFiles/neon_test.dir/depend:
	cd /Users/artem/Documents/диплом/dlib/dlib-19.10/build/temp.macosx-10.13-x86_64-3.6/dlib_build/neon_test_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/artem/Documents/диплом/dlib/dlib-19.10/dlib/cmake_utils/test_for_neon /Users/artem/Documents/диплом/dlib/dlib-19.10/dlib/cmake_utils/test_for_neon /Users/artem/Documents/диплом/dlib/dlib-19.10/build/temp.macosx-10.13-x86_64-3.6/dlib_build/neon_test_build /Users/artem/Documents/диплом/dlib/dlib-19.10/build/temp.macosx-10.13-x86_64-3.6/dlib_build/neon_test_build /Users/artem/Documents/диплом/dlib/dlib-19.10/build/temp.macosx-10.13-x86_64-3.6/dlib_build/neon_test_build/CMakeFiles/neon_test.dir/DependInfo.cmake
.PHONY : CMakeFiles/neon_test.dir/depend

