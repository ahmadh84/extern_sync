# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.6

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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
CMAKE_COMMAND = /lusr/opt/cmake-2.6.4/bin/cmake

# The command to remove a file.
RM = /lusr/opt/cmake-2.6.4/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /lusr/opt/cmake-2.6.4/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /projects/vision/6/jaechul/flann-1.21-src-64/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /projects/vision/6/jaechul/flann-1.21-src-64/tmp

# Include any dependencies generated for this target.
include cpp/tests/CMakeFiles/flann_test.dir/depend.make

# Include the progress variables for this target.
include cpp/tests/CMakeFiles/flann_test.dir/progress.make

# Include the compile flags for this target's objects.
include cpp/tests/CMakeFiles/flann_test.dir/flags.make

cpp/tests/CMakeFiles/flann_test.dir/flann_test.o: cpp/tests/CMakeFiles/flann_test.dir/flags.make
cpp/tests/CMakeFiles/flann_test.dir/flann_test.o: /projects/vision/6/jaechul/flann-1.21-src-64/src/cpp/tests/flann_test.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /projects/vision/6/jaechul/flann-1.21-src-64/tmp/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object cpp/tests/CMakeFiles/flann_test.dir/flann_test.o"
	cd /projects/vision/6/jaechul/flann-1.21-src-64/tmp/cpp/tests && /lusr/opt/gcc-4.0.2/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/flann_test.dir/flann_test.o -c /projects/vision/6/jaechul/flann-1.21-src-64/src/cpp/tests/flann_test.cc

cpp/tests/CMakeFiles/flann_test.dir/flann_test.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/flann_test.dir/flann_test.i"
	cd /projects/vision/6/jaechul/flann-1.21-src-64/tmp/cpp/tests && /lusr/opt/gcc-4.0.2/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /projects/vision/6/jaechul/flann-1.21-src-64/src/cpp/tests/flann_test.cc > CMakeFiles/flann_test.dir/flann_test.i

cpp/tests/CMakeFiles/flann_test.dir/flann_test.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/flann_test.dir/flann_test.s"
	cd /projects/vision/6/jaechul/flann-1.21-src-64/tmp/cpp/tests && /lusr/opt/gcc-4.0.2/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /projects/vision/6/jaechul/flann-1.21-src-64/src/cpp/tests/flann_test.cc -o CMakeFiles/flann_test.dir/flann_test.s

cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.requires:
.PHONY : cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.requires

cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.provides: cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.requires
	$(MAKE) -f cpp/tests/CMakeFiles/flann_test.dir/build.make cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.provides.build
.PHONY : cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.provides

cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.provides.build: cpp/tests/CMakeFiles/flann_test.dir/flann_test.o
.PHONY : cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.provides.build

# Object files for target flann_test
flann_test_OBJECTS = \
"CMakeFiles/flann_test.dir/flann_test.o"

# External object files for target flann_test
flann_test_EXTERNAL_OBJECTS =

cpp/tests/flann_test: cpp/tests/CMakeFiles/flann_test.dir/flann_test.o
cpp/tests/flann_test: cpp/libflann.so
cpp/tests/flann_test: cpp/tests/CMakeFiles/flann_test.dir/build.make
cpp/tests/flann_test: cpp/tests/CMakeFiles/flann_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable flann_test"
	cd /projects/vision/6/jaechul/flann-1.21-src-64/tmp/cpp/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/flann_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
cpp/tests/CMakeFiles/flann_test.dir/build: cpp/tests/flann_test
.PHONY : cpp/tests/CMakeFiles/flann_test.dir/build

cpp/tests/CMakeFiles/flann_test.dir/requires: cpp/tests/CMakeFiles/flann_test.dir/flann_test.o.requires
.PHONY : cpp/tests/CMakeFiles/flann_test.dir/requires

cpp/tests/CMakeFiles/flann_test.dir/clean:
	cd /projects/vision/6/jaechul/flann-1.21-src-64/tmp/cpp/tests && $(CMAKE_COMMAND) -P CMakeFiles/flann_test.dir/cmake_clean.cmake
.PHONY : cpp/tests/CMakeFiles/flann_test.dir/clean

cpp/tests/CMakeFiles/flann_test.dir/depend:
	cd /projects/vision/6/jaechul/flann-1.21-src-64/tmp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /projects/vision/6/jaechul/flann-1.21-src-64/src /projects/vision/6/jaechul/flann-1.21-src-64/src/cpp/tests /projects/vision/6/jaechul/flann-1.21-src-64/tmp /projects/vision/6/jaechul/flann-1.21-src-64/tmp/cpp/tests /projects/vision/6/jaechul/flann-1.21-src-64/tmp/cpp/tests/CMakeFiles/flann_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : cpp/tests/CMakeFiles/flann_test.dir/depend

