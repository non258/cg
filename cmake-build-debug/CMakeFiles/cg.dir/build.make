# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /home/nozomi/.bin/clion-2019.1.4/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/nozomi/.bin/clion-2019.1.4/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nozomi/CLionProjects/cg

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nozomi/CLionProjects/cg/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/cg.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cg.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cg.dir/flags.make

CMakeFiles/cg.dir/main.cpp.o: CMakeFiles/cg.dir/flags.make
CMakeFiles/cg.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nozomi/CLionProjects/cg/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cg.dir/main.cpp.o"
	/usr/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cg.dir/main.cpp.o -c /home/nozomi/CLionProjects/cg/main.cpp

CMakeFiles/cg.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cg.dir/main.cpp.i"
	/usr/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nozomi/CLionProjects/cg/main.cpp > CMakeFiles/cg.dir/main.cpp.i

CMakeFiles/cg.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cg.dir/main.cpp.s"
	/usr/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nozomi/CLionProjects/cg/main.cpp -o CMakeFiles/cg.dir/main.cpp.s

# Object files for target cg
cg_OBJECTS = \
"CMakeFiles/cg.dir/main.cpp.o"

# External object files for target cg
cg_EXTERNAL_OBJECTS =

cg: CMakeFiles/cg.dir/main.cpp.o
cg: CMakeFiles/cg.dir/build.make
cg: CMakeFiles/cg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nozomi/CLionProjects/cg/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cg"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cg.dir/build: cg

.PHONY : CMakeFiles/cg.dir/build

CMakeFiles/cg.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cg.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cg.dir/clean

CMakeFiles/cg.dir/depend:
	cd /home/nozomi/CLionProjects/cg/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nozomi/CLionProjects/cg /home/nozomi/CLionProjects/cg /home/nozomi/CLionProjects/cg/cmake-build-debug /home/nozomi/CLionProjects/cg/cmake-build-debug /home/nozomi/CLionProjects/cg/cmake-build-debug/CMakeFiles/cg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cg.dir/depend
