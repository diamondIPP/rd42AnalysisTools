# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/CrossTalkCorrections.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CrossTalkCorrections.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CrossTalkCorrections.dir/flags.make

CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o: CMakeFiles/CrossTalkCorrections.dir/flags.make
CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o: ../FeedThroughCorrection.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o -c /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/FeedThroughCorrection.cpp

CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/FeedThroughCorrection.cpp > CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.i

CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/FeedThroughCorrection.cpp -o CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.s

CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.requires:

.PHONY : CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.requires

CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.provides: CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.requires
	$(MAKE) -f CMakeFiles/CrossTalkCorrections.dir/build.make CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.provides.build
.PHONY : CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.provides

CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.provides.build: CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o


# Object files for target CrossTalkCorrections
CrossTalkCorrections_OBJECTS = \
"CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o"

# External object files for target CrossTalkCorrections
CrossTalkCorrections_EXTERNAL_OBJECTS =

CrossTalkCorrections: CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o
CrossTalkCorrections: CMakeFiles/CrossTalkCorrections.dir/build.make
CrossTalkCorrections: CMakeFiles/CrossTalkCorrections.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable CrossTalkCorrections"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CrossTalkCorrections.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CrossTalkCorrections.dir/build: CrossTalkCorrections

.PHONY : CMakeFiles/CrossTalkCorrections.dir/build

CMakeFiles/CrossTalkCorrections.dir/requires: CMakeFiles/CrossTalkCorrections.dir/FeedThroughCorrection.cpp.o.requires

.PHONY : CMakeFiles/CrossTalkCorrections.dir/requires

CMakeFiles/CrossTalkCorrections.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CrossTalkCorrections.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CrossTalkCorrections.dir/clean

CMakeFiles/CrossTalkCorrections.dir/depend:
	cd /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/cmake-build-debug /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/cmake-build-debug /Users/diegoalejandro/Desktop/mntSnickers/CrossTalkCorrections/cmake-build-debug/CMakeFiles/CrossTalkCorrections.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CrossTalkCorrections.dir/depend

