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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/junyao/git_OTC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/junyao/git_OTC

# Include any dependencies generated for this target.
include CMakeFiles/TMC_OCTv3.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TMC_OCTv3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TMC_OCTv3.dir/flags.make

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.o: CMakeFiles/TMC_OCTv3.dir/flags.make
CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.o: src/tmcoct_go.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/junyao/git_OTC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.o   -c /home/junyao/git_OTC/src/tmcoct_go.c

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/junyao/git_OTC/src/tmcoct_go.c > CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.i

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/junyao/git_OTC/src/tmcoct_go.c -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.s

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.o: CMakeFiles/TMC_OCTv3.dir/flags.make
CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.o: src/tmcoct_io.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/junyao/git_OTC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.o   -c /home/junyao/git_OTC/src/tmcoct_io.c

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/junyao/git_OTC/src/tmcoct_io.c > CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.i

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/junyao/git_OTC/src/tmcoct_io.c -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.s

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.o: CMakeFiles/TMC_OCTv3.dir/flags.make
CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.o: src/tmcoct_main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/junyao/git_OTC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.o   -c /home/junyao/git_OTC/src/tmcoct_main.c

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/junyao/git_OTC/src/tmcoct_main.c > CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.i

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/junyao/git_OTC/src/tmcoct_main.c -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.s

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.o: CMakeFiles/TMC_OCTv3.dir/flags.make
CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.o: src/tmcoct_mesh.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/junyao/git_OTC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.o   -c /home/junyao/git_OTC/src/tmcoct_mesh.c

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/junyao/git_OTC/src/tmcoct_mesh.c > CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.i

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/junyao/git_OTC/src/tmcoct_mesh.c -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.s

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.o: CMakeFiles/TMC_OCTv3.dir/flags.make
CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.o: src/tmcoct_nr.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/junyao/git_OTC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.o   -c /home/junyao/git_OTC/src/tmcoct_nr.c

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/junyao/git_OTC/src/tmcoct_nr.c > CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.i

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/junyao/git_OTC/src/tmcoct_nr.c -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.s

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.o: CMakeFiles/TMC_OCTv3.dir/flags.make
CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.o: src/tmcoct_bias.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/junyao/git_OTC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.o   -c /home/junyao/git_OTC/src/tmcoct_bias.c

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/junyao/git_OTC/src/tmcoct_bias.c > CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.i

CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/junyao/git_OTC/src/tmcoct_bias.c -o CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.s

# Object files for target TMC_OCTv3
TMC_OCTv3_OBJECTS = \
"CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.o" \
"CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.o" \
"CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.o" \
"CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.o" \
"CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.o" \
"CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.o"

# External object files for target TMC_OCTv3
TMC_OCTv3_EXTERNAL_OBJECTS = \
"/home/junyao/git_OTC/build/tmcoct_go.o" \
"/home/junyao/git_OTC/build/tmcoct_io.o" \
"/home/junyao/git_OTC/build/tmcoct_mesh.o" \
"/home/junyao/git_OTC/build/tmcoct_nr.o" \
"/home/junyao/git_OTC/build/tmcoct_bias.o"

TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/src/tmcoct_go.c.o
TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/src/tmcoct_io.c.o
TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/src/tmcoct_main.c.o
TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/src/tmcoct_mesh.c.o
TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/src/tmcoct_nr.c.o
TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/src/tmcoct_bias.c.o
TMC_OCTv3: build/tmcoct_go.o
TMC_OCTv3: build/tmcoct_io.o
TMC_OCTv3: build/tmcoct_mesh.o
TMC_OCTv3: build/tmcoct_nr.o
TMC_OCTv3: build/tmcoct_bias.o
TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/build.make
TMC_OCTv3: CMakeFiles/TMC_OCTv3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/junyao/git_OTC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking C executable TMC_OCTv3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TMC_OCTv3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TMC_OCTv3.dir/build: TMC_OCTv3

.PHONY : CMakeFiles/TMC_OCTv3.dir/build

CMakeFiles/TMC_OCTv3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TMC_OCTv3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TMC_OCTv3.dir/clean

CMakeFiles/TMC_OCTv3.dir/depend:
	cd /home/junyao/git_OTC && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/junyao/git_OTC /home/junyao/git_OTC /home/junyao/git_OTC /home/junyao/git_OTC /home/junyao/git_OTC/CMakeFiles/TMC_OCTv3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TMC_OCTv3.dir/depend
