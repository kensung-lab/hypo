# SLOG: Simple LOGging library
=======================================================
Slog is a simple C++ library for logging timing and memory (Max RSS) usage.
## Dependencies
- gcc 4.8+ or clang 3.4+
- cmake 3.2+

## Building
CmakeLists is provided in the project root folder.
### Building Library
Run the following commands to build a library `suk.a` in `build/lib` :
```console
  git clone https://github.com/Ritu-Kundu/slog slog
  cd slog
  mkdir build
  cd build
  cmake ..
  make
```
## Usage of the Library
You can use Slog in your project by linking `slog` library (See `Building Library (above)`).
An example usage of the library is as follows:
```console
#include "slog/Monitor.hpp"

int main(int argc, char **argv) {
  // Create monitor instance.
  slog::Monitor monitor;

  // Start the monitor.
  monitor.start();
  
  /* Perform the processing */

  // Print Time spent between the starting of the monitor to this point and Max RSS used so far
  monitor.stop("Message 1");
  // Output will be: "RESOURCES (Message 1): TIME= <elapsed_secs> sec; PEAK RSS (so far)= <peakRSS>MB; CURRENT RSS (so far)= <currentRSS>MB.

  /* Perform the processing */

  // Print Time spent between the previous call to monitor to this point and Max RSS used so far
  monitor.stop("Message 2");

  /* Perform the processing */

  // Print Time spent between the previous call to monitor to this point and Max RSS used so far
  monitor.stop("Message 3");

  /* Perform the processing */

  // Print Time spent between the first call to monitor (to start it) to this point and Max RSS used overall
  monitor.total("Final message");
  // If monitor is to be used again after the call to total, it should be started again.

  return 0;
}
```
### Add Slog into your project via CMake
If you are using CMake for compilation, we recommend adding `slog` as a git submodule with the command git submodule `add https://github.com/Ritu-Kundu/slog external/slog`. Afterwards, modify your top level `CMakeLists.txt` file accordingly:
```console
add_subdirectory(external/slog EXCLUDE_FROM_ALL)
target_link_libraries(your_exe slog)
```
Add `#include "slog/Monitor.hpp"` to access the slog API in your source file.

## Notes
The implementation to compute memory usage makes use of [memory_usage.h](http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use).



