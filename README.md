## Requirements
* ``cmake``
* ``qt5``

### Windows
On Windows you can install the open-source version of the Qt library using the binaries provided by [www.qt.io](https://www.qt.io/download-open-source/). CMake can also be installed via binary packages from [www.cmake.org](https://cmake.org/download/). It is recommended to install the dependencies into the proposed directories; CMake will then find the libraries automatically.

### Linux
Linux users can use their favorite package manager to install the dependencies:

* Ubuntu: ``sudo apt-get install qt5 cmake``
* Arch: ``sudo pacman -S qt5 cmake``
* Fedora: ``sudo yum install qt5 cmake``

### macOS
The easiest way to install the dependencies for macOS is to first install [Homebrew](www.brew.sh). You can then install the dependencies by typing ``brew install qt5 cmake`` into the command line.:w

## Building & Running

~~~Bash
  cd <project-root>
  mkdir build; cd build
  cmake ..
  make
  ./global-illu
~~~

On Windows you can use the graphical UI of CMake to first configure your project and then generate project files for your IDE (for example Visual Studio).
