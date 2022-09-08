# phasedisperser_tilde

A phase dispersion Max object.

The algorithm and parameter names are courtesy of the free and open-source allpassphase VST by
[enummusic](https://github.com/enummusic/allpassphase).


## Building

Required:
- Mac: Clang (e.g. Xcode)
- Windows: MSVC (e.g. Visual Studio)
- CMake
- Ninja (optional)

```sh
git clone https://git.sr.ht/~isabelgk/phasedisperser_tilde --recurse-submodules
mkdir build && cd build
cmake ..
cmake --build .
```