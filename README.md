# LITTLE KISS FFT
A port of [KISS FFT](https://github.com/mborgerding/kissfft) that targets μcontrollers.


## Features

 * Uses int16_t as data type.
 * Unlike KISS FFT, LITTLE KISS FFT will never use doubles, not even for initialization code. At worst, float will be used during setup. (sparingly.)
 * Drop-in code for Arduino IDE, no CMakefile or Makefile.

## License

Like its origin, LITTLE KISS FFT uses BSD-3 License.

