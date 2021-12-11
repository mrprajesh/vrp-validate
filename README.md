# vrp-validate

[![Build](https://github.com/mrprajesh/vrp-validate/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/mrprajesh/vrp-validate/actions/workflows/c-cpp.yml)

### Assumptions

1. We have the input file
2. We have the output file
3. complie `g++ vrp-validate.cpp -o vrp-validate.out`  in gcc 5+

### How to run

```
cat toy.vrp toy.sol | vrp-validate.out
```
