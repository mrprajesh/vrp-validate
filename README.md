# vrp-validate

[![Dev Branch](https://github.com/mrprajesh/vrp-validate/actions/workflows/dev.yml/badge.svg?branch=dev)](https://github.com/mrprajesh/vrp-validate/actions/workflows/dev.yml)
[![Main Branch](https://github.com/mrprajesh/vrp-validate/actions/workflows/main.yml/badge.svg)](https://github.com/mrprajesh/vrp-validate/actions/workflows/main.yml)

### Assumptions

1. We have the input file
2. We have the output file
3. complie `g++ vrp-validate.cpp -o vrp-validate.out`  in gcc 5+

### How to run

```
cat toy.vrp toy.sol | vrp-validate.out
```


