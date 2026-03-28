# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

## Project

ApproxMC is a C++ approximate model counter for CNF formulas with PAC (Probably
Approximately Correct) guarantees. It counts satisfying assignments using
random XOR-based hashing. Provides a C++ library, CLI binary, and Python
bindings (pyapproxmc).

## Build Commands

Out-of-source CMake build (from a `build/` directory inside the repo root):

```bash
cmake -DENABLE_TESTING=ON ..
make -j4
ctest        # run tests
```

Key CMake options:
- `-DENABLE_TESTING=ON` — enables Google Test suite
- `-DSTATICCOMPILE=ON` — static linking
- `-DNOZLIB=ON` — disable zlib/gzip support
- `-DSANITIZE=ON` — enable address/UB sanitizers

Requires C++20. Dependencies: CryptoMiniSat5, Arjun, GMP. Optional: Zlib.

## Architecture

- **`src/approxmc.h` / `src/approxmc.cpp`** — Public API. `AppMC` is the facade
  class exposing `count()`, clause/variable management, and configuration.
- **`src/counter.h` / `src/counter.cpp`** — Core `Counter` class implementing
  the approximate counting algorithm (bounded solution counting, XOR hash
  addition, iteration logic).
- **`src/appmcconfig.h`** — Configuration struct (epsilon, delta, seed,
  verbosity, sparse mode).
- **`src/appmc_constants.h` / `.cpp`** — Pre-computed probability tables for
  sparse XOR hashing.
- **`src/main.cpp`** — CLI binary with argument parsing (uses bundled
  `argparse.hpp`).
- **`python/`** — Python bindings via C extension (`pyapproxmc.cpp`).
- **`tests/simpletest.cpp`** — Google Test suite (gtest is a git submodule in
  `tests/gtest`).

The library wraps CryptoMiniSat as the underlying SAT solver and uses Arjun for
sampling set determination. The counting algorithm adds random XOR constraints
and uses binary search over hash counts to estimate the model count.

## Python Bindings

```bash
pip install pyapproxmc
```

Or build from source — `setup.py` compiles the C++ sources directly into a
Python extension module.
