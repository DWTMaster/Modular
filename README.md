# Modular Arithmetic Library

<p style="text-align: center;">
  <img src="https://img.shields.io/badge/C%2B%2B-20-blue" alt="C++20">
  <img src="https://img.shields.io/badge/Performance-Optimized-brightgreen" alt="Performance">
</p>

## Overview

The **Modular Arithmetic Library** is a highly efficient, flexible, and easy-to-use C++ template-based library for performing modular arithmetic operations. This library supports a wide range of modular arithmetic features, including:

- Basic operations (`+`, `-`, `*`, `/`, `%`) under a modular system.
- Modular exponentiation, binomial coefficients, modular square roots, and primitive root calculations.
- Precomputed factorial and inverse factorial arrays for efficient combinatorics.
- Support for any modulus (prime or non-prime, with runtime checks for non-prime cases).

With its design, the library achieves **3.1x slower performance compared to `unsigned int`**, even when using precomputed factorials, making it an excellent choice for competitive programming and mathematical computations.

---

## Features

### Modular Arithmetic
Supports:
- Addition, subtraction, multiplication, and division under modular arithmetic.
- Modular inverses using Fermat's Little Theorem (for prime moduli) or Extended Euclidean Algorithm.

### Combinatorics
- Efficient computation of binomial coefficients with precomputed factorials and inverse factorials.
- Precomputation is highly optimized for fast combinatorics.

### Advanced Features
- Modular exponentiation (`power(base, exponent)`).
- Modular square root computation using Tonelli–Shanks algorithm.
- Primitive root computation for prime moduli.

### Compile-Time Optimizations
- Supports `constexpr` for compile-time evaluation of certain operations.
- Uses `if constexpr` to enable optimizations when the modulus is known to be prime.

---

## Usage

### Example Code
Below is a usage example of the library:

```cpp
#include <bits/stdc++.h>
#include "modular.cpp"

int main() {
    // MOD = 1e9 + 7
    Mint a = 3, b = 5;

    Mint sum = a + b;
    Mint prod = a * b;
    Mint quot = a / b;
    Mint pow = Mint::power(a, 10);
    Mint bin = Mint::binomial(7, 4);

    std::cout << "Sum: " << sum << std::endl;
    std::cout << "Product: " << prod << std::endl;
    std::cout << "Quotient: " << quot << std::endl;
    std::cout << "Power: " << pow << std::endl;
    std::cout << "Binomial: " << bin << std::endl;
}
```

### Output:
```
Sum: 8
Product: 15
Quotient: 200000002
Power: 59049
Binomial: 35
```

---

## Installation

### Include in Your Project
Simply include the `modular.cpp` file in your project or copy the class to your file and change the constants in the structure Constants.

```cpp
#include "modular.cpp"
```

---

## API Documentation

### Class Template: `Modular<Constants>`
This is the main class that provides modular arithmetic operations.

#### Supported Operators
- Arithmetic: `+`, `-`, `*`, `/`
- Increment/Decrement: `++`, `--`
- Unary: `-` (negation)
- Comparison: `==`, `!=`, `<`

#### Static Methods
- `power(base, exponent)`: Computes `(base^exponent) % MOD`.
- `binomial(n, k)`: Computes the binomial coefficient `C(n, k)` modulo `MOD`.
- `sqrt(x)`: Computes the modular square root of `x`.

#### Instance Methods
- `{object}()`: Returns the raw value of the modular number.
- Explicit conversion to `int`, `uint32_t`, etc.

---

## Performance

The library is designed for competitive programming and high-performance use cases. Below is a breakdown of its performance:

| Feature                               | Time Complexity |
|---------------------------------------|-----------------|
| Addition/Subtraction/Multiplication   | O(1)            |
| Division/Exponentiation               | O(log(MOD))     |
| Precomputing                          | O(k + sqrt(MOD))|
| Binomial/Precomuted division          | O(1)            |   
| Square Root                           | O(log^2(MOD))   |

---
```
