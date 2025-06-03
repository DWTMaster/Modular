template<typename Constants>
class Modular {
 public:
    static constexpr int32_t mod = Constants::MOD;
    static uint32_t fact[Constants::PRECALC_MAX + 1];
    static uint32_t ifact[Constants::PRECALC_MAX + 1];

    /*
     * Pre-calculation
     */

 private:
    static bool constexpr module_is_prime() {
        for (uint32_t i = 2; i * i <= mod; ++i) {
            if (mod % i == 0) { return false; }
        }
        return true;
    }

    static void calculate_factorials() {
        fact[0] = 1;
        for (uint32_t i = 1; i <= Constants::PRECALC_MAX; ++i) {
            fact[i] = static_cast<uint64_t>(fact[i - 1]) * i % mod;
        }
        ifact[Constants::PRECALC_MAX] = inverse(fact[Constants::PRECALC_MAX]);
        for (uint32_t i = Constants::PRECALC_MAX; i > 0; --i) {
            ifact[i - 1] = static_cast<uint64_t>(ifact[i]) * i % mod;
        }
    }

    static void initialize() {
        static_assert(1 < mod && mod <= INT_MAX, "Module must be in range [2; 2^31 - 1]\n");
        if constexpr (!Constants::NEED_PRECALC) { return; }
        static_assert(module_is_prime(), "Precalc is possible only with prime module\n");
        calculate_factorials();
    }

    [[maybe_unused]] static inline const int initializer = (initialize(), 0);


    /*
     * General
     */

    uint32_t value;

    template<typename T>
    static constexpr uint32_t normalize(T number) {
        if (number < static_cast<T>(-mod) || number >= static_cast<T>(mod)) { number %= mod; }
        if (number < 0) { number += mod; }
        return static_cast<uint32_t>(number);
    }

    static constexpr uint32_t inverse(const uint32_t &number) {
        if constexpr (module_is_prime()) {
            return power(number, mod - 2).value;
        }
        int32_t a = static_cast<int32_t>(number), m = mod;
        int32_t u = 0, v = 1;
        while (a != 0) {
            const int32_t t = m / a;
            m -= t * a;
            u -= t * v;
            std::swap(a, m);
            std::swap(u, v);
        }
        if (m != 1) { throw std::runtime_error("Impossible to inverse number\n"); }
        return normalize(u);
    }

 public:
    constexpr Modular() : value() {}

    template<typename T>
    Modular(const T &number) : value(normalize(number)) {}

    /*
     * Static functions
     */

    template<typename T, typename U>
    requires Constants::NEED_PRECALC
    static constexpr Modular binomial(const T &n, const U &k) {
        if (n < k) { return 0; }
        if (n > static_cast<T>(Constants::PRECALC_MAX)) {
            throw std::runtime_error("N is greater than maximum calculated value\n");
        }
        return Modular(fact[n]) * ifact[k] * ifact[n - k];
    }

    template<typename T, typename U>
    static constexpr Modular power(const T &number, U deg) {
        Modular res = 1, x = number;
        if (deg < 0) {
            x = 1 / x;
            deg = -deg;
        }
        while (deg != 0) {
            if (deg & 1) {
                res *= x;
            }
            x *= x;
            deg >>= 1;
        }
        return res;
    }

    template<typename T>
    static constexpr Modular sqrt(const T &number) {
        if (number == 0 || mod == 2) {
            return number;
        }
        if constexpr (!module_is_prime()) {
            throw std::runtime_error("Module is not prime\n");
        }
        if (power(number, (mod - 1) / 2) != 1) {
            throw std::runtime_error("Number is not a quadratic residue\n");
        }
        if constexpr (mod % 4 == 3) {
            return power(number, (mod + 1) / 4);
        }
        uint32_t Q = mod - 1;
        uint32_t S = 0;
        while ((Q & 1) == 0) {
            Q >>= 1;
            ++S;
        }
        uint32_t M = S;
        Modular c, t = power(number, Q), R = power(number, (Q + 1) / 2);
        for (uint32_t z = 2;; ++z) {
            if (power(z, (mod - 1) / 2) == mod - 1) {
                c = power(z, Q);
                break;
            }
        }
        while (t != 1) {
            uint32_t i = 0;
            Modular tmp = t;
            while (tmp != 1) {
                tmp *= tmp;
                ++i;
            }
            Modular b = power(c, 1 << (M - i - 1));
            c = b * b;
            t *= c;
            R *= b;
            M = i;
        }
        return R;
    }

    static constexpr int32_t primitive_root() {
        if constexpr (!module_is_prime()) {
            throw std::runtime_error("Module is not prime\n");
        }
        std::vector<int32_t> primes;
        int32_t number = mod - 1;
        for (int32_t i = 2; i * i <= number; ++i) {
            if (number % i) { continue; }
            primes.push_back(i);
            while (number % i == 0) {
                number /= i;
            }
        }
        if (number != -1) {
            primes.push_back(number);
        }
        for (int r = 2;; ++r) {
            bool ok = true;
            for (int32_t p : primes) {
                if (power(r, (mod - 1) / p) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                return r;
            }
        }
    }

    /*
     * Operators
     */

    const uint32_t &operator ()() const { return value; }

    template<typename T>
    explicit operator T() const { return static_cast<T>(value); }

    Modular &operator +=(const Modular &other) {
        value += other.value;
        if (value >= mod) { value -= mod; }
        return *this;
    }

    template<typename T>
    Modular &operator +=(const T &other) {
        return *this += Modular(other);
    }

    Modular &operator -=(const Modular &other) {
        if (value < other.value) { value += mod; }
        value -= other.value;
        return *this;
    }

    template<typename T>
    Modular &operator -=(const T &other) {
        return *this -= Modular(other);
    }

    Modular &operator *=(const Modular &other) {
        value = static_cast<uint64_t>(value) * other.value % mod;
        return *this;
    }

    template<typename T>
    Modular &operator *=(const T &other) {
        return *this *= Modular(other);
    }

    Modular &operator /=(const Modular &other) {
        if constexpr (Constants::NEED_PRECALC) {
            if (other.value <= Constants::PRECALC_MAX) {
                return *this *= Modular(static_cast<uint64_t>(ifact[other.value]) * fact[other.value - 1]);
            }
        }
        return *this *= inverse(other.value);
    }

    Modular &operator ++() {
        if (++value == mod) { value = 0; }
        return *this;
    }

    Modular &operator --() {
        if (!value) {
            value = mod - 1;
        } else {
            value--;
        }
        return *this;
    }

    Modular operator ++(int) { return Modular(value + 1 < mod ? value + 1 : 0); }

    Modular operator --(int) { return Modular(value ? value - 1 : mod - 1); }

    Modular operator -() { return Modular(value ? mod - value : 0); }

    friend const uint32_t &abs(const Modular &obj) { return obj(); }

    template<typename T>
    friend bool operator ==(const Modular<T> &obj1, const Modular<T> &obj2);

    template<typename T>
    friend bool operator <(const Modular<T> &obj1, const Modular<T> &obj2);

    template<typename U>
    friend std::istream &operator >>(std::istream &in, Modular<U> &obj);
};

template<typename Constants>
uint32_t Modular<Constants>::fact[Constants::PRECALC_MAX + 1];

template<typename Constants>
uint32_t Modular<Constants>::ifact[Constants::PRECALC_MAX + 1];

template<typename T>
bool operator ==(const Modular<T> &obj1, const Modular<T> &obj2) { return obj1.value == obj2.value; }

template<typename T, typename U>
bool operator ==(const Modular<T> &obj, const U &number) { return obj == Modular<T>(number); }

template<typename T, typename U>
bool operator ==(const U &number, const Modular<T> &obj) { return obj == Modular<T>(number); }

template<typename T>
bool operator !=(const Modular<T> &obj1, const Modular<T> &obj2) { return !(obj1 == obj2); }

template<typename T, typename U>
bool operator !=(const Modular<T> &obj, const U &number) { return !(obj == number); }

template<typename T, typename U>
bool operator !=(const U &number, const Modular<T> &obj) { return !(obj == number); }

template<typename T>
bool operator <(const Modular<T> &obj1, const Modular<T> &obj2) { return obj1.value < obj2.value; }

template<typename T>
Modular<T> operator +(const Modular<T> &obj1, const Modular<T> &obj2) { return Modular<T>(obj1) += obj2; }

template<typename T, typename U>
Modular<T> operator +(const Modular<T> &obj, const U &number) { return Modular<T>(obj) += number; }

template<typename T, typename U>
Modular<T> operator +(const U &number, const Modular<T> &obj) { return Modular<T>(obj) += number; }

template<typename T>
Modular<T> operator -(const Modular<T> &obj1, const Modular<T> &obj2) { return Modular<T>(obj1) -= obj2; }

template<typename T, typename U>
Modular<T> operator -(const Modular<T> &obj, const U &number) { return Modular<T>(obj) -= number; }

template<typename T, typename U>
Modular<T> operator -(const U &number, const Modular<T> &obj) { return Modular<T>(number) -= obj; }

template<typename T>
Modular<T> operator *(const Modular<T> &obj1, const Modular<T> &obj2) { return Modular<T>(obj1) *= obj2; }

template<typename T, typename U>
Modular<T> operator *(const Modular<T> &obj, const U &number) { return Modular<T>(obj) *= number; }

template<typename T, typename U>
Modular<T> operator *(const U &number, const Modular<T> &obj) { return Modular<T>(obj) *= number; }

template<typename T>
Modular<T> operator /(const Modular<T> &obj1, const Modular<T> &obj2) { return Modular<T>(obj1) /= obj2; }

template<typename T, typename U>
Modular<T> operator /(const Modular<T> &obj, const U &number) { return Modular<T>(obj) /= number; }

template<typename T, typename U>
Modular<T> operator /(const U &number, const Modular<T> &obj) { return Modular<T>(number) /= obj; }

template<typename T>
std::string to_string(const Modular<T> &number) {
    return to_string(number());
}

template<typename T>
std::istream &operator >>(std::istream &in, Modular<T> &obj) {
    int64_t val;
    in >> val;
    obj.value = Modular<T>::normalize(val);
    return in;
}

template<typename T>
std::ostream &operator <<(std::ostream &out, const Modular<T> &obj) {
    return out << obj();
}


struct Constants {
    static constexpr uint32_t MOD = 1e9 + 7;
    static constexpr bool NEED_PRECALC = true;
    static constexpr uint32_t PRECALC_MAX = 200'000;
};

template
class Modular<Constants>;

using Mint = Modular<Constants>;