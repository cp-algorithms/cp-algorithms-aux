#ifndef CP_ALGO_MATH_LAURENT_HPP
#define CP_ALGO_MATH_LAURENT_HPP
#include "../util/big_alloc.hpp"
#include <memory>
#include <optional>
#include <vector>
#include <cassert>
#include <algorithm>

namespace cp_algo::math {
    // Base provider interface for lazy coefficient evaluation
    template<typename T>
    struct provider {
        mutable big_vector<T> cache;
        mutable int cache_offset = 0; // Index of first cached coefficient
        mutable bool initialized = false;
        mutable bool all_cached = false; // True if all non-zero coeffs are cached
        
        virtual ~provider() = default;
        virtual std::optional<int> valuation() const { return std::nullopt; }
        virtual std::optional<int> degree() const { return std::nullopt; }
        
        // Compute k-th coefficient lazily without caching
        virtual T coeff_lazy(int k) const = 0;
        
        // Double the number of known coefficients (or cache all for finite series)
        // Default: use coeff_lazy, but can be overridden for efficiency (e.g., FFT)
        virtual void double_up() const {
            int old_size = cache.size();
            int new_size = old_size == 0 ? 1 : 2 * old_size;
            
            // Check if we have a known degree
            if(auto deg = degree()) {
                int max_idx = *deg - cache_offset;
                if(max_idx < new_size) {
                    // All coefficients fit, just cache them all
                    new_size = max_idx + 1;
                    all_cached = true;
                }
            }
            
            cache.resize(new_size);
            for(int i = old_size; i < new_size; i++) {
                cache[i] = coeff_lazy(cache_offset + i);
            }
        }
        
        // Get coefficient with caching and doubling (default implementation)
        virtual T coeff(int k) const {
            if(!initialized) {
                auto v = valuation();
                cache_offset = v ? *v : 0;
                initialized = true;
            }
            
            int idx = k - cache_offset;
            if(idx < 0) {
                return T(0); // Below cached range
            }
            
            // If all coeffs are cached and we're beyond cache, return 0
            if(all_cached && idx >= (int)cache.size()) {
                return T(0);
            }
            
            // Extend cache by doubling until we have enough
            while(idx >= (int)cache.size() && !all_cached) {
                double_up();
            }
            
            // Should have the coefficient now
            if(idx < (int)cache.size()) {
                return cache[idx];
            }
            
            return T(0);
        }
        
        // Alias for backwards compatibility
        T get(int k) const {
            return coeff(k);
        }
    };
    
    // Constant provider - returns a single coefficient at position offset
    template<typename T>
    struct constant_provider : provider<T> {
        T value;
        int offset;
        
        constant_provider(T value, int offset = 0) : value(value), offset(offset) {}
        
        std::optional<int> valuation() const override {
            return value != T(0) ? std::optional<int>(offset) : std::nullopt;
        }
        
        std::optional<int> degree() const override {
            return value != T(0) ? std::optional<int>(offset) : std::nullopt;
        }
        
        T coeff_lazy(int k) const override {
            return k == offset ? value : T(0);
        }
        
        T coeff(int k) const override {
            return coeff_lazy(k);
        }
    };
    
    // Polynomial provider - wraps a vector of coefficients
    template<typename T>
    struct polynomial_provider : provider<T> {
        polynomial_provider(big_vector<T> coeffs, int offset = 0) {
            // Find first and last non-zero coefficients
            auto non_zero = [](const T& x) { return x != T(0); };
            auto first = std::ranges::find_if(coeffs, non_zero);
            auto last = std::ranges::find_if(coeffs | std::views::reverse, non_zero);
            
            // Extract non-zero range
            if(first != coeffs.end()) {
                int leading_zeros = first - coeffs.begin();
                int trailing_zeros = last - coeffs.rbegin();
                coeffs = big_vector<T>(first, coeffs.end() - trailing_zeros);
                offset += leading_zeros;
            } else {
                // All zeros
                coeffs.clear();
            }
            
            // Initialize cache directly with the coefficients
            this->cache = std::move(coeffs);
            this->cache_offset = offset;
            this->initialized = true;
            this->all_cached = true;
        }
        
        std::optional<int> valuation() const override {
            if(this->cache.empty()) return std::nullopt;
            return this->cache_offset;
        }
        
        std::optional<int> degree() const override {
            if(this->cache.empty()) return std::nullopt;
            return this->cache_offset + (int)this->cache.size() - 1;
        }
        
        T coeff_lazy(int k) const override {
            int idx = k - this->cache_offset;
            if(idx < 0 || idx >= (int)this->cache.size()) {
                return T(0);
            }
            return this->cache[idx];
        }
        
        T coeff(int k) const override {
            return coeff_lazy(k);
        }
    };
    
    // Base class for unary operations
    template<typename T>
    struct unary_provider : provider<T> {
        std::shared_ptr<provider<T>> operand;
        
        unary_provider(std::shared_ptr<provider<T>> operand)
            : operand(std::move(operand)) {}
        
        virtual T transform(T const& a) const = 0;
        
        std::optional<int> valuation() const override {
            return operand->valuation();
        }
        
        std::optional<int> degree() const override {
            return operand->degree();
        }
        
        T coeff_lazy(int k) const override {
            return transform(operand->coeff_lazy(k));
        }
        
        T coeff(int k) const {
            return transform(operand->coeff(k));
        }
    };
    
    // Base class for binary operations that may have non-trivial bounds
    template<typename T>
    struct binary_provider : provider<T> {
        enum class bound_type { valuation, degree };
        
        std::shared_ptr<provider<T>> lhs, rhs;
        
        binary_provider(std::shared_ptr<provider<T>> lhs, std::shared_ptr<provider<T>> rhs)
            : lhs(std::move(lhs)), rhs(std::move(rhs)) {}
        
        virtual T combine(T const& a, T const& b) const = 0;
        
        std::optional<int> valuation() const override {
            auto lv = lhs->valuation();
            auto rv = rhs->valuation();
            if(!lv || !rv) return std::nullopt;
            
            // If valuations are distinct, result is the minimum
            if(*lv != *rv) {
                return std::min(*lv, *rv);
            }
            
            // Same valuation - check if they cancel
            T val = combine(lhs->coeff_lazy(*lv), rhs->coeff_lazy(*lv));
            if(val != T(0)) {
                return *lv;
            }
            
            // They cancel - undefined
            return std::nullopt;
        }
        
        std::optional<int> degree() const override {
            auto ld = lhs->degree();
            auto rd = rhs->degree();
            if(!ld || !rd) return std::nullopt;
            
            // If degrees are distinct, result is the maximum
            if(*ld != *rd) {
                return std::max(*ld, *rd);
            }
            
            // Same degree - check if they cancel
            T val = combine(lhs->coeff_lazy(*ld), rhs->coeff_lazy(*ld));
            if(val != T(0)) {
                return *ld;
            }
            
            // They cancel - undefined
            return std::nullopt;
        }
        
        T coeff_lazy(int k) const override {
            return combine(lhs->coeff_lazy(k), rhs->coeff_lazy(k));
        }
        
        T coeff(int k) const {
            return combine(lhs->coeff(k), rhs->coeff(k));
        }
    };
    
    // Addition provider
    template<typename T>
    struct add_provider : binary_provider<T> {
        using binary_provider<T>::binary_provider;
        
        T combine(T const& a, T const& b) const override {
            return a + b;
        }
    };
    
    // Subtraction provider
    template<typename T>
    struct subtract_provider : binary_provider<T> {
        using binary_provider<T>::binary_provider;
        
        T combine(T const& a, T const& b) const override {
            return a - b;
        }
    };
    
    // Negation provider
    template<typename T>
    struct negate_provider : unary_provider<T> {
        using unary_provider<T>::unary_provider;
        
        T transform(T const& a) const override {
            return -a;
        }
    };
    
    // Scalar multiplication provider
    template<typename T>
    struct scale_provider : unary_provider<T> {
        T scalar;
        
        scale_provider(std::shared_ptr<provider<T>> operand, T scalar)
            : unary_provider<T>(std::move(operand)), scalar(scalar) {}
        
        T transform(T const& a) const override {
            return a * scalar;
        }
    };
    
    // Multiplication provider (Cauchy product)
    template<typename T>
    struct multiply_provider : provider<T> {
        std::shared_ptr<provider<T>> lhs, rhs;
        
        multiply_provider(std::shared_ptr<provider<T>> lhs, std::shared_ptr<provider<T>> rhs)
            : lhs(std::move(lhs)), rhs(std::move(rhs)) {}
        
        std::optional<int> valuation() const override {
            auto lv = lhs->valuation();
            auto rv = rhs->valuation();
            if(lv && rv) return *lv + *rv;
            return std::nullopt;
        }
        
        std::optional<int> degree() const override {
            auto ld = lhs->degree();
            auto rd = rhs->degree();
            if(ld && rd) return *ld + *rd;
            return std::nullopt;
        }
        
        T coeff_lazy(int k) const override {
            T result = T(0);
            for(int j = 0; j <= k; j++) {
                result += lhs->coeff(j) * rhs->coeff(k - j);
            }
            return result;
        }
    };
    
    // Main Laurent series class
    template<typename T>
    struct laurent {
        std::shared_ptr<provider<T>> impl;
        
        laurent() : impl(std::make_shared<constant_provider<T>>(T(0), 0)) {}
        
        laurent(T value, int offset = 0) 
            : impl(std::make_shared<constant_provider<T>>(value, offset)) {}
        
        laurent(big_vector<T> coeffs, int offset = 0)
            : impl(std::make_shared<polynomial_provider<T>>(std::move(coeffs), offset)) {}
        
        laurent(std::shared_ptr<provider<T>> impl) : impl(std::move(impl)) {}
        
        // Get k-th coefficient (delegates to provider's caching)
        T operator[](int k) const {
            return impl->get(k);
        }
        
        // Arithmetic operations
        laurent operator-() const {
            return std::make_shared<negate_provider<T>>(impl);
        }
        
        laurent operator+(const laurent& other) const {
            return std::make_shared<add_provider<T>>(impl, other.impl);
        }
        
        laurent operator-(const laurent& other) const {
            return std::make_shared<subtract_provider<T>>(impl, other.impl);
        }
        
        laurent operator*(const laurent& other) const {
            return std::make_shared<multiply_provider<T>>(impl, other.impl);
        }
        
        laurent& operator+=(const laurent& other) {
            return *this = *this + other;
        }
        
        laurent& operator-=(const laurent& other) {
            return *this = *this - other;
        }
        
        laurent& operator*=(const laurent& other) {
            return *this = *this * other;
        }
        
        // Scalar multiplication
        laurent operator*(T const& scalar) const {
            return std::make_shared<scale_provider<T>>(impl, scalar);
        }
        
        laurent& operator*=(T const& scalar) {
            return *this = *this * scalar;
        }
    };
    
    template<typename T>
    laurent<T> operator*(T const& scalar, laurent<T> const& series) {
        return series * scalar;
    }
}

#endif // CP_ALGO_MATH_LAURENT_HPP
