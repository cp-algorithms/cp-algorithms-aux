#ifndef CP_ALGO_MATH_LAURENT_HPP
#define CP_ALGO_MATH_LAURENT_HPP
#include "../util/big_alloc.hpp"
#include <memory>
#include <optional>
#include <vector>
#include <cassert>
#include <algorithm>
#include <type_traits>
#include <bit>

#include "cvector.hpp"
#include "convolution.hpp"

namespace cp_algo::math {
    // Base provider interface for lazy coefficient evaluation
    template<typename T>
    struct provider {
        mutable big_vector<T> cache;
        mutable int cache_offset = 0; // Index of first cached coefficient
        mutable bool initialized = false;
        mutable bool all_cached = false; // True if all non-zero coeffs are cached
        
        virtual ~provider() = default;
        virtual int offset() const { return 0; }
        
        // Returns true if this provider requires lazy evaluation (coefficients must be
        // computed in order). False means dependencies can be bulk-cached for FFT.
        // Examples: multiply needs lazy eval, add/subtract/negate/scale don't.
        virtual bool needs_lazy_eval() const { return false; }
        
        // Compute k-th coefficient lazily without caching
        virtual T coeff_lazy(int k) const = 0;
        
        // Double the number of known coefficients (or cache all for finite series)
        // Default: use coeff_lazy, but can be overridden for efficiency (e.g., FFT)
        virtual void double_up() const {
            int old_size = cache.size();
            int new_size = old_size == 0 ? 1 : 2 * old_size;
            
            cache.resize(new_size);
            for(int i = old_size; i < new_size; i++) {
                cache[i] = coeff_lazy(cache_offset + i);
            }
        }
        
        // Get coefficient with caching and doubling (default implementation)
        virtual T coeff(int k) const {
            if(!initialized) {
                cache_offset = offset();
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
            
            if(needs_lazy_eval()) {
                // Sequentially extend cache to the requested index
                while(idx >= (int)cache.size() && !all_cached) {
                    int next_k = cache_offset + (int)cache.size();
                    cache.push_back(coeff_lazy(next_k));
                }
            } else {
                // Extend cache by doubling until we have enough
                while(idx >= (int)cache.size() && !all_cached) {
                    double_up();
                }
            }
            
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
        
        int offset() const override {
            return offset;
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
        
        int offset() const override {
            return this->cache_offset;
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
        
        int offset() const override {
            return operand->offset();
        }
        
        T coeff_lazy(int k) const override {
            return transform(operand->coeff_lazy(k));
        }
        
        T coeff(int k) const {
            return transform(operand->coeff(k));
        }
    };
    
    // Base class for binary operations
    template<typename T>
    struct binary_provider : provider<T> {
        std::shared_ptr<provider<T>> lhs, rhs;
        
        binary_provider(std::shared_ptr<provider<T>> lhs, std::shared_ptr<provider<T>> rhs)
            : lhs(std::move(lhs)), rhs(std::move(rhs)) {}
        
        virtual T combine(T const& a, T const& b) const = 0;
        
        int offset() const override {
            return std::min(lhs->offset(), rhs->offset());
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
        
        int offset() const override {
            return lhs->offset() + rhs->offset();
        }
        
        bool needs_lazy_eval() const override {
            return lhs->needs_lazy_eval() || rhs->needs_lazy_eval();
        }
        
        T coeff_lazy(int k) const override {
            int n = k - offset();
            if(n < 0) return T(0);
            T result = T(0);
            bool lazy_lhs = lhs->needs_lazy_eval();
            bool lazy_rhs = rhs->needs_lazy_eval();
            for(int j = 0; j <= n; j++) {
                int i_l = lhs->offset() + j;
                int i_r = rhs->offset() + (n - j);
                auto a = lazy_lhs ? lhs->coeff(i_l) : lhs->coeff_lazy(i_l);
                auto b = lazy_rhs ? rhs->coeff(i_r) : rhs->coeff_lazy(i_r);
                result += a * b;
            }
            return result;
        }

        void double_up() const override {
            int old_size = this->cache.size();
            int new_size = old_size == 0 ? 1 : 2 * old_size;

            // Lazy path: compute the next coefficient sequentially
            if(needs_lazy_eval()) {
                int k = this->cache_offset + old_size;
                this->cache.push_back(coeff_lazy(k));
                return;
            }

            // Ensure operands have enough cached coefficients for the prefix we need
            int lhs_need = lhs->offset() + new_size - 1;
            int rhs_need = rhs->offset() + new_size - 1;
            lhs->coeff(lhs_need);
            rhs->coeff(rhs_need);

            // Build aligned prefixes starting at operand offsets
            std::vector<T, big_alloc<T>> la(new_size), rb(new_size);
            for(int i = 0; i < new_size; i++) {
                la[i] = lhs->coeff(lhs->offset() + i);
                rb[i] = rhs->coeff(rhs->offset() + i);
            }

            this->cache.resize(new_size);
            convolution_prefix(la, rb, new_size);
            for(int i = old_size; i < new_size && i < (int)la.size(); i++) {
                this->cache[i] = la[i];
            }

            // If both operands are fully cached and we reached their total length, mark as done
            if(lhs->all_cached && rhs->all_cached) {
                size_t total_len = lhs->cache.size() + rhs->cache.size() - 1;
                if((size_t)new_size >= total_len) {
                    this->cache.resize(total_len);
                    this->all_cached = true;
                }
            }
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
