#ifndef CP_ALGO_UTIL_COMPLEX_HPP
#define CP_ALGO_UTIL_COMPLEX_HPP
#include <iostream>
#include <cmath>
namespace cp_algo {
    // Custom implementation, since std::complex is UB on non-floating types
    template<typename T>
    struct complex {
        using value_type = T;
        T x, y;
        constexpr complex() {}
        constexpr complex(T x): x(x), y() {}
        constexpr complex(T x, T y): x(x), y(y) {}
        complex& operator *= (T t) {x *= t; y *= t; return *this;}
        complex& operator /= (T t) {x /= t; y /= t; return *this;}
        complex operator * (T t) const {return complex(*this) *= t;}
        complex operator / (T t) const {return complex(*this) /= t;}
        complex& operator += (complex t) {x += t.x; y += t.y; return *this;}
        complex& operator -= (complex t) {x -= t.x; y -= t.y; return *this;}
        complex operator * (complex t) const {return {x * t.x - y * t.y, x * t.y + y * t.x};}
        complex operator / (complex t) const {return *this * t.conj() / t.norm();}
        complex operator + (complex t) const {return complex(*this) += t;}
        complex operator - (complex t) const {return complex(*this) -= t;}
        complex& operator *= (complex t) {return *this = *this * t;}
        complex& operator /= (complex t) {return *this = *this / t;}
        complex operator - () const {return {-x, -y};}
        complex conj() const {return {x, -y};}
        T norm() const {return x * x + y * y;}
        T abs() const {return std::sqrt(norm());}
        T real() const {return x;}
        T imag() const {return y;}
        T& real() {return x;}
        T& imag() {return y;}
        static constexpr complex polar(T r, T theta) {return {r * cos(theta), r * sin(theta)};}
        auto operator <=> (complex const& t) const = default;
    };
    template<typename T>
    complex<T> operator * (auto x, complex<T> y) {return y *= x;}
    template<typename T> complex<T> conj(complex<T> x) {return x.conj();}
    template<typename T> T norm(complex<T> x) {return x.norm();}
    template<typename T> T abs(complex<T> x) {return x.abs();}
    template<typename T> T& real(complex<T> &x) {return x.real();}
    template<typename T> T& imag(complex<T> &x) {return x.imag();}
    template<typename T> T real(complex<T> const& x) {return x.real();}
    template<typename T> T imag(complex<T> const& x) {return x.imag();}
    template<typename T>
    constexpr complex<T> polar(T r, T theta) {
        return complex<T>::polar(r, theta);
    }
    template<typename T>
    std::ostream& operator << (std::ostream &out, complex<T> x) {
        return out << x.real() << ' ' << x.imag();
    }
}
#endif // CP_ALGO_UTIL_COMPLEX_HPP
