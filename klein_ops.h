#include <klein/klein.hpp>
#include <klein/detail/geometric_product.hpp>

/// Branch cast as rotor
[[nodiscard]] inline kln::rotor KLN_VEC_CALL as_rotor(kln::branch b) noexcept
{
    kln::rotor out;
    out.p1_ = b.p1_;
    return out;
}


/// Rotor cast as branch
[[nodiscard]] inline kln::branch KLN_VEC_CALL as_branch(kln::rotor b) noexcept
{
    kln::branch out;
    out.p1_ = b.p1_;
    return out;
}


/// Line cast as motor
[[nodiscard]] inline kln::motor KLN_VEC_CALL as_motor(kln::line b) noexcept
{
    kln::motor out;
    out.p1_ = b.p1_;
    out.p2_ = b.p2_;
    return out;
}

/// Motor cast as line
[[nodiscard]] inline kln::line KLN_VEC_CALL as_line(kln::motor b) noexcept
{
    kln::line out;
    out.p1_ = b.p1_;
    out.p2_ = b.p2_;
    return out;
}


/// Branch scalar addition
[[nodiscard]] inline kln::rotor KLN_VEC_CALL operator+(float a, kln::branch b) noexcept
{
    kln::rotor out;
    out.p1_ = _mm_add_ss(b.p1_, _mm_set_ss(a));
    return out;
}

/// Motor scalar addition
[[nodiscard]] inline kln::motor KLN_VEC_CALL operator+(float a, kln::motor b) noexcept
{
    kln::motor out;
    out.p1_ = _mm_add_ss(b.p1_, _mm_set_ss(a));
    out.p2_ = b.p2_;
    return out;
}

/// Line scalar addition
[[nodiscard]] inline kln::motor KLN_VEC_CALL operator+(float a, kln::line b) noexcept
{
    return a + as_motor(b);
}

/// Branch times rotor
[[nodiscard]] inline kln::rotor KLN_VEC_CALL operator*(kln::branch a, kln::rotor b) noexcept
{
    kln::rotor out = as_rotor(a);
    return out*b;
}

/// Line times rotor
[[nodiscard]] inline kln::motor KLN_VEC_CALL operator*(kln::line a, kln::rotor b) noexcept
{
    kln::motor out = as_motor(a);
    return out*b;
}

/// Line times motor
[[nodiscard]] inline kln::motor KLN_VEC_CALL operator*(kln::line a, kln::motor b) noexcept
{
    kln::motor out = as_motor(a);
    return out*b;
}
