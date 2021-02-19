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

/// Branch scalar addition
[[nodiscard]] inline kln::rotor KLN_VEC_CALL operator+(float a, kln::branch b) noexcept
{
    kln::rotor out;
    out.p1_ = _mm_add_ss(b.p1_, _mm_set_ss(a));
    return out;
}

/// Branch times rotor
[[nodiscard]] inline kln::rotor KLN_VEC_CALL operator*(kln::branch a, kln::rotor b) noexcept
{
    kln::rotor out = as_rotor(a);
    return out*b;
}
