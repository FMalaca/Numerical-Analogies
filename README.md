# Functions for numerical analogy

**Definition:**
An analogy in $`p`$ holds between four numbers, $`(a, b, c, d)`$, when
the generalized mean in $`p`$ of the extremes
is equal to
the generalized mean in $`p`$ of the means, i.e.,
```math
m_p(a, d) = m_p(b, c)
```
where
```math
m_p(a, d) = \begin{cases}
            \left(\frac{1}{2}(a^p + d^p)\right)^{1/p} \text{, if } p \text{ is non-null real} \\
              a^{1/2}  d^{1/2} \text{, if } p = 0 \\
             \min \{|a|, |d|\} \text{, if } p = -\infty \\
             \max \{|a|, |d|\} \text{, if } p = +\infty,
        \end{cases}
```
and the exponentiations are considered using the principal branch of the complex logarithm, 
which is defined through the principal argument.

**Library purpose:**
The present library computes the analogical powers of a given quadruple of non-null numbers, $`(a, b, c, d)`$.
