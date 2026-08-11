# What is a Fourier Transform?

- Good 3B1B video: https://youtu.be/spUNpyF58BY?si=OCtgf0Zbp54sCz5G
- Fun Veritasium video to get you excited about DFT: https://youtu.be/nmgFG7PUHfo?si=l93nWHY8KPE19dsi

## Fourier Transform vs. Fourier Series

- **Series**: Given a periodic signal, with period $T$ and fundamental frequency $f_0 = 1/T$, you can represent it as a discrete sum of periodic functions (sines and cosines) that have integer multiples of the fundamental frequency.

- **Transform**: Given a continuous, aperiodic signal, you can represent it as a continuous mixture of all possible real frequencies.

A Fourier transform can be understood as the limiting case where the period of the Fourier series approaches infinity.

## Formulae for Fourier Transform & Series

### Continuous Fourier Series (CFS)

For a continuous periodic signal $x(t)$ with period $T$ and fundamental angular frequency $\omega_0 = \frac{2\pi}{T}$:

$$
x(t) = \sum_{k=-\infty}^{\infty} c_k e^{jk\omega_0 t}
$$

Where

$$
c_k = \frac{1}{T} \int_{-T/2}^{T/2} x(t)e^{-jk\omega_0t},dt
$$

### Continuous Fourier Transform (CFT)

For a continuous signal $x(t)$:

$$
X(f) = \int_{-\infty}^{\infty} x(t)e^{-j2\pi ft},dt
$$

$$
x(t) = \int_{-\infty}^{\infty} X(f)e^{j2\pi ft},df
$$

# What is a DISCRETE Fourier Transform?

A Discrete Fourier Transform is a Fourier transform computed from a discrete, finite-length sample of a signal. The DFT can also be interpreted as finding the Fourier series coefficients of the periodic sequence obtained by repeating those $N$ samples indefinitely.

### Formula for DFT

$$
X[k] = \sum_{n=0}^{N-1} x[n]e^{-j\frac{2\pi}{N}kn},
\qquad k=0,1,\ldots,N-1
$$

$$
x[n] = \frac{1}{N}\sum_{k=0}^{N-1}X[k]e^{j\frac{2\pi}{N}kn},
\qquad n=0,1,\ldots,N-1
$$

Whilst the claim that a "series" is a better name may be true, a discrete Fourier transform is still a transform. This becomes clearer when we look at the linear algebraic definition of the DFT. A linear algebraic representation of a DFT is a matrix multiplication of the following form:

$$
\begin{bmatrix}
    X[0]\\
    X[1]\\
    X[2]\\
    \vdots\\
    X[N-1]
\end{bmatrix}
=
\begin{bmatrix}
    1 & 1 & 1 & \dots & 1\\
    1 & \omega_N^1 & \omega_N^2 & \dots & \omega_N^{N-1}\\
    1 & \omega_N^2 & \omega_N^4 & \dots & \omega_N^{2(N-1)}\\
    \vdots & \vdots & \vdots & \ddots & \vdots\\
    1 & \omega_N^{N-1} & \omega_N^{2(N-1)} & \dots & \omega_N^{(N-1)(N-1)}
\end{bmatrix}
\begin{bmatrix}
    x[0]\\
    x[1]\\
    x[2]\\
    \vdots\\
    x[N-1]
\end{bmatrix}
$$

where

$$
\omega_N = e^{-j\frac{2\pi}{N}}
$$

is the principal $N$-th root of unity.

Take a moment to examine the structure of this matrix, and see how the summation notation and indexing was translated.

## Linear Algebra & Parseval's Theorem

A matrix multiplication can be understood as a transformation between two different vector spaces. In the case of the DFT, those two spaces are often interpreted as the time-domain and frequency-domain representations of the same signal.

This matrix $F_N$ is known as a **DFT matrix** or **Fourier matrix**, and is a special type of Vandermonde matrix. A crucial characteristic of $F_N$ is that its columns are mutually orthogonal complex vectors.

Thus, scaling $F_N$ by $1/\sqrt{N}$ produces a unitary matrix. This means that the transformation preserves inner products and therefore preserves the Euclidean norm.

When using the standard unscaled $F_N$, this preservation of length manifests as _Parseval's Theorem_:

$$
\sum_{k=0}^{N-1}|X[k]|^2
=
N\sum_{n=0}^{N-1}|x[n]|^2
$$

or equivalently,

$$
\sum_{n=0}^{N-1}|x[n]|^2
=
\frac{1}{N}\sum_{k=0}^{N-1}|X[k]|^2
$$

Taking inspiration from the fact that the DFT is a linear algebraic operation, we can also understand the DFT in terms of polynomials, like we do with many other linear algebraic operations. This becomes a nice segue to our next topic, the Fast Fourier Transform (FFT).

# FAST Fourier Transform (FFT)

The FFT is a family of algorithms that computes the DFT much faster by exploiting the structure and symmetries of the DFT. Instead of the naive DFT algorithm, which takes $O(N^2)$ time, common FFT algorithms can compute the DFT in $O(N\log N)$ time.

The FFT does **not** produce an approximation of the DFT. It computes the same DFT values (up to numerical precision) while avoiding redundant computation.

## Polynomial Evaluation View

Let's represent the input signal as the coefficients of a degree-$(N-1)$ polynomial:

$$
P(z)
=
x[0] + x[1]z + x[2]z^2 + \dots + x[N-1]z^{N-1}
=
\sum_{n=0}^{N-1}x[n]z^n
$$

Computing the DFT is equivalent to evaluating $P(z)$ at the $N$ distinct $N$-th roots of unity,

$$
z_k = \omega_N^k = e^{-j\frac{2\pi}{N}k}
$$

so that

$$
X[k] = P(\omega_N^k)
$$

These roots of unity sit symmetrically around the unit circle in the complex plane, giving us two key properties.

### Halving / Squaring Property

For even $N$,

$$
(\omega_N^k)^2
=
\left(e^{-j\frac{2\pi}{N}k}\right)^2
=
e^{-j\frac{2\pi}{N/2}k}
=
\omega_{N/2}^k
$$

Squaring an $N$-th root of unity produces an $(N/2)$-th root of unity.

### Negation / Symmetry Property

For even $N$,

$$
\omega_N^{k+N/2}
=
e^{-j\frac{2\pi}{N}(k+N/2)}
=
e^{-j\frac{2\pi k}{N}}e^{-j\pi}
=
-\omega_N^k
$$

Thus, roots opposite each other across the unit circle are exact additive inverses.

## Even / Odd Polynomial Decomposition

We can split $P(z)$ into even- and odd-indexed terms:

$$
P(z)
=
\underbrace{
\left(x[0]+x[2]z^2+\dots\right)
}_{P_{\text{even}}(z^2)}
+
z\cdot
\underbrace{
\left(x[1]+x[3]z^2+\dots\right)
}_{P_{\text{odd}}(z^2)}
$$

Evaluating this at the root $\omega_N^k$ and using the squaring property,

$$
(\omega_N^k)^2=\omega_{N/2}^k
$$

gives

$$
X[k]
=
P_{\text{even}}(\omega_{N/2}^k)
+
\omega_N^kP_{\text{odd}}(\omega_{N/2}^k)
$$

### Reusing Computations (The Butterfly Operation)

Using the negation property,

$$
\omega_N^{k+N/2}=-\omega_N^k
$$

we can compute both $X[k]$ and $X[k+N/2]$ using the exact same sub-evaluations:

$$
\begin{aligned}
X[k]
&=
P_{\text{even}}(\omega_{N/2}^k)
+
\omega_N^kP_{\text{odd}}(\omega_{N/2}^k)
\\
X[k+N/2]
&=
P_{\text{even}}(\omega_{N/2}^k)
-
\omega_N^kP_{\text{odd}}(\omega_{N/2}^k)
\end{aligned}
$$

Notice that

$$
P_{\text{even}}(\omega_{N/2}^k)
$$

and

$$
\omega_N^kP_{\text{odd}}(\omega_{N/2}^k)
$$

only need to be computed once, and are then added or subtracted to yield the two frequency outputs. This pair of operations is called a **butterfly**.

We can recursively apply this even/odd decomposition, eventually reducing the original $N$-point DFT into many small DFTs. When $N$ is a power of two, this produces a binary recursion tree with $\log_2 N$ levels, giving the familiar $O(N\log_2N)$ complexity.

## Example: 8-Point FFT

Let's see what this actually looks like for an $N=8$ FFT.

The naive DFT would calculate each of the 8 outputs using all 8 input samples:

$$
X[k]=\sum_{n=0}^{7}x[n]\omega_8^{kn}
$$

giving $8\times8=64$ complex multiply/add contributions.

The radix-2 FFT instead repeatedly splits the input into even and odd indices:

$$
[x_0,x_1,x_2,x_3,x_4,x_5,x_6,x_7]
$$

becomes

$$
[x_0,x_2,x_4,x_6]
\qquad
[x_1,x_3,x_5,x_7]
$$

and each of those is split again:

$$
[x_0,x_4],\ [x_2,x_6],\ [x_1,x_5],\ [x_3,x_7]
$$

Finally, the 2-point transforms are just additions and subtractions:

$$
\begin{aligned}
A_0 &= x_0+x_4\\
A_1 &= x_0-x_4\\
A_2 &= x_2+x_6\\
A_3 &= x_2-x_6
\end{aligned}
$$

and similarly for the odd-indexed branch.

These results are then combined through three stages of butterflies. Schematically:

$$
8\text{-point DFT}
\rightarrow
2\times4\text{-point DFTs}
\rightarrow
4\times2\text{-point DFTs}
\rightarrow
8\times1\text{-point DFTs}
$$

There are

$$
\log_2(8)=3
$$

stages, with $8/2=4$ butterflies per stage, for a total of

$$
4\times3=12
$$

butterflies.

The important idea is that **the FFT does not calculate every output independently**. It computes smaller transforms once and reuses their results.

This is the fundamental source of the FFT's computational advantage.

## Power of Twos

The radix-2 decomposition requires $N$ to be even at each recursive step. Therefore, it naturally works when $N$ is repeatedly divisible by 2, i.e. when

$$
N=2^m
$$

for some integer $m$.

However, FFTs can also be constructed for many non-power-of-two lengths by using other factorizations.

For example,

$$
N=12=3\times4
$$

can be decomposed using a combination of radix-3 and radix-2/radix-4 stages. Likewise,

$$
N=1000=2^3\times5^3
$$

can be implemented using mixed-radix FFT algorithms involving radix-2 and radix-5 stages.

This gives us several options for non-power-of-two FFT lengths:

- **Zero-padding:** Extend the input to a convenient FFT size, often a power of two.
- **Mixed-radix FFT:** Factor $N$ into smaller integers and use different radix stages.
- **Prime-length algorithms:** For prime $N$, algorithms such as **Rader's algorithm** or **Bluestein's algorithm** can be used.
- **Direct DFT:** For very small or unusual transform sizes, simply computing the DFT directly can sometimes be the simplest hardware solution.

### Why Zero-Padding Is Common in Hardware

If the desired transform length is awkward but approximate frequency resolution is acceptable, zero-padding can be attractive.

For example, suppose we have 1000 samples:

$$
x[0],x[1],\ldots,x[999]
$$

We could append 24 zeros and perform a 1024-point FFT:

$$
1000\rightarrow1024
$$

This allows us to use a highly optimized radix-2/radix-4/radix-8 architecture.

However, **zero-padding does not increase the underlying frequency resolution of the original data**. It simply evaluates the DTFT more densely in frequency. The actual ability to distinguish closely spaced frequencies is still fundamentally determined by the observation duration and sampling characteristics. Watch this MATLAB video if this confuses you: https://youtu.be/QmgJmh2I3Fw?si=6impZFgYffkJ3QVq

Therefore, zero-padding is an implementation convenience or a way to obtain a denser frequency grid, not a way to create additional information.

### Hardware Considerations

In an FPGA, FFT implementations are typically constructed from hardware primitives such as:

- DSP/multiplier blocks for complex multiplication
- BRAM/URAM for storing intermediate data
- adders/subtractors for butterfly operations
- control and address-generation logic
- pipelining and streaming infrastructure

FPGA vendors often provide highly optimized FFT IP that supports a particular set of transform lengths and radix architectures. These implementations may favor powers of two because they allow particularly regular memory-access patterns, butterfly structures, and pipeline architectures.

If an exact non-power-of-two transform is required, the designer can instead implement a mixed-radix architecture, use an algorithm such as Rader or Bluestein, or choose another architecture depending on the required throughput, latency, memory usage, DSP utilization, and numerical precision.

For example, a 1000-point FFT could be implemented as a mixed-radix transform because

$$
1000=2^3\times5^3
$$

rather than padding to 1024. Whether this is actually preferable depends on the hardware constraints: a 1024-point radix-2/4/8 FFT may use more arithmetic operations but have a much simpler and more regular architecture, while a 1000-point mixed-radix FFT avoids computing the extra 24 points at the cost of more complicated radix-5 stages and potentially more complicated control and memory access.

**Now that we've covered the theoretical background of DFTs and FFTs, let's see their implementations in practice and compare!**
