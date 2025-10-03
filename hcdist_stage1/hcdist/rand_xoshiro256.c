/* Modifications by Eric A. Welsh (2018) simply repackages several functions
 * together and renames/reformats them to more easily fit into other codebases,
 * as well as implementing an initialization routine and the recommended
 * uint64 -> double conversion function.  I release this code under the same
 * license as used originally below.
 *
 * I have also added a few functions for generating non-uniform distributions
 *
 *
 *
 * Copyright notice accompanying original xoshiro256 code:
 *
 * Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 *
 * To the extent possible under law, the author has dedicated all copyright
 * and related and neighboring rights to this software to the public domain
 * worldwide. This software is distributed without any warranty.
 *
 * See <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 *
 *
 * Copyright notice accompanying original splitmix64 code:
 *
 * Written in 2015 by Sebastiano Vigna (vigna@acm.org)
 *
 * To the extent possible under law, the author has dedicated all copyright
 * and related and neighboring rights to this software to the public domain
 * worldwide. This software is distributed without any warranty.
 *
 * See <http://creativecommons.org/publicdomain/zero/1.0/>.
 */

#include <stdint.h>
#include <math.h>

static uint64_t states_ss[4]; /* states for xorshiro256** */
static uint64_t states_pp[4]; /* states for xorshiro256++ */
static uint64_t states_p[4];  /* states for xorshiro256+  */
static uint64_t state_sm64;   /* The state can be seeded with any value. */



/* This is a fixed-increment version of Java 8's SplittableRandom generator
 * See http://dx.doi.org/10.1145/2714064.2660195 and 
 *
 * http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html
 *
 * It is a very fast generator passing BigCrush, and it can be useful if
 * for some reason you absolutely want 64 bits of state; otherwise, we
 * rather suggest to use a xoroshiro128+ (for moderately parallel
 * computations) or xorshift1024* (for massively parallel computations)
 * generator.
 */
uint64_t splitmix64(void)
{
    uint64_t z = (state_sm64 += 0x9e3779b97f4a7c15);
    
    z =    (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z =    (z ^ (z >> 27)) * 0x94d049bb133111eb;
    
    return (z ^ (z >> 31));
}


/* rotate function for the xoshiro256 generators */
static inline uint64_t rotl64(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}


/* This is xoshiro256** 1.0, one of our all-purpose, rock-solid generators.
 * It has excellent (sub-ns) speed, a state (256 bits) that is large enough
 * for any parallel application, and it passes all tests we are aware of.
 *
 * For generating just floating-point numbers, xoshiro256+ is even faster.
 *
 * The state must be seeded so that it is not everywhere zero. If you have
 * a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 * output to fill s.
 */
uint64_t xoshiro256_starstar(void)
{
    const uint64_t result = rotl64(states_ss[1] * 5, 7) * 9;

    const uint64_t t = states_ss[1] << 17;

    states_ss[2] ^= states_ss[0];
    states_ss[3] ^= states_ss[1];
    states_ss[1] ^= states_ss[2];
    states_ss[0] ^= states_ss[3];

    states_ss[2] ^= t;

    states_ss[3]  = rotl64(states_ss[3], 45);

    return result;
}


/* This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators.
 * It has excellent (sub-ns) speed, a state (256 bits) that is large enough
 * for any parallel application, and it passes all tests we are aware of.
 *
 * For generating just floating-point numbers, xoshiro256+ is even faster.
 *
 * The state must be seeded so that it is not everywhere zero. If you have
 * a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 * output to fill s.
 */
uint64_t xoshiro256_plusplus(void)
{
    const uint64_t result = rotl64(states_pp[0] + states_pp[3], 23) +
                                   states_pp[0];

    const uint64_t t = states_pp[1] << 17;

    states_pp[2] ^= states_pp[0];
    states_pp[3] ^= states_pp[1];
    states_pp[1] ^= states_pp[2];
    states_pp[0] ^= states_pp[3];

    states_pp[2] ^= t;

    states_pp[3]  = rotl64(states_pp[3], 45);

    return result;
}


/* This is xoshiro256+ 1.0, our best and fastest generator for floating-point
 * numbers. We suggest to use its upper bits for floating-point
 * generation, as it is slightly faster than xoshiro256**. It passes all
 * tests we are aware of except for the lowest three bits, which might
 * fail linearity tests (and just those), so if low linear complexity is
 * not considered an issue (as it is usually the case) it can be used to
 * generate 64-bit outputs, too.
 *
 * We suggest to use a sign test to extract a random Boolean value, and
 * right shifts to extract subsets of bits.
 *
 * The state must be seeded so that it is not everywhere zero. If you have
 * a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 * output to fill s.
 */
uint64_t xoshiro256_plus(void)
{
    const uint64_t result = states_p[0] + states_p[3];
    const uint64_t t = states_p[1] << 17;

    states_p[2] ^= states_p[0];
    states_p[3] ^= states_p[1];
    states_p[1] ^= states_p[2];
    states_p[0] ^= states_p[3];

    states_p[2] ^= t;

    states_p[3] = rotl64(states_p[3], 45);

    return result;
}


/* This is the jump function for the generator. It is equivalent
 * to 2^128 calls to xoshiro256_starstar(); it can be used to generate 2^128
 * non-overlapping subsequences for parallel computations.
 */
void jump_xoshiro256_starstar(void)
{
    static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba,
                                     0xd5a61266f0c9392c,
                                     0xa9582618e03fc9aa,
                                     0x39abdc4529b1661c };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for(i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
    {
        for(b = 0; b < 64; b++)
        {
            if (JUMP[i] & UINT64_C(1) << b)
            {
                s0 ^= states_ss[0];
                s1 ^= states_ss[1];
                s2 ^= states_ss[2];
                s3 ^= states_ss[3];
            }
            xoshiro256_starstar();
        }
    }
        
    states_ss[0] = s0;
    states_ss[1] = s1;
    states_ss[2] = s2;
    states_ss[3] = s3;
}


/* This is the long-jump function for the generator. It is equivalent to
 * 2^192 calls to xoshiro256_starstar(); it can be used to generate 2^64
 * starting points, from each of which jump_xoshiro256_starstar() will
 * generate 2^64 non-overlapping subsequences for parallel distributed
 * computations.
 */
void long_jump_xoshiro256_starstar(void)
{
    static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf,
                                          0xc5004e441c522fb3,
                                          0x77710069854ee241,
                                          0x39109bb02acbe635 };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for(i = 0; i < sizeof(LONG_JUMP) / sizeof(*LONG_JUMP); i++)
    {
        for(b = 0; b < 64; b++)
        {
            if (LONG_JUMP[i] & UINT64_C(1) << b)
            {
                s0 ^= states_ss[0];
                s1 ^= states_ss[1];
                s2 ^= states_ss[2];
                s3 ^= states_ss[3];
            }
            xoshiro256_starstar();
        }
    }
        
    states_ss[0] = s0;
    states_ss[1] = s1;
    states_ss[2] = s2;
    states_ss[3] = s3;
}


/* This is the jump function for the generator. It is equivalent
 * to 2^128 calls to xoshiro256_plusplus(); it can be used to generate 2^128
 * non-overlapping subsequences for parallel computations.
 */
void jump_xoshiro256_plusplus(void)
{
    static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba,
                                     0xd5a61266f0c9392c,
                                     0xa9582618e03fc9aa,
                                     0x39abdc4529b1661c };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for(i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
    {
        for(b = 0; b < 64; b++)
        {
            if (JUMP[i] & UINT64_C(1) << b)
            {
                s0 ^= states_pp[0];
                s1 ^= states_pp[1];
                s2 ^= states_pp[2];
                s3 ^= states_pp[3];
            }
            xoshiro256_plusplus();
        }
    }
        
    states_pp[0] = s0;
    states_pp[1] = s1;
    states_pp[2] = s2;
    states_pp[3] = s3;
}


/* This is the long-jump function for the generator. It is equivalent to
 * 2^192 calls to xoshiro256_plusplus(); it can be used to generate 2^64
 * starting points, from each of which jump_xoshiro256_plusplus() will
 * generate 2^64 non-overlapping subsequences for parallel distributed
 * computations.
 */
void long_jump_xoshiro256_plusplus(void)
{
    static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf,
                                          0xc5004e441c522fb3,
                                          0x77710069854ee241,
                                          0x39109bb02acbe635 };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for(i = 0; i < sizeof(LONG_JUMP) / sizeof(*LONG_JUMP); i++)
    {
        for(b = 0; b < 64; b++)
        {
            if (LONG_JUMP[i] & UINT64_C(1) << b)
            {
                s0 ^= states_pp[0];
                s1 ^= states_pp[1];
                s2 ^= states_pp[2];
                s3 ^= states_pp[3];
            }
            xoshiro256_plusplus();
        }
    }
        
    states_pp[0] = s0;
    states_pp[1] = s1;
    states_pp[2] = s2;
    states_pp[3] = s3;
}


/* This is the jump function for the generator. It is equivalent
 * to 2^128 calls to xoshiro256_plus(); it can be used to generate 2^128
 * non-overlapping subsequences for parallel computations.
 */
void jump_xoshiro256_plus(void)
{
    static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba,
                                     0xd5a61266f0c9392c,
                                     0xa9582618e03fc9aa,
                                     0x39abdc4529b1661c };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for(i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
    {
        for(b = 0; b < 64; b++)
        {
            if (JUMP[i] & UINT64_C(1) << b)
            {
                s0 ^= states_p[0];
                s1 ^= states_p[1];
                s2 ^= states_p[2];
                s3 ^= states_p[3];
            }
            xoshiro256_plus();
        }
    }
        
    states_p[0] = s0;
    states_p[1] = s1;
    states_p[2] = s2;
    states_p[3] = s3;
}


/* This is the long-jump function for the generator. It is equivalent to
 * 2^192 calls to xoshiro256_plus(); it can be used to generate 2^64 starting
 * points, from each of which jump_xoshiro256_plus() will generate 2^64
 * non-overlapping subsequences for parallel distributed computations.
 */
void long_jump_xoshiro256_plus(void)
{
    static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf,
                                          0xc5004e441c522fb3,
                                          0x77710069854ee241,
                                          0x39109bb02acbe635 };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for(i = 0; i < sizeof(LONG_JUMP) / sizeof(*LONG_JUMP); i++)
    {
        for(b = 0; b < 64; b++)
        {
            if (LONG_JUMP[i] & UINT64_C(1) << b)
            {
                s0 ^= states_p[0];
                s1 ^= states_p[1];
                s2 ^= states_p[2];
                s3 ^= states_p[3];
            }
            xoshiro256_plus();
        }
    }
        
    states_p[0] = s0;
    states_p[1] = s1;
    states_p[2] = s2;
    states_p[3] = s3;
}


/* From http://xoshiro.di.unimi.it/:
 *
 * A standard double (64-bit) floating-point number in IEEE floating point
 * format has 52 bits of significand, plus an implicit bit at the left of the
 * significand. Thus, the representation can actually store numbers with 53
 * significant binary digits.
 *
 * This conversion guarantees that all dyadic rationals of the form k / 2^-53
 * will be equally likely. Note that this conversion prefers the high bits of
 * x, but you can alternatively use the lowest bits.
 */
double uint64_to_double(uint64_t value)
{
    return ((value >> 11) * 0x1.0p-53);
}


/* return a double from [0,1) */
double dxoshiro256ss(void)
{
    return uint64_to_double(xoshiro256_starstar());
}

/* return a double from [0,1) */
double dxoshiro256pp(void)
{
    return uint64_to_double(xoshiro256_plusplus());
}

/* return a double from [0,1) */
double dxoshiro256p(void)
{
    return uint64_to_double(xoshiro256_plus());
}


/* initialize states for random number generators, call enough times
 * to escape from zero and reach an average of 50% 1's in the states
 */
void initialize_xoshiro256(uint64_t initial_seed)
{
    int i;
    
    /* initialize splitmix64 state */
    state_sm64 = initial_seed;

    /* Call splitmix64 enough times to escape from zero.
     *
     * I do not know how rapidly splitmix64 escapes from zero.
     * Assume it escapes at least as rapidly as xoshiro256+.
     */
    for (i = 0; i < 100; i++)
        splitmix64();
    
    /* initialize xorshiro256 states */
    states_ss[0] = splitmix64();
    states_ss[1] = splitmix64();
    states_ss[2] = splitmix64();
    states_ss[3] = splitmix64();
    states_pp[0] = splitmix64();
    states_pp[1] = splitmix64();
    states_pp[2] = splitmix64();
    states_pp[3] = splitmix64();
    states_p[0]  = splitmix64();
    states_p[1]  = splitmix64();
    states_p[2]  = splitmix64();
    states_p[3]  = splitmix64();

    /* Call xorshiro256 enough times to escape from zero */
    for (i = 0; i < 100; i++)
        xoshiro256_starstar();
    for (i = 0; i < 100; i++)
        xoshiro256_plusplus();
    for (i = 0; i < 100; i++)
        xoshiro256_plus();
}


/* Wichura 1998 Algorithm AS 241
 *
 * faster than Guirguis 1991
 *
 * 14 mul, 1 dev, (1 sqrt, 1 log) or (2 more mul)
 */
double invcdf_norm_wichura(double p)
{
    double Q, R, ppnd16;

    Q = p - 0.5;

    /* p close to 0.5 */
    if (fabs(Q) <= 0.425)
    {
        R = 0.180625 - Q * Q;

        return  (Q * (((((((2.5090809287301226727E3  * R +
                            3.3430575583588128105E4) * R +
                            6.7265770927008700853E4) * R +
                            4.5921953931549871457E4) * R +
                            1.3731693765509461125E4) * R +
                            1.9715909503065514427E3) * R +
                            1.3314166789178437745E2) * R +
                            3.3871328727963666080E0) /
                     (((((((5.2264952788528545610E3  * R +
                            2.8729085735721942674E4) * R +
                            3.9307895800092710610E4) * R +
                            2.1213794301586595867E4) * R +
                            5.3941960214247511077E3) * R +
                            6.8718700749205790830E2) * R +
                            4.2313330701600911252E1) * R + 1.0));
    }
    else
    {
        if (Q < 0.0)
            R = p;
        else
            R = 1.0 - p;

        if (R <= 0.0)
            return (0.0);
        
        R = sqrt(-log(R));

        /* p neither close to 0.5 nor 0 or 1 */
        if (R <= 5.0)
        {
            R = R - 1.6;

            ppnd16 = (((((((7.74545014278341407640E-4  * R +
                            2.27238449892691845833E-2) * R +
                            2.41780725177450611770E-1) * R +
                            1.27045825245236838258E0)  * R +
                            3.64784832476320460504E0)  * R +
                            5.76949722146069140550E0)  * R +
                            4.63033784615654529590E0)  * R +
                            1.42343711074968357734E0) /
                     (((((((1.05075007164441684324E-9  * R +
                            5.47593808499534494600E-4) * R +
                            1.51986665636164571966E-2) * R +
                            1.48103976427480074590E-1) * R +
                            6.89767334985100004550E-1) * R +
                            1.67638483018380384940E0)  * R +
                            2.05319162663775882187E0)  * R + 1.0);
        }
        /* p near 0 or 1 */
        else
        {
            R = R - 5.0;

            ppnd16 = (((((((2.01033439929228813265E-7  * R +
                            2.71155556874348757815E-5) * R +
                            1.24266094738807843860E-3) * R +
                            2.65321895265761230930E-2) * R +
                            2.96560571828504891230E-1) * R +
                            1.78482653991729133580E0)  * R +
                            5.46378491116411436990E0)  * R +
                            6.65790464350110377720E0) /
                     (((((((2.04426310338993978564E-15 * R +
                            1.42151175831644588870E-7) * R +
                            1.84631831751005468180E-5) * R +
                            7.86869131145613259100E-4) * R +
                            1.48753612908506148525E-2) * R +
                            1.36929880922735805310E-1) * R +
                            5.99832206555887937690E-1) * R + 1.0);
        }
        
        if (Q < 0.0)
            ppnd16 = -ppnd16;
        
        return (ppnd16);
    }
}


/* Use inverse CDF to transform uniform into normal distribution
 *
 * Faster than polar normal if Wichura 1998 AS 241 is used
 */
inline double drand_normal256p_icdf(double mean, double sd)
{
    /* Wichura already handles values > 0.5 */
    return (mean + sd * invcdf_norm_wichura(dxoshiro256p()));
}
