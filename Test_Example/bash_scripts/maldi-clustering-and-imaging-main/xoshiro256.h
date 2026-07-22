extern uint64_t splitmix64(void);    /* period 2^64, radically different method */
extern uint64_t xoshiro256_starstar(void);
extern uint64_t xoshiro256_plusplus(void);
extern uint64_t xoshiro256_plus(void);
extern void     jump_xoshiro256_starstar(void);
extern void     long_jump_xoshiro256_starstar(void);
extern void     jump_xoshiro256_plusplus(void);
extern void     long_jump_xoshiro256_plusplus(void);
extern void     jump_xoshiro256_plus(void);
extern void     long_jump_xoshiro256_plus(void);
extern double   uint64_to_double(uint64_t value);
extern double   dxoshiro256ss(void);
extern double   dxoshiro256pp(void);
extern double   dxoshiro256p(void);
extern void     initialize_xoshiro256(uint64_t initial_seed);
extern double   drand_normal256p_icdf(double mean, double sd);
