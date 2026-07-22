
# holy crap, GCC >= v4.8 diagonstics are *WAY* too verbose !!

# -ffinite-math-only slows some benchmarks down for some reason
gcc -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0 -fno-stack-protector -static -pthread -Wall -Wno-unused-but-set-variable \
  -fno-diagnostics-show-caret -fno-diagnostics-show-option \
  -mfpmath=sse -mavx2 \
  -O3 -funroll-loops -fomit-frame-pointer -momit-leaf-frame-pointer \
  -fno-math-errno -fno-signed-zeros -fno-signaling-nans \
  -fno-trapping-math -mdaz-ftz \
  hcdist.c tree.c text.c rand_xoshiro256.c -o hcdist -lm \
  |& grep -v "In function "

# debug
#gcc -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0 -fno-stack-protector -pthread -Wall -Wno-unused-but-set-variable -fno-diagnostics-show-caret -fno-diagnostics-show-option -mfpmath=sse -mavx2 -O2 -g -funroll-loops -fomit-frame-pointer -momit-leaf-frame-pointer -fno-math-errno -fno-signed-zeros -fno-trapping-math -fno-signaling-nans hcdist.c tree.c text.c rand_xoshiro256.c -o hcdist_test -lm |& grep -v "In function "

# for profiling with gprof
#gcc -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0 -fno-stack-protector -pthread -Wall -Wno-unused-but-set-variable -fno-diagnostics-show-caret -fno-diagnostics-show-option -mfpmath=sse -mavx2 -O2 -g -pg -funroll-loops -fomit-frame-pointer -momit-leaf-frame-pointer -fno-math-errno -fno-signed-zeros -fno-trapping-math -fno-signaling-nans hcdist.c tree.c text.c rand_xoshiro256.c -o hcdist_test -lm |& grep -v "In function "


# --no-merge-notes shuts up nuisance warnings on glibc, but is good bit larger
# so, we'll need to take a sledgehammer to *ALL* warnings instead
strip hcdist 2> /dev/null
