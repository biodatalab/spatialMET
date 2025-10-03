
# -ffinite-math-only slows some benchmarks down for some reason

# modify -I and -L as appropriate to point to your FreeImage installation
gcc -static -pthread -Wall -Wno-unused-but-set-variable \
  -fno-diagnostics-show-caret -fno-diagnostics-show-option \
  -mfpmath=sse -mavx2 -ftree-vectorize \
  -O3 -g -funroll-loops -fno-omit-frame-pointer \
  -I/share/data2/welshea/FreeImage/Dist -L/share/data2/welshea/FreeImage/Dist \
  maldi_image_from_clusters.c tree.c text.c rand_xoshiro256.c cet_colors.c \
  -o maldi_image_from_clusters -lfreeimage -lstdc++ -lm \
  |& grep -v "In function "

# --no-merge-notes shuts up nuisance warnings on glibc, but is good bit larger
# so, we'll need to take a sledgehammer to *ALL* warnings instead
strip maldi_image_from_clusters 2> /dev/null
