

make: 
	# clang test.c /Users/ferxinii/Desktop/TFM/predicates/build/Bin/libpredicates.dylib -o test -Wpedantic -Wextra -Wall -Werror -Wl,-rpath,/Users/ferxinii/Desktop/TFM/predicates/build/Bin/ -O3
	# clang test_lung.c /Users/ferxinii/Desktop/TFM/predicates/build/Bin/libpredicates.dylib -o test_lung -Wpedantic -Wextra -Wall -Werror -Wl,-rpath,/Users/ferxinii/Desktop/TFM/predicates/build/Bin/ -g -O3
	# clang test_dt.c /Users/ferxinii/Desktop/TFM/predicates/build/Bin/libpredicates.dylib -o test_dt -Wpedantic -Wextra -Wall -Werror -Wl,-rpath,/Users/ferxinii/Desktop/TFM/predicates/build/Bin/ -g -O0
	clang test_voronoi.c voronoi.c vd_3d.c bpoly.c algebra.c geometry.c dt_3d_incremental.c simplical_complex.c convhull_3d.c /Users/ferxinii/Desktop/TFM/predicates/build/Bin/libpredicates.dylib -o test_voronoi -Wpedantic -Wextra -Wall -Werror -Wl,-rpath,/Users/ferxinii/Desktop/TFM/predicates/build/Bin/ -O3 -g
