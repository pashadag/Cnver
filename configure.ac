AC_INIT([cnver],[0.8.1],[cnver@cs.toronto.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
CFLAGS='-std=c99'
AC_PROG_CC
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CHECK_LIB([z], [gzread], ,
			     [AC_MSG_ERROR([zlib not found, see http://www.zlib.net])]) 
AC_CONFIG_FILES([
           Makefile
	   src/Makefile
	   src/cluster/Makefile
	   src/high-copy-count-regions/Makefile
	   src/post_analysis/Makefile
	   src/cs2-4.6/Makefile
         ])
AC_OUTPUT
