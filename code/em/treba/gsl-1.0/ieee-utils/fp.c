#include <config.h>

#if defined(HAVE_SPARCLINUX_IEEE_INTERFACE)
#include "fp-sparclinux.c"
#elif defined(HAVE_M68KLINUX_IEEE_INTERFACE)
#include "fp-m68klinux.c"
#elif defined(HAVE_PPCLINUX_IEEE_INTERFACE)
#include "fp-ppclinux.c"
#elif defined(HAVE_X86LINUX_IEEE_INTERFACE)
#include "fp-x86linux.c"
#elif defined(HAVE_HPUX11_IEEE_INTERFACE)
#include "fp-hpux11.c"
#elif defined(HAVE_HPUX_IEEE_INTERFACE)
#include "fp-hpux.c"
#elif defined(HAVE_SUNOS4_IEEE_INTERFACE)
#include "fp-sunos4.c"
#elif defined(HAVE_SOLARIS_IEEE_INTERFACE)
#include "fp-solaris.c"
#elif defined(HAVE_IRIX_IEEE_INTERFACE)
#include "fp-irix.c"
#elif defined(HAVE_AIX_IEEE_INTERFACE)
#include "fp-aix.c"
#elif defined(HAVE_TRU64_IEEE_INTERFACE)
#include "fp-tru64.c"
#elif defined(HAVE_FREEBSD_IEEE_INTERFACE)
#include "fp-freebsd.c"
#elif defined(HAVE_OS2EMX_IEEE_INTERFACE)
#include "fp-os2emx.c"
#elif defined(HAVE_NETBSD_IEEE_INTERFACE)
#include "fp-netbsd.c"
#elif defined(HAVE_OPENBSD_IEEE_INTERFACE)
#include "fp-openbsd.c"
#elif defined(HAVE_DARWIN_IEEE_INTERFACE)
#include "fp-darwin.c"
#else
#include "fp-unknown.c" 
#endif



