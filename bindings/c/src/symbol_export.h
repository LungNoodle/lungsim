
#ifndef AETHER_SYMBOL_EXPORT_H
#define AETHER_SYMBOL_EXPORT_H

#if defined _WIN32 || defined __CYGWIN__
  #ifdef aether_c_EXPORTS
    #ifdef __GNUC__
      #define SHO_PUBLIC __attribute__ ((dllexport))
    #else
      #define SHO_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #else
    #ifdef __GNUC__
      #define SHO_PUBLIC __attribute__ ((dllimport))
    #else
      #define SHO_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #endif
  #define SHO_LOCAL
#else
  #if __GNUC__ >= 4
    #define SHO_PUBLIC __attribute__ ((visibility ("default")))
    #define SHO_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define SHO_PUBLIC
    #define SHO_LOCAL
  #endif
#endif

#endif /* AETHER_SYMBOL_EXPORT_H */

