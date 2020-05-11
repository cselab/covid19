#include "utils.h"

#include <cstdarg>
#include <cstdio>
#include <cstdlib>

void _die [[gnu::format(printf, 3, 4)]] (const char *filename, int line, const char *fmt, ...) {
    fprintf(stderr, "%s:%d ", filename, line);

    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    exit(1);
}
