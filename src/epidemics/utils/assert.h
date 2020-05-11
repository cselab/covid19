#pragma once

void _die [[gnu::format(printf, 3, 4)]] (const char *filename, int line, const char *fmt, ...);

// TODO: assert

#define DIE(...) _die(__FILE__, __LINE__, __VA_ARGS__)
