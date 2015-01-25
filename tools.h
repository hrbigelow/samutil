#ifndef _TOOLS_H
#define _TOOLS_H

FILE *open_if_present(const char *file, const char *mode);

int close_if_present(FILE *fh);

#endif /* _TOOLS_H */
