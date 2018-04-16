#ifndef _SYS_DERP_BYTECOUNTER_H
#define  _SYS_DERP_BYTECOUNTER_H

#include <sys/fs/zfs.h>

uint64_t derp_bytecounter(void *buf, size_t len);
uint64_t derp_to_Q_type(uint64_t bytecount);

#endif /* _SYS_DERP_BYTECOUNTER_H */
