#ifndef _SYS_DERP_AC_H
#define  _SYS_DERP_AC_H

#include <sys/spa.h>
#include <sys/zio.h>
#include <sys/fs/zfs.h>
#include <sys/time.h>


// You can find the variables for bucket size and etc. in zfs.h around line 886.


size_t derp_ac_compress(zio_t *zio, void *dst, size_t s_len, enum zio_compress *compress);
void derp_add_to_rolling_ave(uint64_t new_value, uint64_t *ave, int n);

#define DERP_COMPRESS_RATE_INDEX 0
#define DERP_COMPRESS_RATIO_INDEX 1

// These function don't really mean bytes per second for bps.
#define BYTES_FACTOR (1000*1000*1000)
#define MAX_RATE ((uint64_t)(1000*1000*1000)*10)
#define SECOND (1000*1000*1000)
//*1000
uint64_t derp_calc_bps(uint64_t bytes, hrtime_t time);
uint64_t derp_calc_time_from_bps(uint64_t bps, uint64_t bytes);

// TODO: make separate file containing all cpu monitoring stuff.
//typedef struct derp_cpuload {
//    long unsigned 	last_total_time;
//    long unsigned 	last_idle_time;
//    long unsigned 	cur_total_time;
//    long unsigned 	cur_idle_time;
//    unsigned		idle_time_diff;   // The current idle time diff calculated.
//    unsigned		total_time_diff;
//} d_cpuload;


#endif /* _SYS_DERP_AC_H */
