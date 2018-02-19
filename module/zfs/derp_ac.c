#include <sys/derp_ac.h>
#include <sys/spa_impl.h>
#include <sys/vdev_impl.h>
//#include <sys/zio_impl.h>
//#include <sys/spa.h>
//#include <sys/zfs_context.h>
//#include <sys/fm/fs/zfs.h>
//#include <sys/spa.h>
//#include <sys/spa_impl.h>
//#include <sys/dmu.h>
//#include <sys/dmu_tx.h>
//#include <sys/uberblock_impl.h>
//#include <sys/metaslab.h>
//#include <sys/metaslab_impl.h>
//#include <sys/space_map.h>
//#include <sys/space_reftree.h>
//#include <sys/zio.h>
//#include <sys/zap.h>
//#include <sys/fs/zfs.h>
//#include <sys/arc.h>
//#include <sys/zil.h>
//#include <sys/dsl_scan.h>
//#include <sys/abd.h>
//#include <sys/zvol.h>
//#include <sys/zfs_ratelimit.h>


enum zio_compress compress_funcs[DERP_AC_NUM_OF_COMP_ALGS] = {
		ZIO_COMPRESS_EMPTY,
		ZIO_COMPRESS_LZ4,
		ZIO_COMPRESS_GZIP_1
};


// Taken from: https://en.wikipedia.org/wiki/Linear-feedback_shift_register.
uint16_t lfsr = 0xACE1u;
uint16_t rand_int16(void) {
	uint16_t bit = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 5) ) & 1;
	return lfsr  = (lfsr >> 1) | (bit << 15);
}

uint64_t zio_compress_to_derp_compress(enum zio_compress compress) {
	switch (compress) {
		case ZIO_COMPRESS_INHERIT:
			return -1;
		case ZIO_COMPRESS_ON:
			return -1;
		case ZIO_COMPRESS_OFF:
			return -1;
		case ZIO_COMPRESS_LZJB:
			return -1;
		case ZIO_COMPRESS_EMPTY:
			return 0;
		case ZIO_COMPRESS_GZIP_1:
			return 2;
		case ZIO_COMPRESS_GZIP_2:
			return -1;
		case ZIO_COMPRESS_GZIP_3:
			return -1;
		case ZIO_COMPRESS_GZIP_4:
			return -1;
		case ZIO_COMPRESS_GZIP_5:
			return -1;
		case ZIO_COMPRESS_GZIP_6:
			return -1;
		case ZIO_COMPRESS_GZIP_7:
			return -1;
		case ZIO_COMPRESS_GZIP_8:
			return -1;
		case ZIO_COMPRESS_GZIP_9:
			return -1;
		case ZIO_COMPRESS_ZLE:
			return -1;
		case ZIO_COMPRESS_LZ4:
			return 1;
		case ZIO_COMPRESS_DERP_AC:
			return -1;
		case ZIO_COMPRESS_FUNCTIONS:
			return -1;
	}
	return -1;
}

void derp_add_to_rolling_ave(uint64_t new_value, uint64_t *ave, int n) {
    if (new_value) {  // if the new_value isn't a 0.
        if (*ave) {  // if the average isn't a 0
            uint64_t temp = *ave;
            *ave = (new_value + (temp * (n-1))) / n;
        } else {
            // This is the first value.
            *ave = new_value;
        }
    }
}

uint64_t derp_calc_bps(uint64_t bytes, hrtime_t time) {
	return (bytes*BYTES_FACTOR) / time;
}

uint64_t derp_calc_time_from_bps(uint64_t bps, uint64_t bytes) {
	return (bytes*BYTES_FACTOR) / bps;
}

uint64_t derp_calc_compress_ratio(uint64_t start_bytes, uint64_t end_bytes) {
	return (start_bytes*BYTES_FACTOR) / end_bytes;
}

uint64_t derp_calc_compressed_size(uint64_t start_bytes, uint64_t compress_ratio) {
	return (start_bytes*BYTES_FACTOR) / compress_ratio;
}

void derp_set_default_compress(enum zio_compress *compress) {
	*compress = ZIO_COMPRESS_EMPTY;
}


void derp_update_vdev_stats(vdev_t *vd, uint64_t bc_bucket, uint64_t compress, uint64_t bps, uint64_t ratio)
{
//	uint64_t min_delay = 0;

	if (!vd->vdev_children) { // is leaf
		dprintf("hi vd:%u updating element[%u][%u][%u] with value %u", vd->vdev_id, bc_bucket, compress, COMPRESS_RATE_INDEX, bps);
		derp_add_to_rolling_ave(bps,
				&vd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][compress][COMPRESS_RATE_INDEX],
				1000);
		dprintf("hi updated to %u", vd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][compress][COMPRESS_RATE_INDEX]);
		derp_add_to_rolling_ave(ratio,
				&vd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][compress][COMPRESS_RATIO_INDEX],
				1000);
//		compress_vdev_queue_delay(size, vd);
	} else {
		int i;
		for (i = 0; i < vd->vdev_children; i++) {
			derp_update_vdev_stats(vd->vdev_child[i], bc_bucket, compress, bps, ratio);
//			if (vdev_delay) {
//				if (min_delay == 0) {
//					min_delay = vdev_delay;
//				} else if (vdev_delay < min_delay) {
//					min_delay = vdev_delay;
//				}
//			}
		}
	}
//	return (min_delay);
}

void derp_update_model(vdev_t *rvd, enum zio_compress *compress, hrtime_t compress_delay, uint64_t lsize, uint64_t psize) {
//void derp_update_model(zio_t *pio, enum zio_compress *compress, hrtime_t compress_delay, uint64_t lsize, uint64_t psize) {
	if (compress_delay == 0 || lsize == 0 || psize == 0) return;
	uint64_t derp_compress = zio_compress_to_derp_compress(*compress);
	uint64_t bps = derp_calc_bps(lsize, compress_delay);
	uint64_t ratio = derp_calc_compress_ratio(lsize, psize);
	derp_update_vdev_stats(rvd, 0, derp_compress, bps, ratio);
}

size_t derp_ac_compress(zio_t *zio, void *dst, size_t s_len, enum zio_compress *compress) {
	/* Params:
	 * 	zio:
	 * 	compress:	This is needed so we can change the compression algorithm used which is  set by the calling function.
	 * */
	// If we want to die early and do no compression simply return 0.

	// Get the io device which contains the ac info.
//	zio_t* pio = zio_unique_parent(zio);

//	spa_t *spa = zio->io_spa;
//	vdev_t *rvd = spa->spa_root_vdev;
//	vdev_t *vd = zio->io_vd ? zio->io_vd : rvd;

////	vdev_t *pvd;
////	uint64_t txg = zio->io_txg;
////	vdev_stat_t *vs = &vd->vdev_stat;
//	vdev_stat_ex_t *vsx = &vd->vdev_stat_ex;
////	zio_type_t type = zio->io_type;

	// Decide compress here
//	uint64_t bc_bucket = 0;

	size_t psize;
	vdev_t *rvd = zio->io_spa->spa_root_vdev;
	zio_t *pio = zio_unique_parent(zio);

	derp_set_default_compress(compress);

	if (pio == NULL) {
		psize = zio_compress_data(*compress, zio->io_abd, dst, s_len);
	} else {
		uint16_t i = rand_int16() % DERP_AC_NUM_OF_COMP_ALGS;
		*compress = compress_funcs[i];
		hrtime_t comp_start = gethrtime();
		psize = zio_compress_data(*compress, zio->io_abd, dst, s_len);
		hrtime_t comp_end = gethrtime();
//		rvd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][zio_compress_to_derp_compress(*compress)][0] = comp_end-comp_start;
//		dprintf("Hi there: 1: %u, 2: %u", rvd->vdev_stat_ex.vsx_derp_ac_model[0][0][0], rvd->vdev_stat_ex.vsx_derp_ac_model[1][0][0]);

//		derp_update_model(pio, compress, comp_end-comp_start, zio->io_lsize, psize);
		derp_update_model(rvd, compress, comp_end-comp_start, zio->io_lsize, psize);
	}

//	*compress = ZIO_COMPRESS_LZ4;

//	hrtime_t comp_start = gethrtime();
//	hrtime_t comp_end = gethrtime();

//	vsx->vsx_derp_ac_model[bc_bucket][(uint64_t)compress][COMPRESS_RATE_INDEX] = comp_end-comp_start;

//	uint64_t compression_bps = derp_calc_bps(s_len, comp_end-comp_start);
//	uint64_t compression_ratio = derp_calc_compress_ratio(s_len, psize);
//
//	derp_add_to_rolling_ave(compression_bps,
//			&vsx->vsx_derp_ac_model[bc_bucket][(uint64_t)compress][COMPRESS_RATE_INDEX],
//			1000);
//	derp_add_to_rolling_ave(compression_ratio,
//			&vsx->vsx_derp_ac_model[bc_bucket][(uint64_t)compress][COMPRESS_RATIO_INDEX],
//			1000);

	return psize;

//	return 0;
}
