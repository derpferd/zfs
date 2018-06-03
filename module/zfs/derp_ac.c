#include <sys/derp_ac.h>
#include <sys/spa_impl.h>
#include <sys/vdev_impl.h>
#include <sys/derp_bytecounter.h>
//#include <linux/blkdev_compat.h>


typedef uint64_t derp_compress;

#define derp_for_each_compress(derp_c) \
    for(derp_c = 0; derp_c < DERP_AC_NUM_OF_COMP_ALGS; derp_c++)

enum zio_compress derp_compress_funcs[DERP_AC_NUM_OF_COMP_ALGS] = {
		ZIO_COMPRESS_EMPTY,
		ZIO_COMPRESS_LZ4,
		ZIO_COMPRESS_GZIP_1
};

size_t derp_compress_data(enum zio_compress c, void *src, void *dst, size_t s_len);
static int derp_zeroed_cb(void *data, size_t len, void *private);

typedef struct bps_variables {
	uint64_t bc_bucket;
	derp_compress compress;
//	TODO: add other variables like cpu usage.
} derp_bps_variables_t;

typedef struct ratio_variables {
	uint64_t bc_bucket;
	derp_compress compress;
} derp_ratio_variables_t;

// Taken from: https://en.wikipedia.org/wiki/Linear-feedback_shift_register.
uint16_t lfsr = 0xACE1u;
uint16_t rand_int16(void) {
	uint16_t bit = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 5) ) & 1;
	return lfsr  = (lfsr >> 1) | (bit << 15);
}

derp_compress zio_compress_to_derp_compress(enum zio_compress compress) {
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
		case ZIO_COMPRESS_DERP_AC_TRAIN:
			return -1;
		case ZIO_COMPRESS_AUTO:
			return -1;
		case ZIO_COMPRESS_FUNCTIONS:
			return -1;
	}
	return -1;
}


int derp_load_cur_stats(d_cpuload* stats) {
//	int i;
//	stats->last_idle_time = 80;
//	stats->last_total_time = 100;
//	stats->cur_idle_time = 100;
//	stats->cur_total_time = 0;

#ifdef _SPL_TIME_H
	int i;

//	stats->last_idle_time = 0;
//	stats->last_total_time = 0;
	stats->cur_idle_time = 0;
	stats->cur_total_time = 0;

//	stats->cur_idle_time = kcpustat_cpu(0).cpustat[CPUTIME_USER];
//	stats->cur_total_time = get_idle_time(0);

	for_each_possible_cpu(i) {
		for (int j=0; j < NR_STATS; j++) {
			stats->cur_total_time += kcpustat_cpu(0).cpustat[j];
			if (j == CPUTIME_IDLE || j == CPUTIME_IOWAIT) {
				stats->cur_idle_time += kcpustat_cpu(0).cpustat[j];
			}
		}
//		stats->cur_idle_time = kcpustat_cpu(0).cpustat[CPUTIME_USER];
//		stats->cur_total_time++;
	}
//	stats->cur_idle_time = CPUTIME_USER;
//	get_cpu_idle_time_us(0, &stats->cur_total_time);
#else
	stats->cur_idle_time = 0;
	stats->cur_total_time = 0;
#endif

    return 0;
}

int derp_update_cpuload(d_cpuload* stats) {
	dprintf("Updating cpuload");
    int err;
    // Get cur times from /proc/stat
    err = derp_load_cur_stats(stats);
    if (err) return err;

    // Calc diff
    if (stats->last_total_time > stats->cur_total_time) {  // Time has moved on
    	if (stats->last_idle_time >= stats->cur_idle_time ) {  // The idle time hasn't gone down.
    		stats->idle_time_diff = stats->cur_idle_time - stats->last_idle_time;
    	} else {
    		dprintf("WARNING: This is very bad. The total idle time has gone down!!!");
    	}
    	stats->total_time_diff = stats->cur_total_time - stats->last_total_time;

        // Set last to cur
        stats->last_idle_time = stats->cur_idle_time;
        stats->last_total_time = stats->cur_total_time;
    } else if (stats->last_total_time == stats->cur_total_time) {
    	dprintf("Not updating since nothing changed.");
    } else {
		dprintf("WARNING: This is very bad. We appear to have gone back in time!!!");
    }

    if (stats->total_time_diff != 0) {
    	stats->idle_percent = (stats->idle_time_diff * 100) / stats->total_time_diff;
    } else {
    	stats->idle_percent = 0;
    }

    return 0;
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
	return (end_bytes*BYTES_FACTOR) / start_bytes;
}

uint64_t derp_calc_compressed_size(uint64_t start_bytes, uint64_t compress_ratio) {
	return (start_bytes*BYTES_FACTOR) / compress_ratio;
}

uint64_t derp_apply_compress_ratio(uint64_t bytes, uint64_t ratio) {
	// Hack to overcome intater overflow
	return ((bytes/10000)*(ratio/100000));
//	return (bytes*ratio)/BYTES_FACTOR;
}

uint64_t* derp_get_bps(vdev_stat_ex_t *stats, derp_bps_variables_t *vars) {
	dprintf("derp: getting element[%u][%u][%u]", vars->bc_bucket, vars->compress, DERP_COMPRESS_RATE_INDEX);
	return &stats->vsx_derp_ac_model[vars->bc_bucket][vars->compress][DERP_COMPRESS_RATE_INDEX];
}

uint64_t* derp_get_ratio(vdev_stat_ex_t *stats, derp_ratio_variables_t *vars) {
	dprintf("derp: getting element[%u][%u][%u]", vars->bc_bucket, vars->compress, DERP_COMPRESS_RATIO_INDEX);
	return &stats->vsx_derp_ac_model[vars->bc_bucket][vars->compress][DERP_COMPRESS_RATIO_INDEX];
}

void derp_set_default_compress(enum zio_compress *compress) {
	*compress = ZIO_COMPRESS_EMPTY;
}

vdev_stat_ex_t* derp_get_vdev_ex_stats(vdev_t *vd) {
	if (!vd->vdev_children) {
		return &vd->vdev_stat_ex;
	}
//	else {
	// Just use the first child;
	return derp_get_vdev_ex_stats(vd->vdev_child[0]);
//		for (int i = 0; i < vd->vdev_children; i++) {
//
//		}
//	}
}

vdev_stat_t* derp_get_vdev_stats(vdev_t *vd) {
	if (!vd->vdev_children) {
		return &vd->vdev_stat;
	}
	return derp_get_vdev_stats(vd->vdev_child[0]);
}

uint64_t get_queue_size(vdev_t *vd) {
	if (!vd->vdev_children) {
		return vd->vdev_queue.vq_class[ZIO_PRIORITY_ASYNC_WRITE].vqc_queued_size;
	}
	return get_queue_size(vd->vdev_child[0]);
}

void derp_update_vdev_stats(vdev_t *vd, uint64_t bc_bucket, uint64_t compress, uint64_t bps, uint64_t ratio) {
//	uint64_t min_delay = 0;

	if (!vd->vdev_children) { // is leaf
		dprintf("derp: vd:%u updating element[%u][%u][%u] with value %u", vd->vdev_id, bc_bucket, compress, DERP_COMPRESS_RATE_INDEX, bps);
		derp_add_to_rolling_ave(bps,
				&vd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][compress][DERP_COMPRESS_RATE_INDEX],
				1000);
		dprintf("derp: updated to %u", vd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][compress][DERP_COMPRESS_RATE_INDEX]);
		dprintf("derp: vd:%u updating element[%u][%u][%u] with value %u", vd->vdev_id, bc_bucket, compress, DERP_COMPRESS_RATIO_INDEX, bps);
		derp_add_to_rolling_ave(ratio,
				&vd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][compress][DERP_COMPRESS_RATIO_INDEX],
				1000);
		dprintf("derp: updated to %u", vd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][compress][DERP_COMPRESS_RATIO_INDEX]);
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

void derp_update_model(vdev_t *rvd, enum zio_compress *compress, hrtime_t compress_delay, uint64_t lsize, uint64_t psize, uint64_t bc_bucket) {
//void derp_update_model(zio_t *pio, enum zio_compress *compress, hrtime_t compress_delay, uint64_t lsize, uint64_t psize) {
	if (compress_delay == 0 || lsize == 0 || psize == 0) return;
	uint64_t derp_compress = zio_compress_to_derp_compress(*compress);
	uint64_t bps = derp_calc_bps(lsize, compress_delay);
	uint64_t ratio = derp_calc_compress_ratio(lsize, psize);
	derp_update_vdev_stats(rvd, bc_bucket, derp_compress, bps, ratio);
}

enum zio_compress derp_rand_compress(void) {
	uint16_t i = rand_int16() % DERP_AC_NUM_OF_COMP_ALGS;
	return derp_compress_funcs[i];
//	return ZIO_COMPRESS_EMPTY;
}

enum zio_compress derp_get_best_compress(vdev_stat_ex_t *stats, enum zio_type type, uint64_t bc_bucket) {
	/*
	 * How to get the best compression alg.
	 * Find the compression alg. that maximizes Total Write Rate (TWR).
	 *  TWR = CR(Compression Rate) + DR(Current Disk Rate)*CRatio(Compression Ratio)
	 *  TWT = CR(Compression time) + (DR(Current Disk Rate)*CRatio(Compression Ratio)
	 *
	 */
	dprintf("Yo Dog 1");
	enum zio_compress c;
	derp_compress derp_c;
	derp_bps_variables_t bps_vars;
	derp_ratio_variables_t ratio_vars;
//	uint64_t cur_bps = stats->vsx_derp_disk_bps[type];
//	uint64_t cur_total_bps = stats->vsx_derp_total_bps[type];
//	cur_total_bps = stats->max_bps;
	uint64_t max_disk_bps = stats->max_bps;

	bps_vars.bc_bucket = bc_bucket;
	ratio_vars.bc_bucket = bc_bucket;

	uint64_t max_twr = 0;
	enum zio_compress best_c = ZIO_COMPRESS_EMPTY; // Default to no compression
	dprintf("Picking best...");
	dprintf("\tMax Disk BPS: %lu \tBucket: %u", max_disk_bps, bc_bucket);
	derp_for_each_compress(derp_c) {
		c = derp_compress_funcs[derp_c];
		bps_vars.compress = derp_c;
		ratio_vars.compress = derp_c;

		uint64_t c_rate = *derp_get_bps(stats, &bps_vars);
		uint64_t c_ratio = *derp_get_ratio(stats, &ratio_vars);

		uint64_t twr;
//		uint64_t MAX_RATE = (uint64_t)BYTES_FACTOR * 10;
//		MAX_RATE = MAX_RATE *10;
		if (c_rate > MAX_RATE) {
			twr = derp_apply_compress_ratio(max_disk_bps, c_ratio);
		} else {
			twr = (max_disk_bps*c_rate) / (max_disk_bps + (derp_apply_compress_ratio(c_rate, c_ratio)));
		}
		dprintf("derp:\tc: %s  \trate:%lu   \tratio:%lu  \ttwr:%lu rate: %lu", (&zio_compress_table[c])->ci_name, c_rate, c_ratio, twr, derp_apply_compress_ratio(c_rate, c_ratio));
//		dprintf("derp:　　　　c: %s c_ratio: %u X: %u Y:%u", (&zio_compress_table[c])->ci_name, c_ratio, c_rate*c_ratio, (c_rate*c_ratio)/BYTES_FACTOR);
//		dprintf("derp:　　　　c: %s twr: %u = (%u*%u) / (%u + %u)", (&zio_compress_table[c])->ci_name, twr, cur_bps, c_rate, cur_bps, derp_apply_compress_ratio(c_rate, c_ratio));
//		dprintf("derp:　　　　c: %s twr: %u = (%u) / (%u)", (&zio_compress_table[c])->ci_name, twr, cur_bps*c_rate, cur_bps+derp_apply_compress_ratio(c_rate, c_ratio));
		if (twr == 0) {
			return c;
		}
		if (twr > max_twr) {
			max_twr = twr;
			best_c = c;
		}
	}
	return best_c;
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
//	derp_test();

	size_t psize;
	vdev_t *rvd = zio->io_spa->spa_root_vdev;
	zio_t *pio = zio_unique_parent(zio);
	uint training = 0;
	if (*compress == ZIO_COMPRESS_DERP_AC_TRAIN) {
		training = 1;
	}

	derp_set_default_compress(compress);

	if (pio == NULL) {
		psize = zio_compress_data(*compress, zio->io_abd, dst, s_len);
	} else {

		/*
		 * If the data is all zeroes, we don't even need to allocate
		 * a block for it.  We indicate this by returning zero size.
		 */
		if (abd_iterate_func(zio->io_abd, 0, s_len, derp_zeroed_cb, NULL) == 0)
			return (0);

		// Get stats
		vdev_stat_ex_t* ex_stats = derp_get_vdev_ex_stats(rvd);

//		psize = zio_compress_data(*compress, zio->io_abd, dst, s_len);
		/* No compression algorithms can read from ABDs directly */
		void *tmp = abd_borrow_buf_copy(zio->io_abd, s_len);
		uint64_t bytecount = derp_bytecounter(tmp, s_len);  // TODO: make this read directly from the ABD.
		uint64_t bc_bucket = derp_to_Q_type(bytecount);

		// Commented out the following to debug cpu_stat
//		// Update cpu stats
//		hrtime_t cpu_stat_update_interval = SECOND;
//		if (ex_stats->cpustats.timestamp == 0) {  // Haven't set the stats yet.
//			derp_update_cpuload(&ex_stats->cpustats);
//			ex_stats->cpustats.timestamp = gethrtime();
//		} else if (ex_stats->cpustats.timestamp + cpu_stat_update_interval < gethrtime()) {
//			derp_update_cpuload(&ex_stats->cpustats);
//			ex_stats->cpustats.timestamp = gethrtime();
//		} else {
//			// Skip, it is too soon to update.
//			dprintf("CPU Stat skipping...");
//		}

		// Maybe this is what is breaking it???
		uint64_t vd_queued_size_write = 0;
//		uint64_t vd_queued_size_write = get_queue_size(rvd);

		// TODO: get cpu usage:

		dprintf("derp: bytecount: %u bucket: %u queue: %u", bytecount, bc_bucket, vd_queued_size_write);
//		dprintf("derp: idle time %lu", ex_stats->cpustats.cur_idle_time);
//		dprintf("derp: total time %lu", ex_stats->cpustats.cur_total_time);
//		dprintf("derp: idle %lu%%", ex_stats->cpustats.idle_percent);
		if (training) {
			if (zio->io_type == ZIO_TYPE_WRITE) {
//				vdev_stat_ex_t* ex_stats = kmem_alloc(sizeof(vdev_stat_ex_t), KM_SLEEP);// = derp_get_vdev_ex_stats(rvd);
				vdev_stat_t* vs = kmem_alloc(sizeof(vdev_stat_t), KM_SLEEP);// = derp_get_vdev_stats(rvd);
//				vdev_stat_t test_stats;
				vdev_get_stats(rvd, vs);
				hrtime_t measuring_interval = SECOND;
				//measuring_interval = measuring_interval;
//				if (ex_stats->last_update == 0) {
//					ex_stats->last_update = gethrtime();
//					ex_stats->test = vs->vs_bytes[zio->io_type]
//				} else
				if (ex_stats->old_timestamp == 0) {
					ex_stats->old_timestamp = vs->vs_timestamp;
					ex_stats->old_bytes = vs->vs_bytes[ZIO_TYPE_WRITE];
//					ex_stats->test = vs->vs_bytes[zio->io_type];
					dprintf("derp: Setting up bps counter");
				} else
				if (ex_stats->old_timestamp + measuring_interval < vs->vs_timestamp) {
					uint64_t tmp_bytes = vs->vs_bytes[ZIO_TYPE_WRITE] - ex_stats->old_bytes;
					uint64_t tdelta = vs->vs_timestamp - ex_stats->old_timestamp;
					uint64_t bps = derp_calc_bps(tmp_bytes, tdelta);
//					double scale = (double)SECOND / tdelta;
//					uint64_t tmp_bytes = vs->vs_bytes[zio->io_type];
//					hrtime_t tmp_time = gethrtime();
//					dprintf("derp: temp_max_bps: %lu. %lu bytes in %lu ns", derp_calc_bps(tmp_bytes-ex_stats->test, tmp_time-ex_stats->last_update), tmp_bytes-ex_stats->test, tmp_time-ex_stats->last_update);
					dprintf("derp: bps: %lu. tmp_bytes: %lu in %lu ns", bps, tmp_bytes, tdelta);
//					tmp_bytes = derp_calc_bps(tmp_bytes-ex_stats->test, tmp_time-ex_stats->last_update);
					if (bps > ex_stats->max_bps) {
						ex_stats->max_bps = bps;
					}
					ex_stats->old_timestamp = vs->vs_timestamp;
					ex_stats->old_bytes = vs->vs_bytes[ZIO_TYPE_WRITE];
//					ex_stats->last_update = tmp_time;
//					ex_stats->test = tmp_bytes;
//					dprintf("derp: updated max_bps: %lu", ex_stats->max_bps);
				} else {
//					dprintf("derp: waiting old: %lld new: %lld cur: %lld", ex_stats->old_timestamp, vs->vs_timestamp, gethrtime());
//					dprintf("derp: waiting (%lu + %lu) - %lu more nanosecs", stats->last_update, SECOND, gethrtime());
	//				dprintf("derp: max_bps: %lu", stats->max_bps);
//					dprintf("derp: .");
//					dprintf("derp: test %lld", test_stats.vs_timestamp);
//					dprintf("derp: bytes %lld", test_stats.vs_bytes[ZIO_TYPE_WRITE]);
				}

				kmem_free(vs, sizeof(vdev_stat_t));
			} else {
				dprintf("derp: warning!!!! Training but zio type isn't write.");
			}
			dprintf("derp: Training");
			*compress = derp_rand_compress();
			dprintf("derp: picked rand %u", *compress);
		} else {
			dprintf("derp: Not Training");
			*compress = derp_get_best_compress(ex_stats, zio->io_type, bc_bucket);
			dprintf("derp: picked %u", *compress);
		}
//		dprintf("derp: bytes: %lu", derp_get_vdev_stats(rvd)->vs_bytes[zio->io_type]);
		hrtime_t comp_start = gethrtime();
		psize = derp_compress_data(*compress, tmp, dst, s_len);
		abd_return_buf(zio->io_abd, tmp, s_len);
		hrtime_t comp_end = gethrtime();
//		rvd->vdev_stat_ex.vsx_derp_ac_model[bc_bucket][zio_compress_to_derp_compress(*compress)][0] = comp_end-comp_start;
//		dprintf("Hi there: 1: %u, 2: %u", rvd->vdev_stat_ex.vsx_derp_ac_model[0][0][0], rvd->vdev_stat_ex.vsx_derp_ac_model[1][0][0]);

		if (training) {
			derp_update_model(rvd, compress, comp_end-comp_start, zio->io_lsize, psize, bc_bucket);
		}
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

// This is a rewrite of the zio_compress_data function which allows us to not re-create a buf.
size_t derp_compress_data(enum zio_compress c, void *src, void *dst, size_t s_len)
{
	size_t c_len, d_len;
	zio_compress_info_t *ci = &zio_compress_table[c];

	ASSERT((uint_t)c < ZIO_COMPRESS_FUNCTIONS);
	ASSERT((uint_t)c == ZIO_COMPRESS_EMPTY || ci->ci_compress != NULL);

	if (c == ZIO_COMPRESS_EMPTY)
		return (s_len);

	/* Compress at least 12.5% */
	d_len = s_len - (s_len >> 3);

	c_len = ci->ci_compress(src, dst, s_len, d_len, ci->ci_level);

	if (c_len > d_len)
		return (s_len);

	ASSERT3U(c_len, <=, d_len);
	return (c_len);
}

static int derp_zeroed_cb(void *data, size_t len, void *private)
{
	uint64_t *end = (uint64_t *)((char *)data + len);
	uint64_t *word;

	for (word = data; word < end; word++)
		if (*word != 0)
			return (1);

	return (0);
}
