#include <sys/derp_ac.h>
#include <sys/spa_impl.h>
#include <sys/vdev_impl.h>
#include <sys/derp_bytecounter.h>
//#include <linux/blkdev_compat.h>


typedef uint64_t derp_compress;

#define derp_for_each_compress(derp_c) \
    for(derp_c = 0; derp_c < DERP_AC_NUM_OF_COMP_ALGS; derp_c++)
#define derp_max(a, b) \
	(a > b ? a : b)

#define cpu_stat_update_interval (SECOND/10)
#define bps_measuring_interval SECOND

enum zio_compress derp_compress_funcs[DERP_AC_NUM_OF_COMP_ALGS] = {
		ZIO_COMPRESS_EMPTY,
		ZIO_COMPRESS_LZ4,
		ZIO_COMPRESS_GZIP_1
};

size_t derp_compress_data(enum zio_compress c, void *src, void *dst, size_t s_len);
static int derp_zeroed_cb(void *data, size_t len, void *private);

typedef struct bps_variables {
	uint64_t bc_bucket;
	uint64_t cpu_bucket;
	derp_compress compress;
	uint64_t bps;
} derp_bps_variables_t;

typedef struct ratio_variables {
	uint64_t bc_bucket;
	derp_compress compress;
	uint64_t ratio;
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
#ifdef _SPL_TIME_H
	dprintf("cpu: good");
	int i;

	stats->cur_idle_time = 0;
	stats->cur_total_time = 0;

	for_each_possible_cpu(i) {
		for (int j=0; j < NR_STATS; j++) {
			stats->cur_total_time += kcpustat_cpu(0).cpustat[j];
			if (j == CPUTIME_IDLE || j == CPUTIME_IOWAIT) {
				stats->cur_idle_time += kcpustat_cpu(0).cpustat[j];
			}
		}
	}
#else
	dprintf("cpu: bad");
	stats->cur_idle_time = 0;
	stats->cur_total_time = 0;
#endif

    return 0;
}

int derp_first_cpuload(d_cpuload* stats) {
	dprintf("Setting up cpuload");
    int err;
    // Get cur times from /proc/stat
    err = derp_load_cur_stats(stats);
    if (err) return err;

    // Set last to cur
    stats->last_idle_time = stats->cur_idle_time;
    stats->last_total_time = stats->cur_total_time;

    return 0;
}

int derp_update_cpuload(d_cpuload* stats) {
	dprintf("Updating cpuload");
    int err;
    // Get cur times from /proc/stat
    err = derp_load_cur_stats(stats);
    if (err) return err;

    // Calc diff
    if (stats->last_total_time < stats->cur_total_time) {  // Time has moved on
    	if (stats->last_idle_time <= stats->cur_idle_time ) {  // The idle time hasn't gone down.
    		stats->idle_time_diff = stats->cur_idle_time - stats->last_idle_time;
    	} else {
    		dprintf("WARNING: This is very bad. The total idle time has gone down!!!");
    		return 0;
    	}
    	stats->total_time_diff = stats->cur_total_time - stats->last_total_time;

        // Set last to cur
        stats->last_idle_time = stats->cur_idle_time;
        stats->last_total_time = stats->cur_total_time;
    } else if (stats->last_total_time == stats->cur_total_time) {
    	dprintf("Not updating since nothing changed.");
    } else {
		dprintf("WARNING: This is very bad. We appear to have gone back in time!!!");
		dprintf("Last: %lu Cur: %lu", stats->last_total_time, stats->cur_total_time);

        // Set last to cur in hopes that this will reset the stat so they can be used next time.
        stats->last_idle_time = stats->cur_idle_time;
        stats->last_total_time = stats->cur_total_time;
        return 0;
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
	dprintf("derp: getting bps element[%u][%u][%u]", vars->bc_bucket, vars->cpu_bucket, vars->compress);
	return &stats->vsx_derp_ac_rate_model[vars->bc_bucket][vars->compress][vars->cpu_bucket];
}

uint64_t* derp_get_ratio(vdev_stat_ex_t *stats, derp_ratio_variables_t *vars) {
	dprintf("derp: getting ratio element[%u][%u]", vars->bc_bucket, vars->compress);
	return &stats->vsx_derp_ac_ratio_model[vars->bc_bucket][vars->compress];
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

uint64_t get_queue_time(vdev_t *vd, uint64_t disk_bps) {
	if (!vd->vdev_children) {
		uint64_t size = vd->vdev_queue.vq_class[ZIO_PRIORITY_ASYNC_WRITE].vqc_queued_size;
		return derp_calc_time_from_bps(disk_bps, size);
	}
	return get_queue_time(vd->vdev_child[0], disk_bps);
}

void derp_update_vdev_stats(vdev_t *vd, derp_bps_variables_t* bps_vars, derp_ratio_variables_t* ratio_vars) {
	if (!vd->vdev_children) { // is this the leaf?

		// Old style
		derp_add_to_rolling_ave(bps_vars->bps,
				&vd->vdev_stat_ex.vsx_derp_ac_model[bps_vars->bc_bucket][bps_vars->compress][DERP_COMPRESS_RATE_INDEX],
				1000);
		derp_add_to_rolling_ave(ratio_vars->ratio,
				&vd->vdev_stat_ex.vsx_derp_ac_model[ratio_vars->bc_bucket][ratio_vars->compress][DERP_COMPRESS_RATIO_INDEX],
				1000);

		// New style
		if (bps_vars) {
			dprintf("derp: vd:%u updating rate model[%u][%u][%u] with bps %u", vd->vdev_id, bps_vars->bc_bucket, bps_vars->cpu_bucket, bps_vars->compress, bps_vars->bps);
			derp_add_to_rolling_ave(bps_vars->bps,
					&vd->vdev_stat_ex.vsx_derp_ac_rate_model[bps_vars->bc_bucket][bps_vars->cpu_bucket][bps_vars->compress],
					1000);
			dprintf("derp: new value %lu", vd->vdev_stat_ex.vsx_derp_ac_rate_model[bps_vars->bc_bucket][bps_vars->cpu_bucket][bps_vars->compress]);
		}
		if (ratio_vars) {
			dprintf("derp: vd:%u updating ratio model[%u][%u] with ratio %u", vd->vdev_id, ratio_vars->bc_bucket, ratio_vars->compress, ratio_vars->ratio);
			derp_add_to_rolling_ave(ratio_vars->ratio,
					&vd->vdev_stat_ex.vsx_derp_ac_ratio_model[ratio_vars->bc_bucket][ratio_vars->compress],
					1000);
			dprintf("derp: new value %lu", vd->vdev_stat_ex.vsx_derp_ac_ratio_model[ratio_vars->bc_bucket][ratio_vars->compress]);
		}
	} else {
		int i;
		for (i = 0; i < vd->vdev_children; i++) {
			derp_update_vdev_stats(vd->vdev_child[i], bps_vars, ratio_vars);
		}
	}
}

void derp_update_model(vdev_t *rvd, derp_bps_variables_t* bps_vars, derp_ratio_variables_t* ratio_vars) {
	// A wrapper for the recursive function `derp_update_vdev_stats`
	derp_update_vdev_stats(rvd, bps_vars, ratio_vars);
}

enum zio_compress derp_rand_compress(void) {
	uint16_t i = rand_int16() % DERP_AC_NUM_OF_COMP_ALGS;
	return derp_compress_funcs[i];
//	return ZIO_COMPRESS_EMPTY;
}

enum zio_compress derp_get_best_compress(vdev_t *rvd, uint64_t lsize, derp_bps_variables_t* bps_vars, derp_ratio_variables_t* ratio_vars) {
	/*
	 * How to get the best compression alg.
	 * Find the compression alg. that maximizes Total Write Rate (TWR).
	 *  TWR = CR(Compression Rate) + DR(Current Disk Rate)*CRatio(Compression Ratio)
	 *  TWT = CR(Compression time) + (DR(Current Disk Rate)*CRatio(Compression Ratio)
	 *
	 */
//	dprintf("Yo Dog 1");
	enum zio_compress c;
	derp_compress derp_c;
	vdev_stat_ex_t *stats = derp_get_vdev_ex_stats(rvd);

//	uint64_t cur_bps = stats->vsx_derp_disk_bps[type];
//	uint64_t cur_total_bps = stats->vsx_derp_total_bps[type];
	uint64_t max_disk_bps = stats->max_bps;

	uint64_t min_total_write_time = 0xFFFFFFFFFFFFFFFF;

	uint64_t queue_time = get_queue_time(rvd, max_disk_bps);

	enum zio_compress best_c = ZIO_COMPRESS_EMPTY; // Default to no compression
	derp_set_default_compress(&best_c);

	// Check preconditions
	if (max_disk_bps == 0) {
		dprintf("derp: Warning!!!! no disk bps!");
		return best_c;
	}

	dprintf("Picking best...");
	dprintf("\tMax Disk BPS: %lu", max_disk_bps);
	derp_for_each_compress(derp_c) {
		c = derp_compress_funcs[derp_c];
		bps_vars->compress = derp_c;
		ratio_vars->compress = derp_c;

		uint64_t c_rate = *derp_get_bps(stats, bps_vars);
		uint64_t c_ratio = *derp_get_ratio(stats, ratio_vars);

		uint64_t compress_time;
		if (c_rate > MAX_RATE) {  // Basically the compression rate is so quick that it doesn't have any influence. E.g. No compression.
			compress_time = 0;
		} else if (c_rate == 0) { // We don't know how fast this algorithm is.
			dprintf("derp: Warning!!!! c_rate is 0!");
			compress_time = 0;
		} else {
			compress_time = derp_calc_time_from_bps(c_rate, lsize);
		}

		uint64_t pre_write_time = derp_max(compress_time, queue_time);
		uint64_t estimated_after_compression_bytes = derp_apply_compress_ratio(lsize, c_ratio);
		uint64_t disk_write_time = derp_calc_time_from_bps(max_disk_bps, estimated_after_compression_bytes);
		uint64_t total_write_time = pre_write_time + disk_write_time;

		if (total_write_time < min_total_write_time) {
			min_total_write_time = total_write_time;
			best_c = c;
		}
	}

	// TODO: set real value
	bps_vars->compress = -1;
	ratio_vars->compress = -1;
	return best_c;
}

void update_bps(vdev_stat_ex_t* ex_stats, vdev_stat_t* vs) {
	hrtime_t cur_timestamp = gethrtime();

	if (ex_stats->old_timestamp == 0) {
		ex_stats->old_timestamp = cur_timestamp;
		ex_stats->old_bytes = vs->vs_bytes[ZIO_TYPE_WRITE];
		dprintf("derp: Setting up bps counter");
	} else if (ex_stats->old_timestamp + bps_measuring_interval < cur_timestamp) {
		uint64_t tmp_bytes = vs->vs_bytes[ZIO_TYPE_WRITE] - ex_stats->old_bytes;
		uint64_t tdelta = cur_timestamp - ex_stats->old_timestamp;
		uint64_t bps = derp_calc_bps(tmp_bytes, tdelta);
		dprintf("derp: bps: %lu. tmp_bytes: %lu in %lu ns", bps, tmp_bytes, tdelta);
		if (bps > ex_stats->max_bps) {
			ex_stats->max_bps = bps;
			dprintf("derp new max bps: %lu", bps);
		}
		ex_stats->old_timestamp = cur_timestamp;
		ex_stats->old_bytes = vs->vs_bytes[ZIO_TYPE_WRITE];
	} else {
//		dprintf("derp: skipping bps recalc. old: %lu cur: %lu", ex_stats->old_timestamp, cur_timestamp);
	}
}

void update_cpu(vdev_stat_ex_t* ex_stats) {
	// Update cpu stats
	if (ex_stats->cpustats.timestamp == 0) {  // Haven't set the stats yet.
		derp_first_cpuload(&ex_stats->cpustats);  // TODO: error checking here
		ex_stats->cpustats.timestamp = gethrtime();
	} else if (ex_stats->cpustats.timestamp + cpu_stat_update_interval < gethrtime()) {
		derp_update_cpuload(&ex_stats->cpustats);  // TODO: error checking here
		ex_stats->cpustats.timestamp = gethrtime();
	} else {
	  // Skip, it is too soon to update.
	//      dprintf("CPU Stat skipping...");
	}
}

uint64_t derp_idle_cpu_percent_to_bucket(uint64_t idle_cpu) {
    if (idle_cpu == 0) {
        return 3;
    } else if (idle_cpu < 33) {
		return 2;
	} else if (idle_cpu < 66) {
		return 1;
	} else {
		return 0;
	}
}

size_t derp_ac_compress(zio_t *zio, void *dst, size_t s_len, enum zio_compress *compress) {
	/* Params:
	 * 	zio:
	 * 	compress:	This is needed so we can change the compression algorithm used which is  set by the calling function.
	 * */
	// If we want to die early and do no compression simply return 0.

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

		// Get ex stats
		vdev_stat_ex_t* ex_stats = derp_get_vdev_ex_stats(rvd);
		update_cpu(ex_stats);

		/* No compression algorithms can read from ABDs directly */
		void *tmp = abd_borrow_buf_copy(zio->io_abd, s_len);

		uint64_t bytecount = derp_bytecounter(tmp, s_len);  // TODO: make this read directly from the ABD.
		uint64_t cpu_idle_percent = ex_stats->cpustats.idle_percent;

		uint64_t bc_bucket = derp_to_Q_type(bytecount);
		uint64_t cpu_bucket = derp_idle_cpu_percent_to_bucket(cpu_idle_percent);
		uint64_t vd_queued_size_write = get_queue_size(rvd);

		derp_bps_variables_t bps_vars;
		derp_ratio_variables_t ratio_vars;

		bps_vars.bc_bucket = bc_bucket;
		bps_vars.cpu_bucket = cpu_bucket;
		bps_vars.compress = -1; // if training {Fill this out when we know it} else {this will be set to the best alg}

		ratio_vars.bc_bucket = bc_bucket;
		ratio_vars.compress = -1;

		// Debug print all inputs
		dprintf("derp: bytecount: %u bucket: %u cpu: %u bucket: %u queue: %u lsize: %u", bytecount, bc_bucket, cpu_idle_percent, cpu_bucket, vd_queued_size_write, zio->io_lsize);


		if (training) {
			dprintf("derp: Training");
			if (zio->io_type == ZIO_TYPE_WRITE) {
				update_bps(ex_stats, derp_get_vdev_stats(rvd));
			} else {
				dprintf("derp: warning!!!! Training but zio type isn't write.");
			}

			*compress = derp_rand_compress();
		} else {
			dprintf("derp: Not Training");
			*compress = derp_get_best_compress(rvd, zio->io_lsize, &bps_vars, &ratio_vars);
		}

		dprintf("derp: picked %u", *compress);
		bps_vars.compress = zio_compress_to_derp_compress(*compress);
		ratio_vars.compress = bps_vars.compress;

		hrtime_t comp_start = gethrtime();
		psize = derp_compress_data(*compress, tmp, dst, s_len);
		abd_return_buf(zio->io_abd, tmp, s_len);
		hrtime_t comp_end = gethrtime();

		if (training || DERP_ONLINE) {
			hrtime_t compress_delay = comp_end-comp_start;
			if (compress_delay > 0 && zio->io_lsize > 0 && psize > 0) {  // This makes sure we will not do any division by 0
				// Update structs with the values seen.
				bps_vars.bps = derp_calc_bps(zio->io_lsize, compress_delay);
				ratio_vars.ratio = derp_calc_compress_ratio(zio->io_lsize, psize);

				derp_update_model(rvd, &bps_vars, &ratio_vars);
			} else {
				dprintf("derp: warning!!!!! not updating stats since we would have a division by zero.");
				dprintf("derp: comp_delay: %lu lsize: %lu psize: %lu", compress_delay, zio->io_lsize, psize);
			}
		}
	}

	return psize;
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

