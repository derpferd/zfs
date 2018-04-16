#include <sys/derp_bytecounter.h>


#define SLOTS 256


uint64_t derp_to_Q_type(uint64_t bytecount)
{
//	int ldebug = 0;
    unsigned new = 0;

    /* three types */
    if (bytecount < 33) {
        new = 3;
    } else if (bytecount < 66) {
		new = 2;
	} else if (bytecount < 100) {
		new = 1;
	} else {
		new = 0;
	}

    /* else it remains 0 -- uncompressible */

//    if (ldebug) fprintf(stderr, "Qtype: was: %u now: %u\n", type, new);

    return new;
}


uint64_t derp_bytecounter(void *buf, size_t len)
{
	/* count the appearances of each unique byte between
	 * buf and buf + len.
	 * then, Normalize the values (divide each by len).
	 *
	 * if the value in the bucket is larger than threshold, add one to
	 * the count of unique "present" bytes.
	 *
	 * return the total number of "present bytes"
	 */

	// Maybe put this on the heap so we can use a larger integer type.
	uint16_t counts[SLOTS] = {0};


//	size_t index = 0;
	uint64_t full = 0;
	const unsigned threshold = len/SLOTS;
	uint16_t *count;
	unsigned quit_max = 0;

	for (uint8_t *buf8 = (uint8_t *)buf; buf8 < (uint8_t *)buf+len; buf8++) {

		count = counts + *buf8;

		if (*count < threshold) {

			*count += 1;

			if (*count >= threshold) {
				full++;

				/* we may want to quit out with some answer if the value
				 * reaches a certain level (e.g., 100). If quit_max is set
				 * to zero, we'll optimize this out and not do it.
				 */
				if (quit_max && full >= quit_max) {
					return full;
				}
			}
		}
	}

	return full;

}
