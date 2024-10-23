#ifndef __LIB_CCCS_HPP
#define __LIB_CCCS_HPP

#pragma once

#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <new>

namespace cccs {

	struct CcsResult {
		char *segment_ptr;
		char *ccs_seq_ptr;
	};

	extern "C" {

	CcsResult *ccs(const uint8_t *seq, unsigned int seq_len);

	void free_ccs_result(CcsResult *result);

	} // extern "C"

} // namespace cccs

#endif // __LIB_CCCS_HPP
