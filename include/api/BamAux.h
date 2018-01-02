#include <string.h>

#ifndef BAMAUX_H
#define BAMAUX_H
namespace BamTools {
	struct RefData
	{

		std::string RefName;  //!< name of reference sequence
		int32_t RefLength;    //!< length of reference sequence

		//! constructor
		RefData(const std::string& name = std::string(), const int32_t& length = 0)
			: RefName(name)
			, RefLength(length)
		{}
	};
}
#endif
