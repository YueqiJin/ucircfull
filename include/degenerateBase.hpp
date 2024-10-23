#pragma once

#include <vector>
#include <string>
#include <map>
#include <cstring>

const char __degenerateCharacter[16] = {'-', 'T', 'C', 'Y', 'G', 'K', 'S', 'B', 'A', 'W', 'M', 'H', 'R', 'D', 'V', 'N'};

/// @brief degenerate base class
struct DegenerateBase
{
	unsigned a : 4;
	DegenerateBase() : a(0){};
	DegenerateBase(const char &b);
	inline char toChar() const;
	bool operator&(const DegenerateBase &b) const;
	DegenerateBase &operator=(const int &b);
};
typedef std::vector<DegenerateBase> DegenerateBaseVector;

/// @brief covert base in char to DegenerateBase
/// @param b base in char
/// @return DegenerateBase
DegenerateBase convertBase(const char &b);

/// @brief covert base in string to DegenerateBaseVector
/// @param seq base in string
/// @return DegenerateBaseVector
DegenerateBaseVector convertBaseArray(std::string &&seq);

/// @brief covert base in DegenerateBase to char
/// @return char
inline char DegenerateBase::toChar() const
{
	return __degenerateCharacter[a];
}
