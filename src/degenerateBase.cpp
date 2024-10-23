#include <string>
#include <cstring>
#include <vector>
#include "degenerateBase.hpp"

DegenerateBase::DegenerateBase(const char &b)
{
	switch (b)
	{
	case '-':
		a = 0x0;
		break;
	case 'T':
		a = 0x1;
		break;
	case 'C':
		a = 0x2;
		break;
	case 'Y':
		a = 0x3;
		break;
	case 'G':
		a = 0x4;
		break;
	case 'K':
		a = 0x5;
		break;
	case 'S':
		a = 0x6;
		break;
	case 'B':
		a = 0x7;
		break;
	case 'A':
		a = 0x8;
		break;
	case 'W':
		a = 0x9;
		break;
	case 'M':
		a = 0xA;
		break;
	case 'H':
		a = 0xB;
		break;
	case 'R':
		a = 0xC;
		break;
	case 'D':
		a = 0xD;
		break;
	case 'V':
		a = 0xE;
		break;
	case 'N':
		a = 0xF;
		break;
	}
}

bool DegenerateBase::operator&(const DegenerateBase &b) const
{
	return this->a & b.a;
}

DegenerateBase &DegenerateBase::operator=(const int &b)
{
	this->a = b;
	return *this;
}

DegenerateBase convertBase(const char &b)
{
	DegenerateBase a;
	switch (b)
	{
	case '-':
		a = 0x0;
		break;
	case 'T':
		a = 0x1;
		break;
	case 'C':
		a = 0x2;
		break;
	case 'Y':
		a = 0x3;
		break;
	case 'G':
		a = 0x4;
		break;
	case 'K':
		a = 0x5;
		break;
	case 'S':
		a = 0x6;
		break;
	case 'B':
		a = 0x7;
		break;
	case 'A':
		a = 0x8;
		break;
	case 'W':
		a = 0x9;
		break;
	case 'M':
		a = 0xA;
		break;
	case 'H':
		a = 0xB;
		break;
	case 'R':
		a = 0xC;
		break;
	case 'D':
		a = 0xD;
		break;
	case 'V':
		a = 0xE;
		break;
	case 'N':
		a = 0xF;
		break;
	}
	return a;
}

DegenerateBaseVector convertBaseArray(std::string &&seq)
{
	int n = seq.length();
	DegenerateBaseVector ret(n);
	ret[0] = 0;
	for (int i = 0; i < n; i++)
	{
		ret[i] = convertBase(seq[i]);
	}
	return ret;
}