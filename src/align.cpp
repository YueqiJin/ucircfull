#include <iostream>
#include <string>
#include <cstring>
#include <algorithm>
#include <limits.h>
#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include "degenerateBase.hpp"
#include "align.hpp"
#include "minimap2/minimap.h"
#include "minimap2/ketopt.h"
#include "minimap2/mmpriv.h"
#include "minimap2/bseq.h"

using namespace std;
namespace fs = filesystem;

AlignPatternResult::AlignPatternResult(string seq, int score, int seq1Start, int seq1End, int seq2Start, int seq2End) : mappedString(seq), score(score), seq1Site(seq1Start, seq1End), seq2Site(seq2Start, seq2End) {}
AlignPatternResult::AlignPatternResult(string seq, string qual, int score, int seq1Start, int seq1End, int seq2Start, int seq2End) : mappedString(seq), mappedQual(qual), score(score), seq1Site(seq1Start, seq1End), seq2Site(seq2Start, seq2End) {}
AlignPatternResult::AlignPatternResult() {}

AlignOptions::AlignOptions(int matchScore, int mismatchScore, int openScore, int extendScore) : matchScore(matchScore), mismatchScore(mismatchScore), openScore(openScore), extendScore(extendScore), useScoreMat(false) {}
AlignOptions::AlignOptions(bool useScoreMat, const int *matchScores, const int *mismatchScores, int openScore, int extendScore) : openScore(openScore), extendScore(extendScore), useScoreMat(true)
{
	if (!matchScores)
		scoreMat[1][0] = 6, scoreMat[1][1] = 5, scoreMat[1][2] = 4, scoreMat[1][3] = 3;
	else
		memcpy(scoreMat[1], matchScores, 4 * sizeof(int));
	if (!mismatchScores)
		scoreMat[0][0] = -12, scoreMat[0][1] = -14, scoreMat[0][2] = -16, scoreMat[0][3] = -18;
	else
		memcpy(scoreMat[0], mismatchScores, 4 * sizeof(int));
}

void degenerateAlign(const DegenerateBaseVector &seqNum1, const DegenerateBaseVector &seqNum2, int len1, int len2, int *scoreMat, short *path, bool seq1LeftAlign, bool seq2LeftAlign, int matchScore, int mismatchScore, int openScore, int extendScore)
{
	int mat1 = len1 + 1, mat2 = len2 + 1;
	// clean score matrix
	memset(scoreMat, 0, sizeof(scoreMat));
	// define path for score matrix: {0: current position}, {1: up}, {2: left}, {3: up and left}
	memset(path, 0, sizeof(path));
	// init first row and col
	for (int i = 0; i < mat1; i++)
		*(int *)(scoreMat + i * mat2) = 0, *(int *)(path + i * mat2) = 0;
	for (int j = 0; j < mat2; j++)
		*(int *)(scoreMat + j) = 0, *(int *)(path + j) = 0;

	// calculating score matrix
	int lastScore[4], maxPath, maxScore = 0, maxX = 0, maxY = 0;
	int initzero = (seq1LeftAlign | seq2LeftAlign) ? INT_MIN : 0;
	for (int i = 1; i < mat1; i++)
		for (int j = 1; j < mat2; j++)
		{
			lastScore[0] = initzero, lastScore[1] = *(int *)(scoreMat + i * mat2 - mat2 + j), lastScore[2] = *(int *)(scoreMat + i * mat2 + j - 1), lastScore[3] = *(int *)(scoreMat + i * mat2 - mat2 + j - 1);
			// dual with up
			if (*(short *)(path + i * mat2 - mat2 + j) == 1)
			{
				lastScore[1] += extendScore;
			}
			else
			{
				lastScore[1] += openScore;
			}
			// dual with left
			if (*(short *)(path + i * mat2 + j - 1) == 2)
			{
				lastScore[2] += extendScore;
			}
			else
			{
				lastScore[2] += openScore;
			}
			// dual with match
			lastScore[3] += ((seqNum1[i - 1] & seqNum2[j - 1]) ? matchScore : mismatchScore);
			maxPath = max_element(lastScore, lastScore + 4) - lastScore;
			*(short *)(path + i * mat2 + j) = maxPath;
			*(int *)(scoreMat + i * mat2 + j) = lastScore[maxPath];
		}
}

void degenerateAlign(const DegenerateBaseVector &seqNum1, const DegenerateBaseVector &seqNum2, int len1, int len2, int *scoreMat, short *path, AlignOptions aligner, bool seq1LeftAlign, bool seq2LeftAlign)
{
	int mat1 = len1 + 1, mat2 = len2 + 1;
	// clean score matrix
	memset(scoreMat, 0, sizeof(scoreMat));
	// define path for score matrix: {0: current position}, {1: up}, {2: left}, {3: up and left}
	memset(path, 0, sizeof(path));
	// init first row and col
	for (int i = 0; i < mat1; i++)
		*(int *)(scoreMat + i * mat2) = 0, *(int *)(path + i * mat2) = 0;
	for (int j = 0; j < mat2; j++)
		*(int *)(scoreMat + j) = 0, *(int *)(path + j) = 0;

	// calculating score matrix
	int lastScore[4], maxPath, maxScore = 0, maxX = 0, maxY = 0;
	int initzero = (seq1LeftAlign | seq2LeftAlign) ? INT_MIN : 0;
	for (int i = 1; i < mat1; i++)
		for (int j = 1; j < mat2; j++)
		{
			lastScore[0] = initzero, lastScore[1] = *(int *)(scoreMat + i * mat2 - mat2 + j), lastScore[2] = *(int *)(scoreMat + i * mat2 + j - 1), lastScore[3] = *(int *)(scoreMat + i * mat2 - mat2 + j - 1);
			// dual with up
			if (*(short *)(path + i * mat2 - mat2 + j) == 1)
			{
				lastScore[1] += aligner.extendScore;
			}
			else
			{
				lastScore[1] += aligner.openScore;
			}
			// dual with left
			if (*(short *)(path + i * mat2 + j - 1) == 2)
			{
				lastScore[2] += aligner.extendScore;
			}
			else
			{
				lastScore[2] += aligner.openScore;
			}
			// dual with match
			lastScore[3] += aligner.getAlignScore(seqNum1[i - 1], seqNum2[j - 1]);
			maxPath = max_element(lastScore, lastScore + 4) - lastScore;
			*(short *)(path + i * mat2 + j) = maxPath;
			*(int *)(scoreMat + i * mat2 + j) = lastScore[maxPath];
		}
}

void findPath(int mat1, int mat2, int *scoreMat, short *path, bool seq1RightAlign, bool seq2RightAlign)
{
	int len1 = mat1 + 1, len2 = mat2 + 1;
	int maxScore = INT_MIN, startX, startY;
	for (int i = (seq1RightAlign) ? len1 - 1 : 0; i < len1; i++)
		for (int j = (seq2RightAlign) ? len2 - 1 : 0; j < len2; j++)
		{
			if (*(int *)(scoreMat + i * len2 + j) > maxScore)
				maxScore = *(int *)(scoreMat + i * len2 + j), startX = i, startY = j;
		}
	int nowX = startX, nowY = startY, newX, newY;
	while (*(short *)(path + nowX * len2 + nowY))
	{
		newX = nowX - ((*(short *)(path + nowX * len2 + nowY) == 2) ? 0 : 1), newY = nowY - ((*(short *)(path + nowX * len2 + nowY) == 1) ? 0 : 1);
		nowX = newX, nowY = newY;
	}
	cout << maxScore << " " << nowX + 1 << " " << startX + 1 << " " << nowY + 1 << " " << startY + 1 << endl;
}

AlignPatternResult matchSeq(const DegenerateBaseVector &query, const DegenerateBaseVector &pattern, int *scoreMat, short *path, int len1, int len2, bool seq1RightAlign, bool seq2RightAlign, bool completeN)
{
	len1++, len2++;
	int maxScore = INT_MIN, startX, startY;
	for (int i = (seq1RightAlign) ? len1 - 1 : 0; i < len1; i++)
		for (int j = (seq2RightAlign) ? len2 - 1 : 0; j < len2; j++)
		{
			if (*(int *)(scoreMat + i * len2 + j) > maxScore)
				maxScore = *(int *)(scoreMat + i * len2 + j), startX = i, startY = j;
		}
	int nowX = startX, nowY = startY, newX, newY;
	string ret = "";
	while (*(short *)(path + nowX * len2 + nowY))
	{
		switch (*(short *)(path + nowX * len2 + nowY))
		{
		case 1:
			ret = ret + query[nowX - 1].toChar();
			nowX--;
			break;
		case 2:
			if (completeN)
			{
				ret = ret + "N";
			}
			nowY--;
			break;
		case 3:
			ret = ret + query[nowX - 1].toChar();
			nowX--, nowY--;
			break;
		}
	}
	ret += "\0";
	reverse(ret.begin(), ret.end());
	return AlignPatternResult(ret, maxScore, nowX + 1, startX + 1, nowY + 1, startY + 1);
}

AlignPatternResult matchSeq(const DegenerateBaseVector &query, const string &qual, const DegenerateBaseVector &pattern, int *scoreMat, short *path, int len1, int len2, bool seq1RightAlign, bool seq2RightAlign, bool completeN)
{
	len1++, len2++;
	int maxScore = INT_MIN, startX, startY;
	for (int i = (seq1RightAlign) ? len1 - 1 : 0; i < len1; i++)
		for (int j = (seq2RightAlign) ? len2 - 1 : 0; j < len2; j++)
		{
			if (*(int *)(scoreMat + i * len2 + j) > maxScore)
				maxScore = *(int *)(scoreMat + i * len2 + j), startX = i, startY = j;
		}
	int nowX = startX, nowY = startY, newX, newY;
	string retSeq = "", retQual = "";
	while (*(short *)(path + nowX * len2 + nowY))
	{
		switch (*(short *)(path + nowX * len2 + nowY))
		{
		case 1:
			retSeq = retSeq + query[nowX - 1].toChar();
			retQual = retQual + qual[nowX - 1];
			nowX--;
			break;
		case 2:
			if (completeN)
			{
				retSeq = retSeq + "N";
				retQual = retQual + "!";
			}
			nowY--;
			break;
		case 3:
			retSeq = retSeq + query[nowX - 1].toChar();
			retQual = retQual + qual[nowX - 1];
			nowX--, nowY--;
			break;
		}
	}
	retSeq += "\0";
	retQual += "\0";
	reverse(retSeq.begin(), retSeq.end());
	reverse(retQual.begin(), retQual.end());
	return AlignPatternResult(retSeq, retQual, maxScore, nowX + 1, startX + 1, nowY + 1, startY + 1);
}

int minimap2Align(fs::path fastq, fs::path genome, fs::path sam, string minimap2, int thread, bool strand, int k, float p, bool quiet)
{
	string cmdLine;
	if (strand)
		cmdLine = minimap2 + " -ax splice -uf -k" + to_string(k) + " -p " + to_string(p) + " -t " + to_string(thread) + " " + genome.string() + " " + fastq.string() + " > " + sam.string() + (quiet ? " 2>/dev/null" : "");
	else
		cmdLine = minimap2 + " -ax splice -k" + to_string(k) + " -p " + to_string(p) + " -t " + to_string(thread) + " " + genome.string() + " " + fastq.string() + " > " + sam.string() + (quiet ? " 2>/dev/null" : "");
	return system(cmdLine.c_str());
}

int sam2bam(fs::path sam, fs::path bam, string samtools, int thread)
{
	string cmdLine = samtools + " sort -@ " + to_string(thread) + " " + sam.string() + " -o " + bam.string();
	return system(cmdLine.c_str());
}

int bamIndex(fs::path bam, string samtools)
{
	string cmdLine = samtools + " index " + bam.string();
	return system(cmdLine.c_str());
}

mpAligner::mpAligner(string &seq, int k, bool splice, bool uf)
{
	mm_set_opt(NULL, &idxOpt, &mapOpt);
	if (splice)
	{
		mm_set_opt("splice", &idxOpt, &mapOpt);
	}
	else
	{
		mm_set_opt("map-ont", &idxOpt, &mapOpt);
		idxOpt.batch_size = 0x7fffffffffffffffL;
		mapOpt.mid_occ = 1000;
	}
	idxOpt.k = k;
	mapOpt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
	mapOpt.pri_ratio = 0.8;
	if (uf)
		mapOpt.flag |= MM_F_SPLICE_FOR, mapOpt.flag &= ~MM_F_SPLICE_REV;
	char *s = (char *)calloc(seq.length() + 1, 1);
	memcpy(s, seq.c_str(), seq.length());
	const char *fakeName = "N/A";
	idx = mm_idx_str(idxOpt.w, idxOpt.k, false, idxOpt.bucket_bits, 1, (const char **)&s, (const char **)&fakeName);
	free(s);
	mm_mapopt_update(&mapOpt, idx);
}

mpAligner::~mpAligner()
{
	mm_idx_destroy(idx);
}

mm_reg1_t *mpAligner::map(string &seq, int *nRegs)
{
	mm_tbuf_t *tBuf = mm_tbuf_init();
	mm_reg1_t *ret = mm_map(idx, seq.length(), seq.c_str(), nRegs, tBuf, &mapOpt, 0);
	mm_tbuf_destroy(tBuf);
	return ret;
}

std::string mpAlignerCigar2String(mm_reg1_t *r)
{
	std::string cigar;
	for (int i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
		// printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
		cigar += std::to_string(r->p->cigar[i] >> 4) + MM_CIGAR_STR[r->p->cigar[i] & 0xf];
	return cigar;
}

std::vector<seqan3::cigar> mpAlignerCigar2Seqan3Cigar(mm_reg1_t *r)
{
	std::vector<seqan3::cigar> cigar;
	std::string cigarStr;
	for (int i = 0; i < r->p->n_cigar; ++i)
	{
		seqan3::cigar op;
		cigarStr = std::to_string(r->p->cigar[i] >> 4) + MM_CIGAR_STR[r->p->cigar[i] & 0xf];
		op.assign_string(cigarStr);
		cigar.push_back(op);
	}
	return cigar;
}