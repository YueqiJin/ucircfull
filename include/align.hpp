#pragma once

#include "degenerateBase.hpp"
#include <filesystem>
#include <string>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include "minimap2/minimap.h"

/// @brief alignment result
class AlignPatternResult
{
public:
	std::string mappedString, mappedQual;
	int score;
	std::pair<int, int> seq1Site, seq2Site;

	/// @brief initialize alignment result
	/// @param seq query sequence edited by reference sequence
	/// @param score alignment score
	/// @param seq1Start start site of query sequence
	/// @param seq1End end site of query sequence
	/// @param seq2Start start site of reference sequence
	/// @param seq2End end site of reference sequence
	AlignPatternResult(std::string seq, int score, int seq1Start, int seq1End, int seq2Start, int seq2End);

	/// @brief initialize alignment result with sequence quality
	/// @param seq query sequence edited by reference sequence
	/// @param qual edited query sequence quality
	/// @param score alignment score
	/// @param seq1Start start site of query sequence
	/// @param seq1End end site of query sequence
	/// @param seq2Start start site of reference sequence
	/// @param seq2End end site of reference sequence
	AlignPatternResult(std::string seq, std::string qual, int score, int seq1Start, int seq1End, int seq2Start, int seq2End);

	/// @brief empty alignment result
	AlignPatternResult();
};

/// @brief alignment options
class AlignOptions
{
public:
	int matchScore, mismatchScore, openScore, extendScore;
	bool useScoreMat;
	int scoreMat[2][4]; //[0][] mismatch, [1][] match

	/// @brief set alignment with scoring matrix
	/// @param matchScore score when base matched
	/// @param mismatchScore score when base mismatched
	/// @param openScore score when a open gap
	/// @param extendScore score when a extend gap
	AlignOptions(int matchScore = 3, int mismatchScore = -6, int openScore = -5, int extendScore = -2);

	/// @brief set alignment with scoring matrix, different base when matching or mismatching have different match score
	/// @param useScoreMat true when different base when different match
	/// @param matchScores score matrix when base match
	/// @param mismatchScores score matrix when base mismatch
	/// @param openScore score when a open gap
	/// @param extendScore score when a extend gap
	AlignOptions(bool useScoreMat, const int *matchScores = NULL, const int *mismatchScores = NULL, int openScore = -10, int extendScore = -4);

	/// @brief get match score of base a and b
	/// @param a
	/// @param b
	/// @return score
	inline int getAlignScore(const DegenerateBase &a, const DegenerateBase &b) const;

private:
	/// @brief judge base a and b are matched or mismatched
	/// @param a
	/// @param b
	/// @return true when matched
	inline bool judgeMatch(const DegenerateBase &a, const DegenerateBase &b) const;

	/// @brief get score involve with degenerate base of pattern
	/// @param pattern pattern base
	/// @param match true if match
	/// @return match number
	inline int judgePattern(const DegenerateBase &pattern, const bool &match) const;
};

/// @brief adjust SW algorithm to align degenerate base
/// @param seqNum1 query sequence
/// @param seqNum2 reference sequence
/// @param len1 length of query sequence
/// @param len2 length of reference sequence
/// @param scoreMat score matrix
/// @param path align path matrix to record alignment
/// @param seq1LeftAlign true if left align query sequence
/// @param seq2LeftAlign true if left align reference sequence
/// @param matchScore score when base matched
/// @param mismatchScore score when base mismatched
/// @param openScore score when a open gap
/// @param extendScore score when a extend gap
void degenerateAlign(const DegenerateBaseVector &seqNum1, const DegenerateBaseVector &seqNum2, int len1, int len2, int *scoreMat, short *path, bool seq1LeftAlign = false, bool seq2LeftAlign = false, int matchScore = 3, int mismatchScore = -6, int openScore = -5, int extendScore = -2);

/// @brief adjust SW algorithm to align degenerate base
/// @param seqNum1 query sequence
/// @param seqNum2 reference sequence
/// @param len1 length of query sequence
/// @param len2 length of reference sequence
/// @param scoreMat score matrix
/// @param path align path matrix to record alignment
/// @param aligner alignment options
/// @param seq1LeftAlign true if left align query sequence
/// @param seq2LeftAlign true if left align reference sequence
void degenerateAlign(const DegenerateBaseVector &seqNum1, const DegenerateBaseVector &seqNum2, int len1, int len2, int *scoreMat, short *path, AlignOptions aligner, bool seq1LeftAlign = false, bool seq2LeftAlign = false);

/// @brief get alignment result from score matrix and path matrix and output to the std (debug function)
/// @param mat1 length of query sequence
/// @param mat2 length of reference sequence
/// @param scoreMat score matrix
/// @param path align path matrix
/// @param seq1RightAlign true if right align query sequence
/// @param seq2RightAlign true if right align reference sequence
void findPath(int mat1, int mat2, int *scoreMat, short *path, bool seq1RightAlign = false, bool seq2RightAlign = false);

/// @brief get alignment result from score matrix and path matrix
/// @param query query sequence
/// @param pattern reference sequence
/// @param scoreMat score matrix
/// @param path align path matrix
/// @param len1 length of query sequence
/// @param len2 length of reference sequence
/// @param seq1RightAlign true if right align query sequence
/// @param seq2RightAlign true if right align reference sequence
/// @param completeN true if complete N in query sequence
/// @return alignment result
AlignPatternResult matchSeq(const DegenerateBaseVector &query, const DegenerateBaseVector &pattern, int *scoreMat, short *path, int len1, int len2, bool seq1RightAlign = false, bool seq2RightAlign = false, bool completeN = false);

/// @brief get alignment result from score matrix and path matrix
/// @param query query sequence
/// @param qual query sequence quality
/// @param pattern reference sequence
/// @param scoreMat score matrix
/// @param path align path matrix
/// @param len1 length of query sequence
/// @param len2 length of reference sequence
/// @param seq1RightAlign true if right align query sequence
/// @param seq2RightAlign true if right align reference sequence
/// @param completeN true if complete N in query sequence
/// @return alignment result
AlignPatternResult matchSeq(const DegenerateBaseVector &query, const std::string &qual, const DegenerateBaseVector &pattern, int *scoreMat, short *path, int len1, int len2, bool seq1RightAlign = false, bool seq2RightAlign = false, bool completeN = false);

/// @brief judge if base a and b are matched in degenerate base
/// @param a base a
/// @param b base b
/// @return true if matched
inline bool AlignOptions::judgeMatch(const DegenerateBase &a, const DegenerateBase &b) const
{
	return a.a & b.a;
}

/// @brief get score involve with degenerate base of pattern
/// @param pattern pattern base
/// @param match true if match
/// @return match score
inline int AlignOptions::judgePattern(const DegenerateBase &pattern, const bool &match) const
{
	return (pattern.a & 1) + ((pattern.a >> 1) & 1) + ((pattern.a >> 2) & 1) + ((pattern.a >> 3) & 1) - match;
}

/// @brief get match score of base a and b
/// @param a base a
/// @param b base b
/// @return match score
inline int AlignOptions::getAlignScore(const DegenerateBase &a, const DegenerateBase &b) const
{
	bool tmp = judgeMatch(a, b);
	return useScoreMat ? scoreMat[tmp][judgePattern(b, tmp)] : (tmp ? matchScore : mismatchScore);
}

/// @brief default alignment options
const AlignOptions defaultAlignOption(true, NULL, NULL, -10, -4);

/// @brief align sequence fastq to reference with minimap2
/// @param fastq fastq file
/// @param genome reference genome
/// @param sam output sam file
/// @param minimap2 minimap2 path
/// @param thread number of thread
/// @param strand stranded mode if TRUE, use parameter "-uf"
/// @param k k-mer size
/// @param p minimizer chain drop ratio
/// @param quiet true if quiet mode
/// @return 0 if success
int minimap2Align(std::filesystem::path fastq, std::filesystem::path genome, std::filesystem::path sam, std::string minimap2, int thread, bool strand, int k = 15, float p = 0.5, bool quiet = false);
// int minimap2Realign(seqan3::dna5_vector readSeq, seqan3::dna5_vector refSeq, std::filesystem::path sam, string minimap2, int thread);

/// @brief convert sam to bam
/// @param sam sam file
/// @param bam output bam file
/// @param samtools samtools path
/// @param thread number of thread
/// @return 0 if success
int sam2bam(std::filesystem::path sam, std::filesystem::path bam, std::string samtools, int thread);

/// @brief generate bam index
/// @param bam bam file
/// @param samtools samtools path
/// @return 0 if success
int bamIndex(std::filesystem::path bam, std::string samtools);

/// @brief aligner class for minimap2 api
class mpAligner
{
public:
	mm_idx_t *idx;
	mm_idxopt_t idxOpt;
	mm_mapopt_t mapOpt;

	/// @brief initialize minimap2 aligner
	/// @param seq reference sequence
	/// @param k k-mer size
	/// @param splice true if splice
	/// @param uf true if find GU-AG splice site on forward strand
	mpAligner(std::string &seq, int k = 10, bool splice = false, bool uf = false);
	~mpAligner();

	/// @brief map sequence to reference
	/// @param seq query sequence
	/// @param nRegs number of alignment result (part of return)
	/// @return alignment result
	mm_reg1_t *map(std::string &seq, int *nRegs);
};

/// @brief convert minimap2 mapping result to cigar string
/// @param r minimap2 mapping result
/// @return cigar string
std::string mpAlignerCigar2String(mm_reg1_t *r);

/// @brief extract cigar sequence from minimap2 mapping result
/// @param r minimap2 mapping result
/// @return cigar sequence
std::vector<seqan3::cigar> mpAlignerCigar2Seqan3Cigar(mm_reg1_t *r);