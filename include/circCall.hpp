#pragma once

#include <filesystem>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <boost/icl/interval_map.hpp>
#include "circfullUtils.hpp"

namespace ciri
{
	/// @brief circRNA calling with CIRI-long
	/// @param ciri_long path to CIRI-long
	/// @param inputFile input file
	/// @param outputDir output directory
	/// @param outputPrefix output prefix
	/// @param reference reference genome
	/// @param annotation annotation file
	/// @param nthread number of threads
	/// @return 0 if success, 1 if failed
	int ciriCall(std::string ciri_long, std::filesystem::path inputFile, std::filesystem::path outputDir, std::string outputPrefix, std::filesystem::path reference, std::filesystem::path annotation, int nthread);

	/// @brief circRNA collapse with CIRI-long
	/// @param ciri_long path to CIRI-long
	/// @param inputFile input file
	/// @param outputDir output directory
	/// @param outputPrefix output prefix
	/// @param reference reference genome
	/// @param annotation annotation file
	/// @param nthread number of threads
	/// @return 0 if success, 1 if failed
	int ciriCollapse(std::string ciri_long, std::filesystem::path inputFile, std::filesystem::path outputDir, std::string outputPrefix, std::filesystem::path reference, std::filesystem::path annotation, int nthread);
};

namespace circfull
{

	struct referenceStorage
	{
		std::map<std::string, seqan3::dna5_vector> seqs;
		std::string getSeqByPosition(std::string chr, int start, int end);
	};

	class GtfRecord
	{
	public:
		enum class FeatureType
		{
			gene,
			transcript,
			exon,
			CDS,
			UTR,
			start_codon,
			stop_codon,
			circRNA,
			unknown
		} type;
		std::string chr, geneId, transcriptId, geneName;
		char strand;
		int start, end;
		bool operator==(const GtfRecord &b) const;
	};
	struct GtfCompare
	{
		int operator()(const std::tuple<int, int> &x, const std::tuple<int, int> &k) const
		{
			// if x {start, end} >= y {start, end}, return 1, else, return 0
			auto [xs, xe] = x;
			auto [ys, ye] = k;
			if (xs > ys)
				return 1;
			if (xs < ys)
				return 0;
			if (xe >= ye)
				return 1;
			return 0;
		}
	};
	typedef std::multimap<std::string, GtfRecord> GtfStorage;
	typedef std::multimap<std::string, GtfRecord>::const_iterator GtfIter;

	// typedef std::map<std::string, std::multimap<std::tuple<int, int>, GtfRecord, GtfCompare>> GtfFetcher;
	// std::optional<std::vector<GtfRecord>> fetchGtf(const GtfFetcher & gtfFetcher, std::string chr, int start, int end);

	typedef GtfRecord *GtfRecordPtr;
	class GtfFetcherChrom
	{
	protected:
		typedef boost::icl::closed_interval<int> Interval;
		// use boost interval_map to store GtfRecordPtr
		boost::icl::interval_map<int, std::set<GtfRecordPtr>> fetcher;

	public:
		void addRecord(const GtfRecordPtr &gtfList);
		std::optional<std::set<GtfRecordPtr>> fetch(int start, int end);
	};
	class GtfFetcher
	{
	protected:
		std::map<std::string, GtfFetcherChrom> fetchers;

	public:
		void addGtf(const GtfStorage &gtfList);
		std::optional<std::set<GtfRecordPtr>> fetch(std::string chr, int start, int end);
	};

	referenceStorage readReference(std::filesystem::path const &referencePath);
	GtfStorage readGtf(std::filesystem::path const &gtfPath);

	std::istream &operator>>(std::istream &is, GtfRecord &record);

	/// @brief scan ccs reads from fasta file
	/// @param inFile input fasta file
	/// @param outPrefix output prefix
	/// @param nthread number of threads
	/// @return map<string, tuple<string, string, string>>: readId->(segment, ccs, rawSeq)
	std::map<std::string, std::tuple<std::string, std::string, std::string>> findCcsReads(std::filesystem::path inFile, std::filesystem::path outPrefix, int nthread);
	/// @brief construct consensus sequence from ccs reads with the same UMI
	/// @param ccsSeq map<string, tuple<string, string, string>>: readId->(segment, ccs, rawSeq)
	/// @param umiInfo map<string, string> umiInfo: readId->clustId
	/// @param nthread number of threads
	/// @param outPrefix output prefix
	/// @return map<string, string> consensus sequence: clustId->consensus sequence
	std::map<std::string, std::string> constructCCSConsensus(std::map<std::string, std::tuple<std::string, std::string, std::string>> &ccsSeq, const std::map<std::string, std::string> &umiInfo, int nthread, std::filesystem::path outPrefix);

	namespace RG
	{
		class seqan3Record
		{
		public:
			int refBegin;
			std::string refName, queryId;
			std::vector<seqan3::cigar> cigarSeq;
			seqan3::sam_flag flag;
			seqan3Record(const int &_refBegin, const std::string &_queryId, const std::string &_refName, const std::vector<seqan3::cigar> &_cigarSeq, const seqan3::sam_flag &_flag) : refBegin(_refBegin), queryId(_queryId), refName(_refName), cigarSeq(_cigarSeq), flag(_flag){};
		};

		class MapRecord
		{
		public:
			std::string queryId, chr;
			int refStart, refEnd, queryStart, queryEnd;
			std::string exonStart, exonEnd; //, exonLength;
			unsigned long long int exonLength;
			char strand;
			std::string leftSeq, rightSeq;
			bool dup;
			bool isBest;
			bool operator<(const MapRecord &b) const;
			MapRecord() : dup(false), isBest(false){};
			MapRecord(const std::string &_queryId, const std::string &_chr,
					  const int &_refStart, const int &_refEnd, const int &_queryStart, const int &_queryEnd,
					  const std::string &_exonStart, const std::string &_exonEnd, const unsigned long long int &_exonLength,
					  const char &_strand, const std::string &_leftSeq, const std::string &_rightSeq) : queryId(_queryId), chr(_chr), refStart(_refStart), refEnd(_refEnd),
																										queryStart(_queryStart), queryEnd(_queryEnd),
																										exonStart(_exonStart), exonEnd(_exonEnd), exonLength(_exonLength),
																										strand(_strand), leftSeq(_leftSeq), rightSeq(_rightSeq){};
			MapRecord(const int &_refStart, const int &_refEnd, const int &_queryStart, const int &_queryEnd) : refStart(_refStart), refEnd(_refEnd), queryStart(_queryStart), queryEnd(_queryEnd){};
		};
		typedef std::vector<MapRecord>::iterator MapRecordIter;

		class SpliceMotif
		{
		public:
			char first, second;
			SpliceMotif() = default;
			explicit SpliceMotif(const char *seq);
			explicit SpliceMotif(const std::string &seq);
			explicit SpliceMotif(const std::string &&seq);
			inline explicit operator const std::string() const
			{
				return std::string(1, first) + std::string(1, second);
			};
			std::string str() const;
		};

		class BSRecord
		{
		public:
			std::string queryId, chr;
			// std::vector<char> BSstrand;
			std::string BSstrand;
			std::vector<std::pair<int, int>> BSPos;					  // (BSLeft, BSRight) Back-splicing position upstream and downstream seq position (strand +)
			std::vector<std::pair<SpliceMotif, SpliceMotif>> BSMotif; // (BSLeftMotif, BSRightMotif) BS left and right 2 nucleotides
			// std::string leftSeqMapStatus, rightSeqMapStatus;
			bool leftSeqMapStatus, rightSeqMapStatus;
			bool adjust = false;
			// BSRecord();
			// BSRecord(const BSRecord &b);
			// BSRecord(BSRecord &&b);
			// BSRecord & operator=(const BSRecord &b);
			// BSRecord & operator=(BSRecord &&b);
			//~BSRecord() = default;
		};
		typedef std::vector<BSRecord>::iterator BSRecordIter;

		// typedef std::tuple<std::string, char, int, int, std::string, std::string> CircID;
		class CircID
		{
		public:
			std::string chr;
			char strand;
			int start, end;
			std::string leftSeq, rightSeq;
			bool operator<(const CircID &b) const;
			bool operator==(const CircID &b) const;
			bool operator==(const std::tuple<std::string, char, int, int> &b) const;
			CircID(){};
			CircID(const std::string &_chr, const char &_strand, const int &_start, const int &_end, const std::string &_leftSeq, const std::string &_rightSeq) : chr(_chr), strand(_strand), start(_start), end(_end), leftSeq(_leftSeq), rightSeq(_rightSeq){};
		};

		template <typename KeyType = std::string>
		class ExonIndex
		{
		public:
			std::map<std::string, std::set<KeyType>> exonStartIndex, exonEndIndex;
			std::map<KeyType, char> exonStrandIndex;

			void addExon(const std::string &chr, const int &exonStart, const int &exonEnd, const char &strand, const KeyType &geneId, const int errorLen);
			// void saveIndex(std::filesystem::path outPrefix);
			// void loadIndex(std::filesystem::path outPrefix);
		};

		struct BSseq
		{
			std::string seq, revSeq;
		};

		bool hitBSMotif(const std::string &seq, const BSseq &motif, const char &strand);

		class BSMotifIndex
		{
		public:
			std::vector<BSseq> leftMotif, rightMotif, bothMotif;
			void addIndex(std::string _seq);
			BSMotifIndex(){};
		};

		typedef struct
		{
			int start, end;
			char strand;
			std::string leftSeq, rightSeq;
		} SpliceSite;

		class CircRecord
		{
		public:
			std::string queryId;
			std::string chr;
			char strand;
			int start, end;
			std::vector<std::pair<int, int>> exons;
			// std::vector<std::pair<std::string, std::string>> spliceSeqs;
			unsigned long long int len;

			// update length of circRNA
			void updateLength();
			std::string getCircId();
			std::string getExonsString();
			std::string getTranscriptId();
		};

		std::istream &operator>>(std::istream &in, MapRecord &record);

		std::string to_string(const MapRecord *record);
		std::string to_string(const BSRecord *record);
		std::string to_string(const CircRecord *record);

		ExonIndex<std::string> getExonIndex(const GtfStorage &gtfList, const int errorLen);

		std::vector<MapRecord> parseBam(referenceStorage &genome, std::filesystem::path output, std::filesystem::path bam, int nthread);

		std::tuple<std::vector<MapRecord>, std::vector<MapRecord>, std::vector<MapRecord>> filterCandidateCircReads(/*std::filesystem::path mapInfo*/ std::vector<MapRecord> &records, referenceStorage &genome, std::filesystem::path fastq, std::filesystem::path oPrefix, int nthread);

		std::vector<BSRecord> detectBS(referenceStorage &genome, std::filesystem::path BSListFile, std::vector<MapRecord> &circCandidates, int nthread, std::string oPrefix, std::string minimap2);

		std::vector<BSRecord> filterBS(referenceStorage &genome, std::vector<BSRecord> &BSList, GtfStorage gtfRecords, const ExonIndex<> &exonIndex, std::filesystem::path filteredBSListFile, const std::string &spliceSignal, int nthread);

		std::vector<BSRecord> clusterBS(const std::vector<BSRecord> &BSList, std::filesystem::path BSAdjustFile, int nthread);

		std::vector<CircRecord> constructFullStruct(std::vector<BSRecord> &BSList, std::vector<MapRecord> &mapList, const referenceStorage &genome, const ExonIndex<> &exonIndex, GtfStorage &gtfRecords, std::filesystem::path fullStructFile, int nthread);

		std::vector<CircRecord> clustFullStruct(std::vector<CircRecord> &circList, std::filesystem::path adjFullStructFile, int nthread);

		std::vector<circfull::RG::CircRecord> RG(CircfullOption &opt);
	};
	const int hangLen = 300;
	const int errorBSLen = 40;
	const int errorFSLen = 14;
};

template <typename KeyType>
void circfull::RG::ExonIndex<KeyType>::addExon(const std::string &chr, const int &exonStart, const int &exonEnd, const char &strand, const KeyType &geneId, const int errorLen)
{
	exonStrandIndex[geneId] = strand;
	for (int i = 0; i <= errorLen * 2; i++)
	{
		std::string startKey = chr + "_" + std::to_string(exonStart - errorLen + i);
		std::string endKey = chr + "_" + std::to_string(exonEnd - errorLen + i);
		if (exonStartIndex.find(startKey) == exonStartIndex.end())
			exonStartIndex[startKey] = std::set<KeyType>{geneId};
		else
			exonStartIndex[startKey].insert(geneId);
		if (exonEndIndex.find(endKey) == exonEndIndex.end())
			exonEndIndex[endKey] = std::set<KeyType>{geneId};
		else
			exonEndIndex[endKey].insert(geneId);
	}
}