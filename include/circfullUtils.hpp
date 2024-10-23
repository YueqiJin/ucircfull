#pragma once

#include <filesystem>
#include <string>
#include <set>
#include <map>
#include <seqan3/io/sequence_file/all.hpp>
#include "degenerateBase.hpp"

namespace circfull
{

	struct literal
	{
		const char *str;
		literal(const char *in) : str(in) {}
	};

	/// @brief input stream operator for literal
	/// @param is input stream
	/// @param lit literal
	/// @return input stream
	std::istream &operator>>(std::istream &is, const literal &lit);

	/// @brief get integer from comma separated string
	/// @param str input string
	/// @return array of integers
	std::vector<int> getIntFromCommaSepStr(const std::string &str);

	class CircfullOption
	{
	public:
		// circfull2 run info
		typedef enum
		{
			UMI_EXTRACT,
			UMI_CLUST,
			CIRC_CALL,
			CIRC_FILT,
			CIRC_QUANT,
			CIRCFULL
		} mod;
		mod call;
		size_t nthread;
		bool noUMI;
		bool quiet;
		bool help;
		bool clustSeq;

		// circfull2 software info
		std::string seqkit;
		std::string ciriLong;
		std::string porechop;
		std::string minimap2;
		std::string samtools;

		// circfull2 input info
		std::filesystem::path rawFastq;
		std::filesystem::path referenceFasta;
		std::filesystem::path annotationGTF;

		// circfull2 output info
		std::filesystem::path oPath;
		std::string oPrefix;

		// circFL library info
		DegenerateBaseVector umiSeq;
		std::string anchorx;
		int umiLength;

		// circRNA calling and filtering options
		std::string spliceSignal;
		int circNumThres;

		// circfull2 processes files
		std::filesystem::path umiFasta;
		std::filesystem::path strandFastq;
		bool useBam;
		std::filesystem::path bamFile;
		std::filesystem::path strandNotTrimFastq;
		std::filesystem::path umiClustTab;
		std::filesystem::path consensusFasta;
		std::filesystem::path circCallOutputPath;
		std::filesystem::path circCollapsePath;
		std::filesystem::path rawCircRes;
		std::filesystem::path ccsFasta;
		std::filesystem::path consCircRes;
		std::filesystem::path circCandFile;

		typedef enum
		{
			RG,
			cRG,
			ucRG,
			CIRI
		} CircCallMod;
		CircCallMod modCircCall;
		bool strand;
	};

	/// @brief check if file exist
	/// @param file path
	/// @return true if file existed
	bool checkFile(const std::filesystem::path &file);

	/// @brief return number of rows of the given file
	/// @param filename path to file
	/// @return number of rows
	int getFileNrow(std::filesystem::path file, size_t nthread = 1);

	/// @brief make a directory by force
	/// @param targetDir path to directory, old folder will be deleted if existed
	/// @return 0 when operation success
	int makeDir(const std::filesystem::path &targetDir);

	/// @brief delete a directory
	/// @param targetDir path to directory
	/// @return 0 when operation success
	int delDir(std::filesystem::path targetDir);

	/// @brief check seqkit software existed
	/// @param seqkit path to seqkit
	/// @return 0 when exit succeed
	int checkSeqkit(const std::string &seqkit);

	/// @brief split a fastq file to nSplit files by seqkit
	/// @param seqkit path to seqkit software
	/// @param fastqFile path to fastq
	/// @param outputDir output directory
	/// @param outputPrefix prefix of output files
	/// @param nSplit number of files spliting into
	/// @return 0 when operation success
	int splitFastq(const std::string &seqkit, const std::filesystem::path &fastqFile, const std::filesystem::path &outputDir, const std::string &outputPrefix, int nSplit);
	/*
	/// @brief split a fastq file to nSplit files by seqkit
	/// @param seqkit path to seqkit software
	/// @param fastqFile path to fastq
	/// @param fastaFile path to output
	/// @return 0 when operation success
	int convertFastq2Fasta(std::string seqkit, std::filesystem::path fastqFile, std::filesystem::path fastaFile, int nSplit);
	*/
	/// @brief get first string with delim " " from input
	/// @param text input
	/// @param delim delim char
	/// @return first string
	std::string getFirstString(const std::string &text, const char &delim = ' ');

	/// @brief output a sequence record as fasta format
	/// @param id sequence id
	/// @param seq sequence
	/// @return sequence in fasta format
	std::string outputAsFasta(const std::string &id, const std::string &seq);

	/// @brief conversion vector in seqan3 library to std::string
	/// @tparam T vector type in seqan3 library
	/// @param value vector in seqan3 library
	/// @return string
	template <typename T>
	std::string seqan3ToString(T &&value);

	/// @brief get sequence id from the id line of fastq
	/// @param idString id line of fastq
	/// @return sequence id
	std::string getId(const std::string &idString);

	/// @brief sort seqs by seq_id using seqkit
	/// @param seqkit path to seqkit
	/// @param inputFasta path to input fastq
	/// @param outputFasta path to output fastq
	/// @param nthread number of threads
	/// @return 0 when operation success
	int sortFastqById(const std::string &seqkit, const std::filesystem::path &inputFastq, const std::filesystem::path &outputFastq, size_t nthread);

	/// @brief trim adaptors from fastq string using porechop
	/// @param seqkit path to seqkit
	/// @param inputFasta path to input fastq
	/// @param outputFasta path to output fastq
	/// @param nthread number of threads
	/// @return 0 when operation success
	int trimAdaptor(const std::string &porechop, const std::filesystem::path &fastqFile, const std::filesystem::path &outFile, size_t nthread);

	/// @brief merge multiple files into one file
	/// @param files array of files to be merged
	/// @param outputFile path of output file
	/// @param nFiles number of files
	/// @return 0 when success
	int mergeFile(std::filesystem::path *files, const std::filesystem::path &outputFile, int nFiles);

	/// @brief delete files
	/// @param files path of files
	/// @param nFiles number of files to be removed
	/// @return 0 when success
	int delFile(std::filesystem::path *files, int nFiles = 1);

	/// @brief print time and information
	/// @param text information
	/// @param timeColor color of time
	/// @param textColor color of information
	void printTimeInfo(const std::string &text);

	/// @brief get the length of array
	/// @tparam T
	/// @param array
	/// @return length of the array
	template <class T>
	int getArrayLength(T &array);

	/// @brief count frequencies of items in array
	/// @tparam T type of item in array
	/// @tparam HashType hash type of item in array
	/// @param array input array
	/// @param hashFunc function to hash item
	/// @return map of hashed item and its frequencies
	template <typename T, typename HashType = T>
	std::map<HashType, int> countFrequencies(
		const std::vector<T> &array,
		std::function<HashType(const T &)> hashFunc = [](const T &item)
		{ return static_cast<HashType>(item); });

	/// @brief get most frequent items in array
	/// @tparam T type of item in array
	/// @tparam HashType hash type of item in array
	/// @param array input array
	/// @param hashFunc function to hash item
	/// @return array of most frequent items
	template <typename T, typename HashType = T>
	std::vector<T> getMostFrequent(
		const std::vector<T> &array,
		std::function<HashType(const T &)> hashFunc = [](const T &item)
		{ return static_cast<HashType>(item); });

	/// @brief filter array by a function
	/// @tparam T type of item in array
	/// @param array input array
	/// @param filterFunc function to filter item
	/// @return array of filtered items
	template <typename T>
	std::vector<T> filter(
		const std::vector<T> &array,
		std::function<bool(const T &)> filterFunc);

	/// @brief filter array by a function and return a new array
	/// @tparam vectorType type of array
	/// @tparam vectorItemType type of item in array
	/// @tparam resultItem type of item in result array
	/// @param array input array
	/// @param filterFunc function to filter item
	/// @param resultFunc function to convert item to result item
	/// @return array of filtered items as result array
	template <typename vectorType, typename vectorItemType, typename resultItem>
	std::vector<resultItem> filter(
		const vectorType &array,
		std::function<bool(const vectorItemType &)> filterFunc,
		std::function<resultItem(const vectorItemType &)> resultFunc = [](const vectorItemType &item)
		{ return item; });

	/// @brief split array by a function
	/// @tparam T type of item in array
	/// @param array input array
	/// @param filterFunc function returning a bool value to split item
	/// @return tuple of two arrays, first array is true, second array is false
	template <typename T>
	std::tuple<std::vector<T>, std::vector<T>> splitByBool(
		const std::vector<T> &array,
		std::function<bool(const T &)> filterFunc);

	/// @brief mutate array by a function
	/// @tparam FromType type of item in input array
	/// @tparam ToType type of item in result array
	/// @param array input array
	/// @param mutateFunc function to mutate item
	/// @return array of mutated items
	template <typename FromType, typename ToType>
	std::vector<ToType> mutate(
		const std::vector<FromType> &array,
		std::function<ToType(const FromType &)> mutateFunc);

	/// @brief get max value of array by a function
	/// @tparam T type of item in array
	/// @tparam ToType type of value to compare
	/// @param array input array
	/// @param valueFunc function to get value of item
	/// @return max value
	template <typename T, typename ToType>
	std::vector<T> maxElements(const std::vector<T> &array, std::function<ToType(const T &)> valueFunc);

	/// @brief parse command line args to CircfullOption
	/// @param opt CircfullOption to save all args
	/// @param argc arg count from main function
	/// @param argv arg values from main function
	/// @return return value of program running
	void parseArgs(CircfullOption &opt, int argc, char **argv);

	/// @brief split string by delim
	/// @param str input string
	/// @param delim delim char
	/// @return array of integers
	std::vector<int> split(const std::string &str, const char &delim = ' ');

	extern bool debug;

};

template <typename T>
std::string circfull::seqan3ToString(T &&value)
{
	std::ostringstream oss;
	for (auto &&c : value)
	{
		oss << (seqan3::to_char(c));
	}
	return oss.str();
}

template <typename T, typename HashType>
std::map<HashType, int> circfull::countFrequencies(
	const std::vector<T> &array,
	std::function<HashType(const T &)> hashFunc)
{

	std::map<HashType, int> freq;
	for (const auto &item : array)
	{
		HashType hashValue = hashFunc(item);
		auto it = freq.find(hashValue);
		if (it != freq.end())
		{
			it->second++;
		}
		else
		{
			freq[hashValue] = 1;
		}
	}
	return freq;
}

template <typename T, typename HashType>
std::vector<T> circfull::getMostFrequent(
	const std::vector<T> &array,
	std::function<HashType(const T &)> hashFunc)
{

	std::map<HashType, int> freq = countFrequencies<T, HashType>(array, hashFunc);
	std::vector<T> mostFrequent;
	int maxFreq = 0;
	std::set<HashType> maxFreqItem;
	for (const auto &item : freq)
		if (item.second > maxFreq)
			maxFreq = item.second;
	for (const auto &item : freq)
		if (item.second == maxFreq)
			maxFreqItem.insert(item.first);
	for (const auto &item : array)
		if (maxFreqItem.find(hashFunc(item)) != maxFreqItem.end())
			mostFrequent.push_back(item);
	return mostFrequent;
}

template <typename T>
std::vector<T> circfull::filter(
	const std::vector<T> &array,
	std::function<bool(const T &)> filterFunc)
{

	std::vector<T> filtered;
	for (const auto &item : array)
		if (filterFunc(item))
			filtered.push_back(item);
	return filtered;
}

template <typename vectorType, typename vectorItemType, typename resultItem>
std::vector<resultItem> circfull::filter(
	const vectorType &array,
	std::function<bool(const vectorItemType &)> filterFunc,
	std::function<resultItem(const vectorItemType &)> resultFunc)
{

	std::vector<resultItem> filtered;
	for (const auto &item : array)
		if (filterFunc(item))
			filtered.push_back(resultFunc(item));
	return filtered;
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> circfull::splitByBool(
	const std::vector<T> &array,
	std::function<bool(const T &)> filterFunc)
{

	std::vector<T> res1, res2;
	res1.reserve(array.size());
	res2.reserve(array.size());
	for (const auto &item : array)
	{
		if (filterFunc(item))
			res1.push_back(item);
		else
			res2.push_back(item);
	}
	return {res1, res2};
}

template <typename FromType, typename ToType>
std::vector<ToType> mutate(
	const std::vector<FromType> &array,
	std::function<ToType(const FromType &)> mutateFunc)
{

	std::vector<ToType> ret;
	for (FromType &i : array)
		ret.push_back(mutateFunc(i));
	return ret;
}

template <typename T, typename ToType>
ToType max(const std::vector<T> &array, std::function<ToType(const T &)> valueFunc)
{
	ToType maxValue = valueFunc(array[0]);
	for (const T &item : array)
	{
		ToType value = valueFunc(item);
		if (value > maxValue)
			maxValue = value;
	}
	return maxValue;
}

template <typename T, typename ToType>
std::vector<T> circfull::maxElements(const std::vector<T> &array, std::function<ToType(const T &)> valueFunc)
{
	ToType maxValue = max<T, ToType>(array, valueFunc);
	std::vector<T> maxElements = filter<T>(array, [valueFunc, &maxValue](const T &item)
										   { return valueFunc(item) == maxValue; });
	return maxElements;
}