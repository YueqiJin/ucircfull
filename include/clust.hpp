#pragma once

#include <string>
#include <map>
#include <filesystem>

/// @brief clust UMI sequence with cdhit api
/// @param fastaFile input fasta file
/// @param oPrefix output prefix
/// @param wordSize word size
/// @param nTread number of thread
/// @param identity identity threshold
/// @param clustInfo cluster information (part of results)
/// @return number of cluster
int clustUMISeq(const std::filesystem::path &fastaFile, const std::filesystem::path &oPrefix, int wordSize, int nTread, double identity, std::map<std::string, int> &clustInfo);

/// @brief generate consensus sequence from multiple sequence
/// @param seqs sequences
/// @param quals quality of sequences
/// @return consensus sequence
std::string generateConsensusSeq(const std::vector<std::string> &seqs, const std::vector<std::string> &quals);

/// @brief generate consensus sequence from multiple sequence
/// @param seqs sequences
/// @return consensus sequence
std::string generateConsensusSeq(const std::vector<std::string> &seqs);
