#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include "circCall.hpp"

namespace circfull
{

	std::string referenceStorage::getSeqByPosition(std::string chr, int start, int end)
	{
		// check if seqs[chr] position: start : end is not empty
		if (seqs.find(chr) == seqs.end())
		{
			std::cout << "Error: chromosome " << chr << " not found in reference genome" << std::endl;
			return "";
		}
		if (seqs[chr].size() < end)
			return "NA";
		return seqan3ToString(seqs[chr] | seqan3::views::slice(start, end));
	}

	std::istream &operator>>(std::istream &is, GtfRecord &record)
	{
		std::string tmp, tmpLine;
		is >> record.chr >> tmpLine;
		is >> tmp; // read in feature type
		if (tmp == "gene")
			record.type = GtfRecord::FeatureType::gene;
		else if (tmp == "transcript")
			record.type = GtfRecord::FeatureType::transcript;
		else if (tmp == "exon")
			record.type = GtfRecord::FeatureType::exon;
		else if (tmp == "CDS")
			record.type = GtfRecord::FeatureType::CDS;
		else if (tmp == "UTR")
			record.type = GtfRecord::FeatureType::UTR;
		else if (tmp == "start_codon")
			record.type = GtfRecord::FeatureType::start_codon;
		else if (tmp == "stop_codon")
			record.type = GtfRecord::FeatureType::stop_codon;
		else
			record.type = GtfRecord::FeatureType::unknown;
		is >> record.start >> record.end >> tmpLine;
		is >> tmp; // read in strand
		record.strand = tmp[0];
		is >> tmpLine;
		std::getline(is, tmp); // get additional informations splited by ';', format: key "value";
		std::istringstream iss{tmp};
		std::string informationLine, valueLine;
		while (std::getline(iss, informationLine, ';'))
		{
			while (informationLine[0] == ' ' || informationLine[0] == '\t')
				informationLine.erase(0, 1);
			std::istringstream iss2{informationLine};
			std::getline(iss2, valueLine, ' ');
			std::getline(iss2, tmpLine, '"');
			if (valueLine == "gene_id")
				std::getline(iss2, record.geneId, '"');
			else if (valueLine == "transcript_id")
				std::getline(iss2, record.transcriptId, '"');
			else if (valueLine == "gene_name")
				std::getline(iss2, record.geneName, '"');
		}
		return is;
	}

	referenceStorage readReference(std::filesystem::path const &referencePath)
	{
		referenceStorage storage{};
		seqan3::sequence_file_input reference_in{referencePath};
		for (auto &&record : reference_in)
		{
			storage.seqs[getId(record.id())] = record.sequence();
		}
		return storage;
	}

	GtfStorage readGtf(std::filesystem::path const &gtfPath)
	{
		GtfStorage storage{};
		std::ifstream gtfFile{gtfPath};
		std::string line;
		GtfRecord record;
		while (std::getline(gtfFile, line))
		{
			if (line[0] == '#')
				continue;
			std::istringstream iss{line};
			iss >> record;
			// if (record.type == GtfRecord::FeatureType::exon)
			storage.insert({record.geneId, record});
		}
		return storage;
	}

	namespace RG
	{

		SpliceMotif::SpliceMotif(const char *seq)
		{
			if (strlen(seq) != 2)
				throw std::invalid_argument("SpliceMotif: seq length must be 2");
			first = seq[0];
			second = seq[1];
		}

		SpliceMotif::SpliceMotif(const std::string &seq)
		{
			if (seq.length() != 2)
				throw std::invalid_argument("SpliceMotif: seq length must be 2");
			first = seq[0];
			second = seq[1];
		}

		SpliceMotif::SpliceMotif(const std::string &&seq)
		{
			if (seq.length() != 2)
				throw std::invalid_argument("SpliceMotif: seq length must be 2");
			first = seq[0];
			second = seq[1];
		}
		std::string SpliceMotif::str() const
		{
			return std::string(1, first) + std::string(1, second);
		}

		/*BSRecord::BSRecord() : adjust(false), leftSeqMapStatus(false), rightSeqMapStatus(false), queryId(""), chr(""), BSMotif(), BSstrand(), BSPos(){
			// reserve enough space for vectors
			BSMotif.reserve(10);
			BSstrand.reserve(10);
			BSPos.reserve(10);
			// reserve enough space for string
			queryId.reserve(40);
			chr.reserve(10);
		}

		BSRecord::BSRecord(const BSRecord &b) : BSRecord::BSRecord() {
			adjust = b.adjust;
			leftSeqMapStatus = b.leftSeqMapStatus;
			rightSeqMapStatus = b.rightSeqMapStatus;
			BSMotif = b.BSMotif;
			BSstrand = b.BSstrand;
			BSPos = b.BSPos;
			queryId = b.queryId;
			chr = b.chr;
		}

		BSRecord::BSRecord(BSRecord &&b) : BSRecord::BSRecord() {
			adjust = b.adjust;
			leftSeqMapStatus = b.leftSeqMapStatus;
			rightSeqMapStatus = b.rightSeqMapStatus;
			BSMotif = std::move(b.BSMotif);
			BSstrand = std::move(b.BSstrand);
			BSPos = std::move(b.BSPos);
			queryId = std::move(b.queryId);
			chr = std::move(b.chr);
		}

		BSRecord & BSRecord::operator=(const BSRecord &b) {
			adjust = b.adjust;
			leftSeqMapStatus = b.leftSeqMapStatus;
			rightSeqMapStatus = b.rightSeqMapStatus;
			BSMotif = b.BSMotif;
			BSstrand = b.BSstrand;
			BSPos = b.BSPos;
			queryId = b.queryId;
			chr = b.chr;
			return *this;
		}

		BSRecord & BSRecord::operator=(BSRecord &&b) {
			adjust = b.adjust;
			leftSeqMapStatus = b.leftSeqMapStatus;
			rightSeqMapStatus = b.rightSeqMapStatus;
			BSMotif = std::move(b.BSMotif);
			BSstrand = std::move(b.BSstrand);
			BSPos = std::move(b.BSPos);
			queryId = std::move(b.queryId);
			chr = std::move(b.chr);
			return *this;
		}*/

		bool MapRecord::operator<(const MapRecord &b) const
		{
			if (queryId == b.queryId)
				return queryStart < b.queryStart;
			return queryId < b.queryId;
		}

		std::istream &operator>>(std::istream &in, MapRecord &record)
		{
			std::string tmp;
			in >> record.queryId >>
				record.chr >> record.refStart >> record.refEnd >>
				record.queryStart >> record.queryEnd >>
				record.exonStart >> record.exonEnd >> record.exonLength >>
				tmp >> record.leftSeq >> record.rightSeq;
			record.strand = tmp[0];
			return in;
		}

		bool hitBSMotif(const std::string &seq, const BSseq &motif, const char &strand)
		{
			if (strand == '+')
				return seq == motif.seq;
			if (strand == '-')
				return seq == motif.revSeq;
			return (seq == motif.seq || seq == motif.revSeq);
		}

		std::string to_string(const MapRecord *record)
		{
			std::string out;
			if (record == NULL)
				out = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
			else
				out = record->queryId + "\t" + record->chr + "\t" +
					  std::to_string(record->refStart) + "\t" + std::to_string(record->refEnd) + "\t" +
					  std::to_string(record->queryStart) + "\t" + std::to_string(record->queryEnd) + "\t" +
					  record->exonStart + "\t" + record->exonEnd + "\t" + std::to_string(record->exonLength) + "\t" +
					  record->strand + "\t" + record->leftSeq + "\t" + record->rightSeq + "\n";
			return out;
		}

		std::string to_string(const BSRecord *record)
		{
			std::string out;
			if (record == NULL)
				return "";
			else
			{
				std::string BSLeft = std::to_string(record->BSPos.begin()->first);
				std::string BSRight = std::to_string(record->BSPos.begin()->second);
				for (auto iter = record->BSPos.begin() + 1; iter != record->BSPos.end(); iter++)
				{
					BSLeft += "," + std::to_string(iter->first);
					BSRight += "," + std::to_string(iter->second);
				}
				std::string BSMotifLeft = record->BSMotif.begin()->first.str();
				std::string BSMotifRight = record->BSMotif.begin()->second.str();
				for (auto iter = record->BSMotif.begin() + 1; iter != record->BSMotif.end(); iter++)
				{
					BSMotifLeft += "," + iter->first.str();
					BSMotifRight += "," + iter->second.str();
				}
				std::string BSstrand = std::string(1, *(record->BSstrand.begin()));
				for (auto iter = record->BSstrand.begin() + 1; iter != record->BSstrand.end(); iter++)
					BSstrand += "," + std::string(1, *iter);
				out = record->queryId + "\t" + record->chr + "\t" + BSLeft + "\t" + BSRight + "\t" + BSstrand + "\t" + BSMotifLeft + "\t" + BSMotifRight;
			}
			if (record->adjust)
				out += std::string("\t") + std::string(record->leftSeqMapStatus ? "true" : "false") + "\t" + std::string(record->rightSeqMapStatus ? "true" : "false");
			out += "\n";
			return out;
		}

		std::string to_string(const CircRecord *record)
		{
			std::string out;
			if (record == NULL)
				return "";
			else
			{
				if (record->queryId != "")
					out += record->queryId + "\t";
				out += record->chr + ":" + std::to_string(record->start) + "|" + std::to_string(record->end) +
					   "\t" + record->strand + "\t";
				std::string exonStart = std::to_string(record->exons.begin()->first);
				std::string exonEnd = std::to_string(record->exons.begin()->second);
				for (auto iter = record->exons.begin() + 1; iter != record->exons.end(); iter++)
				{
					exonStart += "," + std::to_string(iter->first);
					exonEnd += "," + std::to_string(iter->second);
				}
				out += exonStart + "\t" + exonEnd + "\t" + std::to_string(record->len) + "\t";
				/*exonStart = record->spliceSeqs.begin()->first;
				exonEnd = record->spliceSeqs.begin()->second;
				for (auto iter = record->spliceSeqs.begin() + 1; iter != record->spliceSeqs.end(); iter++) {
					exonStart += "," + iter->first;
					exonEnd += "," + iter->second;
				}
				out += exonStart + "\t" + exonEnd + "\t";*/
				out += "\n";
			}
			return out;
		}

		ExonIndex<std::string> getExonIndex(const GtfStorage &gtfList, const int errorLen)
		{
			// filter exons
			ExonIndex exonIndex;
			for (auto &[geneId, record] : gtfList)
				if (record.type == GtfRecord::FeatureType::exon)
					exonIndex.addExon(record.chr, record.start, record.end, record.strand, record.geneId, errorLen);
			return exonIndex;
		}

		void BSMotifIndex::addIndex(std::string _seq)
		{
			std::string leftSeq = _seq.substr(0, 2);
			std::string rightSeq = _seq.substr(2, 2);
			std::string revLeftSeq = seqan3ToString(rightSeq | seqan3::views::char_to<seqan3::dna5> | seqan3::views::complement | std::views::reverse);
			std::string revRightSeq = seqan3ToString(leftSeq | seqan3::views::char_to<seqan3::dna5> | seqan3::views::complement | std::views::reverse);
			leftMotif.push_back({leftSeq, revLeftSeq});
			rightMotif.push_back({rightSeq, revRightSeq});
			bothMotif.push_back({leftSeq + rightSeq, revLeftSeq + revRightSeq});
		}

		bool CircID::operator<(const CircID &b) const
		{
			return chr < b.chr || (chr == b.chr && start < b.start) || (chr == b.chr && start == b.start && end < b.end);
		}

		bool CircID::operator==(const CircID &b) const
		{
			return chr == b.chr && start == b.start && end == b.end && strand == b.strand;
		}

		bool CircID::operator==(const std::tuple<std::string, char, int, int> &b) const
		{
			return chr == std::get<0>(b) && strand == std::get<1>(b) && start == std::get<2>(b) && end == std::get<3>(b);
		}

		void CircRecord::updateLength()
		{
			len = 0;
			for (auto &exon : exons)
				len += exon.second - exon.first + 1;
		}

		std::string CircRecord::getCircId()
		{
			return chr + ":" + std::to_string(start) + "-" + std::to_string(end) + ":" + strand;
		}

		std::string CircRecord::getExonsString()
		{
			std::string ret = std::to_string(exons.begin()->first) + "-" + std::to_string(exons.begin()->second);
			for (auto iter = exons.begin() + 1; iter != exons.end(); iter++)
				ret += "," + std::to_string(iter->first) + "-" + std::to_string(iter->second);
			return ret;
		}

		std::string CircRecord::getTranscriptId()
		{
			return getCircId() + "|" + getExonsString();
		}
	}

	std::optional<std::set<GtfRecordPtr>> GtfFetcher::fetch(std::string chr, int start, int end)
	{
		std::vector<GtfRecord> ret;
		if (fetchers.find(chr) == fetchers.end())
			return std::nullopt;
		return fetchers[chr].fetch(start, end);
	}

	std::optional<std::set<GtfRecordPtr>> GtfFetcherChrom::fetch(int start, int end)
	{
		auto [lowerIter, upperIter]{fetcher.equal_range({start, end})};
		if (lowerIter == upperIter)
			return std::nullopt;
		std::set<GtfRecordPtr> ret;
		for (auto iter = lowerIter; iter != upperIter; iter++)
		{
			ret.insert(iter->second.begin(), iter->second.end());
		}
		return ret;
	}
}