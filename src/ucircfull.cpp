#include "circfullUtils.hpp"
#include "ucircfull.hpp"

int circfull_pipeline(circfull::CircfullOption &opt)
{
	return 0;
};

int main(int argc, char **argv)
{
	circfull::CircfullOption opt;
	circfull::parseArgs(opt, argc, argv);
	std::cout << "              _           ____      ____" << std::endl
			  << "  __  _______(_)_________/ __/_  __/ / /" << std::endl
			  << " / / / / ___/ / ___/ ___/ /_/ / / / / / " << std::endl
			  << "/ /_/ / /__/ / /  / /__/ __/ /_/ / / /  " << std::endl
			  << "\\__,_/\\___/_/_/   \\___/_/  \\__,_/_/_/   " << std::endl
			  << std::endl;
	switch (opt.call)
	{
	case circfull::CircfullOption::mod::UMI_EXTRACT:
		return circfull::circfull_umi_extract(opt);
		break;
	case circfull::CircfullOption::mod::UMI_CLUST:
		return circfull::circfull_umi_clust(opt);
		break;
	case circfull::CircfullOption::mod::CIRC_CALL:
		return circfull::circ_call_cpp(opt);
		break;
	case circfull::CircfullOption::mod::CIRC_FILT:
		return circfull::filterCircCand(opt);
		break;
	}
	return 0;
}