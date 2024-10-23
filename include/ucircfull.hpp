/************************************************************************************
 *  
 *               _           ____      ____
 *   __  _______(_)_________/ __/_  __/ / /   circfull2 - Identify circRNAs from
 *  / / / / ___/ / ___/ ___/ /_/ / / / / /                full-length circRNA
 * / /_/ / /__/ / /  / /__/ __/ /_/ / / /                 sequencing data
 * \__,_/\___/_/_/   \___/_/  \__,_/_/_/      https://github.com/yangence/circfull2
 *
 * Copyright (C) 2024 Yangence Lab at Peking University
 *
 ************************************************************************************/

#pragma once

#include "circfullUtils.hpp"

namespace circfull
{
	int circfull_umi_extract(CircfullOption &opt);
	int circfull_umi_clust(CircfullOption &opt);
	int circ_call(CircfullOption &opt);
	int circ_call_cpp(CircfullOption &opt);
	int circ_call_test(CircfullOption &opt);
	int filterCircCand(CircfullOption &opt);
}
