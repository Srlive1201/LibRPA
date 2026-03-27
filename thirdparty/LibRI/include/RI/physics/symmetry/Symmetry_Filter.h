#pragma once

#include "../../global/Array_Operator.h"

#include <array>
#include <map>
#include <set>

#define NO_SEC_RETURN_TRUE if(this->irreducible_sector_.empty()) return true;

namespace RI
{
	template<typename TA, typename TC, typename Tdata>
	class Symmetry_Filter
	{
		using TAC = std::pair<TA, TC>;
		using Tab = std::pair<TA, TA>;
		using TabR = std::pair<Tab, TC>;
		using Tsec = std::map<Tab, std::set<TC>>;

	  public:
		Symmetry_Filter(const TC& period_in, const Tsec& irsec)
			:period(period_in), irreducible_sector_(irsec) {}
		bool in_irreducible_sector(const TA& Aa, const TAC& Ab) const
		{
			NO_SEC_RETURN_TRUE;
			using namespace Array_Operator;
			const Tab& ap = { Aa, Ab.first };
			const auto ptr = this->irreducible_sector_.find(ap);
			if ( ptr!= this->irreducible_sector_.end())
				if (ptr->second.find(Ab.second % this->period) != ptr->second.end())
					return true;
			return false;
		}
		bool in_irreducible_sector(const TAC& Aa, const TAC& Ab) const
		{
			NO_SEC_RETURN_TRUE;
			using namespace Array_Operator;
			const TC dR = (Ab.second - Aa.second) % this->period;
			const std::pair<TA, TA> ap = { Aa.first, Ab.first };
			const auto ptr = this->irreducible_sector_.find(ap);
			if (ptr != this->irreducible_sector_.end())
				if (ptr->second.find(dR) != ptr->second.end())
					return true;
			return false;
		}
		bool is_Aa_in_irreducible_sector(const TA& Aa) const
		{
			NO_SEC_RETURN_TRUE;
			for (const auto& apRs : this->irreducible_sector_)
				if (apRs.first.first == Aa)
					return true;
			return false;
		}
		bool is_Ab_in_irreducible_sector(const TA& Ab) const
		{
			NO_SEC_RETURN_TRUE;
			for (const auto& apRs : this->irreducible_sector_)
				if (apRs.first.second == Ab)
					return true;
			return false;
		}
		TabR get_abR(const TA& Aa, const TAC& Ab) const
		{
			return { {Aa,Ab.first}, Ab.second % this->period };
		}
		TabR get_abR(const TAC& Aa, const TAC& Ab) const
		{
			return { {Aa.first,Ab.first}, (Ab.second - Aa.second) % this->period };
		}

	  public:	// private:
		const TC& period;
		const Tsec& irreducible_sector_;
	};

}

#undef NO_SEC_RETURN_TRUE