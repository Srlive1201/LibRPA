// ===================
//  Author: Peize Lin
//  date: 2022.06.22
// ===================

#pragma once

#include <vector>

namespace Comm
{

namespace Communicate_Vector
{
	template<typename Tkey>
	void traverse_keys(
		const std::vector<Tkey> &keys,
		std::function<void(const Tkey&)> &func)
	{
		for(const Tkey &key : keys)
			func(key);
	}
}

}