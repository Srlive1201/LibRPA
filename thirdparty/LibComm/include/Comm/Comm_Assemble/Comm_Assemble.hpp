//=======================
// AUTHOR : Peize Lin
// DATE :   2022-07-06
//=======================

#pragma once

#include "Comm_Assemble.h"

namespace Comm
{

template<typename Tkey, typename Tvalue, typename Tdatas_provide, typename Tkeys_require, typename Tdatas_require>
Comm_Assemble<Tkey,Tvalue,Tdatas_provide,Tkeys_require,Tdatas_require>::Comm_Assemble(const MPI_Comm &mpi_comm_in)
	:traverse_keys_provide(comm_keys.traverse_keys_provide),
	 set_value_require(comm_trans.set_value_recv),
	 flag_lock_set_value(comm_trans.flag_lock_set_value),
	 init_datas_local(comm_trans.init_datas_local),
	 add_datas(comm_trans.add_datas),
	 comm_keys(mpi_comm_in),
	 comm_trans(mpi_comm_in){}

template<typename Tkey, typename Tvalue, typename Tdatas_provide, typename Tkeys_require, typename Tdatas_require>
void Comm_Assemble<Tkey,Tvalue,Tdatas_provide,Tkeys_require,Tdatas_require>::communicate(
	const Tdatas_provide &datas_provide,
	const Tkeys_require &keys_require,
	Tdatas_require &datas_require)
{
	const std::vector<std::vector<Tkey>> keys_trans = comm_keys.trans( datas_provide, keys_require );
	comm_trans.traverse_isend = [&](
		const Tdatas_provide &datas_provide,
		const int rank_isend,
		std::function<void(const Tkey&, const Tvalue&)> &func)
	{
		for(const Tkey &key : keys_trans[rank_isend])
			func(key, this->get_value_provide(key, datas_provide));
	};
	comm_trans.communicate(datas_provide, datas_require);
}

}