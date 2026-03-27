//=======================
// AUTHOR : Peize Lin
// DATE :   2022-07-06
//=======================

#pragma once

#include "../Comm_Keys/Comm_Keys_31-gather.h"
#include "../Comm_Trans/Comm_Trans.h"
#include "../Comm_Tools.h"

#include <mpi.h>
#include <functional>

namespace Comm
{

template<typename Tkey, typename Tvalue, typename Tdatas_provide, typename Tkeys_require, typename Tdatas_require>
class Comm_Assemble
{
public:
	Comm_Assemble(const MPI_Comm &mpi_comm_in);

	std::function<
		void(
			const Tdatas_provide &keys_provide_mine,
			std::function<void(const Tkey&)> &func)>
		&traverse_keys_provide;
	std::function<
		const Tvalue&(
			const Tkey &key,
			const Tdatas_provide &datas_provide)>
		get_value_provide;
	std::function<
		void(
			Tkey &&key,
			Tvalue &&value,
			Tdatas_require &datas_require)>
		&set_value_require;

	Comm_Tools::Lock_Type &flag_lock_set_value;
	std::function<
		Tdatas_require(
			const int rank_recv)>
		&init_datas_local;
	std::function<
		void(
			Tdatas_require &&datas_local,
			Tdatas_require &datas_recv)>
		&add_datas;

	void communicate(
		const Tdatas_provide &datas_provide,
		const Tkeys_require &keys_require,
		Tdatas_require &datas_require);

private:
	Comm_Keys_31_SenderTraversal<Tkey,Tdatas_provide,Tkeys_require> comm_keys;
	Comm_Trans<Tkey,Tvalue,Tdatas_provide,Tdatas_require> comm_trans;
};

}

#include "Comm_Assemble.hpp"