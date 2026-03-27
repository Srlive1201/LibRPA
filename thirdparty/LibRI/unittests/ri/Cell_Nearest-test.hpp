// ===================
//  Author: Peize Lin
//  date: 2022.12.30
// ===================

#pragma once

#include "../print_stl.h"

#include "RI/ri/Cell_Nearest.h"
#include "../global/Tensor-test.h"

namespace Cell_Nearest_Test
{
	static void main()
	{
		using TA = std::string;
		using Tcell = int;
		constexpr int Ndim = 2;
		using Tpos = double;
		constexpr int Npos = 2;
		using Tatom_pos = std::array<Tpos,Npos>;

		std::map<TA,Tatom_pos> atoms_pos;
		std::array<Tatom_pos,Ndim> latvec;
		std::array<Tcell,Ndim> period;

		auto test = [&atoms_pos, &latvec, &period]()
		{
			RI::Cell_Nearest<TA, Tcell, Ndim, Tpos, Npos> stress;
			stress.init(atoms_pos, latvec, period);
			std::cout<<stress.cells_nearest_continuous<<std::endl;

			for(Tcell idim0=0; idim0<period[0]; ++idim0)
				for(Tcell idim1=0; idim1<period[1]; ++idim1)
					std::cout<<idim0<<"\t"<<idim1<<"\t\t"<<stress.get_cell_nearest_discrete("H", "C", {idim0,idim1})<<std::endl;
		};

		{
			atoms_pos["H"] = {0, 0};
			atoms_pos["C"] = {-8, 2};
			latvec[0]={10,0};
			latvec[1]={0,10};
			period = {2,1};
			test();
			/*
				C	C	0		0
					H	-0.8	0.2
				H	C	0.8		-0.2
					H	0		0

				0	0		0	0
				1	0		1	0
			*/
		}

		{
			atoms_pos["H"] = {0, 0};
			atoms_pos["C"] = {8, 2};
			latvec[0]={10,0};
			latvec[1]={0,10};
			period = {2,1};
			test();
			/*
				C	C	0		0
					H	0.8		0.2
				H	C	-0.8	-0.2
					H	0		0

				0	0		0	0
				1	0		-1	0
			*/
		}

		{
			atoms_pos["H"] = {0, 0};
			atoms_pos["C"] = {18, 2};
			latvec[0]={10,0};
			latvec[1]={0,10};
			period = {2,1};
			test();
			/*
				C	C	0		0
					H	1.8		0.2

				H	C	-1.8	-0.2
					H	0		0


				0	0		-2	0
				1	0		-1	0
			*/
		}
	}
}