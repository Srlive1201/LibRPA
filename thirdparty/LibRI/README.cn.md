# LibRI

其他语言版本：[English](./README.md)，[中文](./README.cn.md)


## RI 方法

常规密度泛函理论（DFT）方法不足以正确描述某些相互作用。
通常需要使用更高层次的方法，如杂化密度泛函、无规相近似（RPA）、GW方法等，以进行精确计算。
然而，这些方法所需计算资源极大。

对于局域基组，resolution-of-the-identity（RI）显著降低了内存需求与计算成本。


## LibRI 计算需求

LibRI 为只包含头文件的 C++ 库，用以计算 RI 形式下的高阶方法。

以下为 LibRI 所需：

- C++编译器，需支持 C++14 标准，且需支持 OpenMP 线程并行。
- MPI 库，用于进程间数据通讯。
- BLAS 与 LAPACK 库，用于加速张量运算。
  > 若 BLAS 与 LAPACK 库使用 Math Kernel Library (MKL)，则建议在 include 任意 LibRI 头文件前定义宏 `__MKL_RI`。LibRI 中部分函数将在编译时自动替换为 MKL 中的函数。
- [cereal](https://uscilab.github.io/cereal/) 库，用于数据序列化与反序列化，为纯头文件库。
- [LibComm](https://github.com/abacusmodeling/LibComm.git)  库，用于进程间数据传输，为纯头文件库。
