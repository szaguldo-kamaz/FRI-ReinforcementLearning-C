# FRI-ReinforcementLearning-C
*Fuzzy Rule Interpolation-based Reinforcement Learning* (FRIRL - ANSI C + AVX port)

For details about FRIRL itself, please see the main [FRI-ReinforcementLearning repository](https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning/).

This C + AVX port (AVX2 with inline assembly) is a bit behind the MATLAB implementation regarding features, but the main engine is fully functional, and indeed it is significantly faster.

About the AVX optimization parts, you can find more information in this paper:

* D. Vincze, "Parallelization by Vectorization in Fuzzy Rule Interpolation Adapted to FRIQ-Learning" in Proc. of the 2018 World Symposium on Digital Intelligence for Systems and Machines (DISA), 2018, pp. 131-136, [doi: 10.1109/DISA.2018.8490614](https://doi.org/10.1109%2FDISA.2018.8490614).
