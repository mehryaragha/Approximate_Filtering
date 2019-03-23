# Approximate Computing for an Extended Kalman Filtering Tracker

An implementation of online approximator for an extended Kalman filter tracker implemented over FPGA is presented.

## Getting started
Please run Demo_detailed_visualisation.m and Demo.m to simulate Figure 4 of the paper, which contains dynamic approximation via four different Kullback–Leibler (KL) thresholds.

Also, please run Demo_Exact_Computation.m for exact computation, i.e target tracking with an extended Kalman filter without any approximation.

## Paper

```
@inproceedings{Emambakhsh:2017,
  title={Learning to approximate computing at run-time},
  author={P. Garcia, M. Emambakhsh and A. Wallace},
  booktitle={IET 3rd International Conference on Intelligent Signal Processing (ISP)},
  pages={1--8},
  year={2017}
}
```
