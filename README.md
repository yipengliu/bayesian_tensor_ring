# demo for Beyesian tensor ring

The codes are for the paper:

Zhen Long, Ce Zhu, Jiani Liu, Yipeng Liu, "Bayesian Low Rank Tensor Ring for Image Recovery,"  IEEE Transactions on Image Processing, 2021. DOI: 10.1109/TIP.2021.3062195.

For experiments on color image completion: 
demo_image_recovery.m 

For setting up, we can add path by the following sentence:
addpath('mylib');  addpath('TestImage');


Paper abstract:

Low rank tensor ring based data recovery can recover missing image entries in signal acquisition and transformation. The recently proposed tensor ring (TR) based completion algorithms generally solve the low rank optimization problem by alternating least squares method with predefined ranks, which may easily lead to overfitting when the unknown ranks are set too large and only a few measurements are available. In this paper, we present a Bayesian low rank tensor ring completion method for image recovery by automatically learning the low-rank structure of data. A multiplicative interaction model is developed for low rank tensor ring approximation, where sparsity-inducing hierarchical prior is placed over horizontal and frontal slices of core factors. Compared with most of the existing methods, the proposed one is free of parameter-tuning, and the TR ranks can be obtained by Bayesian inference. Numerical experiments, including synthetic data, real-world color images and YaleFace dataset, show that the proposed method outperforms state-of-the-art ones, especially in terms of recovery accuracy.
