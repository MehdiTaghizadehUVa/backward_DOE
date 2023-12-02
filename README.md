# Improving Accuracy and Computational Efficiency of Optimal Design of Experiments via Greedy Backward Approach

## Overview
This repository contains the implementation of the algorithms and methodologies presented in the paper "Improving Accuracy and Computational Efficiency of Optimal Design of Experiments via Greedy Backward Approach," published in the International Journal for Uncertainty Quantification in June 2024.

## Paper Abstract
Our research introduces a novel non-intrusive polynomial chaos expansion (PCE) technique, focusing on improving the efficiency and accuracy of surrogate model construction for uncertainty quantification. The proposed backward greedy algorithm deviates from traditional sample point selection strategies in Design of Experiments (DoEs), utilizing a coherence-optimal sample pool and selectively removing less influential samples. This method shows significant improvements in computational efficiency and robustness of DoEs, with various numerical examples demonstrating its efficacy.

## Repository Content

- `Algorithms/`: Implementation of the backward greedy algorithm as outlined in the paper.
- `Examples/`: Sample scripts and datasets used for the numerical experiments in the paper.
- `Documentation/`: Additional notes, algorithm descriptions, and usage guides.
- `LICENSE`: Details of the licensing and usage rights.

## Key Features

- Implementation of the backward greedy algorithm for optimal DoEs.
- Techniques for minimizing the curse of dimensionality in PCE.
- Coherence-optimal sampling strategies.
- Alphabetic optimality criteria (A-, D-, and S-optimality) implementations.
- Comprehensive examples with different complexity and dimensionality.

## Citation

If you use this code or our methodology in your research, please cite our paper. The paper can be accessed [here](https://www.researchgate.net/publication/371340182_IMPROVING_ACCURACY_AND_COMPUTATIONAL_EFFICIENCY_OF_OPTIMAL_DESIGN_OF_EXPERIMENTS_VIA_GREEDY_BACKWARD_APPROACH):

Taghizadeh, M., Xiu, D., & Alemazkoor, N. (2024). Improving Accuracy and Computational Efficiency of Optimal Design of Experiments via Greedy Backward Approach. International Journal for Uncertainty Quantification, 4(6), 204-220. DOI: 10.1615/Int.J.UncertaintyQuantification.2023046204


## Contributors

- Mehdi Taghizadeh (University of Virginia)
- Dongbin Xiu (The Ohio State University)
- Negin Alemazkoor (University of Virginia)
