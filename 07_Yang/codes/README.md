<div align="center">

# DeepHAM: A global solution method for heterogeneous agent models with aggregate shocks

Jiequn Han, Yucheng Yang, Weinan E

[![arXiv](https://img.shields.io/badge/arXiv-2112.14377-b31b1b.svg)](https://arxiv.org/abs/2112.14377)
[![SSRN](https://img.shields.io/badge/SSRN-3990409-133a6f.svg)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3990409)
[![PDF](https://img.shields.io/badge/PDF-8A2BE2)](https://yangycpku.github.io/files/DeepHAM_paper.pdf)

Link to code repository: https://github.com/frankhan91/DeepHAM

</div>


## Dependencies
* Quick installation of conda environment for Python: ``conda env create -f environment.yml``

## Running
### Quick start for the Krusell-Smith (KS) model under default configs:
To use DeepHAM to solve the competitive equilibrium of the KS model, run
```
python train_KS.py
```
To evaluate the Bellman error of the solution of the KS model, run
```
python validate_KS.py
```

Sample scripts for solving the KS model in the Slurm system are provided in the folder ``src/slurm_scripts``

### Solve the model in Fernandez-Villaverde, Hurtado, and Nuno (2019):
```
python train_JFV.py
```
```
python validate_JFV.py
```
Details on the model setup and algorithm can be found in our paper.

## Citation
If you find this work helpful, please consider starring this repo and citing our paper using the following Bibtex.
```bibtex
@article{HanYangE2021deepham,
  title={Deep{HAM}: A global solution method for heterogeneous agent models with aggregate shocks},
  author={Han, Jiequn and Yang, Yucheng and E, Weinan},
  journal={arXiv preprint arXiv:2112.14377},
  year={2021}
}
```

## Contact
Please contact us at jiequnhan@gmail.com and yucheng.yang@uzh.ch if you have any questions.
