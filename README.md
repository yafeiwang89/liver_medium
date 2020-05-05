## Source code of nanoHUB GUI-liver tissue mechanobiology  
This repository contains code and data for the nanoHUB app: https://nanohub.org/tools/pc4livermedium 

### Citing
If you use the code in your work, please cite [Impact of tumor-parenchyma biomechanics on liver metastatic progression: a multi-model approach](https://www.biorxiv.org/content/10.1101/2020.05.04.074989v1).
```text
@article {Wang2020.05.04.074989,
	author = {Wang, Yafei and Brodin, Erik and Nishii, Kenichiro and Frieboes, Hermann B and Mumenthaler, Shannon and Sparks, Jessica L and Macklin, Paul},
	title = {Impact of tumor-parenchyma biomechanics on liver metastatic progression: a multi-model approach},
	elocation-id = {2020.05.04.074989},
	year = {2020},
	doi = {10.1101/2020.05.04.074989},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2020/05/05/2020.05.04.074989},
	eprint = {https://www.biorxiv.org/content/early/2020/05/05/2020.05.04.074989.full.pdf},
	journal = {bioRxiv}
}

```

* * * 

## Release summary:

### 1.3:
* This version is actually same as 2.0. It's released in order to match with nanoHUB gui version

### 2.0: 
* Add an animate tap
* Fix bugs in `src/custom_modules/liver_functions.cpp`
```c++
void HCT116_phenotype_update_function( Cell* pCell , Phenotype& phenotype, double dt )
{
    static int start_phase_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_negative ); // 0
    static int end_phase_index = Ki67_advanced.find_phase_index( PhysiCell_constants::Ki67_positive_premitotic ); // 1
}
```

### 1.0:
This is initial release. 
