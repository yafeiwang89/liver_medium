## Source code of nanoHUB GUI-liver tissue mechanobiology  
This repository contains code and data for the nanoHUB app: https://nanohub.org/tools/pc4livermedium 

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
