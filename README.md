# Delineating urban areas using building density (21.07.08)

**Reference**: De Bellefon M.P., Combes P.P., Duranton G., Gobillon L. and Gorin C., 2019. Delineating urban areas using building density, *Journal of Urban Economics*.

**Contact**: Clément Gorin, gorinclem@gmail.com

**Description**: This note contains the guidelines to run the scripts `delineation_raster.R` and `delineation_building.R` that compute urban delineations as described in the article referenced above. The first script computes delineation using an aggregated raster of densities as input, while the second uses a disaggregated dataset of densities (e.g. buildings). The scripts are R executables that can be used either interactively or from the command line. The scripts produce the following delineations.

1. Urban areas
2. Urban cores (multiple)
3. Urban areas with cores (multiple)

## Delineations

**Inputs**: 

The script `delineation_raster.R` takes two raster layers as inputs. The raster `density` has a single layer containing the values used to compute the delineations (e.g. building volume, population). The raster `unlivable` contains a binary layer where the pixels that are not eligible for the bootstrap take the value `1`, and `0` otherwise (e.g. rivers, mountains). Each raster should be stored in Geotiff format (i.e. `[file].tif`) with the same extent, resolution and projection. Besides, pixels located outside the area of interest (e.g. country or region borders) must have missing values, as in the example provided.

The script `delineation_building.R` takes a dataset and two raster layers as inputs. The dataset `density` is a csv file with two columns named "density" and "pixel". The first column contains the disaggregated densities (e.g. volume of individual buildings) and the second contains the integer identifier of the raster cell in which the entity is located. The raster `reference` has a single layer containing those integer identifiers. The raster `unlivable` and the constrains on the raster format mentioned above remain the same.

**Outputs**: The scripts produce several rasters whose name include the main parameters. Those names are structured as `[density]d[bandwidth]b[bootstraps]_[type].tif`. Regarding the `type`, delineations for urban areas are noted `*_ur.tif`, urban cores, `*_co.tif`, and urban areas with cores, `*_cc.tif`. Pixels belonging to the same delineation receive the same identifier, while those that are not part of a delineation have missing values. The identifiers are attributed according to the relative size of the delineations (e.g. pixels of the largest delineation take the value `1`). The script also outputs the threshold rasters used to compute urban areas (`*_ut.tif`) and cores (`*_ct.tif`).

Parameters | Description | Type | Default | Script
---|---|---|---|---
`density` | Path to the raster or the csv file containing the values for the criterion used to compute the delineations. See details above. | Character | None | All
`reference` | Path to the raster containing the pixel identifiers. See details above. | Character | None | Building
`unlivable` | Path to the binary raster indicating the pixels that are not eligible for the bootstrap. | Character | None | All
`outdir` | Path to the directory in which output rasters are written. The folder is created if it does not exist. | Character | None | All
`tmpdir` | Path to the directory in which temporary files are written. The folder is created if it does not exist. | Character | None | All
`nboots` | Number of bootstrapped counterfactual densities. If the memory is insufficient, the densities are written to the temporary folder. | Integer | `100` | All
`niter` | Number of iterations to produce urban cores with increasing degrees of density. | Integer | `1` | All
`bandwidth` | Bandwidth (diameter) of the bi-squared smoothing kernel expressed in pixels. This value should depend on the resolution of the input rasters. | Integer | `15` | All
`quantile` | Cut-off percentile for the counterfacual distributions above which the density is considered to be significantly above randomness. | Integer | `95` | All
`replace` | Controls whether the bootstrap for counterfactual densities is performed with (`1`) or without (`0`) replacement. | Integer | `1` | All
`joinall` | Attributes the same identifier to distant delineations that are up to `[parameter value]` pixels apart using 8 connectivity. | Integer | `0` | All
`joinunl` | Attributes the same identifier to distant delineations that are up to `[parameter value]` unlivable pixels apart using 8 connectivity. | Integer | `0` | All
`filter` | Removes the delineations whose number of pixels is smaller than or equal to the parameter value. | Integer | `0` | All
`workers` | Number of cores used for parallel computing. The value `-1` indicates all cores save one. Values exceeding the available cores are collapsed to `-1`. | Integer | `-1` | All
`memory` | Amount of memory used for the computations. The value `-1` indicates all available memory. Values exceeding the available memory are collapsed to `-1`. | Integer |  `-1` | All
`seed` | Bootstrap seed for reproducibility.  L’Ecuyer seeds are used to ensure that the bootstraps performed in parallel are independent. | Integer | `1` | All

**Computations**: When multiple cores are available (i.e. `workers`), computations are performed asynchronously. Forking is used on Unix machines and PSOCK on Windows, in which case the globals are identified and exported to the worker processes. The recommended amount of memory available for computations is 2 GB per worker (i.e. `memory`). Depending on the size of the input data, the number of bootstraps and the operation system, the computations may exceed the available memory. In this case, temporary files are written to the temporary folder (i.e. `tmpdir`), which should be located on a SSD drive to speed up read and write operations. For computational efficiency, all functions are pre-compiled.

**Requirements**: The scripts were written using 4.0.3. The required packages are automatically installed upon execution when they are missing. You may need to run the script with administrator rights the first time.

Packages | Version used | Details
--- | --- | ---
compiler | 4.0.3 | Byte code compiler
data.table | 1.13.2 | Two-dimensional structures 
fst | 0.9.4 | Input-output operations
imager | 0.42.3 | Image processing
raster | 3.4.5 | Raster processing
future.apply | 1.6.0 | Asynchronous processing
memuse | 4.1.0 | Memory management
pacman | 0.5.1 | Package management
optparse | 1.6.6 | Argument management
tictoc | 1.0.0 | Measures execution time

## Usage

The scripts `delineation_raster.R` and `delineation_building.R` can be executed from the Linux or Mac Terminals or the Windows Command Prompt using the commands below. Alternatively, do-file `delineation_wrappers.do` is a wrapper to allows to execute the scripts using Stata. As a example, the `input` folder of the project contains the inputs for the two types of delineation for the Paris area. The instructions below compute the delineations using this input. The different output rasters can be found in the `output` folder of the project.

**Linux or Mac Terminal**

```
# Changes the current directory (adapt path)
cd [...]/delineation

# Delineation from raster
Rscript delineation_raster.R --help
Rscript delineation_raster.R --density=input/raster/density.tif --unlivable=input/raster/unlivable.tif --outdir=output/raster --tmpdir=temporary --nboots=1000 --bandwidth=15 --quantile=95

# Delineation from building
Rscript delineation_building.R --help
Rscript delineation_building.R --density=input/building/density.csv --reference=input/building/reference.tif --unlivable=input/building/unlivable.tif --outdir=output/building --tmpdir=temporary --nboots=100 --bandwidth=15 --quantile=95
```

**Windows Command Prompt**

```
:: Changes the current directory (adapt path)
cd [...]\delineation

:: Declares Rscript as a local or global variable (adapt R version). The Rscript.exe environment variable is not defined by default on Windows
set Rscript="C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe"

:: Delineation from raster
%Rscript% delineation_raster.R --help
%Rscript% delineation_raster.R --density=input\raster\density.tif --unlivable=input\raster\unlivable.tif --outdir=output\raster --tmpdir=temporary --nboots=1000 --bandwidth=15 --quantile=95

:: Delineation from building
%Rscript% delineation_building.R --help
%Rscript% delineation_building.R --density=input\building\density.csv --reference=input\building\reference.tif --unlivable=input\building\unlivable.tif --outdir=output\building --tmpdir=temporary --nboots=100 --bandwidth=15 --quantile=95
```

**Stata wrappers**

```
* Declares Rscript as a global variable (adapt R version)
global Rscript "/usr/local/bin/Rscript"                         // Linux or Mac
global Rscript "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe" // Windows

* Changes the current directory (adapt path)
cd [...]/delineation

* Sources wrappers
quietly do delineation_wrappers.do

* Delineation from raster
delineation_raster, help
delineation_raster, density("input/raster/density.tif") unlivable("input/raster/unlivable.tif") outdir("output/raster") tmpdir("temporary") nboots(1000) bandwidth(15) quantile(95)

* Delineation from buildings
delineation_building, help
delineation_building, density("input/building/density.csv") reference("input/building/reference.tif") unlivable("input/building/unlivable.tif") outdir("output/building") tmpdir("temporary") nboots(100) bandwidth(15) quantile(95)
```

**Windows Power Shell**

```
# Changes the current directory (adapt path)
cd [...]\delineation

# Delineation from raster (adapt R version)
& "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe" delineation_raster.R --help
& "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe" delineation_raster.R --density=input\raster\density.tif --unlivable=input\raster\unlivable.tif --outdir=output\raster --tmpdir=temporary --nboots=1000 --bandwidth=15 --quantile=95

* Delineation from buildings (adapt R version)
& "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe" delineation_building.R --help
& "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe" delineation_building.R --density=input\building\density.csv --reference=input\building\reference.tif --unlivable=input\building\unlivable.tif --outdir=output\building --tmpdir=temporary --nboots=100 --bandwidth=15 --quantile=95
```