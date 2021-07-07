# Delineating urban areas using building density (21.06.28)

**Reference**: De Bellefon M.P., Combes P.P., Duranton G., Gobillon L. and Gorin C., 2019. Delineating urban areas using building density, *Journal of Urban Economics*.

**Contact**: Clément Gorin, gorinclem@gmail.com

**Description**: This note contains the guidelines to run the script `delineation.R` that computes urban delineations as described in the article referenced above. The script is a R executable that can be used either interactively or from the command line. The script produces the following delineations.

1. Urban areas
2. Urban cores
3. Urban areas with cores

## Delineations

**Inputs**: The script takes two raster layers as inputs. The raster `density` has a single layer containing the values used to compute the delineations (e.g. building volume, population). The raster `unlivable` contains a binary layer where the pixels that are not eligible for the bootstrap take the value `1`, and `0` otherwise (e.g. rivers, mountains). Each raster should be stored in Geotiff format (i.e. `[file].tif`) with the same extent, resolution and projection. Besides, pixels located outside the area of interest (e.g. country or region borders) must have missing values, as in the example provided.

**Outputs**: The script produces several rasters whose name include the main parameters. Those names are structured as `[density]d[bandwidth]b[bootstraps]_[type].tif`. Regarding the `type`, delineations for urban areas are noted `*_ur.tif`, urban cores, `*_co.tif`, and urban areas with cores, `*_cc.tif`. Pixels belonging to the same delineation receive the same identifier, while those that are not part of a delineation have missing values. The identifiers are attributed according to the relative size of the delineations (e.g. pixels of the largest delineation take the value `1`). The script also outputs the threshold rasters used to compute urban areas (`*_ut.tif`) and cores (`*_ct.tif`).

Parameters | Description | Type | Default
---|---|---|---
`density` | Path to the raster containing the values for the criterion used to compute the delineations. | Character | None
`unlivable` | Path to the binary raster indicating the pixels that are not eligible for the bootstrap. | Character | None
`outdir` | Path to the directory in which output rasters are written. The folder is created if it does not exist. | Character | None
`tmpdir` | Path to the directory in which temporary files are written. The folder is created if it does not exist. <span style="color:red">**If the folder exists, its content is cleared upon exectution**</span>. | Character | None
`nboots` | Number of bootstrapped counterfactual densities. If the memory is insufficient, the densities are written to the temporary folder. | Integer | `100`
`niter` | Number of iterations to produce urban cores with increasing degrees of density. | Integer | `1`
`bandwidth` | Bandwidth (diameter) of the bi-squared smoothing kernel expressed in pixels. This value should depend on the resolution of the input rasters. | Integer | `15`
`quantile` | Cut-off percentile for the counterfacual distributions above which the density is considered to be significantly above randomness. | Integer | `95`
`replace` | Controls whether the bootstrap for counterfactual densities is performed with (`1`) or without (`0`) replacement. | Integer | `1`
`joinall` | Attributes the same identifier to distant delineations that are up to `[parameter value]` pixels apart using 8 connectivity. | Integer | `0`
`joinunl` | Attributes the same identifier to distant delineations that are up to `[parameter value]` unlivable pixels apart using 8 connectivity. | Integer | `0`
`filter` | Removes the delineations whose number of pixels is smaller than or equal to the parameter value. | Integer | `0`
`workers` | Number of cores used for parallel computing. The value `-1` indicates all cores save one. Values exceeding the available cores are collapsed to `-1`. | Integer | `-1`
`memory` | Amount of memory used for the computations. The value `-1` indicates all available memory. Values exceeding the available memory are collapsed to `-1`. | Integer |  `-1`
`seed` | Bootstrap seed for reproducibility.  L’Ecuyer seeds are used to ensure that the bootstraps performed in parallel are independent. | Integer | `1`

**Computations**: When multiple cores are available (i.e. `workers`), computations are performed asynchronously. Forking is used on Unix machines and PSOCK on Windows, in which case the globals are identified and exported to the worker processes. The recommended amount of memory available for computations is 2 GB per worker (i.e. `memory`). Depending on the size of the input data, the number of bootstraps and the operation system, the computations may exceed the available memory. In this case, temporary files are written to the temporary folder (i.e. `tmpdir`), which should be located on a SSD drive to speed up read and write operations. For computational efficiency, all functions are pre-compiled.

**Requirements**: The required packages are automatically installed upon execution when they are missing. You may need to run the script with administrator rights the first time.

Packages | Version used | Details
--- | --- | ---
R | 4.0.3 | R language
compiler | 4.0.3 | Byte code compiler
fst | 0.9.4 | Input-output operations
imager | 0.42.3 | Image processing
raster | 3.4.5 | Raster processing
future.apply | 1.6.0 | Asynchronous processing
memuse | 4.1.0 | Memory management
pacman | 0.5.1 | Package management
optparse | 1.6.6 | Argument management
tictoc | 1.0.0 | Measures execution time

## Usage

The script `delineation.R` can be executed from the Linux or Mac Terminals or the Windows Command Prompt using the commands below. Alternatively, do-file `delineation_wrapper.do` is a wrapper to allows to execute the script `delineation.R` from Stata. As a example, the `input` folder of the project contains the input rasters (i.e. `density.tif` and `unlivable.tif`) for the Paris area. The instructions below computes delineations using this input data. The different output rasters can be found in the `output` folder of the project.

**Linux or Mac Terminal**

```
# Changes the current directory (adapt path)
cd [...]/delineation

# Displays help
Rscript delineation.R --help

# Executes script
Rscript delineation.R --density=input/density.tif --unlivable=input/unlivable.tif --outdir=output --tmpdir=temporary --nboots=100 --bandwidth=15 --quantile=95
```

**Windows Command Prompt**

```
:: Changes the current directory (adapt path)
cd [...]\delineation

:: Declares Rscript as a local or global variable (adapt R version). The Rscript.exe environment variable is not defined by default on Windows
set Rscript="C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe"

:: Displays help
%Rscript% delineation.R --help

:: Executes script
%Rscript% delineation.R --density=input\density.tif --unlivable=input\unlivable.tif --outdir=output --tmpdir=temporary --nboots=100 --bandwidth=15 --quantile=95
```

**Windows Power Shell**

```
# Changes the current directory (adapt path)
cd [...]\delineation

# Displays help (adapt R version)
& "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe" delineation.R --help

# Executes script (adapt R version)
& "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe" delineation.R --density=input\density.tif --unlivable=input\unlivable.tif --outdir=output --tmpdir=temporary --nboots=100 --bandwidth=15 --quantile=95
```

**Stata wrapper**

```
* Declares Rscript as a global variable (adapt R version)
global Rscript "C:\Program Files\R\R-4.0.3\bin\x64\Rscript.exe"

* Changes the current directory (adapt path)
cd [...]/delineation

* Sources wrapper
quietly do delineation_wrapper.do

* Displays help
delineation, help

* Executes script
delineation, density("input/density.tif") unlivable("input/unlivable.tif") outdir("output") tmpdir("temporary") nboots(100) bandwidth(15) quantile(95)
```