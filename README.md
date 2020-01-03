# CLEAR-SKY DETECTION METHODOLOGY LIBRARY

## What is a clear-sky detection methodology?
A clear-sky detection (CSD) methodology is one that analysies time-series of irradiance data and determines one of two things, (1) whether or not the sky at the time of measurement was free of clouds, or (2) whether or not there was a clear line of sight to the sun at the time of measurement.

## Introduction
Welcome to the CSD methodology library. 
This is a collection of as many CSD methodologies as was found in literature and subsequently coded into Matlab. 
Note that many of the models required some form of interpretation or conversion. 
Each model (found in the 'models' directory) state the citation of the literature article or conference proceedings from where the model was interpreted. 
In the case when the methodology was provided in code, the disclaimer reads "Coded and Converted by..." as opposed to just "Coded by..."

Whilst the CSD methodologies from literature do state their express intentions between type (1) and type (2) methodologies, our investigation determines that they are often not adept at their intended use case.
As such, it is recommended to consider the different methodologies (start by using the Example.m script) in order to get a feel for their strengths and weaknesses. 

![Time series of clear-sky detection from every CSD model, the same graphic as produced by Example.m](https://github.com/JamieMBright/csd-library/raw/master/Example_output.PNG "Example_output.PNG")

## Download
To download or clone this repository, click the above link or [click here to visit the csd-library Github Repository](https://github.com/JamieMBright/csd-library).
Or simply [download the whole master repository](https://github.com/JamieMBright/csd-library/archive/master.zip) as a `.zip` file.

## Input Data
The different types of input into each methodology are defined with the assumption that all data is at 1-minute resolution, and that all data is a single column vector, and also that each of the same index across different variables corresponds to the same time-step.
Note that no methodology needs all of these, however, all are needed should one wish to run all methods. The fewest inputs is two, the most is six. Some need clear-sky irradiance, others do not, some it doesn't matter even it if is a poor clear-sky estimate due to optimisation. 

  `ghi` = Global horizontal irradiance (Wm-2), column vector. 
 
  `ghics` = Clear-sky global horizontal irradiance (Wm-2), column vector. 

  `dni` = Direct normal irradiance (Wm-2), column vector.
 
  `dnics` = Clear-sky direct normal irradiance (Wm-2), column vector. 

  `dif` = Diffuse horizontal irradiance (Wm-2), column vector. 

  `zen` = Solar zenith angle (degrees), column vector.

  `exth` = Horizontal projection of extraterrestrial irradiance (Wm-2), column vector.

  `aod` = Total aerosol extinction at 550nm (dimensionless), column vector.
  
  `LST` =  the local solar time, datevector.

  `plot_figure` = if this variable is defined, a figure is plotted illustrating the outcome of the select CSD method, e.g. 1;
               
               
A summary of the different model's inputs are detailed in the below table. Note that the top 16 models are cloudless-sky (type 1) methodologies, and the bottom 4 are line of sight detection methodologies (type 2). 

| Author | Abbreviation	| `zen`	| `exth`	| `ghi`	| `dni`	| `dif` | `ghics`	| `dnics`	| `difcs`	| `aod` |  `LST`
| ---------------------| ------ |:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|Alia-Martinez et al. (2016) | AliaMartinez | • |  | • |  |  | • |  | | | •
|Batlles | Batlles | • | • | • |  | • |  |  | | |
|Bright et al. (2020) | BrightSunCSDc | • | • | • |  | • | • |  | • | | • 
|Bright et al. (2020) | BrightSunCSDs | • | • | • |  | • | • |  | • | | • 
|Ellis et al. (2018) | Ellis |  |  | • |  |  | • | | | |
|Garcia et al. (2014) | Garcia | • |  | • |  | • |  | | | •| •
|Ineichen (2006) | Ineichen06 | • | • |  | • |  |  |  | | |
|Ineichen et al. (2009) | Ineichen09 | • | • | • |  |  | | |  | |
|Ineichen (2016) | Ineichen16 | • | • | • | • | • |  | | | •|
|Inman et al. (2015) | Inman |  |  | • | • |  | • | • | | |
|Lefevre (2013) | Lefevre | • | • | • |  | • |  |  | | |
|Long and Ackerman (2000) | Long | • |  | • |  | • |  | | | | •
|Perez et al. (1990) | Perez | • |  |  | • | • |  |  | | |
|Polo et al. (2009) | Polo | • |  | • |  |  | • |  | | | •
|Quesada-Ruiz et al. (2015) | Quesada |  |  | • |  |  | • | | | | 
|Reno and Hansen (2016) | Reno |  |  | • |  |  | • |  | | |
|Xie and Liu (2013) | Xie | • | • | • | • |  | • | • | | |
|Zhang et al. (2018) | Zhang |  |  | • |  |  | • |  | | |
|||||||||
|Gueymard (2013)|	Gueymard	|•|	      |	    |•    |	    |	      |•      |	    | |
|Larraneta et al. (2017)|Larraneta|•|	|	    |•    |	    |	      |•	    |     | |
|Ruiz-Arias et al. (2018)|RuizArias|•||•	  |     |•    |       |       |   |  | •	
|Shen et al. (2018)	|Shen	|	|         |•    |•    |     |       |•	   | |•    |
|Zhao et al. (2018)	|Zhao	||	|	    |•    |	    |	      |•	    |     | |

## Output Data
Each CSD model offers a single output which is of identical dimensions to the input data. Pay particular attention to the fact that 1 indicates cloud/not-clear and 0 indicates clear. This can, of course be swapped. 

  csd = clear-sky detection. A value of 1 means that clouds were detected whereas a flag of 0 means that the hour is clear.
  
  **1 = cloud**
  
  **0 = clear**
       
## Example
An example of all the methodologies in action is provided in `Example.m`. 
This script should work first time with no inputs and produce a nice figure to see every model's outputs. 
In the "Example.m" script, sample data is provided from Adelaide Airport irradiance station operated by the Bureau of Meteorology, Australia. 
It is situated at -34.9524 degrees of latitude and 138.5196 degrees of longitude (East of prime meridian) and at an elevation of 2m above sea level. 
The time stamps are not provided. 
You can use this to decide which is the most suitable method for your needs, as some have better applications than others.
For eaxmple, the Zhang model is by far the best at isolating only entirely clear days, something not offered by other models, however, it is hyper-conservative in that it rejects the vast majority of days where clear-sky conditions are also observed.
Contrastingly, the Ellis/Reno/Inman are the most adept in identifying any possible period where it could be considered clear-sky, though one could argue that they are better type 2 methodologies than those advertised as so.

To run all the models:
1. Download and unzip the csd-library
2. Open `Example.m` with Matlab
3. Click `Run` or hit `F5`

You should be able to see the convention for all the input data, and how each model is called in turn.
Once you have decided which model you wish to use, for example if you were to choose the _RuizArias_ methodology, one must only call `csd =RuizArias2018CSD(ghi,dif,zen)` in order to obtain clear-sky detection flags by replacing `ghi`, `dif` and `zen` with your own data.

## Contributors
The literature review in order to identify all the different methodologies of CSD was carried out by Christian A. Gueymard.
The model interpretation and coding was carried out by Jamie M. Bright with support, proofing and testing by David Lingfors. 
We pay particular thanks to Gueymard, Jose Ruiz-Arias, Yu Xie, Javier Antonanzas, Ben Ellis and Wenqi 'Flora' Zhang for providing coded versions or support for their models. 

## Future
This repository will be managed by Jamie M. Bright. Do not hesitate to get in touch should you have new models or corrections to existing ones that you wish to be amended.       

## Usage rights
All models within this library have been interpreted from literature, converted to Matlab from publically available code, or done so with permission of the author.
These methods are freely available to be downloaded and adopted, though check the license before use! A **strict requirement** for using any of the models that may result in publication is appropriate recognition through citation of the original source as stated in each model [1], this repository as the location for obtaining the models [2] and furthermore reference the following work [3] which was the purpose for the conception of this library:

[1] Found at the top of the script of each model script.

[2] Bright, Jamie M. 2018. Clear-sky Detection Library. Github repository, url <https://github.com/JamieMBright/csd-library>.

[3] Gueymard, Christian. A., Bright, Jamie M. and Lingfors, David. 2018. A posteriori clear-sky identification methods in solar irradiance time series: Review and validation. Solar Energy. In Review.
