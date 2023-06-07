This is a github repository that holds all the code required to reproduce the results in the Integrative Modeling papers that assess the injury on the cetacean populations of the Gulf of Mexico DeepWater Horizon explosion in 2010. 

These papers were created within the CARMMHA project, where a team of marine mammal health scientists conducted cross-discipline research that included veterinary assessments of managed animals, field assessments with wild populations, and integrative statistical modeling to understand how the DWH oil spill affected Gulf of Mexico marine mammal health.

There are two main papers:

1. A paper focused on bottlenose dolphins in Barataria Bay, where all the animals are assumed to have been exposed to oil. The published version of the paper is at: https://conbio.onlinelibrary.wiley.com/doi/10.1111/cobi.13878
2. A paper focusing on 15 pelagic taxonomic units, including both dolphins :dolphin: and whales :whale:, where the different taxa were exposed to different extents to oil. The published version of the paper now in press at Marine Ecology Progress Series will be added here as soon as it is published online.

Here we describe how you can navigate through this material. We note upfront that while there are two different papers, the underlying population dynamics simulation code that allows us to simulate populations of animals in the presence of oil and in the absence of oil is common to both papers. 

Hence, the results for both papers can  be obtained by running the same code just by using different arguments to the key functions involved in running the population dynamics model and that evaluate the injury metrics.

There are two master files that facilitate access to the material provided:

* :boom: :fire: :arrow_right: :dolphin: [Bottlenose dolphin Barataria Bay paper master file](https://htmlpreview.github.io/?https://github.com/TiagoAMarques/CARMMHApapersSI/blob/master/FolderArchitecture2runCode/BND_ElectronicSupplements.html);

* :boom: :fire: :arrow_right: :whale: [Pelagic paper master file](https://htmlpreview.github.io/?https://github.com/TiagoAMarques/CARMMHApapersSI/blob/master/FolderArchitecture2runCode/Supplement_S1.html). Note that these files were submitted as supplementary material to MEPS (files Supplement S1 to Supplement S8) and they live as such, cast in stone, in the folder [CastInStoneMEPS](https://htmlpreview.github.io/?https://github.com/TiagoAMarques/CARMMHApapersSI/CastInStoneMEPS).

A short paper submitted to Marine Mammal Science, describing the effects on Bottlenose dolphin Barataria Bay population dynamics of a proposed Mississippi river diversion project (Mid-Barataria Sediment Diversion) also uses a slight adaptation of the code required for the main papers. That paper is published here: https://onlinelibrary.wiley.com/doi/10.1111/mms.12930 

The master file to navigate those files is here:

* :droplet: :arrow_right: :dolphin: [Bottlenose dolphin Barataria Bay river diversion paper master file](https://htmlpreview.github.io/?https://github.com/TiagoAMarques/CARMMHApapersSI/blob/master/FolderArchitecture2runCode/Diversion_ElectronicSupplements.html)

Feel free to drop us any comments at tiago.marques at st-andrews.ac.uk
