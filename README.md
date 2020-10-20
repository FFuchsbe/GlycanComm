# GycanComm

These are the scripts used in 

Information transfer in mammalian glycan-based communication

Felix F. Fuchsberger1, Dongyoon Kim1, Marten Kagelmacher1, Robert Wawrzinek1, Christoph Rademacher1
1Max-Planck-Institute of Colloids and Interfaces, Department of Biomolecular Systems, Am Mühlenberg 1 14424 Potsdam, 

Installation time: about 20 min
Software used: Anaconda Navigator, comes with Spider and python version 3.8. Also tested with this version.
Teste on Windows 7/64 bit

How to use the scripts (step by step instruction):

1) install the anaconda navigator: https://www.anaconda.com/products/individual

Once installed you are presented with a bunch of applications in the anaconda navigator.
For the use of the code you will be using "Spyder (Scientific PYthon DEvelopment enviRoment)"

2) Install FlowCytometryTools
Next a package not integrated within anaconda needs to be installed, called FlowCytometryTools more info here:
https://eyurtsev.github.io/FlowCytometryTools/install.html

To install FlowCytometryTools, in the anaconda navigator, launch the first application, the command promt "CMD.exe Promt"
in the application run the following command:
pip install flowcytometrytools

(you run the command by simply typing or copy-pasing the above sentence, then hit enter)

3) Use the code provided here in Github
We provide an example dataset to try the code. Save that dataset anywhere on your computer.

Next launch Spyder and download or copy paste one of the three python scripts of this paper into spyder.
Now save the script you just copied to the location(directory) of your data.
The script  can access data in subdirectories of that location but this has to be specified.
IF you use our example data set, save the scripts into the Example data folder. This has the sub directories: "Channel capacity" and "Correlation"
By default the scripts can access the data in those directories.
To run the scripts hit the green play button "run" or F5.

You should see the calculations in the console on the right, which should take no more than 5 min with the example dataset.
The dataset provided dataset of TNF-a should give you a channel capacity of 1.345 bit and 1.373 bit for discrete and continuous input respectively.

The scripts themselves contain information on how they are to be used on datasets.
  
In case of errors you can contact the first author of the paper.
