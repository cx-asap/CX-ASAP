# Welcome to CX-ASAP!

## Introduction

CX-ASAP is a collection of modules and pipelines designed to automate crystallographic analysis. It is particularly useful for automatically analysing large datasets where experiments have been performed on the same crystal. Consider using this software if you regularly perform experiments such as:

* Variable Temperature Studies
* Variable Pressure Studies
* Positional Mapping Studies 

## Main branch or Dev branch?

The Main branch contains the first official release of CX-ASAP! It is focused on post data-reduction only. While developing this software, however, a suite of other tools were written which you may find useful. These can be found on the Dev branch. Please note, however, that the Dev branch has not been fully tested and requires software that is incompatable with Windows (ie XDS). 

The Main branch is fully compatable with Windows and has been tested more rigorously :)

## Installation

You have two options to install CX-ASAP: either straight onto your computer, or into a virtual environment. It is highly recommended that you choose the virtual environment option if you are familiar with using them (or are happy to learn how to :) ) 

First, navigate in the command line into the CX-ASAP folder (where Makefile is). 

If you want to install CX-ASAP straight onto your computer, execute the below instruction:

`make install-quick`

If you want to install CX-ASAP into a virtual environment (remembering you will need to activate it manually), execute the below instruction:

`make install-venv`

This should also be done every time you re-download the code from github! 

## Additional Installation Options (Linux and Mac only) 

If you would like command-line completion (ie hitting the tab key to auto-complete options), then you will also need to enter the below code which updates your bashrc file (Linux) or bash_profile (Mac):

`make cxasap-complete` 

If you would like an alias to easily edit the conf.yaml file, there is an option to add this into your bashrc file (Linux) or bash_profile (Mac). If you run this installation option, you will be able to type 'cxasap_yaml' into the terminal to easily open the conf.yaml file. This will be easier than finding it and opening it manually. To do this, enter the below command:

`make yaml-alias`

Note that these steps are not necessary, but will make your life easier :) 

## Running the Code

Open your terminal (in your virtual environment if you installed it there) and execute the below command:

`cxasap`

From here, the commandline will guide you through usage of the code.

## Citations

If you include any data or analysis output from CX-ASAP in your publications, please include a citation to the code (LINK TO PAPER HERE)

## Contributing

We welcome contributions from the crystallographic community :) 

See the files included in the 'documentation' folder for information on how to contribute to this package. 

## Troubleshooting

If you receive an error upon executing the code that reads like this...

`Click will abort further execution because Python 3 was configured to use ASCII as encoding for the environment.`

... you will need to change the encoding of your computer. To see what options are available, open your terminal and run the below command:

`locale -a`

Choose an option from the output (we will call it X). Ideally, it should match the country of your computer. Ie, for Australian users, it will be:

`en_AU.utf8`

Note that there are variations and some computers may express it as 'UTF-8'. 

Once you have chosen your encoding, run the below two commands (where X is your encoding):

`export LC_ALL=X`
`export LANG=X`

Try executing CX-ASAP once again and the error should have gone. 

For more information, you can read the below link:

`https://click.palletsprojects.com/en/7.x/python3/` 

If you receive errors while running a module or pipeline, you can view the error log which may give an indication of what went wrong. You can view it in two ways:

1. execute the command 'cxasap errors' - this will print it out into the terminal
2. navigate to cx_asap/error_logs/error_output.txt - you will be able to open this file in a text editor

## Requirements

Required Crystallography Programs (note that these must be in your path and commandline executable):

* Platon
* shredCIF
* SHELXL
* xprep (dev branch only)
* XDS (dev branch only)

Python requirements are listed in requirements.txt and will automatically be installed upon installation of CX-ASAP

## All the different CX-ASAPs explained

* CX-ASAP = name of the software package
* cxasap = what you type into the commandline to execute the code 
* cx_asap = directory containing the conf.yaml file
* cx-asap = name of the repo on GitHub

## Authors
* Amy Thompson
* Dr Kate Smith
* Dr Daniel Eriksson
* Dr Jack Clegg
* Dr Jason Price
