CWD := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

ifeq ($(OS),Windows_NT)     # is Windows_NT on XP, 2000, 7, Vista, 10...
    detected_OS := Windows
else
    detected_OS := $(shell uname)  # same as "uname -s"
endif

VIRTUAL_ENV := $(CWD)/cxasap_venv
ifeq ($(OS),Windows_NT)
PYTHON ?= python
else
PYTHON ?= python3
endif

# installs packages onto the base system
install-quick:
	$(PYTHON) -m pip install --editable .
	$(PYTHON) $(CWD)/cx_asap/system_files/setup_sys_yaml.py

# uses the base system - must explicitly use python3
test-quick:
	$(PYTHON) $(CWD)/cx_asap/system_files/test_installation.py

# create python3 venv, activate it and use the venv python path (python3)
install-venv:
	$(PYTHON) -m venv $(VIRTUAL_ENV)
ifeq ($(OS),Windows_NT)
	$(VIRTUAL_ENV)/Scripts/python.exe -m pip install --editable . && $(VIRTUAL_ENV)/Scripts/python.exe $(CWD)/cx_asap/system_files/setup_sys_yaml.py
	@echo virtual environment successfully created
	@echo please activate by copying below command into your command prompt:
	@echo $(VIRTUAL_ENV)/Scripts/activate.bat
else
	$(VIRTUAL_ENV)/bin/python -m pip install --editable . && $(VIRTUAL_ENV)/bin/python $(CWD)/cx_asap/system_files/setup_sys_yaml.py
	@echo virtual environment successfully created
	@echo please activate by copying below command into your terminal:
	@echo . $(VIRTUAL_ENV)/bin/activate
endif
	
# run this test ONLY if venv is activated - fails otherwise
test-venv:
ifeq ($(OS),Windows_NT)
	$(VIRTUAL_ENV)/Scripts/python.exe $(CWD)/cx_asap/system_files/test_installation.py
else
	$(VIRTUAL_ENV)/bin/python $(CWD)/cx_asap/system_files/test_installation.py
endif

cxasap-complete:
ifeq ($(strip $(detected_OS)),Darwin) # MacOS
	@echo . $(CWD)/cx_asap/system_files/cxasap-complete.sh >> ~/.bash_profile
	@echo cxasap completion successfully added to ~/.bash_profile
else ifeq ($(strip $(detected_OS)),Linux)
	@echo . $(CWD)/cx_asap/system_files/cxasap-complete.sh >> ~/.bashrc
	@echo cxasap completion successfully added to ~/.bashrc
else
	@echo your operating system is not yet supported for tab completion
endif

test-complete:
	cxasap test

yaml-alias:
ifeq ($(strip $(detected_OS)),Darwin) # MacOS
	@echo alias cxasap_yaml="'"edit $(CWD)/cx_asap/conf.yaml"'" >> ~/.bash_profile
	@echo yaml alias "'"cxasap_yaml"'" successfully added to ~/.bash_profile
else ifeq ($(strip $(detected_OS)),Linux)
	@echo alias cxasap_yaml="'"xed $(CWD)/cx_asap/conf.yaml"'" >> ~/.bashrc
	@echo yaml alias "'"cxasap_yaml"'" successfully added to ~/.bashrc
else
	@echo your operating system is not yet supported for creating a yaml alias
endif