CWD := $(abspath $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST))))))

ifeq ($(OS),Windows_NT)     # is Windows_NT on XP, 2000, 7, Vista, 10...
    detected_OS := Windows
else
    detected_OS := $(shell uname)  # same as "uname -s"
endif

VIRTUAL_ENV = $(CWD)/cxasap_venv

$(VIRTUAL_ENV):
	python3 -m venv $(VIRTUAL_ENV)

install-quick:
	pip install -e . 
	python $(CWD)/cx_asap/system_files/setup_sys_yaml.py

install-venv: 
	python3 -m venv $(VIRTUAL_ENV)
	. $(VIRTUAL_ENV)/bin/activate && pip install -e .
	python $(CWD)/cx_asap/system_files/setup_sys_yaml.py
test-install:
	python $(CWD)/cx_asap/system_files/test_installation.py

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

test-complete:
	cxasap test
