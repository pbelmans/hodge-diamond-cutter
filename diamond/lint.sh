#!/bin/sh

mkinit --black __init__.py > __init__.py

black *.py
isort --profile black *.py
flake8 --max-line-length 88 *.py
ruff check --ignore=E501,E741 *.py
