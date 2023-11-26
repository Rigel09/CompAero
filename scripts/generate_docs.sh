#!/bin/bash

if [[ ! -d "CompAero" ]]; then
	echo "Error: This script must be ran from the top level project directory"
	exit
fi

rm -rf docs/*
pdoc --html CompAero
mv html/CompAero/* docs/
rm -rf html/

echo "Documentation Generation Complete"
