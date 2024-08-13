#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 <xml_file_path>"
    echo
    echo "This script runs a Docker container for the htgts_pipeline_lts with the specified XML file."
    echo "It processes the XML file to generate libseqInfo.txt and libseqInfo2.txt within the same directory."
    echo
    echo "Arguments:"
    echo "  xml_file_path   Path to the XML file to process."
}

# Check if at least one argument is provided
if [ $# -ne 1 ]; then
    echo "Error: Missing XML file path."
    show_help
    exit 1
fi

# Check for help argument
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    show_help
    exit 0
fi

xmlFile="$1"
data_folder=$(dirname "$xmlFile")
data_folder_name=$(basename "$xmlFile")

# Check if the XML file exists
if [ ! -f "$xmlFile" ]; then
    echo "Error: The specified XML file does not exist."
    exit 1
fi

# Run the Docker container
docker run -v "${data_folder}:/Data" repbioinfo/htgts_pipeline_lts python3 /Algorithm/sample_sheetTolibInfo.py "/Data/$data_folder_name" /Data/libseqInfo.txt /Data/libseqInfo2.txt

if [ $? -eq 0 ]; then
    echo "The process has completed successfully."
else
    echo "Error: The process encountered an error."
    exit 1
fi
