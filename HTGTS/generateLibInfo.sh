#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 <xml_file_path> <configType>"
    echo
    echo "This script runs a Docker container for the htgts_pipeline_lts with the specified XML file."
    echo "It processes the XML file to generate libseqInfo.txt and libseqInfo2.txt within the same directory."
    echo
    echo "Arguments:"
    echo "  xml_file_path     Path to the XML (Excel) file to process."
    echo "  configType        Configuration type required by the script."
}

# Check if exactly two arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Missing arguments."
    show_help
    exit 1
fi

# Check for help argument
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    show_help
    exit 0
fi

xmlFile="$1"
configType="$2"
data_folder=$(dirname "$xmlFile")
data_folder_name=$(basename "$xmlFile")

# Check if the XML file exists
if [ ! -f "$xmlFile" ]; then
    echo "Error: The specified XML file does not exist."
    exit 1
fi

# Prepare the Docker command
docker_command="docker run -v \"${data_folder}:/Data\" repbioinfo/htgts_pipeline_lts_v16 python3 /Algorithm/sample_sheetTolibInfo.py /Data/$data_folder_name /Data/libseqInfo.txt /Data/libseqInfo2.txt $configType"

# Echo the command to be executed
echo "Running the following command:"
echo $docker_command

# Execute the Docker command
eval $docker_command

if [ $? -eq 0 ]; then
    echo "The process has completed successfully."
else
    echo "Error: The process encountered an error."
    exit 1
fi
