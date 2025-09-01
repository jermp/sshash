#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status

# Default values
RELEASE_NAME="sshash-latest"

# Process command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --release-name)
      if [[ -n "$2" ]]; then
        RELEASE_NAME="$2"
        shift # past argument
        shift # past value
      else
        echo "Error: --release-name requires a value"
        exit 1
      fi
      ;;
    *)
      echo "Error: Unknown parameter passed: $1"
      exit 1
      ;;
  esac
done

ZIP_OUTPUT_FILE="$RELEASE_NAME.zip"
TAR_GZ_OUTPUT_FILE="$RELEASE_NAME.tar.gz"

git submodule init
git submodule update --recursive

echo "Creating zip archive..."
git archive --format=zip --output="$ZIP_OUTPUT_FILE" --prefix="$RELEASE_NAME/" HEAD

echo "Creating tar.gz archive..."
git archive --format=tar.gz --output="$TAR_GZ_OUTPUT_FILE" --prefix="$RELEASE_NAME/" HEAD

echo "Archives created: $ZIP_OUTPUT_FILE and $TAR_GZ_OUTPUT_FILE"

echo "Computing sha256 of '$TAR_GZ_OUTPUT_FILE'..."
shasum -a 256 $TAR_GZ_OUTPUT_FILE

exit 0
